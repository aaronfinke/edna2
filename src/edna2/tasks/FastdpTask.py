# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ["A. Finke"]
__license__ = "MIT"
__date__ = "20/01/2023"

import os
import shutil
import tempfile
import numpy as np
from pathlib import Path
import math
import gzip
import time
import re
import json
from datetime import datetime
import socket
STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

# for the os.chmod
from stat import *

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsCCTBX
from edna2.utils import UtilsImage

logger = UtilsLogging.getLogger()

from edna2.tasks.XDSTasks import XDSTask
from edna2.tasks.CCP4Tasks import AimlessTask
from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults
from edna2.tasks.WaitFileTask import WaitFileTask

class FastdpTask(AbstractTask):

    def setFailure(self):
        self._dictInOut["isFailure"] = True
        if self.doUploadIspyb:
            if self.integrationId is not None and self.programId is not None:
                ISPyBStoreAutoProcResults.setIspybToFailed(
                    dataCollectionId=self.dataCollectionId,
                    autoProcProgramId=self.programId, 
                    autoProcIntegrationId=self.integrationId, 
                    processingCommandLine=self.processingCommandLine, 
                    processingPrograms=self.processingPrograms, 
                    isAnom=False, 
                    timeStart=self.startDateTime, 
                    timeEnd=datetime.now().isoformat(timespec='seconds')
                )
    def getInDataSchema(self):
        return {
            "type": "object",
            "properties":
            {
                "dataCollectionId": {"type":["integer","null"]},
                "onlineAutoProcessing": {"type":["boolean","null"]},
                "waitForFiles": {"type":["boolean","null"]},
                "doUploadIspyb": {"type":["boolean","null"]},
                "masterFilePath": {"type":["string","null"]},
                "spaceGroup": {"type":["integer","string"]},
                "unitCell": {"type":["string","null"]},
                "residues": {"type":["integer","null"]},
                "anomalous": {"type":["boolean","null"]},
                "imageNoStart": {"type":["integer","null"]},
                "imageNoEnd": {"type":["integer","null"]},
                "workingDirectory": {"type":["string","null"]},
                
            }
        }

    def run(self, inData):
        UtilsLogging.addLocalFileHandler(logger, self.getWorkingDirectory()/"EDNA_fastdp.log")

        self.timeStart = time.perf_counter()
        self.tmpdir = None
        self.startDateTime =  datetime.now().isoformat(timespec='seconds')
        self.processingPrograms="fast_dp"
        self.processingCommandLine = ""
        self.lowRes = 50
        self.setLogFileName('fastDp.log')
        self.tmpdir = None
        self.imageDirectory = None
        self.fileTemplate = None
        pathToStartImage = None
        pathToEndImage = None
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.anomalous = inData.get("anomalous",False)
        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        self.imageNoStart = inData.get("imageNoStart",None)
        self.imageNoEnd = inData.get("imageNoEnd",None)
        self.masterFilePath = inData.get("masterFilePath",None)
        self.doUploadIspyb = inData.get("doUploadIspyb",False)
        self.waitForFiles = inData.get("waitForFiles",True)
        self.integrationId, self.programId = None, None
        outData = {}



        logger.info("fastdp auto processing started")
        logger.info(f"Running on {socket.gethostname()}")

        try:
            logger.info(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        #set up SG and unit cell
        self.spaceGroupNumber, self.spaceGroupString = UtilsCCTBX.parseSpaceGroup(self.spaceGroup)

        # set up unit cell
        if self.unitCell is not None:
            self.unitCell = UtilsCCTBX.parseUnitCell(self.unitCell)
        else:
            logger.info("No unit cell supplied")

        # get masterfile name
        if self.masterFilePath is None:
            if self.dataCollectionId:
                self.masterFilePath = UtilsIspyb.getXDSMasterFilePath(self.dataCollectionId)
                if self.masterFilePath is None or not self.masterFilePath.exists():
                    logger.error("dataCollectionId could not return master file path, exiting.")
                    self.setFailure()
                    return

            else:
                logger.error("No dataCollectionId or masterfile, exiting.")
                self.setFailure()
                return 

        # now we have masterfile name, need number of images and first/last file
        dataCollectionWS3VO = None
        if self.imageNoStart is None or self.imageNoEnd is None:
            if self.dataCollectionId:
                try:
                    dataCollectionWS3VO = UtilsIspyb.findDataCollection(self.dataCollectionId)
                    self.imageNoStart = dataCollectionWS3VO.startImageNumber
                    numImages = dataCollectionWS3VO.numberOfImages
                    self.imageNoEnd = numImages - self.imageNoStart + 1
                except:
                    logger.error("Could not access number of images from ISPyB")
                    self.imageNoStart = 1
                    numImages = UtilsImage.getNumberOfImages(self.masterFilePath)
                    self.imageNoEnd = numImages - self.imageNoStart + 1
            else:
                self.imageNoStart = 1
                numImages = UtilsImage.getNumberOfImages(self.masterFilePath)
                self.imageNoEnd = numImages - self.imageNoStart + 1
        
        if self.imageNoEnd - self.imageNoStart < 8:
            #if self.imageNoEnd - self.imageNoStart < -1:
            logger.error("There are fewer than 8 images, aborting")
            self.setFailure()
            return
        
        dataH5ImageList = UtilsImage.generateDataFileListFromH5Master(self.masterFilePath)
        pathToStartImage = dataH5ImageList[0]
        pathToEndImage = dataH5ImageList[-1]
                
        listPrefix = dataCollectionWS3VO.fileTemplate.split("_") if dataCollectionWS3VO else Path(self.masterFilePath).name.split("_")

        #generate pyarch prefix
        if UtilsConfig.isALBA():
            self.pyarchPrefix = "ap_{0}_{1}".format("_".join(listPrefix[:-2]),
                                                       listPrefix[-2])
        elif UtilsConfig.isMAXIV():
            self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])
        else:
            if len(listPrefix) > 2:
                self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])
            elif len(listPrefix) > 1:
                self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[:-2], listPrefix[-2])
            else:
                self.pyarchPrefix = "ap_{0}_run".format(listPrefix[0])


        if self.waitForFiles:
            logger.info("Waiting for start image: {0}".format(pathToStartImage))
            waitFileFirst = WaitFileTask(inData= {
                "file":pathToStartImage,
                "expectedSize": 100000
            })
            waitFileFirst.execute()
            if waitFileFirst.outData["timedOut"]:
                logger.warning("Timeout after {0:d} seconds waiting for the first image {1}!".format(waitFileFirst.outData["timeOut"], pathToStartImage))
            
            logger.info("Waiting for end image: {0}".format(pathToEndImage))
            waitFileLast = WaitFileTask(inData= {
                "file":pathToEndImage,
                "expectedSize": 100000
            })
            waitFileLast.execute()
            if waitFileLast.outData["timedOut"]:
                logger.warning("Timeout after {0:d} seconds waiting for the last image {1}!".format(waitFileLast.outData["timeOut"], pathToEndImage))

        self.resultsDirectory = Path(self.getWorkingDirectory() / "results")
        self.resultsDirectory.mkdir(parents=True, exist_ok=True)

        self.timeStart = datetime.now().isoformat(timespec='seconds')

        
        if self.doUploadIspyb:
            #set ISPyB to running
            self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine = self.processingCommandLine,
                processingPrograms = self.processingPrograms,
                isAnom = self.anomalous,
                timeStart = self.timeStart)
        

        #set up command line
        fastdpSetup = UtilsConfig.get(self,"fastdpSetup", None)
        fastdpExecutable = UtilsConfig.get(self,"fastdpExecutable", "fast_dp")
        maxNumJobs = UtilsConfig.get(self, "maxNumJobs", None)
        numJobs = UtilsConfig.get(self, "numJobs", None)
        numCores = UtilsConfig.get(self, "numCores", None)
        pathToNeggiaPlugin = UtilsConfig.get(self, "pathToNeggiaPlugin", None)
        highResolutionLimit = inData.get("highResolutionLimit", None)
        atom = inData.get("atom", None)
        if atom is None and self.anomalous:
            atom = "X"
        beamX = inData.get("beamX",None)
        beamY = inData.get("beamY",None)

        if fastdpSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + fastdpSetup + '\n'
        commandLine += " {0}".format(fastdpExecutable)
        #add flags, if present
        commandLine += " -J {0}".format(maxNumJobs) if maxNumJobs else ""
        commandLine += " -j {0}".format(numJobs) if numJobs else ""
        commandLine += " -k {0}".format(numCores) if numCores else ""
        commandLine += " -R {0}".format(self.lowRes) if self.lowRes else ""
        commandLine += " -r {0}".format(highResolutionLimit) if highResolutionLimit else ""
        commandLine += " -s {0}".format(self.spaceGroupNumber) if self.spaceGroupNumber else ""
        commandLine += " -c \"{cell_a},{cell_b},{cell_c},{cell_alpha},{cell_beta},{cell_gamma}\"".format(**self.unitCell) if self.unitCell else ""
        commandLine += " -a {0}".format(atom) if atom else ""
        commandLine += " -b {0},{1}".format(beamX, beamY) if beamX and beamY else ""
        commandLine += " -l {0}".format(pathToNeggiaPlugin) if pathToNeggiaPlugin else ""
        commandLine += " {0}".format(self.masterFilePath)

        logger.info("fastdp command is {}".format(commandLine))

        if self.onlineAutoProcessing:
            returnCode = self.submitCommandLine(commandLine, partition=UtilsConfig.get(self,"slurm_partition",None),ignoreErrors=False,jobName="EDNA2_fastdp")
            if returnCode != 0:
                self.setFailure()
                return
        else:
            try:
                self.runCommandLine(commandLine, listCommand=[])
            except RuntimeError:
                self.setFailure()
                return
        
        self.endDateTime = datetime.now().isoformat(timespec='seconds')
        
        self.fastDpResultFiles = {
            "correctLp": self.getWorkingDirectory() / "CORRECT.LP",
            "gxParmXds": self.getWorkingDirectory() / "GXPARM.XDS",
            "xdsAsciiHkl" : self.getWorkingDirectory() / "XDS_ASCII.HKL",
            "aimlessLog" : self.getWorkingDirectory() / "aimless.log",
            "fastDpMtz" : self.getWorkingDirectory() / "fast_dp.mtz",
            "fastDpLog" : self.getLogFileName(),
            "fastDpJson" : self.getWorkingDirectory() / "fast_dp.json"
        }
        fastDpJson = self.getWorkingDirectory() / "fast_dp.json"

        # Add aimless.log if present
        pathToAimlessLog = self.fastDpResultFiles.get("aimlessLog")
        if pathToAimlessLog.exists():
            pyarchAimlessLog = self.pyarchPrefix + "_aimless.log"
            shutil.copy(pathToAimlessLog, self.resultsDirectory / pyarchAimlessLog)
            self.aimlessData = AimlessTask.extractAimlessResults(pathToAimlessLog)
            # logger.debug(f"aimlessData = {self.aimlessData}")
        #extract ISa...
        correctLpResults = XDSTask.parseCorrectLp(inData={"correctLp":self.getWorkingDirectory() / "CORRECT.LP",
                                                          "gxParmXds":self.getWorkingDirectory() / "GXPARM.XDS" })
        # logger.debug(f"correctLpResults: {correctLpResults}")
        self.ISa = correctLpResults.get("ISa")
        logger.debug(f"Isa = {self.ISa}")

        with open(fastDpJson,"r") as fp:
            self.fastDpResults = json.load(fp)

        # Add XDS_ASCII.HKL if present and gzip it
        pathToXdsAsciiHkl = self.fastDpResultFiles.get("xdsAsciiHkl")
        if pathToXdsAsciiHkl.exists():
            pyarchXdsAsciiHkl = self.pyarchPrefix + "_XDS_ASCII.HKL.gz"
            with open(pathToXdsAsciiHkl, 'rb') as f_in:
                with gzip.open(os.path.join(self.resultsDirectory, pyarchXdsAsciiHkl), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        # Add fast_dp.mtz if present and gzip it
        pathToFastDpMtz = self.fastDpResultFiles.get("fastDpMtz")
        if pathToFastDpMtz.exists():
            pyarchFastDpMtz = self.pyarchPrefix + "_fast_dp.mtz"
            shutil.copy(pathToFastDpMtz, self.resultsDirectory /  pyarchFastDpMtz)
        
        # add fast_dp.log to results/pyarch directory
        pyarchFastDpLog = self.pyarchPrefix + "_fast_dp.log"
        shutil.copy(self.getWorkingDirectory() / "fast_dp.log", self.resultsDirectory / pyarchFastDpLog)

        self.resultFilePaths = list(self.resultsDirectory.iterdir())

        if inData.get("test",False):
            self.tmpdir = tempfile.TemporaryDirectory() 
            pyarchDirectory = Path(self.tmpdir.name)
            self.pyarchDirectory = self.storeDataOnPyarch(pyarchDirectory=pyarchDirectory)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch()

        autoProcResults = self.generateAutoProcResultsContainer(self.programId, self.integrationId, isAnom=self.anomalous)        
        if self.doUploadIspyb:
            with open(self.resultsDirectory / "fast_dp_ispyb.json","w") as fp:
                json.dump(autoProcResults, fp, indent=2,default=lambda o:str(o))
            ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=autoProcResults, workingDirectorySuffix="uploadFinal")
            ispybStoreAutoProcResults.execute()
            
        if self.tmpdir is not None:
            self.tmpdir.cleanup()
        
        outData = autoProcResults
        outData["anomalous"] = self.if_anomalous_signal(pathToAimlessLog, threshold=1.0)
        if outData["anomalous"]:
            logger.info("Significant anomalous signal found.")

        return outData


        
    def generateAutoProcResultsContainer(self, programId, integrationId, isAnom):
        autoProcResultsContainer = {
            "dataCollectionId": self.dataCollectionId
        }

        autoProcProgramContainer = {
            "autoProcProgramId" : programId,
            "processingCommandLine" : self.processingCommandLine,
            "processingPrograms" : self.processingPrograms,
            "processingStatus" : "SUCCESS",
            "processingStartTime" : self.startDateTime,
            "processingEndTime" : self.endDateTime,
        }
        autoProcResultsContainer["autoProcProgram"] = autoProcProgramContainer

        autoProcContainer = {
            "autoProcProgramId" : programId,
            "spaceGroup" : self.fastDpResults.get("spacegroup"),
            "refinedCellA" : self.fastDpResults.get("unit_cell")[0],
            "refinedCellB" : self.fastDpResults.get("unit_cell")[1],
            "refinedCellC" : self.fastDpResults.get("unit_cell")[2],
            "refinedCellAlpha" : self.fastDpResults.get("unit_cell")[3],
            "refinedCellBeta" : self.fastDpResults.get("unit_cell")[4],
            "refinedCellGamma" : self.fastDpResults.get("unit_cell")[5],
        }
        autoProcResultsContainer["autoProc"] = autoProcContainer

        autoProcAttachmentContainerList = []
        for file in self.pyarchDirectory.iterdir():
            attachmentContainer = {
                "file" : file,
            }
            autoProcAttachmentContainerList.append(attachmentContainer)

        autoProcResultsContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList

        autoProcScalingStatisticsContainer = self.fastDpJsonToISPyBScalingStatistics(self.fastDpResults, aimlessResults=self.aimlessData, isAnom=False)
        autoProcResultsContainer["autoProcScalingStatistics"] = autoProcScalingStatisticsContainer

        if self.fastDpResultFiles["gxParmXds"].exists():
            xdsRerun = XDSTask.parseCorrectLp(inData=self.fastDpResultFiles)

            autoProcIntegrationContainer = {
                "autoProcIntegrationId" : integrationId,
                "autoProcProgramId" : programId,
                "startImageNumber" : self.imageNoStart,
                "endImageNumber" : self.imageNoEnd,
                "refinedDetectorDistance" : xdsRerun.get("refinedDiffractionParams").get("crystal_to_detector_distance"),
                "refinedXbeam" : xdsRerun.get("refinedDiffractionParams").get("direct_beam_detector_coordinates")[0],
                "refinedYbeam" : xdsRerun.get("refinedDiffractionParams").get("direct_beam_detector_coordinates")[1],
                "rotationAxisX" : xdsRerun.get("gxparmData").get("rot")[0],
                "rotationAxisY" : xdsRerun.get("gxparmData").get("rot")[1],
                "rotationAxisZ" : xdsRerun.get("gxparmData").get("rot")[2],
                "beamVectorX" : xdsRerun.get("gxparmData").get("beam")[0],
                "beamVectorY" : xdsRerun.get("gxparmData").get("beam")[1],
                "beamVectorZ" : xdsRerun.get("gxparmData").get("beam")[2],
                "cellA" : xdsRerun.get("refinedDiffractionParams").get("cell_a"),
                "cellB" : xdsRerun.get("refinedDiffractionParams").get("cell_b"),
                "cellC" : xdsRerun.get("refinedDiffractionParams").get("cell_c"),
                "cellAlpha" : xdsRerun.get("refinedDiffractionParams").get("cell_alpha"),
                "cellBeta" : xdsRerun.get("refinedDiffractionParams").get("cell_beta"),
                "cellGamma" : xdsRerun.get("refinedDiffractionParams").get("cell_gamma"),
                "anomalous" : isAnom,
                "dataCollectionId" : self.dataCollectionId,
            }
            autoProcResultsContainer["autoProcIntegration"] = autoProcIntegrationContainer


        return autoProcResultsContainer

    def storeDataOnPyarch(self, pyarchDirectory=None):
        # create paths on Pyarch
        if pyarchDirectory is None:
            pyarchDirectory = UtilsPath.createPyarchFilePath(
                self.resultFilePaths[0]
            ).parent
            if not pyarchDirectory.exists():
                pyarchDirectory.mkdir(parents=True, exist_ok=True, mode=0o755)
                logger.debug(f"pyarchDirectory: {pyarchDirectory}")
            for resultFile in [f for f in self.resultFilePaths if f.exists()]:
                resultFilePyarchPath = UtilsPath.createPyarchFilePath(resultFile)
                try:
                    logger.info(f"Copying {resultFile} to pyarch directory")
                    shutil.copy(resultFile, resultFilePyarchPath)
                except Exception as e:
                    logger.warning(
                        f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}"
                    )
                    logger.warning(e)
        else:
            for resultFile in [f for f in self.resultFilePaths if f.exists()]:
                try:
                    logger.info(f"Copying {resultFile} to pyarch directory")
                    resultFilePyarchPath = pyarchDirectory / Path(resultFile).name
                    shutil.copy(resultFile, resultFilePyarchPath)
                except Exception as e:
                    logger.warning(
                        f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}"
                    )
                    logger.warning(e)
                
        return pyarchDirectory


    def fastDpJsonToISPyBScalingStatistics(self, fastDpResults, aimlessResults=None, isAnom=False):
        autoProcScalingContainer = []
        for shell,result in fastDpResults["scaling_statistics"].items():
            resultsShell = {}
            resultsShell["scalingStatisticsType"] = shell
            resultsShell["resolutionLimitLow"] = result.get("res_lim_low" ,None)
            resultsShell["resolutionLimitHigh"] = result.get("res_lim_high" ,None)
            resultsShell["rmerge"] = result.get("r_merge",None) * 100
            resultsShell["rmeasAllIplusIminus"] = result.get("r_meas_all_iplusi_minus",None) * 100
            resultsShell["nTotalObservations"] = result.get("n_tot_obs" ,None)
            resultsShell["nTotalUniqueObservations"] = result.get("n_tot_unique_obs" ,None)
            resultsShell["meanIoverSigI"] = result.get("mean_i_sig_i" ,None)
            resultsShell["completeness"] = result.get("completeness" ,None)
            resultsShell["multiplicity"] = result.get("multiplicity" ,None)
            resultsShell["anomalousCompleteness"] = result.get("anom_completeness" ,None)
            resultsShell["anomalousMultiplicity"] = result.get("anom_multiplicity" ,None)
            resultsShell["anomalous"] = isAnom
            resultsShell["ccHalf"] = result.get("cc_half",None)
            resultsShell["ccAno"] = result.get("cc_anom",None) * 100
            if self.ISa and shell == "overall":
                resultsShell["isa"] = self.ISa

            autoProcScalingContainer.append(resultsShell)
        if aimlessResults is not None:
            for shell,result in aimlessResults.items():
                f = next(item for item in autoProcScalingContainer if item["scalingStatisticsType"] == shell)
                f["rpimWithinIplusIminus"] = result.get("rpimWithinIplusIminus") * 100  if result.get("rpimWithinIplusIminus") else None
                f["rpimAllIplusIminus"] = result.get("rpimAllIplusIminus") * 100 if result.get("rpimAllIplusIminus") else None
                f["sigAno"] = result.get("sigAno") 
        # sorting these so EXI will put sigAno/ISa in the right position...
        try:        
            autoProcScalingContainer = [next(item for item in autoProcScalingContainer if item['scalingStatisticsType'] == 'overall'),
                                        next(item for item in autoProcScalingContainer if item['scalingStatisticsType'] == 'innerShell'),
                                        next(item for item in autoProcScalingContainer if item['scalingStatisticsType'] == 'outerShell')]
        except:
            logger.error("autoProcScalingContainer could not be sorted")
            pass
        return autoProcScalingContainer

    def if_anomalous_signal(self, aimless_log, threshold = 1.0):
        """Grab the anomalous CC RCR value and see if it is 
        sufficiently large to run fast_ep. Generally, a value 
        greater than 1 indicates a significant anomalous signal."""
        cc_rcr = 0.0
        try:
            with open(aimless_log,'r') as fp:
                for line in fp:
                    if "$TABLE:  Correlations CC(1/2) within dataset, XDSdataset" in line:
                        while "Overall" not in line: 
                            line = next(fp)
                        cc_rcr = float(line.split()[3])
        except:
            pass
        if cc_rcr >= threshold:
            return True
        else:
            return False
