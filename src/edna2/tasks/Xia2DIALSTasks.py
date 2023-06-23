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
__date__ = "25/04/2023"

import os
import shutil
import tempfile
import numpy as np
from pathlib import Path
import gzip
import time
import re
import distutils
import json
import socket
import traceback

from cctbx import sgtbx
from datetime import datetime

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsXML


logger = UtilsLogging.getLogger()

from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, ISPyBStoreAutoProcStatus, createIntegrationId
from edna2.tasks.WaitFileTask import WaitFileTask


STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

class Xia2DialsTask(AbstractTask):
    def setFailure(self):
        self._dictInOut["isFailure"] = True
        if self.dataCollectionId:
            if self.integrationId is not None and self.programId is not None:
                ISPyBStoreAutoProcResults.setIspybToFailed(
                    dataCollectionId=self.dataCollectionId,
                    autoProcProgramId=self.programId, 
                    autoProcIntegrationId=self.integrationId, 
                    processingCommandLine=self.processingCommandLine, 
                    processingPrograms=self.processingPrograms, 
                    isAnom=False, 
                    timeStart=self.startDateTime, 
                    timeEnd=datetime.now().isoformat(timespec="seconds")
                )
                self.logToIspyb(self.integrationId,
                    'Indexing', 'Failed', 'AutoPROC ended')

    def run(self, inData):
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.processingPrograms="xia2DIALS"
        self.processingCommandLine = ""

        self.setLogFileName(f"xia2DIALS_{self.startDateTimeFormatted}.log")
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.tmpdir = None
        directory = None
        template = None
        pathToStartImage = None
        pathToEndImage = None

        self.doAnom = inData.get("doAnom",False)
        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.lowResLimit = inData.get("lowResolutionLimit",None)
        self.highResLimit = inData.get("highResolutionLimit",None)

        self.proteinAcronym = "AUTOMATIC"
        self.sampleName = "DEFAULT"

        logger.info("Xia2DIALS processing started")
        logger.info(f"Running on {socket.gethostname()}")
        try:
            logger.info(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass


        if self.spaceGroup != 0:
            try:
                spaceGroupInfo = sgtbx.space_group_info(self.spaceGroup).symbol_and_number()
                self.spaceGroupString = spaceGroupInfo.split("No. ")[0][:-2]
                self.spaceGroupNumber = int(spaceGroupInfo.split("No. ")[1][:-1])
                logger.info("Supplied space group is {}, number {}".format(self.spaceGroupString, self.spaceGroupNumber))
            except:
                logger.debug("Could not parse space group")
                self.spaceGroupNumber = 0
        else:
            self.spaceGroupNumber = 0
            self.spaceGroupString = ""            
            logger.info("No space group supplied")

        # need both SG and unit cell
        if self.spaceGroup != 0 and self.unitCell is not None:
            try:
                unitCellList = [float(x) for x in self.unitCell.split(",")]
                #if there are zeroes parsed in, need to deal with it
                if 0.0 in unitCellList:
                    raise Exception
                self.unitCell = {
                            "cell_a": unitCellList[0],
                            "cell_b": unitCellList[1],
                            "cell_c": unitCellList[2],
                            "cell_alpha": unitCellList[3],
                            "cell_beta": unitCellList[4],
                            "cell_gamma": unitCellList[5]
                            }
                logger.info("Supplied unit cell is {cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**self.unitCell))
            except:
                logger.debug("could not parse unit cell")
                self.unitCell = None
        else:
            logger.info("No unit cell supplied")


        if self.dataCollectionId is not None:
            self.integrationId = None
            self.programId = None
            dataCollectionWS3VO = UtilsIspyb.findDataCollection(self.dataCollectionId)
            if dataCollectionWS3VO is not None:
                ispybDataCollection = dict(dataCollectionWS3VO)
                logger.debug("ispybDataCollection: {}".format(ispybDataCollection))
                directory = ispybDataCollection.get("imageDirectory")
                if UtilsConfig.isEMBL():
                    template = ispybDataCollection["fileTemplate"].replace("%05d", "#" * 5)
                elif UtilsConfig.isMAXIV():
                    template = ispybDataCollection["fileTemplate"]
                else:
                    template = ispybDataCollection["fileTemplate"].replace("%04d", "####")
                self.imageNoStart = inData.get("imageNoStart", ispybDataCollection["startImageNumber"])
                numImages = ispybDataCollection["numberOfImages"]
                self.imageNoEnd = inData.get("imageNoEnd", (numImages - self.imageNoStart + 1))  
                pathToStartImage = os.path.join(directory, template % self.imageNoStart)
                pathToEndImage = os.path.join(directory, template % self.imageNoEnd)
            else:
                directory = self.dataInput.dirN.value
                template = self.dataInput.templateN.value
                self.imageNoStart = self.dataInput.fromN.value
                self.imageNoEnd = self.dataInput.toN.value
                if UtilsConfig.isEMBL():
                    fileTemplate = template.replace("#####", "%05d")
                else:
                    fileTemplate = template.replace("####", "%04d")

                pathToStartImage = os.path.join(directory, fileTemplate % self.imageNoStart)
                pathToEndImage = os.path.join(directory, fileTemplate % self.imageNoEnd)

            # Determine pyarch prefix
            if UtilsConfig.isALBA():
                listPrefix = template.split("_")
                self.pyarchPrefix = "di_{0}_{1}".format("_".join(listPrefix[:-2]),
                                                        listPrefix[-2])
            else:
                listPrefix = template.split("_")
                self.pyarchPrefix = "di_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])

            if self.imageNoEnd - self.imageNoStart < 8:
                #if self.imageNoEnd - self.imageNoStart < -1:
                logger.error("There are fewer than 8 images, aborting")
                self.setFailure()
                return
            
            logger.info(f"Starting:ending image numbers: {self.imageNoStart}:{self.imageNoEnd}")
            logger.info(f"dataCollectionId: {self.dataCollectionId}")

            proteinAcronym, sampleName = UtilsIspyb.getProteinAcronymAndSampleNameFromDataCollectionId(self.dataCollectionId)
            if proteinAcronym is not None and sampleName is not None:
                # only alphanumerics and underscores are allowed
                proteinAcronym_corrected = re.sub(r"\W", '_', proteinAcronym)
                sampleName_corrected = re.sub(r"\W", '_', sampleName)
                self.proteinAcronym = proteinAcronym_corrected
                self.sampleName = sampleName_corrected
            logger.info(f"Protein Acronym:{self.proteinAcronym}, sample name:{self.sampleName}")


        #make results directory
        self.resultsDirectory = self.getWorkingDirectory() / "results"
        self.resultsDirectory.mkdir(exist_ok=True, parents=True, mode=0o755)

        #make pyarch directory 
        if inData.get("test",False):
            self.tmpdir = tempfile.TemporaryDirectory() 
            self.pyarchDirectory = Path(self.tmpdir.name)
        else:
            reg = re.compile(r"(?:/gpfs/offline1/visitors/biomax/|/data/visitors/biomax/)")
            pyarchDirectory = re.sub(reg, "/data/staff/ispybstorage/visitors/biomax/", str(self.resultsDirectory))
            self.pyarchDirectory = Path(pyarchDirectory)
            try:
                self.pyarchDirectory.mkdir(exist_ok=True,parents=True, mode=0o755)
                logger.info(f"Created pyarch directory: {self.pyarchDirectory}")
            except OSError as e:
                logger.error(f"Error when creating pyarch_dir: {e}")
                self.tmpdir = tempfile.TemporaryDirectory() 
                self.pyarchDirectory = Path(self.tmpdir.name)
        
        isH5 = False
        if any(beamline in pathToStartImage for beamline in ["id23eh1", "id29"]):
            minSizeFirst = 6000000
            minSizeLast = 6000000
        elif any(beamline in pathToStartImage for beamline in ["id23eh2", "id30a1"]):
            minSizeFirst = 2000000
            minSizeLast = 2000000
        elif any(beamline in pathToStartImage for beamline in ["id30a3"]):
            minSizeFirst = 100000
            minSizeLast = 100000
            pathToStartImage = os.path.join(directory,
                                            self.eiger_template_to_image(template, self.imageNoStart))
            pathToEndImage = os.path.join(directory,
                                          self.eiger_template_to_image(template, self.imageNoEnd))
            isH5 = True
        elif UtilsConfig.isMAXIV():
            minSizeFirst = 100000
            minSizeLast = 100000
            pathToStartImage = os.path.join(directory,
                                            self.eiger_template_to_image(template, self.imageNoStart))
            pathToEndImage = os.path.join(directory,
                                          self.eiger_template_to_image(template, self.imageNoEnd))
            isH5 = True
        else:
            minSizeFirst = 1000000
            minSizeLast = 1000000        

        logger.info("Waiting for start image: {0}".format(pathToStartImage))
        waitFileFirst = WaitFileTask(inData= {
            "file":pathToStartImage,
            "expectedSize": minSizeFirst
        })
        waitFileFirst.execute()
        if waitFileFirst.outData["timedOut"]:
            logger.warning("Timeout after {0:d} seconds waiting for the first image {1}!".format(waitFileFirst.outData["timeOut"], pathToStartImage))
        
        logger.info("Waiting for end image: {0}".format(pathToEndImage))
        waitFileLast = WaitFileTask(inData= {
            "file":pathToEndImage,
            "expectedSize": minSizeLast
        })
        waitFileLast.execute()
        if waitFileLast.outData["timedOut"]:
            logger.warning("Timeout after {0:d} seconds waiting for the last image {1}!".format(waitFileLast.outData["timeOut"], pathToEndImage))

        self.timeStart = datetime.now().isoformat(timespec="seconds")
        

        if isH5:
            masterFilePath = os.path.join(directory,
                                self.eiger_template_to_master(template))
        else:
            logger.error("Only supporing HDF5 data at this time. Stopping.")
            self.setFailure()
            return

        xia2DIALSExecinData = {
            "masterFilePath" : masterFilePath,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd" : self.imageNoEnd,
            "proteinAcronym" : self.proteinAcronym,
            "sampleName" : self.sampleName,
            "spaceGroupNumber" : self.spaceGroupNumber,
            "spaceGroupString" : self.spaceGroupString,
            "unitCell" : self.unitCell,
            "lowResLimit" : self.lowResLimit,
            "highResLimit" : self.highResLimit,
            "doAnom" : self.doAnom
        }
        try: 
            self.integrationId,self.programId = createIntegrationId(
                                    self, 
                                    "Creating integration ID", 
                                    isAnom=self.doAnom)
        except Exception as e:
            logger.error("Could not get integration ID: \n{0}".format(traceback.format_exc(e)))
        logger.info(f"integrationID: {self.integrationId}, programId: {self.programId}")

        self.logToIspyb(self.integrationId,
                            'Indexing', 'Launched', 'Xia2 started')

        timeOut = inData.get("timeOut",None)
        if timeOut is None:
            timeOut = UtilsConfig.get(self,"timeOut",3600)
        xia2DIALSExec = Xia2DialsExecTask(inData=xia2DIALSExecinData, workingDirectorySuffix="0")
        xia2DIALSExec.setTimeout(timeOut)
        xia2DIALSExec.execute()
        if xia2DIALSExec.isFailure():
            if xia2DIALSExec["timeoutExit"] == True:
                logger.error(f"Operation timed out after {timeOut} s.")
            self.setFailure()
            return
        self.timeEnd = datetime.now().isoformat(timespec="seconds")
        self.logToIspyb(self.integrationId,
            'Indexing', 'Successful', 'Xia2Dials finished')

        xia2DIALSExecDir = Path(xia2DIALSExec.outData["workingDirectory"])
        logger.debug(f"Working directory is {xia2DIALSExecDir}")

        for resultFile in Path(xia2DIALSExecDir / "DataFiles").glob("*"):
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{resultFile.name}"
            UtilsPath.systemCopyFile(resultFile,targetFile)
        
        for logFile in Path(xia2DIALSExecDir / "LogFiles").glob("*"):
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{logFile.name}"
            UtilsPath.systemCopyFile(logFile,targetFile)

        xia2Html = xia2DIALSExecDir /  "xia2.html"
        if xia2Html.exists():
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{xia2Html.name}"
            UtilsPath.systemCopyFile(xia2Html,targetFile)
        
        xia2Txt = xia2DIALSExecDir /  "xia2.txt"
        if xia2Txt.exists():
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{xia2Txt.name}"
            UtilsPath.systemCopyFile(xia2Txt,targetFile)

        for file in self.resultsDirectory.glob("*"):
            targetFile = self.pyarchDirectory / file.name
            UtilsPath.systemCopyFile(file,targetFile)


        # run xia2.ispyb_json
        xia2JsonIspybTask = Xia2JsonIspybTask(inData={"xia2DialsExecDir":str(xia2DIALSExecDir)}, workingDirectorySuffix="final")
        xia2JsonIspybTask.execute()
        xia2JsonFile = xia2JsonIspybTask.outData.get("ispyb_json",None)

        if xia2JsonFile is not None:
            logger.info("ispyb.json successfully created")
            xia2AutoProcContainer = self.loadAndFixJsonOutput(xia2JsonFile)

            if self.dataCollectionId is not None:
                ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=xia2AutoProcContainer, workingDirectorySuffix="uploadFinal")
                ispybStoreAutoProcResults.execute()

        if inData.get("test",False):
            self.tmpdir.cleanup()
            


    def eiger_template_to_master(self, fmt):
        if UtilsConfig.isMAXIV():
            fmt_string = fmt.replace("%06d", "master")
        else:
            fmt_string = fmt.replace("####", "1_master")
        return fmt_string

    def eiger_template_to_image(self, fmt, num):
        import math
        fileNumber = int(math.ceil(num / 100.0))
        if UtilsConfig.isMAXIV():
            fmt_string = fmt.replace("%06d", "data_%06d" % fileNumber)
        else:
            fmt_string = fmt.replace("####", "1_data_%06d" % fileNumber)
        return fmt_string.format(num)


    def logToIspyb(self, integrationId, step, status, comments=""):
        # hack in the event we could not create an integration ID
        if integrationId is None:
            logger.error('could not log to ispyb: no integration id')
            return
        
        statusInput = {
            "dataCollectionId": self.dataCollectionId,
            "autoProcIntegration" : {
                "autoProcIntegrationId": integrationId,
            },
            "autoProcProgram": {
            },
            "autoProcStatus": {
                "autoProcIntegrationId": integrationId,
                "step":  step,
                "status": status,
                "comments": comments,
                "bltimeStamp": datetime.now().isoformat(timespec='seconds'),
            }
        }
        autoprocStatus = ISPyBStoreAutoProcStatus(inData=statusInput, workingDirectorySuffix="")
        autoprocStatus.execute()

    def loadAndFixJsonOutput(self,jsonFile, trunc_len=256):
        """fixes some of the output from the xia2 JSON output
        for submission to ISPyB."""
        autoProcContainer = {
            "dataCollectionId" : self.dataCollectionId
        }
        with open(jsonFile,'r') as fp:
            jsonFile = json.load(fp)

        autoProcContainer["autoProcProgram"] = jsonFile["AutoProcProgramContainer"]["AutoProcProgram"]
        autoProcContainer["autoProcProgram"]["autoProcProgramId"] = self.programId

        autoProcContainer["autoProcProgram"]["processingPrograms"] = autoProcContainer["autoProcProgram"].pop("processingProgram")
        autoProcContainer["autoProcProgram"]["processingCommandLine"] = ""
        autoProcContainer["autoProcProgram"]["processingStatus"] = "SUCCESS"
        autoProcContainer["autoProcProgram"]["processingPrograms"] = "xia2DIALS"
        autoProcContainer["autoProcProgram"]["processingStartTime"] = self.timeStart
        autoProcContainer["autoProcProgram"]["processingEndTime"] = self.timeEnd

        autoProcContainer["autoProc"] = jsonFile["AutoProc"]
        autoProcContainer["autoProc"]["autoProcProgramId"] = self.programId

        autoProcContainer["autoProc"]["refinedCellA"] = autoProcContainer["autoProc"].pop("refinedCell_a")
        autoProcContainer["autoProc"]["refinedCellB"] = autoProcContainer["autoProc"].pop("refinedCell_b")
        autoProcContainer["autoProc"]["refinedCellC"] = autoProcContainer["autoProc"].pop("refinedCell_c")
        autoProcContainer["autoProc"]["refinedCellAlpha"] = autoProcContainer["autoProc"].pop("refinedCell_alpha")
        autoProcContainer["autoProc"]["refinedCellBeta"] = autoProcContainer["autoProc"].pop("refinedCell_beta")
        autoProcContainer["autoProc"]["refinedCellGamma"] = autoProcContainer["autoProc"].pop("refinedCell_gamma")

        autoProcContainer["autoProcIntegration"] = jsonFile["AutoProcScalingContainer"]["AutoProcIntegrationContainer"][0]["AutoProcIntegration"]
        autoProcContainer["autoProcIntegration"]["autoProcIntegrationId"] = self.integrationId
        autoProcContainer["autoProcIntegration"]["autoProcProgramId"] = self.programId
        autoProcContainer["autoProcIntegration"]["cellA"] = autoProcContainer["autoProcIntegration"].pop("cell_a")
        autoProcContainer["autoProcIntegration"]["cellB"] = autoProcContainer["autoProcIntegration"].pop("cell_b")
        autoProcContainer["autoProcIntegration"]["cellC"] = autoProcContainer["autoProcIntegration"].pop("cell_c")
        autoProcContainer["autoProcIntegration"]["cellAlpha"] = autoProcContainer["autoProcIntegration"].pop("cell_alpha")
        autoProcContainer["autoProcIntegration"]["cellBeta"] = autoProcContainer["autoProcIntegration"].pop("cell_beta")
        autoProcContainer["autoProcIntegration"]["cellGamma"] = autoProcContainer["autoProcIntegration"].pop("cell_gamma")
        autoProcContainer["autoProcIntegration"]["refinedXbeam"] = autoProcContainer["autoProcIntegration"].pop("refinedXBeam")
        autoProcContainer["autoProcIntegration"]["refinedYbeam"] = autoProcContainer["autoProcIntegration"].pop("refinedYBeam")
        autoProcContainer["autoProcIntegration"]["anomalous"] = self.doAnom
        autoProcContainer["autoProcIntegration"]["dataCollectionId"] = self.dataCollectionId


        # round up some values...
        for k,v in autoProcContainer["autoProcIntegration"].items():
            if isinstance(v,float):
                if any(x for x in ["cellA","cellB","cellC"] if x in k):
                    autoProcContainer["autoProcIntegration"][k] = round(v,3)
                else:
                    autoProcContainer["autoProcIntegration"][k] = round(v,2)

        autoProcContainer["autoProcScalingHasInt"] = {
            "autoProcIntegrationId" : self.integrationId
        }
        autoProcContainer["autoProcScalingStatistics"] = jsonFile["AutoProcScalingContainer"]["AutoProcScalingStatistics"]

        try:
            Isa = self.getIsa()
        except:
            logger.error("Could not get Isa value from logs.")
            Isa = None


        for shell in autoProcContainer["autoProcScalingStatistics"]:
            shell["rmerge"] = shell.pop("rMerge") * 100
            shell["ccAno"] = shell.pop("ccAnomalous")
            shell["meanIoverSigI"] = shell.pop("meanIOverSigI")
            shell["rmeasAllIplusIminus"] = shell.pop("rMeasAllIPlusIMinus") * 100
            shell["rmeasWithinIplusIminus"] = shell.pop("rMeasWithinIPlusIMinus") * 100
            shell["rpimWithinIplusIminus"] = shell.pop("rPimWithinIPlusIMinus") * 100
            shell["rpimAllIplusIminus"] = shell.pop("rPimAllIPlusIMinus") * 100
            if shell["scalingStatisticsType"] == "overall":
                shell["isa"] = round(float(Isa),2)

            for k,v in shell.items():
                if isinstance(v,float):
                    shell[k] = round(v,2)


        autoProcAttachmentContainerList = []
        for file in self.pyarchDirectory.iterdir():
            attachmentContainer = {
                "file" : file,
            }
            autoProcAttachmentContainerList.append(attachmentContainer)

        autoProcContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList


        return autoProcContainer
    
    def getIsa(self):
        """extract the last estimated Isa value from _SCALE.log file."""
        Isa = None
        for scaleLogFile in self.resultsDirectory.glob("*SCALE.log"):
            with open(scaleLogFile,"r") as fp:
                for line in fp:
                    if "estimated I/sigma asymptotic limit" in line:
                        Isa = line.split()[-1]
        return Isa

class Xia2DialsExecTask(AbstractTask):
    def run(self, inData):
        outData = {}
        logger.debug(f"working directory is {self.getWorkingDirectory()}")
        self.masterFilePath = inData["masterFilePath"]
        self.imageNoStart = inData["imageNoStart"]
        self.imageNoEnd = inData["imageNoEnd"]
        self.proteinAcronym = inData["proteinAcronym"]
        self.sampleName = inData["sampleName"]

        self.spaceGroupNumber = inData.get("spaceGroupNumber",0)
        self.spaceGroupString = inData.get("spaceGroupString","")

        self.unitCell = inData.get("unitCell",None)
        self.lowResLimit = inData.get("lowResolutionLimit",None)
        self.highResLimit = inData.get("highResolutionLimit",None)
        self.doAnom = inData.get("doAnom",False)

        outData["workingDirectory"] = str(self.getWorkingDirectory())

        xia2DialsSetup = UtilsConfig.get("Xia2DialsTask","xia2DialsSetup", None)
        logger.debug(f"xia2DialsSetup: {xia2DialsSetup}")
        xia2DialsExecutable = UtilsConfig.get("Xia2DialsTask","xia2DialsExecutable", "xia2")
        maxNoProcessors = UtilsConfig.get("Xia2DialsTask", "maxNoProcessors", os.cpu_count())
        xia2DialsFastMode = distutils.util.strtobool(UtilsConfig.get("Xia2DialsTask","xia2DialsFastMode"))

        #prepare nproc, njobs for dials.integrate
        dialsIntegratePhil = self.getWorkingDirectory() / "dials_integrate.phil"
        cpusPerJob = 8
        numJobs = 4
        with open(dialsIntegratePhil,'w') as fp:
            fp.write(f"""integration {{
    block {{
        size = Auto
        units = *degrees radians frames
    }}
    mp {{
        nproc={cpusPerJob}
        njobs={numJobs}
    }}
}}""")
        
        #set up command line
        if xia2DialsSetup:
            commandLine = f". {xia2DialsSetup} \n"
        else:
            commandLine = ""
        commandLine += f" {xia2DialsExecutable}"
        #add flags, if present
        commandLine += " pipeline=dials"
        if self.doAnom:
            commandLine += " atom=X"
        commandLine += f" image={self.masterFilePath}:{self.imageNoStart}:{self.imageNoEnd}"
        if maxNoProcessors:
            commandLine += f" nproc={int(maxNoProcessors)}"

        commandLine += f" project={self.proteinAcronym} crystal={self.sampleName}"

        if xia2DialsFastMode:
            commandLine += " dials.fast_mode=True"
        if self.spaceGroupNumber != 0:
            commandLine += f" space_group={self.spaceGroupString}"
            commandLine += " unit_cell={cell_a},{cell_b},{cell_c},{cell_alpha},{cell_beta},{cell_gamma}".format(**self.unitCell)

        if self.lowResLimit is not None or self.highResLimit is not None:
            low = self.lowResLimit if self.lowResLimit else 1000.0
            high = self.highResLimit if self.highResLimit else 0.1
            commandLine += f" resolution.d_min={low}"
            commandLine += f" resolution.d_max={high}"
        commandLine += f" integrate.mosaic=new dials.integrate.phil_file={dialsIntegratePhil}"
        logger.info("xia2Dials command is {}".format(commandLine))

        try:
            self.runCommandLine(commandLine, listCommand=[])
        except RuntimeError:
            self.setFailure()
            return
        return outData

        

class Xia2JsonIspybTask(AbstractTask):
    def run(self, inData):
        xia2DialsExecDir = inData["xia2DialsExecDir"]
        xia2DialsSetup = UtilsConfig.get("Xia2DialsTask","xia2DialsSetup", None)

        if xia2DialsSetup:
            commandLine = f". {xia2DialsSetup} \n"
        else:
            commandLine = ""

        commandLine += f" cd {xia2DialsExecDir}; \n"
        commandLine += f" xia2.ispyb_json"

        logger.info("xia2.ispyb_json command is {}".format(commandLine))
        outData = {
            "ispyb_json" : None
        }
        try:
            self.runCommandLine(commandLine, listCommand=[])
        except RuntimeError:
            self.setFailure()
            return outData
        out_file = xia2DialsExecDir + "/ispyb.json"
        outData["ispyb_json"] = str(out_file)
        return outData
    

        
        
        






