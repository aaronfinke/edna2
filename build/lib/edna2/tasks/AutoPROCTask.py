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
import json
from datetime import datetime
STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

# for the os.chmod
from stat import *

from cctbx import sgtbx

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsXML


logger = UtilsLogging.getLogger()

from edna2.tasks.XDSTasks import XDSTask
from edna2.tasks.CCP4Tasks import AimlessTask
from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults
from edna2.tasks.WaitFileTask import WaitFileTask

class AutoPROCTask(AbstractTask):

    def setFailure(self):
        self._dictInOut["isFailure"] = True
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
    
    def run(self, inData):
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec='seconds')
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.processingPrograms="edna2autoPROC"
        self.processingCommandLine = ""

        self.setLogFileName(f'autoPROC_{self.startDateTimeFormatted}.log')
        self.dataCollectionId = inData.get("dataCollectionId")
        self.tmpdir = None
        directory = None
        template = None
        self.imageNoStart = None
        self.imageNoEnd = None
        pathToStartImage = None
        pathToEndImage = None

        self.doAnom = inData.get("doAnom",False)

        logger.debug("Working directory is {0}".format(self.getWorkingDirectory()))

        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.lowResLimit = inData.get("lowResolutionLimit",None)
        self.highResLimit = inData.get("highResolutionLimit",None)

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
            logger.info("No space group supplied")

        # need both SG and unit cell
        if self.spaceGroup != 0 and self.unitCell is not None:
            try:
                unitCellList = [float(x) for x in self.unitCell.split(',')]
                #if there are zeroes parsed in, need to deal with it
                if 0.0 in unitCellList:
                    raise Exception
                self.unitCell = {"cell_a": unitCellList[0],
                            "cell_b": unitCellList[1],
                            "cell_c": unitCellList[2],
                            "cell_alpha": unitCellList[3],
                            "cell_beta": unitCellList[4],
                            "cell_gamma": unitCellList[5]}
                logger.info("Supplied unit cell is {cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**self.unitCell))
            except:
                logger.debug("could not parse unit cell")
                self.unitCell = None
        else:
            logger.info("No unit cell supplied")


        if self.dataCollectionId is not None:
            identifier = str(self.dataCollectionId)
            dataCollectionWS3VO = UtilsIspyb.findDataCollection(self.dataCollectionId)
            ispybDataCollection = dict(dataCollectionWS3VO)
            logger.debug("ispybDataCollection: {}".format(ispybDataCollection))
            if ispybDataCollection is not None:
                directory = ispybDataCollection.get("imageDirectory")
                if UtilsConfig.isEMBL():
                    template = ispybDataCollection["fileTemplate"].replace("%05d", "#" * 5)
                elif UtilsConfig.isMAXIV():
                    template = ispybDataCollection["fileTemplate"]
                else:
                    template = ispybDataCollection["fileTemplate"].replace("%04d", "####")
                self.imageNoStart = ispybDataCollection["startImageNumber"]
                numImages = ispybDataCollection["numberOfImages"]
                self.imageNoEnd = numImages - self.imageNoStart + 1
                pathToStartImage = os.path.join(directory, template % self.imageNoStart)
                pathToEndImage = os.path.join(directory, template % self.imageNoEnd)
            else:
                identifier = str(int(time.time()))
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

            if self.imageNoEnd - self.imageNoStart < 8:
                #if self.imageNoEnd - self.imageNoStart < -1:
                logger.error("There are fewer than 8 images, aborting")
                self.setFailure()
                return
        
        #make results directory
        self.resultsDirectory = self.getWorkingDirectory() / "results"
        self.resultsDirectory.mkdir(exist_ok=True, parents=True, mode=0o755)

        #make pyarch directory 
        reg = re.compile(r"(?:/gpfs/offline1/visitors/biomax/|/data/visitors/biomax/)")
        pyarchDirectory = re.sub(reg, "/data/staff/ispybstorage/visitors/biomax/", str(self.resultsDirectory))
        self.pyarchDirectory = Path(pyarchDirectory)
        try:
            self.pyarchDirectory.mkdir(exist_ok=True,parents=True, mode=0o755)
        except OSError as e:
            logger.error(f"Error when creating pyarch_dir: {e}")
            self.tmpdir = tempfile.TemporaryDirectory() 
            self.pyarchDirectory = Path(self.tmpdir.name)

        for file in self.resultsDirectory.iterdir():
            pyarchFile = UtilsPath.createPyarchFilePath(file)
            shutil.copy(file,pyarchFile)
        
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

        self.timeStart = datetime.now().isoformat(timespec='seconds')

        if inData.get("dataCollectionId") is not None:
            #set ISPyB to running
            self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine = self.processingCommandLine,
                processingPrograms = self.processingPrograms,
                isAnom = self.doAnom,
                timeStart = self.timeStart)
        
        # Determine pyarch prefix
        if UtilsConfig.isALBA():
            listPrefix = template.split("_")
            self.pyarchPrefix = "ap_{0}_{1}".format("_".join(listPrefix[:-2]),
                                                       listPrefix[-2])
        else:
            listPrefix = template.split("_")
            self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])

        if isH5:
            masterFilePath = os.path.join(directory,
                                self.eiger_template_to_master(template))
        else:
            logger.error("Only supporing HDF5 data at this time. Stopping.")
            self.setFailure()
            return

        #set up command line
        autoPROCSetup = UtilsConfig.get(self,"autoPROCSetup", None)
        autoPROCExecutable = UtilsConfig.get(self,"autoPROCExecutable", "process")
        maxNoProcessors = UtilsConfig.get(self, "maxNoProcessors", None)
        pathToNeggiaPlugin = UtilsConfig.get(self, "pathToNeggiaPlugin", None)
        autoPROCmacro = UtilsConfig.get(self, "macro", None)

        if autoPROCSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + autoPROCSetup + '\n'
        commandLine += " {0}".format(autoPROCExecutable)
        #add flags, if present
        commandLine += " -B -xml -nthreads {0}".format(maxNoProcessors)
        
        commandLine += " -h5 {0}".format(masterFilePath)
        commandLine += " autoPROC_XdsKeyword_LIB={0}".format(pathToNeggiaPlugin) if pathToNeggiaPlugin else ""
        commandLine += " -ANO" if self.doAnom else ""

        if autoPROCmacro is not None:
            for macro in autoPROCmacro.split():
                commandLine += " -M {0}".format(macro)

        if self.spaceGroup != 0 and self.unitCell is not None:
            commandLine += " symm=\"{0}\"".format(self.spaceGroup)
            commandLine += " cell=\"{cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}\"".format(**self.unitCell)

        if self.lowResLimit is not None or self.highResLimit is not None:
            low = self.lowResLimit if self.lowResLimit else 1000.0
            high = self.highResLimit if self.highResLimit else 0.1
            commandLine += " -R {0} {1}".format(low,high)

        config = UtilsConfig.getConfig()
        config.optionxform = str
        aP_config = config["AutoPROCTask"]
        logger.debug(f"{aP_config}")
        for k,v in aP_config.items():
            if k.startswith("autoPROC_"):
                logger.debug(f"autoPROC option: {k}={v}")
                commandLine += " {0}={1}".format(k,v)

        logger.info("autoPROC command is {}".format(commandLine))
        
        try:
            self.runCommandLine(commandLine, listCommand=[])
        except RuntimeError:
            self.setFailure()
            return
        self.endDateTime = datetime.now().isoformat(timespec='seconds')

        ispybXML = self.getWorkingDirectory() / "autoPROC.xml"
        if ispybXML.is_file():
            self.outData["ispybXML"] = str(ispybXML)
        ispybXML_staraniso = self.getWorkingDirectory() / "autoPROC_staraniso.xml"
        if ispybXML_staraniso.is_file():
            self.outData["ispybXML_staraniso"] = str(ispybXML_staraniso)
        
        autoProcContainer = UtilsXML.dictfromXML(ispybXML)
        


        return
        
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

        #no fast_dp.json, task fails
        if not fastDpJson.exists():
            logger.error("fast_dp.json not found- error log below")
            errorLog = self.getErrorLog()
            logger.error(errorLog)
            self.setFailure()
            return

        # Add aimless.log if present
        pathToAimlessLog = self.fastDpResultFiles.get("aimlessLog")
        if pathToAimlessLog.exists():
            pyarchAimlessLog = self.pyarchPrefix + "_aimless.log"
            shutil.copy(pathToAimlessLog, self.resultsDirectory /  pyarchAimlessLog)
            shutil.copy(pathToAimlessLog, self.pyarchDirectory / pyarchAimlessLog)
            self.aimlessData = AimlessTask.extractAimlessResults(pathToAimlessLog)
            logger.debug(f"aimlessData = {self.aimlessData}")
        #extract ISa...
        correctLpResults = XDSTask.parseCorrectLp(inData={"correctLp":self.getWorkingDirectory() / "CORRECT.LP",
                                                          "gxParmXds":self.getWorkingDirectory() / "GXPARM.XDS" })
        logger.debug(f"correctLpResults: {correctLpResults}")
        self.ISa = correctLpResults.get("ISa")
        logger.debug(f"Isa = {self.ISa}")

        with open(fastDpJson,"r") as fp:
            self.fastDpResults = json.load(fp)

        
        # Add XDS_ASCII.HKL if present and gzip it
        pathToXdsAsciiHkl = self.fastDpResultFiles.get("xdsAsciiHkl")
        if pathToXdsAsciiHkl.exists():
            pyarchXdsAsciiHkl = self.pyarchPrefix + "_XDS_ASCII.HKL.gz"
            with open(pathToXdsAsciiHkl, 'rb') as f_in:
                with gzip.open(os.path.join(self.pyarchDirectory, pyarchXdsAsciiHkl), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            if self.resultsDirectory:
                shutil.copy(self.pyarchDirectory / pyarchXdsAsciiHkl, self.resultsDirectory /  pyarchXdsAsciiHkl)
        
        # Add fast_dp.mtz if present and gzip it
        pathToFastDpMtz = self.fastDpResultFiles.get("fastDpMtz")
        if pathToFastDpMtz.exists():
            pyarchFastDpMtz = self.pyarchPrefix + "_fast_dp.mtz"
            shutil.copy(pathToFastDpMtz, self.resultsDirectory /  pyarchFastDpMtz)
            shutil.copy(pathToFastDpMtz, self.pyarchDirectory / pyarchFastDpMtz)
        
        # add fast_dp.log to results/pyarch directory
        shutil.copy(self.getWorkingDirectory() / "fast_dp.log", self.resultsDirectory / "fast_dp.log")
        shutil.copy(self.getWorkingDirectory() / "fast_dp.log", self.pyarchDirectory / "fast_dp.log")

        autoProcResults = self.generateAutoProcResultsContainer(self.programId, self.integrationId, isAnom=False)

        ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=autoProcResults, workingDirectorySuffix="uploadFinal")
        ispybStoreAutoProcResults.execute()
        if self.tmpdir is not None:
            self.tmpdir.cleanup()


        
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


