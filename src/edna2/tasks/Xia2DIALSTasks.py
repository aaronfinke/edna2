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

from cctbx import sgtbx
from datetime import datetime

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsXML


logger = UtilsLogging.getLogger()

from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, ISPyBStoreAutoProcStatus
from edna2.tasks.WaitFileTask import WaitFileTask


STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

class Xia2DialsTask(AbstractTask):
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
        self.dataCollectionId = inData.get("dataCollectionId")
        self.tmpdir = None
        directory = None
        template = None
        pathToStartImage = None
        pathToEndImage = None

        self.doAnom = inData.get("doAnom",False)

        logger.debug("Working directory is {0}".format(self.getWorkingDirectory()))

        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.lowResLimit = inData.get("lowResolutionLimit",None)
        self.highResLimit = inData.get("highResolutionLimit",None)

        self.proteinAcronym = "AUTOMATIC"
        self.sampleName = "DEFAULT"

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
            identifier = str(self.dataCollectionId)
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
            logger.info(f"Starting:ending image numbers: {self.imageNoStart}:{self.imageNoEnd}")

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

        # if inData.get("dataCollectionId") is not None:
        #     #set ISPyB to running
        #     self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
        #         dataCollectionId=self.dataCollectionId,
        #         processingCommandLine = self.processingCommandLine,
        #         processingPrograms = self.processingPrograms,
        #         isAnom = self.doAnom,
        #         timeStart = self.timeStart)
        
        # Determine pyarch prefix
        if UtilsConfig.isALBA():
            listPrefix = template.split("_")
            self.pyarchPrefix = "di_{0}_{1}".format("_".join(listPrefix[:-2]),
                                                       listPrefix[-2])
        else:
            listPrefix = template.split("_")
            self.pyarchPrefix = "di_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])

        if isH5:
            masterFilePath = os.path.join(directory,
                                self.eiger_template_to_master(template))
        else:
            logger.error("Only supporing HDF5 data at this time. Stopping.")
            self.setFailure()
            return

        self.xia2DIALSExecDir: Path = self.getWorkingDirectory() / "Xia2DialsExec_0"
        inc_x = 1
        while self.xia2DIALSExecDir.is_dir():
            self.xia2DIALSExecDir = self.getWorkingDirectory() / "Xia2DialsExec_{0}".format(inc_x)
            inc_x += 1
        self.xia2DIALSExecDir.mkdir(exist_ok=True,parents=True)


        xia2DialsSetup = UtilsConfig.get(self,"xia2DialsSetup", None)
        xia2DialsExecutable = UtilsConfig.get(self,"xia2DialsExecutable", "xia2")
        maxNoProcessors = UtilsConfig.get(self, "maxNoProcessors", os.cpu_count())
        xia2DialsFastMode = distutils.util.strtobool(UtilsConfig.get(self,"xia2DialsFastMode"))

        #prepare nproc, njobs for dials.integrate
        dialsIntegratePhil = self.xia2DIALSExecDir / "dials_integrate.phil"
        cpusPerJob = int(maxNoProcessors) // 8
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
        commandLine += f" cd {self.xia2DIALSExecDir};"
        commandLine += f" {xia2DialsExecutable}"
        #add flags, if present
        commandLine += " pipeline=dials"
        # commandLine += f" working_directory={self.xia2DIALSExecDir}"
        if self.doAnom:
            commandLine += " atom=X"
        commandLine += f" image={masterFilePath}:{self.imageNoStart}:{self.imageNoEnd}"
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
        
        # self.logToIspyb(self.integrationId,
        #             'Indexing', 'Launched', 'Xia2Dials started')
        
        logger.info("xia2Dials command is {}".format(commandLine))

        #try this timeout
        pathToFinished = str(self.xia2DIALSExecDir / "DataFiles/xia2.cif")
        logger.info("Waiting for end: {0}".format(pathToFinished))

        waitFileFinished = WaitFileTask(inData= {
            "file":pathToFinished,
            "expectedSize": 100
        })
        waitFileFinished.run()
        try:
            self.runCommandLine(commandLine, listCommand=[])
        except RuntimeError:
            self.setFailure()
            return
        
        self.endDateTime = datetime.now().isoformat(timespec="seconds")

        return

        for resultFile in Path(self.xia2DIALSExecDir / "DataFiles").glob("*"):
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{resultFile.name}"
            UtilsPath.systemCopyFile(resultFile,targetFile)


        # run xia2.ispyb_json
        xia2JsonIspybTask = Xia2JsonIspybTask(inData={"xia2DialsExecDir":str(self.xia2DIALSExecDir)}, workingDirectorySuffix="final")
        xia2JsonIspybTask.execute()

        xia2JsonFile = xia2JsonIspybTask.outData.get("ispyb_json",None)

        if xia2JsonFile is not None:
            logger.info("ispyb.json successfully created")
        
        xia2AutoProcContainer = self.loadAndFixJsonOutput(xia2JsonFile)

        xia2AutoProcContainer["dataCollectionId"] = self.dataCollectionId
        xia2AutoProcContainer["autoProcProgram"]["autoProcProgramId"] = self.programId
        xia2AutoProcContainer["autoProc"]["autoProcProgramId"] = self.programId
        xia2AutoProcContainer["autoProcIntegration"]["autoProcIntegrationId"] = self.integrationId
        xia2AutoProcContainer["autoProcIntegration"]["autoProcProgramId"] = self.programId
        xia2AutoProcContainer["autoProcScalingHasInt"] = {
            "autoProcIntegrationId" : self.integrationId
        }

        self.logToIspyb(self.integrationId,
                    'Indexing', 'Successful', 'Xia2Dials finished')


        # ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=xia2AutoProcContainer, workingDirectorySuffix="uploadFinal")
        # ispybStoreAutoProcResults.execute()

            

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
    def loadAndFixJsonOutput(self,jsonFile):
        """fixes some of the output from the xia2 JSON output."""
        autoProcContainer = {
            "autoProcScalingStatistics" : []
        }
        with open(jsonFile,'r') as fp:
            jsonFile = json.load(fp)
        
        autoProcContainer["autoProcProgram"] = jsonFile["AutoProcProgramContainer"]["AutoProcProgram"]
        autoProcContainer["autoProcProgram"]["processingPrograms"] = autoProcContainer["autoProcProgram"].pop("processingProgram")
        autoProcContainer["autoProc"] = jsonFile["AutoProc"]
        autoProcContainer["autoProc"]["refinedCellA"] = autoProcContainer["autoProc"].pop("refinedCell_a")
        autoProcContainer["autoProc"]["refinedCellB"] = autoProcContainer["autoProc"].pop("refinedCell_b")
        autoProcContainer["autoProc"]["refinedCellC"] = autoProcContainer["autoProc"].pop("refinedCell_c")
        autoProcContainer["autoProc"]["refinedCellAlpha"] = autoProcContainer["autoProc"].pop("refinedCell_alpha")
        autoProcContainer["autoProc"]["refinedCellBeta"] = autoProcContainer["autoProc"].pop("refinedCell_beta")
        autoProcContainer["autoProc"]["refinedCellGamma"] = autoProcContainer["autoProc"].pop("refinedCell_gamma")

        autoProcContainer["autoProcIntegration"] = jsonFile["AutoProcScalingContainer"]["AutoProcIntegrationContainer"][0]["AutoProcIntegration"]
        autoProcContainer["autoProcIntegration"]["cellA"] = autoProcContainer["autoProcIntegration"].pop("cell_a")
        autoProcContainer["autoProcIntegration"]["cellB"] = autoProcContainer["autoProcIntegration"].pop("cell_b")
        autoProcContainer["autoProcIntegration"]["cellC"] = autoProcContainer["autoProcIntegration"].pop("cell_c")
        autoProcContainer["autoProcIntegration"]["cellAlpha"] = autoProcContainer["autoProcIntegration"].pop("cell_alpha")
        autoProcContainer["autoProcIntegration"]["cellBeta"] = autoProcContainer["autoProcIntegration"].pop("cell_beta")
        autoProcContainer["autoProcIntegration"]["cellGamma"] = autoProcContainer["autoProcIntegration"].pop("cell_gamma")
        autoProcContainer["autoProcIntegration"]["refinedXbeam"] = autoProcContainer["autoProcIntegration"].pop("refinedXBeam")
        autoProcContainer["autoProcIntegration"]["refinedYbeam"] = autoProcContainer["autoProcIntegration"].pop("refinedYBeam")
        

        autoProcContainer["autoProcScalingStatistics"] = jsonFile["AutoProcScalingContainer"]["AutoProcScalingStatistics"]

        autoProcContainer["autoProcScaling"] = jsonFile["AutoProcScalingContainer"]["AutoProcScaling"]


        return autoProcContainer



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
    

        
        
        






