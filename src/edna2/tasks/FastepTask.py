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
import gzip
import time
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

class FastSADPhasingTask(AbstractTask):

    def setFailure(self):
        self._dictInOut["isFailure"] = True
        # if self.doUploadIspyb:
            # if self.integrationId is not None and self.programId is not None:
            #     ISPyBStoreAutoProcResults.setIspybToFailed(
            #         dataCollectionId=self.dataCollectionId,
            #         autoProcProgramId=self.programId, 
            #         autoProcIntegrationId=self.integrationId, 
            #         processingCommandLine=self.processingCommandLine, 
            #         processingPrograms=self.processingPrograms, 
            #         isAnom=False, 
            #         timeStart=self.startDateTime, 
            #         timeEnd=datetime.now().isoformat(timespec='seconds')
            #     )
    def getInDataSchema(self):
        return {
            "type": "object",
            "properties":
            {
                "dataCollectionId": {"type":["integer","null"]},
                "onlineAutoProcessing": {"type":["boolean","null"]},
                "doUploadIspyb": {"type":["boolean","null"]},
                "workingDirectory": {"type":["string","null"]},
                "mtzFile": {"type":["string","null"]}
            }
        }

    def run(self, inData):
        UtilsLogging.addLocalFileHandler(logger, self.getWorkingDirectory()/"EDNA_fastdp.log")

        self.timeStart = time.perf_counter()
        self.tmpdir = None
        self.startDateTime =  datetime.now().isoformat(timespec='seconds')
        self.processingPrograms="fast_ep"
        self.processingCommandLine = ""
        self.setLogFileName('fastEp.log')
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        self.mtzFile = inData.get("mtzFile",None)
        self.doUploadIspyb = inData.get("doUploadIspyb",False)
        self.waitForFiles = inData.get("waitForFiles",True)
        self.integrationId, self.programId = None, None
        outData = {}



        logger.info("fast_ep auto phasing started")
        logger.info(f"Running on {socket.gethostname()}")

        try:
            logger.info(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        self.resultsDirectory = Path(self.getWorkingDirectory() / "results")
        self.resultsDirectory.mkdir(parents=True, exist_ok=True)

        self.timeStart = datetime.now().isoformat(timespec='seconds')

        
        if self.doUploadIspyb:
            #set ISPyB to running
            self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine = self.processingCommandLine,
                processingPrograms = self.processingPrograms,
                timeStart = self.timeStart)
        

        #set up command line
        fastEpSetup = UtilsConfig.get(self,"fastEpSetup", None)
        fastEpExecutable = UtilsConfig.get(self,"fastEpExecutable", "fast_ep")
        ncpu = UtilsConfig.get(self, "ncpu", None)
        nmachine = UtilsConfig.get(self, "nmachines", None)

        if fastEpSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + fastEpSetup + '\n'
        commandLine += " {0}".format(fastEpExecutable)
        commandLine += " data={0}".format(self.mtzFile)
        commandLine += " cpu={0}".format(ncpu) if ncpu else ""
        commandLine += " machines={0}".format(nmachine) if nmachine else ""
        #add flags, if present

        logger.info("fastdp command is {}".format(commandLine))

        if self.onlineAutoProcessing:
            returnCode = self.submitCommandLine(commandLine, partition=UtilsConfig.get(self,"slurm_partition",None),ignoreErrors=False,jobName="EDNA2_fastep")
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
            "resultMtz": self.getWorkingDirectory() / "sad.mtz",
            "fastEpLog": self.getWorkingDirectory() / "fast_ep.log",
            "fastEpReport": self.getWorkingDirectory() / "fastep_report.html",
            "cootScript": self.getWorkingDirectory() / "coot.sh",
            "shelxcOutput" : self.getWorkingDirectory() / "shelxc.log",
            "shelxdOutput" : self.getWorkingDirectory() / "sad_fa.lst",
            "shelxdHeavyAtoms": self.getWorkingDirectory() / "sad_fa.pdb",
            "shelxeOutput": self.getWorkingDirectory() / "sad.lst",
            "fastEpJson": self.getWorkingDirectory() / "fast_ep_data.json"
        }
        fastEpJson = self.getWorkingDirectory() / "fast_dp.json"

        #copy to results directory
        for k,v in self.fastDpResultFiles.items():
            resultDirPath = self.resultsDirectory / v.name
            try:
                shutil.copy(v,resultDirPath)
            except Exception as e:
                logger.error(f"Couldn't copy file {v}: {e}")

        self.resultFilePaths = list(self.resultsDirectory.iterdir())



        if inData.get("test",False):
            self.tmpdir = tempfile.TemporaryDirectory() 
            self.pyarchDirectory = Path(self.tmpdir.name)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch()

        # autoProcResults = self.generateAutoProcResultsContainer(self.programId, self.integrationId, isAnom=self.anomalous)        
        # if self.doUploadIspyb:
        #     with open(self.resultsDirectory / "fast_dp_ispyb.json","w") as fp:
        #         json.dump(autoProcResults, fp, indent=2,default=lambda o:str(o))
        #     ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=autoProcResults, workingDirectorySuffix="uploadFinal")
        #     ispybStoreAutoProcResults.execute()
            
        outData = self.fastDpResultFiles
        
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
        #create paths on Pyarch
        if pyarchDirectory is None:
            pyarchDirectory = UtilsPath.createPyarchFilePath(self.resultFilePaths[0]).parent
            if not pyarchDirectory.exists():
                pyarchDirectory.mkdir(parents=True, exist_ok=True, mode=0o755)
                logger.debug(f"pyarchDirectory: {pyarchDirectory}")
        for resultFile in [f for f in self.resultFilePaths if f.exists()]:
            resultFilePyarchPath = UtilsPath.createPyarchFilePath(resultFile)
            try:
                logger.info(f"Copying {resultFile} to pyarch directory")
                shutil.copy(resultFile,resultFilePyarchPath)
            except Exception as e:
                logger.warning(f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}")
                logger.warning(e)
        return pyarchDirectory


