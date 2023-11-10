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
from iotbx import mtz
from cctbx.miller import build_set
from cctbx.crystal import symmetry as crystal_symmetry

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
                "fast_dpMtzFile": {"type":["string","null"]},
                "checkDataFirst": {"type":["boolean","null"]},
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
        self.mtzFile = inData.get("fast_dpMtzFile",None)
        self.doUploadIspyb = inData.get("doUploadIspyb",False)
        self.waitForFiles = inData.get("waitForFiles",True)
        self.checkDataFirst = inData.get("checkDataFirst",False)
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

        if self.checkDataFirst:
            checkData = self.checkForPhasingDataQuality(mtzFile=self.mtzFile)
            if not checkData:
                logger.error("Data quality insufficient for phasing. Exiting")
                return
            else:
                logger.info("Data quality check passed.")

        # if self.doUploadIspyb:
        #     #set ISPyB to running
        #     self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
        #         dataCollectionId=self.dataCollectionId,
        #         processingCommandLine = self.processingCommandLine,
        #         processingPrograms = self.processingPrograms,
        #         timeStart = self.timeStart)
        

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

        logger.info("fast_ep command is {}".format(commandLine))

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
            pyarchDirectory = Path(self.tmpdir.name)
            self.pyarchDirectory = self.storeDataOnPyarch(pyarchDirectory=pyarchDirectory)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch()

        # autoProcResults = self.generateAutoProcResultsContainer(self.programId, self.integrationId, isAnom=self.anomalous)        
        # if self.doUploadIspyb:
        #     with open(self.resultsDirectory / "fast_dp_ispyb.json","w") as fp:
        #         json.dump(autoProcResults, fp, indent=2,default=lambda o:str(o))
        #     ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=autoProcResults, workingDirectorySuffix="uploadFinal")
        #     ispybStoreAutoProcResults.execute()
            
        outData = self.fastDpResultFiles
        try:
            fastEpData = json.load(open(self.fastDpResultFiles["fastEpJson"],'r'))
            for k,v in fastEpData.items():
                outData[k] = v
        except Exception as e:
            logger.error(f"could not parse fast_ep json output:{e}")

        
        return outData

    @staticmethod
    def checkForPhasingDataQuality(mtzFile) -> bool:
        '''Decide whether to run fast_ep or no based on the actual data based on
        the following criteria:

        completeness > 80% to dmin or 2.0 whichever the lower
        dI / s(dI) > 1.0 if resolution lower than 2.0, > 0.8 if better than
        1.5, smoothly varying in between.'''
        m = mtz.object(mtzFile)

        mas = m.as_miller_arrays()
        data = None

        for ma in mas:
            if not ma.anomalous_flag():
                continue
            data = ma
            break

        if not data:
            logger.error("No anomalous data found.")
            return False

        d_min, d_max = sorted(data.resolution_range())

        if d_min > 2.0:
            logger.info(f"high resolution is > 2.0 Å")
            differences = data.anomalous_differences()
            signal_to_noise = sum(abs(differences.data())) / \
                sum(differences.sigmas())
            completeness = data.completeness()

            if completeness < 0.8:
                logger.warning("Completeness of mtz is less than 0.8")
                return False
            if signal_to_noise < 1.0:
                logger.warning("Overall S/N is less than 1")
                return False
            return True

        else:
            logger.info(f"high resolution is < 2.0 Å")
            data2 = data.resolution_filter(d_min = 2.0)
            differences = data2.anomalous_differences()
            signal_to_noise = sum(abs(differences.data())) / \
                sum(differences.sigmas())
            completeness = data2.completeness()

            if completeness < 0.8:
                logger.warning("Completeness of mtz is less than 0.8")
                return False
            if signal_to_noise < 0.8:
                logger.warning("Overall S/N is less than 1")
                return False
            return True

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




