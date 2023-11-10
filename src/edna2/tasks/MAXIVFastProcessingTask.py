#
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

import json
import os
import socket
import time
from datetime import datetime

STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

# for the os.chmod

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsCCTBX
from textwrap import dedent
import subprocess


logger = UtilsLogging.getLogger()

from edna2.tasks.WaitFileTask import WaitFileTask
from edna2.tasks.Edna2ProcTask import Edna2ProcTask
from edna2.tasks.FastdpTask import FastdpTask
from edna2.tasks.FastSADPhasingTask import FastSADPhasingTask


class MAXIVFastProcessingTask(AbstractTask):
    """
    Runs four autoprocessing pipelines.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "required": ["dataCollectionId"],
            "properties": {
                "dataCollectionId": {"type": ["integer", "null"]},
                "masterFilePath": {"type": ["string", "null"]},
                "imageNoStart": {"type": ["integer", "null"]},
                "imageNoEnd": {"type": ["integer", "null"]},
                "numImages": {"type": ["integer", "null"]},
                "spaceGroup": {"type": ["integer", "string"]},
                "unitCell": {"type": ["string", "null"]},
                "residues": {"type": ["integer", "null"]},
                "anomalous": {"type": ["boolean", "null"]},
                "workingDirectory": {"type": ["string", "null"]},
                "pdb": {"type": ["string", "null"]},
                "test": {"type": ["boolean", "null"]},
                "waitForFiles": {"type": ["boolean", "null"]},
                "doUploadIspyb": {"type": ["boolean", "null"]},
            },
        }

    def run(self, inData):
        self.timeStart = time.perf_counter()
        UtilsLogging.addLocalFileHandler(
            logger, self.getWorkingDirectory() / "MAXIVAutoProcessing.log"
        )
        logger.info("MAX IV Autoprocessing started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        self.anomFlag = False
        outData = {}
        self.startDateTime = datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.tmpdir = None
        self.imageNoStart = inData.get("imageNoStart", None)
        self.imageNoEnd = inData.get("imageNoEnd", None)
        self.numImages = inData.get("numImages", None)
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.masterFilePath = inData.get("masterFilePath", None)
        self.anomalous = inData.get("anomalous", False)
        self.spaceGroup = inData.get("spaceGroup", 0)
        self.unitCell = inData.get("unitCell", None)
        self.residues = inData.get("residues", None)
        self.workingDirectory = inData.get(
            "workingDirectory", self.getWorkingDirectory()
        )
        self.doUploadIspyb = inData.get("doUploadIspyb", True)
        self.waitForFiles = inData.get("waitForFiles", True)

        self.pdb = inData.get("pdb", None)
        self.test = inData.get("test", False)

        self.proteinAcronym = "AUTOMATIC"
        self.sampleName = "DEFAULT"

        try:
            logger.debug(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        # set up SG and unit cell
        self.spaceGroupNumber, self.spaceGroupString = UtilsCCTBX.parseSpaceGroup(
            self.spaceGroup
        )

        # set up unit cell
        if self.unitCell is not None:
            self.unitCell = UtilsCCTBX.parseUnitCell_str(self.unitCell)
        else:
            logger.info("No unit cell supplied")

        # get masterfile name
        if self.masterFilePath is None:
            if self.dataCollectionId:
                self.masterFilePath = UtilsIspyb.getXDSMasterFilePath(
                    self.dataCollectionId
                )
                if self.masterFilePath is None or not self.masterFilePath.exists():
                    logger.error(
                        "dataCollectionId could not return master file path, exiting."
                    )
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
                    dataCollectionWS3VO = UtilsIspyb.findDataCollection(
                        self.dataCollectionId
                    )
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
        elif self.imageNoStart and self.imageNoEnd:
            numImages = self.imageNoEnd - self.imageNoStart + 1
        else:
            self.imageNoStart = 1
            numImages = UtilsImage.getNumberOfImages(self.masterFilePath)
            self.imageNoEnd = numImages - self.imageNoStart + 1

        if numImages < 8:
            # if self.imageNoEnd - self.imageNoStart < -1:
            logger.error("There are fewer than 8 images, aborting")
            self.setFailure()
            return

        if self.waitForFiles:
            dataH5ImageList = UtilsImage.generateDataFileListFromH5Master(
                self.masterFilePath
            )
            pathToStartImage = dataH5ImageList[0]
            pathToEndImage = dataH5ImageList[-1]

            logger.info("Waiting for start image: {0}".format(pathToStartImage))
            waitFileFirst = WaitFileTask(
                inData={"file": pathToStartImage, "expectedSize": 5_000_000},
                workingDirectorySuffix="firstImage",
            )
            waitFileFirst.execute()
            if waitFileFirst.outData["timedOut"]:
                logger.warning(
                    "Timeout after {0:d} seconds waiting for the first image {1}!".format(
                        waitFileFirst.outData["timeOut"], pathToStartImage
                    )
                )

            logger.info("Waiting for end image: {0}".format(pathToEndImage))
            waitFileLast = WaitFileTask(
                inData={"file": pathToEndImage, "expectedSize": 5_000_000},
                workingDirectorySuffix="lastImage",
            )
            waitFileLast.execute()
            if waitFileLast.outData["timedOut"]:
                logger.warning(
                    "Timeout after {0:d} seconds waiting for the last image {1}!".format(
                        waitFileLast.outData["timeOut"], pathToEndImage
                    )
                )

        edna2ProcTask = Edna2ProcTask(
            inData={
                "onlineAutoProcessing": False,
                "dataCollectionId": self.dataCollectionId,
                "masterFilePath": str(self.masterFilePath)
                if self.masterFilePath
                else None,
                "unitCell": self.unitCell,
                "spaceGroup": self.spaceGroup,
                "imageNoStart": self.imageNoStart,
                "imageNoEnd": self.imageNoEnd,
                "anomalous": self.anomalous,
                "waitForFiles": False,
                "doUploadIspyb": True,
                "test": self.test,
                "timeOut": 1800,
            },
            workingDirectorySuffix="0",
        )

        fastDpTask = FastdpTask(
            inData={
                "onlineAutoProcessing": False,
                "dataCollectionId": self.dataCollectionId,
                "masterFilePath": str(self.masterFilePath)
                if self.masterFilePath
                else None,
                "unitCell": self.unitCell,
                "spaceGroup": self.spaceGroup,
                "masterFilePath": self.masterFilePath,
                "imageNoStart": self.imageNoStart,
                "imageNoEnd": self.imageNoEnd,
                "waitForFiles": False,
                "doUploadIspyb": True,
                "anomalous": False,
                "test": self.test,
                "timeOut": 1800,
            },
            workingDirectorySuffix="0",
        )

        # imgQualityDozor.start()
        edna2ProcTask.start()
        fastDpTask.start()

        # imgQualityDozor.join()
        fastDpTask.join()
        edna2ProcTask.join()

        # set the logic for anomalous processing: if both fastdp and edna2proc say it's anomalous,
        # then set it to anomalous. otherwise don't (unless one or the other fails)
        logger.debug("Checking for anomalous signal...")

        if edna2ProcTask.isSuccess() and fastDpTask.isSuccess():
            logger.debug("EDNA2Proc and fastdp successful")
            outData = {
                "edna2Proc": edna2ProcTask.outData,
                "fastDp": fastDpTask.outData,
            }
            self.anomalous = edna2ProcTask.outData.get(
                "anomalous", False
            ) and fastDpTask.outData.get("anomalous", False)
            logger.debug(f"self.anomalous = {self.anomalous}")
        elif edna2ProcTask.isSuccess():
            outData = {
                "edna2Proc": edna2ProcTask.outData,
            }
            self.anomalous = edna2ProcTask.outData.get("anomalous", False)
        elif fastDpTask.isSuccess():
            outData = {
                "fastDp": fastDpTask.outData,
            }
            self.anomalous = fastDpTask.outData.get("anomalous", False)

        if self.anomalous:
            logger.debug("Anomalous flag switched on.")

        comments = ""
        # add comments to ISPyB entry
        if self.anomalous:
            comments += "Strong anomalous signal detected. "
        if edna2ProcTask.isSuccess():
            if edna2ProcTask.outData.get("reindex", False):
                comments += "Data could not be processed by EDNA2Proc with supplied unit cell and space group! "
            if edna2ProcTask.outData.get("twinning", False):
                comments += "Twinning detected by phenix.xtriage! "
            if edna2ProcTask.outData.get("pseudotranslation", False):
                comments += "Pseudotranslation detected by phenix.xtriage! "
            if comments:
                UtilsIspyb.updateDataCollectionGroupComments(
                    self.dataCollectionId, comments
                )

        autoPROCTaskinData = {
            "onlineAutoProcessing": False,
            "workingDirectory": str(self.getWorkingDirectory()),
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath": str(self.masterFilePath) if self.masterFilePath else None,
            "anomalous": self.anomalous,
            "test": self.test,
            "doUploadIspyb": self.doUploadIspyb,
            "waitForFiles": False,
        }

        xia2DialsTaskinData = {
            "onlineAutoProcessing": False,
            "workingDirectory": str(self.getWorkingDirectory()),
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath": str(self.masterFilePath) if self.masterFilePath else None,
            "anomalous": self.anomalous,
            "test": self.test,
            "doUploadIspyb": self.doUploadIspyb,
            "waitForFiles": False,
        }

        doFastSADPhasing = False
        if self.anomalous and fastDpTask.isSuccess():
            mtzFile = fastDpTask.outData.get("mtzFileForFastPhasing", None)
            if mtzFile:
                checkData = FastSADPhasingTask.checkForPhasingDataQuality(
                    mtzFile=mtzFile
                )
                if not checkData:
                    logger.error("Data quality insufficient for phasing. Exiting")
                else:
                    logger.info("Data quality check passed.")
                    doFastSADPhasing = True
            else:
                logger.error("Could not get MTZ file for fast phasing.")

        if doFastSADPhasing:
            logger.info("Starting Fast SAD Phasing...")
            fastSADPhasingTask = FastSADPhasingTask(
                inData={
                    "test": self.test,
                    "dataCollectionId": self.dataCollectionId,
                    "fast_dpMtzFile": mtzFile,
                    "onlineAutoProcessing": False,
                    "checkDataFirst": False,
                },
                workingDirectorySuffix="0",
            )
            fastSADPhasingTask.start()

        autoProcSlurminDataJson = self.getWorkingDirectory() / "inDataAutoPROC.json"
        try:
            with open(autoProcSlurminDataJson, "w+") as fp:
                json.dump(autoPROCTaskinData, fp, indent=4)
        except Exception as e:
            logger.error(f"generating autoPROC json failed: {e}")
        if autoProcSlurminDataJson.is_file():
            autoProcSlurm = f"""\
            #!/bin/bash
            #SBATCH --exclusive
            #SBATCH -t 02:00:00
            #SBATCH --mem=0
            #SBATCH --partition=fujitsu
            #SBATCH -J "EDNA2_aP"
            #SBATCH --output EDNA2job_%j.out
            #SBATCH --chdir {self.getWorkingDirectory()}
            source /mxn/groups/sw/mxsw/env_setup/edna2_proc.sh
            run_edna2.py --inDataFile {autoProcSlurminDataJson} AutoPROCTask
            """
            autoProcSlurm = dedent(autoProcSlurm)

            aPinputstring = autoProcSlurm.encode("ascii")
            out = subprocess.run(
                "sbatch",
                input=aPinputstring,
                cwd=self.getWorkingDirectory(),
                capture_output=True,
            )
            aPJobId = out.stdout.decode("ascii").strip("\n").split()[-1]
            logger.info(f"AutoPROCJob submitted to Slurm with jobId {aPJobId}")

        xia2DialsSlurminDataJson = self.getWorkingDirectory() / "inDataXia2.json"
        try:
            with open(xia2DialsSlurminDataJson, "w+") as fp:
                json.dump(xia2DialsTaskinData, fp, indent=4)
        except Exception as e:
            logger.error(f"generating Xia2 json failed: {e}")
        if xia2DialsSlurminDataJson.is_file():
            xia2DIALSSlurm = f"""\
            #!/bin/bash
            #SBATCH --exclusive
            #SBATCH -t 02:00:00
            #SBATCH --mem=0
            #SBATCH --partition=fujitsu
            #SBATCH -J "EDNA2_x2d"
            #SBATCH --output EDNA2job_%j.out
            #SBATCH --chdir {self.getWorkingDirectory()}
            source /mxn/groups/sw/mxsw/env_setup/edna2_proc.sh
            run_edna2.py --inDataFile {xia2DialsSlurminDataJson} Xia2DIALSTask
            """
            xia2DIALSSlurm = dedent(xia2DIALSSlurm)

            x2DinpuString = xia2DIALSSlurm.encode("ascii")
            out = subprocess.run(
                "sbatch",
                input=x2DinpuString,
                cwd=self.getWorkingDirectory(),
                capture_output=True,
            )
            x2DJobId = out.stdout.decode("ascii").strip("\n").split()[-1]
            logger.info(f"Xia2DIALS submitted to Slurm with jobId {x2DJobId}")

        if doFastSADPhasing:
            fastSADPhasingTask.join()
            if fastSADPhasingTask.isSuccess():
                logger.info("fast SAD phasing completed.")
            else:
                logger.error("Fast SAD phasing failed.")

        return outData
