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

import os
import math
import shutil
import tempfile
import traceback
import numpy as np
from pathlib import Path
import re
import h5py
import socket
import time
from datetime import datetime
import json

STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

# for the os.chmod
from stat import *

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsCCTBX

import logging

logger = UtilsLogging.getLogger()

from edna2.tasks.XDSTasks import XDSIndexing, XDSIntegration, XDSRerunCorrect, XdsstatTask
from edna2.tasks.SubWedgeAssembly import SubWedgeAssembly
from edna2.tasks.CCP4Tasks import (
    PointlessTask,
    AimlessTask,
    CTruncateTask,
    UniqueifyTask,
)
from edna2.tasks.PhenixTasks import MmtbxXtriageTask
from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, ISPyBStoreAutoProcStatus
from edna2.tasks.WaitFileTask import WaitFileTask
from edna2.tasks.Edna2ReportTask import Edna2ReportTask


class Edna2ProcTask(AbstractTask):
    """
    Runs XDS both with and without anomalous
    treament of reflections, given the path
    to a master file.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "dataCollectionId": {"type": ["integer", "null"]},
                "onlineAutoProcessing": {"type": ["boolean", "null"]},
                "waitForFiles": {"type": ["boolean", "null"]},
                "doUploadIspyb": {"type": ["boolean", "null"]},
                "masterFilePath": {"type": ["string", "null"]},
                "spaceGroup": {"type": ["integer", "string"]},
                "unitCell": {"type": ["string", "null"]},
                "residues": {"type": ["integer", "null"]},
                "anomalous": {"type": ["boolean", "null"]},
                "imageNoStart": {"type": ["integer", "null"]},
                "imageNoEnd": {"type": ["integer", "null"]},
                "workingDirectory": {"type": ["string", "null"]},
            },
        }

    # sets ISPyB to FAILED if it's already logged
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
                    isAnom=self.anomalous,
                    timeStart=self.startDateTime,
                    timeEnd=datetime.now().isoformat(timespec="seconds"),
                )

    def run(self, inData):
        UtilsLogging.addLocalFileHandler(logger, self.getWorkingDirectory() / "EDNA2Proc.log")
        logger.info("EDNA2Proc started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        self.tmpdir = None
        self.timeStart = time.perf_counter()
        self.startDateTime = datetime.now().isoformat(timespec="seconds")
        self.processingPrograms = "EDNA2_proc"
        self.processingCommandLine = ""
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.anomalous = inData.get("anomalous", False)
        self.spaceGroup = inData.get("spaceGroup", 0)
        self.unitCell = inData.get("unitCell", None)
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing", False)
        self.imageNoStart = inData.get("imageNoStart", None)
        self.imageNoEnd = inData.get("imageNoEnd", None)
        self.masterFilePath = inData.get("masterFilePath", None)
        self.waitForFiles = inData.get("waitForFiles", True)
        self.doUploadIspyb = inData.get("doUploadIspyb", False)
        self.reindex = False
        self.reintegrate = False
        outData = {}
        self.resultFilePaths = []
        self.pseudoTranslation = False
        self.twinning = False

        try:
            logger.debug(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        if self.dataCollectionId:
            (
                proteinAcronym,
                sampleName,
            ) = UtilsIspyb.getProteinAcronymAndSampleNameFromDataCollectionId(
                self.dataCollectionId
            )
            if proteinAcronym is not None and sampleName is not None:
                # only alphanumerics and underscores are allowed
                try:
                    proteinAcronym = re.sub(r"\W", "_", proteinAcronym)
                    proteinAcronym = re.sub(r"^([_0-9])",r"x\1", proteinAcronym)
                    sampleName = re.sub(r"\W", "_", sampleName)
                    sampleName = re.sub(r"^([_0-9])",r"x\1", sampleName)
                    self.proteinAcronym = proteinAcronym
                    self.sampleName = sampleName
                except Exception as e:
                    logger.error(f"Error parsing proteinAcronym and sampleName: {e}")
                    self.proteinAcronym = "AUTOMATIC"
                    self.sampleName = "DEFAULT"
            else:
                self.proteinAcronym = "AUTOMATIC"
                self.sampleName = "DEFAULT"
        else:
            self.proteinAcronym = "AUTOMATIC"
            self.sampleName = "DEFAULT"

        # set up SG and unit cell
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
            # if self.imageNoEnd - self.imageNoStart < -1:
            logger.error("There are fewer than 8 images, aborting")
            self.setFailure()
            return

        dataH5ImageList = UtilsImage.generateDataFileListFromH5Master(self.masterFilePath)
        pathToStartImage = dataH5ImageList[0]
        pathToEndImage = dataH5ImageList[-1]

        listPrefix = (
            dataCollectionWS3VO.fileTemplate.split("_")
            if dataCollectionWS3VO
            else Path(self.masterFilePath).name.split("_")
        )

        # generate pyarch prefix
        if UtilsConfig.isALBA():
            self.pyarchPrefix = "ap_{0}_{1}".format("_".join(listPrefix[:-2]), listPrefix[-2])
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
            waitFileFirst = WaitFileTask(inData={"file": pathToStartImage, "expectedSize": 100000})
            waitFileFirst.execute()
            if waitFileFirst.outData["timedOut"]:
                logger.warning(
                    "Timeout after {0:d} seconds waiting for the first image {1}!".format(
                        waitFileFirst.outData["timeOut"], pathToStartImage
                    )
                )

            logger.info("Waiting for end image: {0}".format(pathToEndImage))
            waitFileLast = WaitFileTask(inData={"file": pathToEndImage, "expectedSize": 100000})
            waitFileLast.execute()
            if waitFileLast.outData["timedOut"]:
                logger.warning(
                    "Timeout after {0:d} seconds waiting for the last image {1}!".format(
                        waitFileLast.outData["timeOut"], pathToEndImage
                    )
                )

        # get integrationIDs and programIDs, set them to running
        self.integrationId = None
        self.programId = None
        if self.doUploadIspyb:
            try:
                self.integrationId, self.programId = self.createIntegrationId(
                    "Creating integration ID", isAnom=self.anomalous
                )
                logger.info(f"integrationID: {self.integrationId}, programId: {self.programId}")
            except Exception as e:
                logger.error("Could not get integration ID: \n{0}".format(traceback.format_exc(e)))

        # set working directory, results directory, log file
        workingDirectory = self.getWorkingDirectory()
        self.resultsDirectory = Path(workingDirectory / "results")
        self.resultsDirectory.mkdir(parents=True, exist_ok=True)
        self.setLogFileName(f"edna2_proc.log")

        # XXX todo: I hate this
        if not inData.get("subWedge"):
            logger.info("Generating subwedge...")
            # num_of_images, image_list = Edna2ProcTask.generateImageListFromH5Master(inData=inData)
            (
                self.imgNumLow,
                self.imgNumHigh,
                self.imageList,
            ) = self.generateImageListFromH5Master_fast(masterFilePath=self.masterFilePath)
            self.subWedgeAssembly = SubWedgeAssembly(inData=self.imageList, workingDirectorySuffix="0")
            self.subWedgeAssembly.execute()
            self.xdsIndexingInData = self.subWedgeAssembly.outData
        else:
            self.xdsIndexingInData = inData

        self.xdsIndexingInData["unitCell"] = self.unitCell
        self.xdsIndexingInData["spaceGroupNumber"] = self.spaceGroupNumber
        self.xdsIndexingInData["isAnom"] = self.anomalous
        self.xdsIndexingInData["onlineAutoProcessing"] = self.onlineAutoProcessing

        self.indexing = XDSIndexing(inData=self.xdsIndexingInData, workingDirectorySuffix="init")

        logger.info("XDS Indexing started")

        if self.doUploadIspyb:
            self.logToIspyb(self.integrationId, "Indexing", "Launched", "XDS started")

        self.indexing.execute()

        time1 = time.perf_counter()
        self.timeXdsIndexing = time1 - self.timeStart

        if self.indexing.isFailure() and self.unitCell is not None:
            logger.info("Indexing Failed. Rerunning indexing with no unit cell and no space group...")
            self.reindex = True
            self.xdsIndexingInDataRound2 = self.xdsIndexingInData
            self.xdsIndexingInDataRound2["unitCell"] = None
            self.xdsIndexingInDataRound2["spaceGroupNumber"] = 0
            self.indexingRound2 = XDSIndexing(inData=self.xdsIndexingInDataRound2, workingDirectorySuffix="round2")
            logger.info("Starting reindexing")
            self.indexingRound2.execute()

            time1 = time.perf_counter()
            self.timeXdsIndexing = time1 - self.timeStart

            if self.indexingRound2.isFailure():
                logger.error("Rerunning indexing failed. Exiting")
                if self.doUploadIspyb:
                    self.logToIspyb(
                        self.integrationId,
                        "Indexing",
                        "Failed",
                        "XDS failed after {0:.1f}s".format(self.timeXdsIndexing),
                    )
                self.setFailure()
                return
            else:
                self.indexing = self.indexingRound2
        elif self.indexing.isFailure():
            logger.error("Indexing Failed. Exiting")
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Indexing",
                    "Failed",
                    "XDS failed after {0:.1f}s".format(self.timeXdsIndexing),
                )

            self.setFailure()
            return
        else:
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Indexing",
                    "Successful",
                    "XDS finished after {0:.1f}s".format(self.timeXdsIndexing),
                )
            logger.info(f"XDS indexing time took {self.timeXdsIndexing:0.1f} seconds")
            logger.info(
                "Indexing successful. a= {cell_a}, b= {cell_b}, c= {cell_c}, al = {cell_alpha}, be = {cell_beta}, ga = {cell_gamma}".format(
                    **self.indexing.outData["idxref"]["unitCell"]
                )
            )

        # Now set up integration
        self.integrationInData = self.indexing.outData
        del self.integrationInData["workingDirectory"]
        self.integrationInData["subWedge"] = self.indexing.inData["subWedge"]
        self.integrationInData["onlineAutoProcessing"] = self.onlineAutoProcessing
        self.integration = XDSIntegration(inData=self.integrationInData, workingDirectorySuffix="init")
        logger.info("Integration started")

        if self.doUploadIspyb:
            self.logToIspyb(self.integrationId, "Integration", "Launched", "XDS started")

        self.integration.execute()

        time2 = time.perf_counter()
        self.timeXdsIntegration = time2 - time1

        resCutoffFlag = False
        if self.integration.isSuccess():
            # calculate resolution cutoff. Reintegrate if no CC above 30%
            completenessEntries = self.integration.outData.get("completenessEntries", None)
            firstResCutoff = self.getResCutoff(completenessEntries)
            if firstResCutoff is None:
                resCutoffFlag = True

        if (self.integration.isFailure() and self.unitCell is not None and not self.reindex) or (
            resCutoffFlag and self.unitCell is not None
        ):
            logger.info(
                "First round of integration failed. Rerunning indexing and integration with no unit cell and no space group..."
            )
            self.reintegrate = True
            self.xdsIndexingInDataReintRound2 = self.xdsIndexingInData
            self.xdsIndexingInDataReintRound2["unitCell"] = None
            self.xdsIndexingInDataReintRound2["spaceGroupNumber"] = 0
            self.indexingReintRound2 = XDSIndexing(
                inData=self.xdsIndexingInDataReintRound2,
                workingDirectorySuffix="reInt_round2",
            )
            logger.info("Starting Reindexing")
            self.indexingReintRound2.execute()

            time1 = time.perf_counter()
            self.timeXdsIndexing = time1 - self.timeStart

            if self.indexingReintRound2.isFailure():
                logger.error("Rerunning indexing failed. Exiting")
                if self.doUploadIspyb:
                    self.logToIspyb(
                        self.integrationId,
                        "Indexing",
                        "Failed",
                        "XDS failed after {0:.1f}s".format(self.timeXdsIndexing),
                    )
                self.setFailure()
                return
            else:
                logger.info("Reindexing Successful.")
                self.indexing = self.indexingReintRound2
            # Now set up reintegration
            self.reintegrationInData = self.indexing.outData
            del self.reintegrationInData["workingDirectory"]
            self.reintegrationInData["subWedge"] = self.indexing.inData["subWedge"]
            self.reintegrationInData["onlineAutoProcessing"] = self.onlineAutoProcessing
            self.reintegration = XDSIntegration(inData=self.reintegrationInData, workingDirectorySuffix="reInt_round2")
            logger.info("Starting Reintegration")

            self.reintegration.execute()

            if self.reintegration.isFailure():
                logger.error("Error at integration step. Stopping.")
                if self.doUploadIspyb:
                    self.logToIspyb(
                        self.integrationId,
                        "Integration",
                        "Failed",
                        "XDS failed after {0:.1f}s".format(self.timeXdsIntegration),
                    )
                self.setFailure()
                return
            else:
                logger.info("Reintegration Successful.")
                self.integration = self.reintegration
        elif self.integration.isFailure():
            logger.error("Error at integration step. Stopping.")
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Integration",
                    "Failed",
                    "XDS failed after {0:.1f}s".format(self.timeXdsIntegration),
                )
            self.setFailure()
            return
        else:
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Integration",
                    "Successful",
                    "XDS finished after {0:.1f}s".format(self.timeXdsIntegration),
                )
            logger.info(f"XDS integration time took {self.timeXdsIntegration:0.1f} seconds")
            logger.info("Integration Successful.")

        logger.info("Starting first resolution cutoff...")
        self.completenessEntries = self.integration.outData["completenessEntries"]
        self.firstResCutoff = self.getResCutoff(self.completenessEntries)
        if self.firstResCutoff is None:
            logger.error("No bins with CC1/2 greater than 30%")
            logger.error("Something could be wrong, or the completeness could be too low!")
            logger.error("bravais lattice/SG could be incorrect or something more insidious like")
            logger.error("incorrect parameters in XDS.INP like distance, X beam, Y beam, etc.")
            logger.error("Stopping")
            self.setFailure()
            return
        logger.info(f"Resolution cutoff is {self.firstResCutoff}")

        # copy the XDS.INP file from the successful run into the results directory.
        xds_INP_result_path = self.resultsDirectory / f"{self.pyarchPrefix}_successful_XDS.INP"
        integrateLp_path = self.resultsDirectory / f"{self.pyarchPrefix}_INTEGRATE.LP"
        correctLp_path = self.resultsDirectory / f"{self.pyarchPrefix}_CORRECT.LP"
        integrateHkl_path = self.resultsDirectory / f"{self.pyarchPrefix}_INTEGRATE.HKL"
        xdsAsciiHkl_path = self.resultsDirectory / f"{self.pyarchPrefix}_XDS_ASCII.HKL"
        UtilsPath.systemCopyFile(Path(self.integration.outData["integrateHkl"]), integrateHkl_path)
        UtilsPath.systemCopyFile(Path(self.integration.outData["xdsInp"]), xds_INP_result_path)
        UtilsPath.systemCopyFile(Path(self.integration.outData["integrateLp"]), integrateLp_path)
        UtilsPath.systemCopyFile(Path(self.integration.outData["correctLp"]), correctLp_path)
        UtilsPath.systemCopyFile(Path(self.integration.outData["xdsAsciiHkl"]), xdsAsciiHkl_path)

        self.resultFilePaths.extend(
            [
                xds_INP_result_path,
                integrateLp_path,
                correctLp_path,
                integrateHkl_path,
                xdsAsciiHkl_path,
            ]
        )

        # run pointless
        pointlessTaskinData = {
            "input_file": self.integration.outData["xdsAsciiHkl"],
            "output_file": "ep_pointless_unmerged.mtz",
        }
        logger.info("Starting pointless task...")
        self.pointlessTask = PointlessTask(inData=pointlessTaskinData, workingDirectorySuffix="init")
        self.pointlessTask.execute()

        # logger.debug(f"Pointless output: {self.pointlessTask.outData}")

        # now rerun CORRECT with corrected parameters
        rerunCor_data = {
            "xdsInp": self.integration.outData["xdsInp"],
            "spotXds": self.integration.inData["spotXds"],
            "gxParmXds": self.integration.outData["gxParmXds"],
            "xParmXds": self.integration.outData["xParmXds"],
            "gainCbf": self.integration.inData["gainCbf"],
            "blankCbf": self.integration.inData["blankCbf"],
            "xCorrectionsCbf": self.integration.inData["xCorrectionsCbf"],
            "yCorrectionsCbf": self.integration.inData["yCorrectionsCbf"],
            "bkginitCbf": self.integration.inData["bkginitCbf"],
            "integrateLp": self.integration.outData["integrateLp"],
            "integrateHkl": self.integration.outData["integrateHkl"],
            "sg_nr_from_pointless": self.pointlessTask.outData["sgnumber"],
            "cell_from_pointless": self.pointlessTask.outData["cell"],
            "subWedge": self.integration.inData["subWedge"],
            "resCutoff": self.firstResCutoff,
            "onlineAutoProcessing": False,
            "isAnom": self.anomalous,
        }

        logger.info("Rerunning CORRECT with the unit cell/SG from POINTLESS...")
        self.xdsRerun = XDSRerunCorrect(inData=rerunCor_data, workingDirectorySuffix="first")

        self.xdsRerun.execute()

        time3 = time.perf_counter()
        self.timeRerunCorrect = time3 - time2

        if self.xdsRerun.isFailure():
            logger.error("Rerun of CORRECT failed")
            self.setFailure()
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Scaling",
                    "Failed",
                    "Scaling failed after {0:.1}s".format(self.timeRerunCorrect),
                )
            return
        else:
            logger.info(f"Rerun of CORRECT finished.")
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Scaling",
                    "Successful",
                    "Scaling finished in {0:.1f}s".format(self.timeRerunCorrect),
                )

        logger.info("Starting second resolution cutoff...")
        self.completenessEntries = self.xdsRerun.outData["completenessEntries"]

        self.resCutoff = self.getResCutoff(self.completenessEntries)
        if self.resCutoff is None:
            logger.error("Error in determining resolution after CORRECT rerun.")
            logger.error("No bins with CC1/2 greater than 30%")
            logger.error("Something could be wrong, or the completeness could be too low!")
            logger.error("bravais lattice/SG could be incorrect or something more insidious like")
            logger.error("incorrect parameters in XDS.INP like distance, X beam, Y beam, etc.")
            logger.error("Stopping")
            self.setFailure()
            if self.doUploadIspyb:
                self.logToIspyb(self.integrationId, "Scaling", "Failed", "resolution cutoffs failed")
            self.setFailure()
            return

        if self.doUploadIspyb:
            self.logToIspyb(
                self.integrationId,
                "Scaling",
                "Successful",
                "Resolution cutoffs finished",
            )

        self.bins = [
            x["res"] for x in self.xdsRerun.outData["completenessEntries"] if x["include_res_based_on_cc"] is True
        ]

        self.pointlessTaskReruninData = {
            "input_file": self.xdsRerun.outData["xdsAsciiHkl"],
            "output_file": "ep__pointless_unmerged.mtz",
        }
        logger.info("Starting pointless tasks...")
        self.pointlessTaskRerun = PointlessTask(inData=self.pointlessTaskReruninData, workingDirectorySuffix="rerun")

        self.pointlessTaskRerun.execute()
        if self.pointlessTaskRerun.isFailure():
            logger.error("Pointless task failed.")
            self.setFailure()
            return

        self.aimlessTaskinData = {
            "input_file": self.pointlessTaskRerun.outData["pointlessUnmergedMtz"],
            "output_file": "ep__aimless.mtz",
            "start_image": self.indexing.outData["start_image"],
            "end_image": self.indexing.outData["end_image"],
            "dataCollectionId": self.dataCollectionId,
            "res": self.resCutoff,
            "anomalous": True,
        }
        logger.info("Starting aimless...")
        self.aimlessTask = AimlessTask(inData=self.aimlessTaskinData, workingDirectorySuffix="first")
        self.aimlessTask.execute()

        logger.info(f"Aimless finished.")

        highAnomSignalFound = self.if_anomalous_signal(self.aimlessTask.outData["aimlessLog"], threshold=1.0)

        if highAnomSignalFound and not self.anomalous:
            logger.info("Rerunning CORRECT/AIMLESS with anomalous flags due to significant anomalous signal found.")
            self.anomalous = True
            rerunCor_Anomdata = {
                "xdsInp": self.integration.outData["xdsInp"],
                "spotXds": self.integration.inData["spotXds"],
                "gxParmXds": self.integration.outData["gxParmXds"],
                "gainCbf": self.integration.inData["gainCbf"],
                "blankCbf": self.integration.inData["blankCbf"],
                "xCorrectionsCbf": self.integration.inData["xCorrectionsCbf"],
                "yCorrectionsCbf": self.integration.inData["yCorrectionsCbf"],
                "bkginitCbf": self.integration.inData["bkginitCbf"],
                "integrateLp": self.integration.outData["integrateLp"],
                "integrateHkl": self.integration.outData["integrateHkl"],
                "sg_nr_from_pointless": self.pointlessTask.outData["sgnumber"],
                "cell_from_pointless": self.pointlessTask.outData["cell"],
                "subWedge": self.integration.inData["subWedge"],
                "resCutoff": self.firstResCutoff,
                "onlineAutoProcessing": False,
                "isAnom": True,
            }

            logger.info("Rerunning CORRECT with the unit cell/SG from POINTLESS and with anomalous flag on...")
            self.xdsRerunAnom = XDSRerunCorrect(inData=rerunCor_Anomdata, workingDirectorySuffix="rerun_anom")

            self.xdsRerunAnom.execute()

            if self.xdsRerunAnom.isFailure():
                logger.error("Rerun of CORRECT failed")
                self.setFailure()
                if self.doUploadIspyb:
                    self.logToIspyb(
                        self.integrationId,
                        "Scaling",
                        "Failed",
                        "Scaling failed after {0:.1}s".format(self.timeRerunCorrect),
                    )
                return
            else:
                logger.info(f"Rerun of CORRECT finished.")
                if self.doUploadIspyb:
                    self.logToIspyb(
                        self.integrationId,
                        "Scaling",
                        "Successful",
                        "Scaling finished in {0:.1f}s".format(self.timeRerunCorrect),
                    )

            logger.info("Starting third resolution cutoff...")
            self.completenessEntries = self.xdsRerun.outData["completenessEntries"]

            self.resCutoff = self.getResCutoff(self.completenessEntries)
            if self.resCutoff is None:
                logger.error("Error in determining resolution after CORRECT rerun.")
                logger.error("No bins with CC1/2 greater than 30%")
                logger.error("Something could be wrong, or the completeness could be too low!")
                logger.error("bravais lattice/SG could be incorrect or something more insidious like")
                logger.error("incorrect parameters in XDS.INP like distance, X beam, Y beam, etc.")
                logger.error("Stopping")
                self.setFailure()
                if self.doUploadIspyb:
                    self.logToIspyb(
                        self.integrationId,
                        "Scaling",
                        "Failed",
                        "resolution cutoffs failed",
                    )
                self.setFailure()
                return
            logger.info(f"Resolution cutoff is {self.resCutoff}")
            if self.doUploadIspyb:
                self.logToIspyb(
                    self.integrationId,
                    "Scaling",
                    "Successful",
                    "Resolution cutoffs finished",
                )

            self.bins = [
                x["res"]
                for x in self.xdsRerunAnom.outData["completenessEntries"]
                if x["include_res_based_on_cc"] is True
            ]

            self.pointlessTaskRerunAnominData = {
                "input_file": self.xdsRerun.outData["xdsAsciiHkl"],
                "output_file": "ep__pointless_unmerged.mtz",
            }
            logger.info("Starting pointless tasks...")
            self.pointlessTaskRerunAnom = PointlessTask(
                inData=self.pointlessTaskRerunAnominData,
                workingDirectorySuffix="rerunAnom",
            )

            self.pointlessTaskRerunAnom.execute()

            self.aimlessTaskInDataAnom = {
                "input_file": self.pointlessTaskRerunAnom.outData["pointlessUnmergedMtz"],
                "output_file": "ep__aimless.mtz",
                "start_image": self.indexing.outData["start_image"],
                "end_image": self.indexing.outData["end_image"],
                "dataCollectionId": self.dataCollectionId,
                "res": self.resCutoff,
                "anomalous": True,
            }
            logger.info("Starting aimless with anomalous flag on...")
            self.aimlessTaskAnom = AimlessTask(inData=self.aimlessTaskInDataAnom, workingDirectorySuffix="anom")
            self.aimlessTaskAnom.execute()

            logger.info(f"Aimless with anomalous flag finished.")

            self.xdsRerun = self.xdsRerunAnom
            self.pointlessTaskRerun = self.pointlessTaskRerunAnom
            self.aimlessTask = self.aimlessTaskAnom

        xdsStatinData = {"xdsAsciiHkl": self.xdsRerun.outData["xdsAsciiHkl"]}

        self.xdsStat = XdsstatTask(inData=xdsStatinData)
        self.xdsStat.start()
        
        logger.info("Start xdsstat...")

        pointlessUnmergedMtzPath = self.resultsDirectory / (f"{self.pyarchPrefix}_ep__pointless_unmerged.mtz")
        self.resultFilePaths.append(pointlessUnmergedMtzPath)
        self.pointlessLog = self.pointlessTaskRerun.outData["pointlessLog"]
        UtilsPath.systemCopyFile(
            Path(self.pointlessTaskRerun.outData["pointlessUnmergedMtz"]),
            pointlessUnmergedMtzPath,
        )

        aimlessMergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless.mtz"
        aimlessUnmergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless_unmerged.mtz.gz"
        aimlessLogPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless.log"

        self.resultFilePaths.extend([aimlessMergedMtzPath, aimlessUnmergedMtzPath, aimlessLogPath])

        UtilsPath.systemCopyFile(Path(self.aimlessTask.outData["aimlessMergedMtz"]), aimlessMergedMtzPath)
        UtilsPath.systemCopyFile(
            Path(self.aimlessTask.outData["aimlessUnmergedMtz"]),
            aimlessUnmergedMtzPath,
        )
        UtilsPath.systemCopyFile(Path(self.aimlessTask.outData["aimlessLog"]), aimlessLogPath)

        logger.info("Start mmtbx.xtriage run...")
        self.phenixXTriageTaskData = {
            "input_file": aimlessUnmergedMtzPath,
            "workingDirectory": self.getWorkingDirectory(),
        }
        self.phenixXTriageTask = MmtbxXtriageTask(inData=self.phenixXTriageTaskData, workingDirectorySuffix="final")
        self.phenixXTriageTask.start()

        # now run truncate/unique
        tempFile = tempfile.NamedTemporaryFile(
            suffix=".mtz",
            prefix="tmp2-",
            dir=self.aimlessTask.getWorkingDirectory(),
            delete=False,
        )
        truncateOut = tempFile.name
        tempFile.close()
        os.chmod(truncateOut, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
        shutil.chown(truncateOut, group=self.getWorkingDirectory().group())

        logger.info("Start ccp4/ctruncate...")
        self.truncate = CTruncateTask(
            inData={
                "inputFile": aimlessMergedMtzPath,
                "outputFile": truncateOut,
                "isAnom": True,
                "colin": "IMEAN,SIGIMEAN",
                "colano": "I(+),SIGI(+),I(-),SIGI(-)",
            },
            workingDirectorySuffix="final",
        )
        self.truncate.execute()
        if self.truncate.isFailure():
            logger.error("Error running truncate. As this is a rare occurrence, it usually means")
            logger.error("something is seriously wrong with the data. Stopping here.")

        uniqueifyData = {
            "inputFile": self.truncate.outData["cTruncateOutputMtz"],
            "outputFile": "truncate_unique.mtz",
        }
        self.uniqueify = UniqueifyTask(inData=uniqueifyData, workingDirectorySuffix="final")

        self.uniqueify.start()

        self.xdsStat.join()
        self.phenixXTriageTask.join()
        self.uniqueify.join()
        if self.uniqueify.isFailure():
            logger.error("Error running uniqueify. As this is a rare occurrence, it usually means")
            logger.error("something is seriously wrong with the data. Stopping here.")
        if self.phenixXTriageTask.isSuccess():
            logger.info("mmtbx.xtriage finished.")

        truncateLog = self.resultsDirectory / f"{self.pyarchPrefix}_ctruncate.log"
        uniqueMtz = self.resultsDirectory / f"{self.pyarchPrefix}_final.mtz"
        phenixXTriageTaskLog = self.resultsDirectory / f"{self.pyarchPrefix}_xtriage.log"
        
        xdsstatLpFile = self.xdsStat.outData["logFile"]

        UtilsPath.systemCopyFile(Path(self.truncate.outData["cTruncateLogPath"]), truncateLog)
        UtilsPath.systemCopyFile(Path(self.uniqueify.outData["uniqueifyOutputMtz"]), uniqueMtz)
        UtilsPath.systemCopyFile(Path(self.phenixXTriageTask.outData["logFile"]), phenixXTriageTaskLog)

        self.endDateTime = datetime.now().isoformat(timespec="seconds")

        headers = ["File Type", "File name"]
        self.edna2ReportFileList = {
                "mtz_files": [    
                    headers,
                    [
                    "Aimless Merged Data",
                    aimlessMergedMtzPath.name,
                    ],
                    [
                    "Aimless Unmerged Data",
                    aimlessUnmergedMtzPath.name,
                    ],
                    [
                    "Truncate Merged Data (with free column)",
                    uniqueMtz.name,
                    ],
                ],
                "unmerged_files": [
                    headers,
                    [
                    "XDS_ASCII.HKL Output",
                    xdsAsciiHkl_path.name,
                    ],
                    [
                    "INTEGRATE.HKL from XDS",
                    integrateHkl_path.name
                    ]
                ],
                "other_files": [
                    headers,
                    [
                    "XDS.INP file",
                    xds_INP_result_path.name
                    ]

                ],
                "log_files": [
                    [
                    "Aimless Log",
                    aimlessLogPath.name,
                    ],
                    [
                    "Truncate Log",
                    truncateLog.name,
                    ],
                    [
                    "INTEGRATE.LP from XDS",
                    integrateLp_path.name,
                    ],
                    [
                    "CORRECT.LP from XDS",
                    correctLp_path.name
                    ]
                ]
            }
        
        #generate report
        self.Edna2ProcReportInData = {
            "datasetName": self.sampleName,
            "aimlessXml": self.aimlessTask.outData["aimlessXml"],
            "aimlessLog": self.aimlessTask.outData["aimlessLog"],
            "pointlessLog": self.pointlessLog,
            "cTruncateLog": self.truncate.outData["cTruncateLogPath"],
            "integrateLp": self.integration.outData["integrateLp"],
            "xtriageOutput": self.phenixXTriageTask.outData,
            "xdsstatLp": self.xdsStat.outData["logFile"],
            "XdsIndexingLp": self.indexing.outData["idxrefLp"],
            "XdsCorrectLp": self.xdsRerun.outData["correctLp"],
            "pointlessOutput": self.pointlessTask.outData,
            "files": self.edna2ReportFileList,
        }
        logger.info("Generating report...")
        self.edna2Report = Edna2ReportTask(inData=self.Edna2ProcReportInData)
        self.edna2Report.execute()
        edna2Report = self.resultsDirectory / "edna2proc_summary.html"
        UtilsPath.systemCopyFile(Path(self.edna2Report.outData["htmlFile"]), edna2Report)

        if inData.get("test", False):
            self.tmpdir = tempfile.TemporaryDirectory()
            pyarchDirectory = Path(self.tmpdir.name)
            self.pyarchDirectory = self.storeDataOnPyarch(pyarchDirectory=pyarchDirectory)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch()

        # Let's get results into a container for ispyb
        self.autoProcResultsContainer = self.generateAutoProcScalingResultsContainer(
            programId=self.programId,
            integrationId=self.integrationId,
            isAnom=self.anomalous,
        )

        # now send it to ISPyB
        if self.doUploadIspyb:
            logger.info("Sending data to ISPyB...")
            self.ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(
                inData=self.autoProcResultsContainer, workingDirectorySuffix="final"
            )
            self.ispybStoreAutoProcResults.execute()
            if self.ispybStoreAutoProcResults.isFailure():
                logger.error("ISPyB Store autoproc results failed.")
                # self.setFailure()
                # return
            if self.phenixXTriageTask.isSuccess():
                if self.phenixXTriageTask.outData.get("hasTwinning"):
                    self.twinning = True
                if self.phenixXTriageTask.outData.get("hasPseudotranslation"):
                    self.pseudoTranslation = True

        outData = self.autoProcResultsContainer
        outData["anomalous"] = self.anomalous
        outData["reindex"] = self.reindex or self.reintegrate
        outData["twinning"] = self.twinning
        outData["pseudotranslation"] = self.pseudoTranslation


        self.timeEnd = time.perf_counter()
        logger.info(f"Time to process was {self.timeEnd-self.timeStart:0.4f} seconds")
        if self.tmpdir is not None:
            self.tmpdir.cleanup()
        self.timeEnd = time.perf_counter()
        logger.info(f"EDNA2Proc Completed. Process time: {self.timeEnd-self.timeStart:.1f} seconds")
        outData["processTime"] = self.timeEnd - self.timeStart
        return outData

    def storeDataOnPyarch(self, pyarchDirectory=None):
        # create paths on Pyarch
        if pyarchDirectory is None:
            pyarchDirectory = UtilsPath.createPyarchFilePath(self.resultFilePaths[0]).parent
            if not pyarchDirectory.exists():
                pyarchDirectory.mkdir(parents=True, exist_ok=True, mode=0o755)
                logger.debug(f"pyarchDirectory: {pyarchDirectory}")
            for resultFile in [f for f in self.resultFilePaths if f.exists()]:
                resultFilePyarchPath = UtilsPath.createPyarchFilePath(resultFile)
                try:
                    logger.info(f"Copying {resultFile} to pyarch directory")
                    shutil.copy2(resultFile, resultFilePyarchPath)
                except Exception as e:
                    logger.warning(f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}")
                    logger.warning(e)
        else:
            for resultFile in [f for f in self.resultFilePaths if f.exists()]:
                try:
                    logger.info(f"Copying {resultFile} to pyarch directory")
                    resultFilePyarchPath = pyarchDirectory / Path(resultFile).name
                    shutil.copy2(resultFile, resultFilePyarchPath)
                except Exception as e:
                    logger.warning(f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}")
                    logger.warning(e)

        return pyarchDirectory

    def generateAutoProcScalingResultsContainer(self, programId, integrationId, isAnom):
        autoProcResultsContainer = {"dataCollectionId": self.dataCollectionId}

        autoProcProgramContainer = {
            "autoProcProgramId": programId,
            "processingCommandLine": self.processingCommandLine,
            "processingPrograms": self.processingPrograms,
            "processingStatus": "SUCCESS",
            "processingStartTime": self.startDateTime,
            "processingEndTime": self.endDateTime,
        }
        autoProcResultsContainer["autoProcProgram"] = autoProcProgramContainer

        pointlessTask = self.pointlessTaskRerun
        autoProcContainer = {
            "autoProcProgramId": programId,
            "spaceGroup": pointlessTask.outData["sgstr"],
            "refinedCellA": pointlessTask.outData["cell"]["length_a"],
            "refinedCellB": pointlessTask.outData["cell"]["length_b"],
            "refinedCellC": pointlessTask.outData["cell"]["length_c"],
            "refinedCellAlpha": pointlessTask.outData["cell"]["angle_alpha"],
            "refinedCellBeta": pointlessTask.outData["cell"]["angle_beta"],
            "refinedCellGamma": pointlessTask.outData["cell"]["angle_gamma"],
        }
        autoProcResultsContainer["autoProc"] = autoProcContainer

        autoProcAttachmentContainerList = []
        for file in self.pyarchDirectory.iterdir():
            attachmentContainer = {
                "file": file,
            }
            autoProcAttachmentContainerList.append(attachmentContainer)

        autoProcResultsContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList
        xdsRerun = self.xdsRerun.outData

        autoProcIntegrationContainer = {
            "autoProcIntegrationId": integrationId,
            "autoProcProgramId": programId,
            "startImageNumber": self.imgNumLow,
            "endImageNumber": self.imgNumHigh,
            "refinedDetectorDistance": xdsRerun.get("refinedDiffractionParams").get("crystal_to_detector_distance"),
            "refinedXbeam": xdsRerun.get("refinedDiffractionParams").get("direct_beam_detector_coordinates")[0],
            "refinedYbeam": xdsRerun.get("refinedDiffractionParams").get("direct_beam_detector_coordinates")[1],
            "rotationAxisX": xdsRerun.get("gxparmData").get("rot")[0],
            "rotationAxisY": xdsRerun.get("gxparmData").get("rot")[1],
            "rotationAxisZ": xdsRerun.get("gxparmData").get("rot")[2],
            "beamVectorX": xdsRerun.get("gxparmData").get("beam")[0],
            "beamVectorY": xdsRerun.get("gxparmData").get("beam")[1],
            "beamVectorZ": xdsRerun.get("gxparmData").get("beam")[2],
            "cellA": xdsRerun.get("refinedDiffractionParams").get("cell_a"),
            "cellB": xdsRerun.get("refinedDiffractionParams").get("cell_b"),
            "cellC": xdsRerun.get("refinedDiffractionParams").get("cell_c"),
            "cellAlpha": xdsRerun.get("refinedDiffractionParams").get("cell_alpha"),
            "cellBeta": xdsRerun.get("refinedDiffractionParams").get("cell_beta"),
            "cellGamma": xdsRerun.get("refinedDiffractionParams").get("cell_gamma"),
            "anomalous": isAnom,
            "dataCollectionId": self.dataCollectionId,
        }
        autoProcResultsContainer["autoProcIntegration"] = autoProcIntegrationContainer

        autoProcScalingHasIntContainer = {
            "autoProcIntegrationId": integrationId,
        }
        autoProcResultsContainer["autoProcScalingHasInt"] = autoProcScalingHasIntContainer

        # add the scaling statistics to a list of containers...
        autoProcScalingStatisticsContainerList = []
        aimlessResults = self.aimlessTask.outData.get("aimlessResults")

        for shell, result in aimlessResults.items():
            autoProcScalingStatisticsContainer = {}
            autoProcScalingStatisticsContainer["scalingStatisticsType"] = shell
            for k, v in result.items():
                autoProcScalingStatisticsContainer[k] = v
            if shell == "overall":
                autoProcScalingStatisticsContainer["isa"] = xdsRerun.get("ISa", 0.0)
            # need to make a few adjustments for ISPyB...
            autoProcScalingStatisticsContainer["rmerge"] *= 100
            autoProcScalingStatisticsContainer["ccAno"] *= 100
            autoProcScalingStatisticsContainer["rmeasWithinIplusIminus"] *= 100
            autoProcScalingStatisticsContainer["rmeasAllIplusIminus"] *= 100
            autoProcScalingStatisticsContainer["rpimWithinIplusIminus"] *= 100
            autoProcScalingStatisticsContainer["rpimAllIplusIminus"] *= 100
            autoProcScalingStatisticsContainerList.append(autoProcScalingStatisticsContainer)

        autoProcResultsContainer["autoProcScalingStatistics"] = autoProcScalingStatisticsContainerList

        return autoProcResultsContainer

    @staticmethod
    def generateImageListFromH5Master(inData):
        """Given an h5 master file, generate an image list for SubWedgeAssembly."""
        image_path = Path(inData["imagePath"])
        m = re.search(r"\S+_\d{1,2}(?=_master.h5)", image_path.name)
        image_list_stem = m.group(0)

        image_list = []
        master_file = h5py.File(inData["imagePath"], "r")
        for data_file in master_file["/entry/data"].keys():
            image_nr_high = int(master_file["/entry/data"][data_file].attrs["image_nr_high"])
            image_nr_low = int(master_file["/entry/data"][data_file].attrs["image_nr_low"])
            for i in range(image_nr_low, image_nr_high + 1):
                image_list.append(f"{str(image_path.parent)}/{image_list_stem}_{i:06}.h5")
        master_file.close()
        return len(image_list), {"imagePath": image_list}

    @staticmethod
    def generateImageListFromH5Master_fast(masterFilePath):
        """Given an h5 master file, generate an image list for SubWedgeAssembly."""
        masterFilePath = Path(masterFilePath)
        m = re.search(r"\S+_\d{1,2}(?=_master.h5)", masterFilePath.name)
        image_list_stem = m.group(0)

        image_list = []
        with h5py.File(masterFilePath, "r") as master_file:
            data_file_low = list(master_file["/entry/data"].keys())[0]
            data_file_high = list(master_file["/entry/data"].keys())[-1]
            image_nr_high = int(master_file["/entry/data"][data_file_high].attrs["image_nr_high"])
            image_nr_low = int(master_file["/entry/data"][data_file_low].attrs["image_nr_low"])
            image_list.append(f"{str(masterFilePath.parent)}/{image_list_stem}_{image_nr_low:06}.h5")
        return image_nr_low, image_nr_high, {"imagePath": image_list}

    def getResCutoff(self, completeness_entries):
        """
        get resolution cutoff based on CORRECT.LP
        suggestion.
        """
        if completeness_entries is None:
            return None
        return min(
            [x["res"] for x in completeness_entries if x["include_res_based_on_cc"]],
            default=None,
        )

    # Proxy since the API changed and we can now log to several ids
    def logToIspyb(self, integrationId, step, status, comments=""):
        if integrationId is not None:
            if type(integrationId) is list:
                for item in integrationId:
                    self.logToIspybImpl(item, step, status, comments)
            else:
                self.logToIspybImpl(integrationId, step, status, comments)
                # if status == "Failed":
                #     for strErrorMessage in self.getListOfErrorMessages():
                #         self.logToIspybImpl(integrationId, step, status, strErrorMessage)

    def logToIspybImpl(self, integrationId, step, status, comments=""):
        # hack in the event we could not create an integration ID
        if integrationId is None:
            logger.error("could not log to ispyb: no integration id")
            return

        statusInput = {
            "dataCollectionId": self.dataCollectionId,
            "autoProcIntegration": {
                "autoProcIntegrationId": integrationId,
            },
            "autoProcProgram": {},
            "autoProcStatus": {
                "autoProcIntegrationId": integrationId,
                "step": step,
                "status": status,
                "comments": comments,
                "bltimeStamp": datetime.now().isoformat(timespec="seconds"),
            },
        }

        autoprocStatus = ISPyBStoreAutoProcStatus(inData=statusInput, workingDirectorySuffix="")
        autoprocStatus.execute()

    def createIntegrationId(self, comments, isAnom=False):
        """
        gets integrationID and programID,
        sets processing status to RUNNING.
        """
        statusInput = {
            "dataCollectionId": self.dataCollectionId,
            "autoProcIntegration": {
                "anomalous": isAnom,
            },
            "autoProcProgram": {
                "processingCommandLine": self.processingCommandLine,
                "processingPrograms": self.processingPrograms,
                "processingStatus": "RUNNING",
                "processingStartTime": self.startDateTime,
            },
            "autoProcStatus": {
                "step": "Indexing",
                "status": "Launched",
                "comments": comments,
                "bltimeStamp": datetime.now().isoformat(timespec="seconds"),
            },
        }
        autoprocStatus = ISPyBStoreAutoProcStatus(inData=statusInput, workingDirectorySuffix="createIntegrationId")

        # get our EDNAproc status id
        autoprocStatus.execute()
        return (
            autoprocStatus.outData["autoProcIntegrationId"],
            autoprocStatus.outData["autoProcProgramId"],
        )

    def if_anomalous_signal(self, aimless_log, threshold=1.0):
        """Grab the anomalous CC RCR value and see if it is
        sufficiently large to run fast_ep. Generally, a value
        greater than 1 indicates a significant anomalous signal."""
        cc_rcr = 0.0
        summary_switch = False
        try:
            with open(aimless_log, "r") as fp:
                for line in fp:
                    if "$TABLE:  Correlations CC(1/2) within dataset" in line:
                        while "Overall" not in line:
                            line = next(fp)
                        cc_rcr = float(line.split()[3])
                        if not cc_rcr >= threshold:
                            return False
                    if "<!--SUMMARY_BEGIN-->" in line:
                        summary_switch = True
                    if summary_switch and "the anomalous signal is weak" in line:
                        return False
                    if "<!--SUMMARY_END-->" in line:
                        summary_switch = False
        except:
            return False
        return True
