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
import socket

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
from edna2.utils import UtilsCCTBX
from edna2.utils import UtilsImage


logger = UtilsLogging.getLogger()

from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, ISPyBStoreAutoProcStatus
from edna2.tasks.WaitFileTask import WaitFileTask


class AutoPROCTask(AbstractTask):
    def getInDataSchema(self):
        return {
            "type": "object",
            "required": ["dataCollectionId"],
            "properties": {
                "onlineAutoProcessing": {"type": ["boolean", "null"]},
                "dataCollectionId": {"type": ["integer", "null"]},
                "masterFilePath": {"type": ["string", "null"]},
                "spaceGroup": {"type": ["integer", "string"]},
                "unitCell": {"type": ["string", "null"]},
                "residues": {"type": ["integer", "null"]},
                "anomalous": {"type": ["boolean", "null"]},
                "imageNoStart": {"type": ["integer", "null"]},
                "imageNoEnd": {"type": ["integer", "null"]},
                "workingDirectory": {"type": ["string", "null"]},
                "waitForFiles": {"type": ["boolean", "null"]},
                "doUploadIspyb": {"type": ["boolean", "null"]},
            },
        }

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
                    timeEnd=datetime.now().isoformat(timespec="seconds"),
                )
                self.logToIspyb(
                    self.integrationId, "Indexing", "Failed", "AutoPROC ended"
                )

            if (
                self.integrationIdStaraniso is not None
                and self.programIdStaraniso is not None
            ):
                ISPyBStoreAutoProcResults.setIspybToFailed(
                    dataCollectionId=self.dataCollectionId,
                    autoProcProgramId=self.programIdStaraniso,
                    autoProcIntegrationId=self.integrationIdStaraniso,
                    processingCommandLine=self.processingCommandLine,
                    processingPrograms=self.processingPrograms,
                    isAnom=False,
                    timeStart=self.startDateTime,
                    timeEnd=datetime.now().isoformat(timespec="seconds"),
                )
                self.logToIspyb(
                    self.integrationIdStaraniso, "Indexing", "Failed", "AutoPROC ended"
                )

    def run(self, inData):
        UtilsLogging.addLocalFileHandler(
            logger, self.getWorkingDirectory() / "EDNA_autoPROC.log"
        )
        logger.info("AutoPROC processing started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        self.timeStart = time.perf_counter()
        self.startDateTime = datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.processingPrograms = "autoPROC"
        self.processingProgramsStaraniso = "autoPROC_staraniso"
        self.processingCommandLine = ""
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing", False)
        # self.setLogFileName(f"autoPROC_{self.startDateTimeFormatted}.log")
        self.dataCollectionId = inData.get("dataCollectionId")
        self.tmpdir = None
        self.imageNoStart = inData.get("imageNoStart", None)
        self.imageNoEnd = inData.get("imageNoEnd", None)
        self.doUploadIspyb = inData.get("doUploadIspyb", False)
        self.waitForFiles = inData.get("waitForFiles", True)
        self.masterFilePath = inData.get("masterFilePath", None)
        self.integrationId, self.programId = None, None
        self.integrationIdStaraniso, self.programIdStaraniso = None, None

        try:
            logger.debug(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        pathToStartImage = None
        pathToEndImage = None

        self.anomalous = inData.get("anomalous", False)

        logger.debug("Working directory is {0}".format(self.getWorkingDirectory()))

        self.spaceGroup = inData.get("spaceGroup", 0)
        self.unitCell = inData.get("unitCell", None)
        self.lowResLimit = inData.get("lowResolutionLimit", None)
        self.highResLimit = inData.get("highResolutionLimit", None)

        # set up SG and unit cell
        self.spaceGroupNumber, self.spaceGroupString = UtilsCCTBX.parseSpaceGroup(
            self.spaceGroup
        )

        # set up unit cell
        if self.unitCell is not None:
            self.unitCell = UtilsCCTBX.parseUnitCell(self.unitCell)

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

        if self.imageNoEnd - self.imageNoStart < 8:
            # if self.imageNoEnd - self.imageNoStart < -1:
            logger.error("There are fewer than 8 images, aborting")
            self.setFailure()
            return

        dataH5ImageList = UtilsImage.generateDataFileListFromH5Master(
            self.masterFilePath
        )
        pathToStartImage = dataH5ImageList[0]
        pathToEndImage = dataH5ImageList[-1]

        listPrefix = (
            dataCollectionWS3VO.fileTemplate.split("_")
            if dataCollectionWS3VO
            else Path(self.masterFilePath).name.split("_")
        )

        # generate pyarch prefix
        if UtilsConfig.isALBA():
            self.pyarchPrefix = "ap_{0}_{1}".format(
                "_".join(listPrefix[:-2]), listPrefix[-2]
            )
        elif UtilsConfig.isMAXIV():
            self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])
        else:
            if len(listPrefix) > 2:
                self.pyarchPrefix = "ap_{0}_run{1}".format(
                    listPrefix[-3], listPrefix[-2]
                )
            elif len(listPrefix) > 1:
                self.pyarchPrefix = "ap_{0}_run{1}".format(
                    listPrefix[:-2], listPrefix[-2]
                )
            else:
                self.pyarchPrefix = "ap_{0}_run".format(listPrefix[0])

        # make results directory
        self.resultsDirectory = self.getWorkingDirectory() / "results"
        self.resultsDirectory.mkdir(exist_ok=True, parents=True, mode=0o755)

        if self.waitForFiles:
            logger.info("Waiting for start image: {0}".format(pathToStartImage))
            waitFileFirst = WaitFileTask(
                inData={"file": pathToStartImage, "expectedSize": 100000}
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
                inData={"file": pathToEndImage, "expectedSize": 100000}
            )
            waitFileLast.execute()
            if waitFileLast.outData["timedOut"]:
                logger.warning(
                    "Timeout after {0:d} seconds waiting for the last image {1}!".format(
                        waitFileLast.outData["timeOut"], pathToEndImage
                    )
                )

        self.timeStart = datetime.now().isoformat(timespec="seconds")

        if self.doUploadIspyb:
            # set ISPyB to running
            (
                self.integrationId,
                self.programId,
            ) = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine=self.processingCommandLine,
                processingPrograms=self.processingPrograms,
                isAnom=self.anomalous,
                timeStart=self.timeStart,
            )
            (
                self.integrationIdStaraniso,
                self.programIdStaraniso,
            ) = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine=self.processingCommandLine,
                processingPrograms=self.processingProgramsStaraniso,
                isAnom=self.anomalous,
                timeStart=self.timeStart,
            )

        autoPROCExecinData = {
            "dataCollectionId": self.dataCollectionId,
            "onlineAutoProcessing": self.onlineAutoProcessing,
            "masterFilePath": self.masterFilePath,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd": self.imageNoEnd,
            "spaceGroupNumber": self.spaceGroupNumber,
            "spaceGroupString": self.spaceGroupString,
            "unitCell": self.unitCell,
            "lowResLimit": self.lowResLimit,
            "highResLimit": self.highResLimit,
            "anomalous": self.anomalous,
        }
        if self.doUploadIspyb:
            self.logToIspyb(
                self.integrationId, "Indexing", "Launched", "AutoPROC started"
            )
            self.logToIspyb(
                self.integrationIdStaraniso, "Indexing", "Launched", "AutoPROC started"
            )

        autoPROCExec = AutoPROCExecTask(
            inData=autoPROCExecinData, workingDirectorySuffix="0"
        )
        autoPROCExec.execute()
        if autoPROCExec.isFailure():
            self.setFailure()
            return

        self.autoPROCExecDir = Path(autoPROCExec.outData["workingDirectory"])
        self.endDateTime = datetime.now().isoformat(timespec="seconds")

        ispybXml = self.autoPROCExecDir / "autoPROC.xml"
        if ispybXml.is_file():
            self.outData["ispybXml"] = str(ispybXml)
            autoProcContainer = self.autoPROCXMLtoISPyBdict(
                xml_path=ispybXml,
                data_collection_id=self.dataCollectionId,
                program_id=self.programId,
                integration_id=self.integrationId,
                processing_programs=self.processingPrograms,
                anomalous=self.anomalous,
            )

        ispybXmlStaraniso = self.autoPROCExecDir / "autoPROC_staraniso.xml"
        if ispybXmlStaraniso.is_file():
            self.outData["ispybXml_staraniso"] = str(ispybXmlStaraniso)
            autoProcContainerStaraniso = self.autoPROCXMLtoISPyBdict(
                xml_path=ispybXmlStaraniso,
                data_collection_id=self.dataCollectionId,
                program_id=self.programIdStaraniso,
                integration_id=self.integrationIdStaraniso,
                processing_programs=self.processingProgramsStaraniso,
                anomalous=self.anomalous,
            )

        # get CIF Files and gzip them
        autoPROCStaranisoAllCif = (
            self.autoPROCExecDir / "Data_1_autoPROC_STARANISO_all.cif"
        )
        autoPROCStaranisoAllCifGz = (
            self.resultsDirectory / f"{self.pyarchPrefix}_autoPROC_STARANISO_all.cif.gz"
        )
        autoPROCTruncateAllCif = (
            self.autoPROCExecDir / "Data_2_autoPROC_TRUNCATE_all.cif"
        )
        autoPROCTruncateAllCifGz = (
            self.resultsDirectory / f"{self.pyarchPrefix}_autoPROC_TRUNCATE_all.cif.gz"
        )
        autoPROCXdsAsciiHkl = self.autoPROCExecDir / "XDS_ASCII.HKL"
        autoPROCXdsAsciiHklGz = (
            self.resultsDirectory / f"{self.pyarchPrefix}_XDS_ASCII.HKL.gz"
        )

        try:
            logger.debug(f"gzip'ing {autoPROCStaranisoAllCif}")
            with open(autoPROCStaranisoAllCif, "rb") as fp_in:
                with gzip.open(autoPROCStaranisoAllCifGz, "wb") as fp_out:
                    shutil.copyfileobj(fp_in, fp_out)
        except:
            logger.error(f"gzip'ing {autoPROCStaranisoAllCif} failed.")
        try:
            logger.debug(f"gzip'ing {autoPROCTruncateAllCif}")
            with open(autoPROCTruncateAllCif, "rb") as fp_in:
                with gzip.open(autoPROCTruncateAllCifGz, "wb") as fp_out:
                    shutil.copyfileobj(fp_in, fp_out)
        except:
            logger.error(f"gzip'ing {autoPROCTruncateAllCif} failed.")
        try:
            logger.debug(f"gzip'ing {autoPROCXdsAsciiHkl}")
            with open(autoPROCXdsAsciiHkl, "rb") as fp_in:
                with gzip.open(autoPROCXdsAsciiHklGz, "wb") as fp_out:
                    shutil.copyfileobj(fp_in, fp_out)
        except:
            logger.error(f"gzip'ing {autoPROCXdsAsciiHkl} failed.")

        # copy files to results directory
        autoPROCLogFile = autoPROCExec.outData.get("logFile")
        autoPROCReportPdf = self.autoPROCExecDir / "report.pdf"
        autoPROCStaranisoReportPdf = self.autoPROCExecDir / "report_staraniso.pdf"
        autoPROCStaranisoAllDataUniqueMtz = (
            self.autoPROCExecDir / "staraniso_alldata-unique.mtz"
        )
        autoPROCStaranisoAllDataUniqueStats = (
            self.autoPROCExecDir / "staraniso_alldata-unique.stats"
        )
        autoPROCStaranisoAllDataUniqueTable1 = (
            self.autoPROCExecDir / "staraniso_alldata-unique.table1"
        )
        autoPROCSummaryInlinedHtml = self.autoPROCExecDir / "summary_inlined.html"
        autoPROCSummaryTarGz = self.autoPROCExecDir / "summary.tar.gz"
        autoPROCTruncateUniqueMtz = self.autoPROCExecDir / "truncate-unique.mtz"
        autoPROCTruncateUniqueStats = self.autoPROCExecDir / "truncate-unique.stats"
        autoPROCTruncateUniqueTable1 = self.autoPROCExecDir / "truncate-unique.table1"

        autoPROCLogFile_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_autoPROC.log"
        )
        autoPROCReportPdf_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_report.pdf"
        )
        autoPROCStaranisoReportPdf_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_report_staraniso.pdf"
        )
        autoPROCStaranisoAllDataUniqueMtz_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_staraniso_alldata-unique.mtz"
        )
        autoPROCStaranisoAllDataUniqueStats_resultsDir = (
            self.resultsDirectory
            / f"{self.pyarchPrefix}_staraniso_alldata-unique.stats"
        )
        autoPROCStaranisoAllDataUniqueTable1_resultsDir = (
            self.resultsDirectory
            / f"{self.pyarchPrefix}_staraniso_alldata-unique.table1"
        )
        autoPROCSummaryInlinedHtml_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_summary_inlined.html"
        )
        autoPROCSummaryTarGz_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_summary.tar.gz"
        )
        autoPROCTruncateUniqueMtz_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_truncate-unique.mtz"
        )
        autoPROCTruncateUniqueStats_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_truncate-unique.stats"
        )
        autoPROCTruncateUniqueTable1_resultsDir = (
            self.resultsDirectory / f"{self.pyarchPrefix}_truncate-unique.table1"
        )

        autoProcAttachmentContainerList = []
        autoProcAttachmentContainerStaranisoList = []

        for files in [
            (autoPROCLogFile, autoPROCLogFile_resultsDir),
            (autoPROCReportPdf, autoPROCReportPdf_resultsDir),
            (autoPROCStaranisoReportPdf, autoPROCStaranisoReportPdf_resultsDir),
            (
                autoPROCStaranisoAllDataUniqueMtz,
                autoPROCStaranisoAllDataUniqueMtz_resultsDir,
            ),
            (
                autoPROCStaranisoAllDataUniqueStats,
                autoPROCStaranisoAllDataUniqueStats_resultsDir,
            ),
            (
                autoPROCStaranisoAllDataUniqueTable1,
                autoPROCStaranisoAllDataUniqueTable1_resultsDir,
            ),
            (autoPROCSummaryInlinedHtml, autoPROCSummaryInlinedHtml_resultsDir),
            (autoPROCSummaryTarGz, autoPROCSummaryTarGz_resultsDir),
            (autoPROCTruncateUniqueMtz, autoPROCTruncateUniqueMtz_resultsDir),
            (autoPROCTruncateUniqueStats, autoPROCTruncateUniqueStats_resultsDir),
            (autoPROCTruncateUniqueTable1, autoPROCTruncateUniqueTable1_resultsDir),
        ]:
            UtilsPath.systemCopyFile(files[0], files[1])

            pyarchFile = (
                UtilsPath.createPyarchFilePath(files[1])
                if self.doUploadIspyb
                else files[1]
            )

            attachmentContainer = {
                "file": pyarchFile,
            }
            if "staraniso" in str(files[1]):
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
            elif "truncate-unique" in str(files[1]):
                autoProcAttachmentContainerList.append(attachmentContainer)
            else:
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
                autoProcAttachmentContainerList.append(attachmentContainer)

        for file in [
            autoPROCStaranisoAllCifGz,
            autoPROCTruncateAllCifGz,
            autoPROCXdsAsciiHklGz,
        ]:
            pyarchFile = (
                UtilsPath.createPyarchFilePath(file) if self.doUploadIspyb else file
            )
            attachmentContainer = {
                "file": pyarchFile,
            }
            if "STARANISO" in str(file):
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
            elif "TRUNCATE" in str(file):
                autoProcAttachmentContainerList.append(attachmentContainer)
            else:
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
                autoProcAttachmentContainerList.append(attachmentContainer)

        autoProcContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList
        autoProcContainerStaraniso[
            "autoProcProgramAttachment"
        ] = autoProcAttachmentContainerStaranisoList

        # save as json
        autoProcContainerJson = self.resultsDirectory / "autoPROC.json"
        autoProcContainerStaranisoJson = (
            self.resultsDirectory / "autoPROC_staraniso.json"
        )

        with open(autoProcContainerJson, "w") as fp:
            json.dump(autoProcContainer, fp, indent=2, default=lambda o: str(o))

        with open(autoProcContainerStaranisoJson, "w") as fp:
            json.dump(
                autoProcContainerStaraniso, fp, indent=2, default=lambda o: str(o)
            )

        self.resultFilePaths = list(self.resultsDirectory.iterdir())
        if inData.get("test", False):
            self.tmpdir = tempfile.TemporaryDirectory()
            self.pyarchDirectory = Path(self.tmpdir.name)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch()

        if self.doUploadIspyb:
            self.logToIspyb(
                self.integrationId, "Indexing", "Successful", "AutoPROC finished"
            )
            self.logToIspyb(
                self.integrationIdStaraniso,
                "Indexing",
                "Successful",
                "AutoPROC finished",
            )

            ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(
                inData=autoProcContainer, workingDirectorySuffix="uploadFinal"
            )
            ispybStoreAutoProcResults.execute()
            ispybStoreAutoProcResultsStaraniso = ISPyBStoreAutoProcResults(
                inData=autoProcContainerStaraniso,
                workingDirectorySuffix="uploadFinal_staraniso",
            )
            ispybStoreAutoProcResultsStaraniso.execute()

        outData = {
            "autoPROC": autoProcContainer,
            "autoPROC_staraniso": autoProcContainerStaraniso,
        }

        if self.tmpdir is not None:
            self.tmpdir.cleanup()
        return outData

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
                logger.info(
                    f"Copying {resultFile} to pyarch directory {resultFilePyarchPath}"
                )
                shutil.copy(resultFile, resultFilePyarchPath)
            except Exception as e:
                logger.warning(
                    f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}"
                )
                logger.warning(e)
        return pyarchDirectory

    @staticmethod
    def autoPROCXMLtoISPyBdict(
        xml_path,
        data_collection_id=None,
        program_id=None,
        integration_id=None,
        processing_programs=None,
        anomalous=False,
        trunc_len=256,
    ):
        dict_data = UtilsXML.dictfromXML(xml_path)
        autoProcXMLContainer = dict_data["AutoProcContainer"]
        autoProcContainer = {"dataCollectionId": data_collection_id}
        autoProcContainer["autoProcProgram"] = {}
        for k, v in autoProcXMLContainer["AutoProcProgramContainer"][
            "AutoProcProgram"
        ].items():
            if isinstance(v, str) and len(v) > trunc_len:
                autoProcContainer["autoProcProgram"][k] = v[: trunc_len - 1]
                logger.warning(
                    f'string {k} truncated for loading to ISPyB to {trunc_len} characters: "{v}"'
                )
            else:
                autoProcContainer["autoProcProgram"][k] = v
        # fix some entries in autoProcProgram
        autoProcContainer["autoProcProgram"]["processingStatus"] = "SUCCESS"
        autoProcContainer["autoProcProgram"]["processingPrograms"] = processing_programs
        autoProcContainer["autoProcProgram"].pop("processingEnvironment")
        autoProcContainer["autoProcProgram"].pop("processingMessage")

        timeStart = autoProcContainer["autoProcProgram"]["processingStartTime"]
        timeEnd = autoProcContainer["autoProcProgram"]["processingEndTime"]
        timeStart_structtime = time.strptime(timeStart, "%a %b %d %H:%M:%S %Z %Y")
        timeEnd_structtime = time.strptime(timeEnd, "%a %b %d %H:%M:%S %Z %Y")
        autoProcContainer["autoProcProgram"]["processingStartTime"] = time.strftime(
            "%Y-%m-%dT%H:%M:%S", timeStart_structtime
        )
        autoProcContainer["autoProcProgram"]["processingEndTime"] = time.strftime(
            "%Y-%m-%dT%H:%M:%S", timeEnd_structtime
        )

        autoProcContainer["autoProcProgram"]["autoProcProgramId"] = program_id
        autoProcContainer["autoProc"] = autoProcXMLContainer["AutoProc"]
        # fix some entries in autoProc
        autoProcContainer["autoProc"]["refinedCellA"] = autoProcContainer[
            "autoProc"
        ].pop("refinedCell_a")
        autoProcContainer["autoProc"]["refinedCellB"] = autoProcContainer[
            "autoProc"
        ].pop("refinedCell_b")
        autoProcContainer["autoProc"]["refinedCellC"] = autoProcContainer[
            "autoProc"
        ].pop("refinedCell_c")
        autoProcContainer["autoProc"]["refinedCellAlpha"] = autoProcContainer[
            "autoProc"
        ].pop("refinedCell_alpha")
        autoProcContainer["autoProc"]["refinedCellBeta"] = autoProcContainer[
            "autoProc"
        ].pop("refinedCell_beta")
        autoProcContainer["autoProc"]["refinedCellGamma"] = autoProcContainer[
            "autoProc"
        ].pop("refinedCell_gamma")
        autoProcContainer["autoProc"]["autoProcProgramId"] = program_id
        del autoProcContainer["autoProc"]["wavelength"]

        autoProcContainer["autoProcScaling"] = autoProcXMLContainer[
            "AutoProcScalingContainer"
        ]["AutoProcScaling"]
        autoProcContainer["autoProcScalingStatistics"] = autoProcXMLContainer[
            "AutoProcScalingContainer"
        ]["AutoProcScalingStatistics"]

        autoProcContainer["autoProcIntegration"] = autoProcXMLContainer[
            "AutoProcScalingContainer"
        ]["AutoProcIntegrationContainer"]["AutoProcIntegration"]
        autoProcContainer["autoProcIntegration"]["autoProcProgramId"] = program_id
        autoProcContainer["autoProcIntegration"][
            "autoProcIntegrationId"
        ] = integration_id
        # fix some entries in AutoProcScalingIntegration
        autoProcContainer["autoProcIntegration"]["cellA"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("cell_a")
        autoProcContainer["autoProcIntegration"]["cellB"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("cell_b")
        autoProcContainer["autoProcIntegration"]["cellC"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("cell_c")
        autoProcContainer["autoProcIntegration"]["cellAlpha"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("cell_alpha")
        autoProcContainer["autoProcIntegration"]["cellBeta"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("cell_beta")
        autoProcContainer["autoProcIntegration"]["cellGamma"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("cell_gamma")
        autoProcContainer["autoProcIntegration"]["refinedXbeam"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("refinedXBeam")
        autoProcContainer["autoProcIntegration"]["refinedYbeam"] = autoProcContainer[
            "autoProcIntegration"
        ].pop("refinedYBeam")
        autoProcContainer["autoProcIntegration"]["anomalous"] = anomalous

        autoProcScalingHasIntContainer = {
            "autoProcIntegrationId": integration_id,
        }
        autoProcContainer["autoProcScalingHasInt"] = autoProcScalingHasIntContainer

        for k, v in autoProcContainer["autoProc"].items():
            autoProcContainer["autoProc"][k] = AutoPROCTask.convertStrToIntOrFloat(v)

        for k, v in autoProcContainer["autoProcIntegration"].items():
            autoProcContainer["autoProcIntegration"][
                k
            ] = AutoPROCTask.convertStrToIntOrFloat(v)

        for k, v in autoProcContainer["autoProcScaling"].items():
            autoProcContainer["autoProcScaling"][
                k
            ] = AutoPROCTask.convertStrToIntOrFloat(v)

        for shell in autoProcContainer["autoProcScalingStatistics"]:
            for k, v in shell.items():
                # should they be ints, floats, or strings? I don't know,
                # but seems like they shouldn't be strings...
                shell[k] = AutoPROCTask.convertStrToIntOrFloat(v)

                # rMeas, rPim, and rMerge need to be multiplied by 100
                if any(x for x in ["rMerge", "rMeas", "rPim", "ccAno"] if x in k):
                    shell[k] *= 100
            shell["rmerge"] = shell.pop("rMerge")
            shell["rmeasWithinIplusIminus"] = shell.pop("rMeasWithinIPlusIMinus")
            shell["rmeasAllIplusIminus"] = shell.pop("rMeasAllIPlusIMinus")
            shell["rpimWithinIplusIminus"] = shell.pop("rPimWithinIPlusIMinus")
            shell["rpimAllIplusIminus"] = shell.pop("rPimAllIPlusIMinus")
            shell["meanIoverSigI"] = shell.pop("meanIOverSigI")
            shell["ccAno"] = shell.pop("ccAnomalous")
            shell["sigAno"] = shell.pop("DanoOverSigDano")

        return autoProcContainer

    @staticmethod
    def convertStrToIntOrFloat(v: str):
        """
        Tries to convert a string to an int first, then a float.
        If it doesn't work, returns the string.
        """
        if isinstance(v, str):
            try:
                v = int(v)
            except:
                try:
                    v = float(v)
                except:
                    pass
        return v

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

        autoprocStatus = ISPyBStoreAutoProcStatus(
            inData=statusInput, workingDirectorySuffix=""
        )
        autoprocStatus.execute()


class AutoPROCExecTask(AbstractTask):
    def run(self, inData):
        logger.info("AutoPROC execution started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        outData = {}
        logger.debug(f"working directory is {self.getWorkingDirectory()}")
        self.onlineAutoProcessing = inData["onlineAutoProcessing"]
        self.masterFilePath = inData["masterFilePath"]
        self.imageNoStart = inData["imageNoStart"]
        self.imageNoEnd = inData["imageNoEnd"]

        self.spaceGroupNumber = inData.get("spaceGroupNumber", 0)
        self.spaceGroupString = inData.get("spaceGroupString", "")

        self.unitCell = inData.get("unitCell", None)
        self.lowResLimit = inData.get("lowResolutionLimit", None)
        self.highResLimit = inData.get("highResolutionLimit", None)
        self.anomalous = inData.get("anomalous", False)

        # set up command line
        self.autoPROCExecDir = self.getWorkingDirectory() / "AutoPROCExec_0"
        inc_x = 1
        while self.autoPROCExecDir.is_dir():
            self.autoPROCExecDir = (
                self.getWorkingDirectory() / "AutoPROCExec_{0}".format(inc_x)
            )
            inc_x += 1

        outData["workingDirectory"] = str(self.autoPROCExecDir)

        autoPROCSetup = UtilsConfig.get("AutoPROCTask", "autoPROCSetup", None)
        autoPROCExecutable = UtilsConfig.get(
            "AutoPROCTask", "autoPROCExecutable", "process"
        )
        maxNoProcessors = UtilsConfig.get("AutoPROCTask", "maxNoProcessors", None)
        autoPROCmacros = UtilsConfig.get("AutoPROCTask", "macros", None)

        if autoPROCSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + autoPROCSetup + "\n"
        commandLine += " {0}".format(autoPROCExecutable)
        # add flags, if present
        commandLine += " -B"
        commandLine += " -d {0}".format(str(self.autoPROCExecDir))
        commandLine += " -nthreads {0}".format(maxNoProcessors)

        commandLine += " -h5 {0}".format(self.masterFilePath)
        commandLine += " -ANO" if self.anomalous else ""

        if autoPROCmacros is not None:
            for macro in autoPROCmacros.split():
                commandLine += " -M {0}".format(macro)

        if self.spaceGroupString != "" and self.unitCell is not None:
            commandLine += ' symm="{0}"'.format(self.spaceGroupString)
            commandLine += ' cell="{cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}"'.format(
                **self.unitCell
            )

        if self.lowResLimit is not None or self.highResLimit is not None:
            low = self.lowResLimit if self.lowResLimit else 1000.0
            high = self.highResLimit if self.highResLimit else 0.1
            commandLine += " -R {0} {1}".format(low, high)

        config = UtilsConfig.getConfig()
        config.optionxform = str
        aP_config = config["AutoPROCTask"]
        # logger.debug(f"{aP_config}")
        for k, v in aP_config.items():
            if k.startswith("autoPROC_"):
                logger.info(f"autoPROC option: {k}={v}")
                commandLine += " {0}={1}".format(k, v)

        logger.info("autoPROC command is {}".format(commandLine))

        if self.onlineAutoProcessing:
            returncode = self.submitCommandLine(
                commandLine,
                jobName="EDNA2_aP",
                timeout=UtilsConfig.get("AutoPROCTask","slurm_timeout"),
                partition=UtilsConfig.get("AutoPROCTask","slurm_partition"),
                ignoreErrors=False,
            )
            outData["logFile"] = self.getSlurmLogPath()
            if returncode != 0:
                self.setFailure()
                return
        else:
            try:
                self.runCommandLine(commandLine, listCommand=[])
                outData["logFile"] = self.getLogPath()
            except RuntimeError:
                self.setFailure()
                return

        return outData
