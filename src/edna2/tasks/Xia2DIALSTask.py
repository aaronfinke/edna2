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

from datetime import datetime

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsImage
from edna2.utils import UtilsCCTBX


logger = UtilsLogging.getLogger()

from edna2.tasks.ISPyBTasks import (
    ISPyBStoreAutoProcResults,
    ISPyBStoreAutoProcStatus,
    createIntegrationId,
)
from edna2.tasks.WaitFileTask import WaitFileTask


STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"


class Xia2DIALSTask(AbstractTask):
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
            self.logToIspyb(self.integrationId, "Indexing", "Failed", "Xia2DIALS ended")

    def run(self, inData):
        UtilsLogging.addLocalFileHandler(
            logger, self.getWorkingDirectory() / "EDNA_xia2DIALS.log"
        )
        logger.info("Xia2DIALS processing started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        self.timeStart = time.perf_counter()
        self.startDateTime = datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.processingPrograms = "XIA2_DIALS"
        self.processingCommandLine = ""

        self.setLogFileName(f"xia2DIALS_{self.startDateTimeFormatted}.log")
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.tmpdir = None
        pathToStartImage = None
        pathToEndImage = None

        self.anomalous = inData.get("anomalous", False)
        self.spaceGroup = inData.get("spaceGroup", 0)
        self.unitCell = inData.get("unitCell", None)
        self.lowResLimit = inData.get("lowResolutionLimit", None)
        self.highResLimit = inData.get("highResolutionLimit", None)
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing", False)
        self.masterFilePath = inData.get("masterFilePath", None)
        self.doUploadIspyb = inData.get("doUploadIspyb", False)
        self.waitForFiles = inData.get("waitForFiles", True)
        self.imageNoStart = inData.get("imageNoStart", None)
        self.imageNoEnd = inData.get("imageNoEnd", None)
        self.pyarchDirectory = None
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

        (
            proteinAcronym,
            sampleName,
        ) = UtilsIspyb.getProteinAcronymAndSampleNameFromDataCollectionId(
            self.dataCollectionId
        )
        if proteinAcronym is not None and sampleName is not None:
            # only alphanumerics and underscores are allowed
            proteinAcronym_corrected = re.sub(r"\W", "_", proteinAcronym)
            sampleName_corrected = re.sub(r"\W", "_", sampleName)
            self.proteinAcronym = proteinAcronym_corrected
            self.sampleName = sampleName_corrected
        else:
            self.proteinAcronym = "AUTOMATIC"
            self.sampleName = "DEFAULT"

        logger.info(
            f"Protein Acronym:{self.proteinAcronym}, sample name:{self.sampleName}"
        )

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

        xia2DIALSExecinData = {
            "onlineAutoProcessing": self.onlineAutoProcessing,
            "masterFilePath": self.masterFilePath,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd": self.imageNoEnd,
            "proteinAcronym": self.proteinAcronym,
            "sampleName": self.sampleName,
            "spaceGroupNumber": self.spaceGroupNumber,
            "spaceGroupString": self.spaceGroupString,
            "unitCell": self.unitCell,
            "lowResLimit": self.lowResLimit,
            "highResLimit": self.highResLimit,
            "anomalous": self.anomalous,
        }
        if self.doUploadIspyb:
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

        if self.doUploadIspyb:
            self.logToIspyb(
                self.integrationId, "Indexing", "Launched", "Xia2Dials Launched"
            )

        xia2DIALSExec = Xia2DialsExecTask(
            inData=xia2DIALSExecinData, workingDirectorySuffix="0"
        )
        xia2DIALSExec.execute()
        if xia2DIALSExec.isFailure():
            self.setFailure()
            return
        self.timeEnd = datetime.now().isoformat(timespec="seconds")

        xia2DIALSExecDir = Path(xia2DIALSExec.outData["workingDirectory"])
        logger.debug(f"Working directory is {xia2DIALSExecDir}")

        for resultFile in Path(xia2DIALSExecDir / "DataFiles").glob("*"):
            targetFile = (
                self.resultsDirectory / f"{self.pyarchPrefix}_{resultFile.name}"
            )
            UtilsPath.systemCopyFile(resultFile, targetFile)

        for logFile in Path(xia2DIALSExecDir / "LogFiles").glob("*"):
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{logFile.name}"
            UtilsPath.systemCopyFile(logFile, targetFile)

        xia2Html = xia2DIALSExecDir / "xia2.html"
        if xia2Html.exists():
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{xia2Html.name}"
            UtilsPath.systemCopyFile(xia2Html, targetFile)

        xia2Txt = xia2DIALSExecDir / "xia2.txt"
        if xia2Txt.exists():
            targetFile = self.resultsDirectory / f"{self.pyarchPrefix}_{xia2Txt.name}"
            UtilsPath.systemCopyFile(xia2Txt, targetFile)

        outData = {}

        # run xia2.ispyb_json
        logger.info("Running xia2.ispyb_json... ")
        xia2JsonIspybTask = Xia2JsonIspybTask(
            inData={"xia2DialsExecDir": str(xia2DIALSExecDir)},
            workingDirectorySuffix="final",
        )
        xia2JsonIspybTask.execute()
        xia2JsonFile = xia2JsonIspybTask.outData.get("ispyb_json", None)

        if xia2JsonFile is None:
            logger.info("xia2.ispyb_json failed. Exiting")
            return outData

        logger.info("ispyb.json successfully created")

        self.resultFilePaths = list(self.resultsDirectory.iterdir())
        if inData.get("test", False):
            self.tmpdir = tempfile.TemporaryDirectory()
            pyarchDirectory = Path(self.tmpdir.name)
            self.pyarchDirectory = self.storeDataOnPyarch(pyarchDirectory=pyarchDirectory)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch()

        xia2AutoProcContainer = self.loadAndFixJsonOutput(xia2JsonFile)
        outData = xia2AutoProcContainer

        if self.doUploadIspyb:
            self.logToIspyb(
                self.integrationId, "Indexing", "Successful", "Xia2Dials finished"
            )

        if self.doUploadIspyb:
            ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(
                inData=xia2AutoProcContainer, workingDirectorySuffix="uploadFinal"
            )
            ispybStoreAutoProcResults.execute()

        if inData.get("test", False):
            self.tmpdir.cleanup()
            
        logger.info("Xia2DIALS Completed.")

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

    def loadAndFixJsonOutput(self, jsonFile, trunc_len=256):
        """fixes some of the output from the xia2 JSON output
        for submission to ISPyB."""
        autoProcContainer = {"dataCollectionId": self.dataCollectionId}
        with open(jsonFile, "r") as fp:
            jsonFile = json.load(fp)

        autoProcContainer["autoProcProgram"] = jsonFile["AutoProcProgramContainer"][
            "AutoProcProgram"
        ]
        autoProcContainer["autoProcProgram"]["autoProcProgramId"] = self.programId

        autoProcContainer["autoProcProgram"]["processingPrograms"] = autoProcContainer[
            "autoProcProgram"
        ].pop("processingProgram")
        autoProcContainer["autoProcProgram"]["processingCommandLine"] = ""
        autoProcContainer["autoProcProgram"]["processingStatus"] = "SUCCESS"
        autoProcContainer["autoProcProgram"]["processingPrograms"] = "xia2DIALS"
        autoProcContainer["autoProcProgram"]["processingStartTime"] = self.timeStart
        autoProcContainer["autoProcProgram"]["processingEndTime"] = self.timeEnd

        autoProcContainer["autoProc"] = jsonFile["AutoProc"]
        autoProcContainer["autoProc"]["autoProcProgramId"] = self.programId

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

        autoProcContainer["autoProcIntegration"] = jsonFile["AutoProcScalingContainer"][
            "AutoProcIntegrationContainer"
        ][0]["AutoProcIntegration"]
        autoProcContainer["autoProcIntegration"][
            "autoProcIntegrationId"
        ] = self.integrationId
        autoProcContainer["autoProcIntegration"]["autoProcProgramId"] = self.programId
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
        autoProcContainer["autoProcIntegration"]["anomalous"] = self.anomalous
        autoProcContainer["autoProcIntegration"][
            "dataCollectionId"
        ] = self.dataCollectionId

        # round up some values...
        for k, v in autoProcContainer["autoProcIntegration"].items():
            if isinstance(v, float):
                if any(x for x in ["cellA", "cellB", "cellC"] if x in k):
                    autoProcContainer["autoProcIntegration"][k] = round(v, 3)
                else:
                    autoProcContainer["autoProcIntegration"][k] = round(v, 2)

        autoProcContainer["autoProcScalingHasInt"] = {
            "autoProcIntegrationId": self.integrationId
        }
        autoProcContainer["autoProcScalingStatistics"] = jsonFile[
            "AutoProcScalingContainer"
        ]["AutoProcScalingStatistics"]

        try:
            Isa = self.getIsa()
        except:
            logger.error("Could not get Isa value from logs.")
            Isa = None

        for shell in autoProcContainer["autoProcScalingStatistics"]:
            shell["rmerge"] = shell.pop("rMerge") * 100
            shell["ccAno"] = shell.pop("ccAnomalous") * 100
            shell["meanIoverSigI"] = shell.pop("meanIOverSigI")
            shell["rmeasAllIplusIminus"] = shell.pop("rMeasAllIPlusIMinus") * 100
            shell["rmeasWithinIplusIminus"] = shell.pop("rMeasWithinIPlusIMinus") * 100
            shell["rpimWithinIplusIminus"] = shell.pop("rPimWithinIPlusIMinus") * 100
            shell["rpimAllIplusIminus"] = shell.pop("rPimAllIPlusIMinus") * 100
            if shell["scalingStatisticsType"] == "overall":
                shell["isa"] = round(float(Isa), 2)

            for k, v in shell.items():
                if isinstance(v, float):
                    shell[k] = round(v, 2)

        autoProcAttachmentContainerList = []
        for file in self.pyarchDirectory.iterdir():
            attachmentContainer = {
                "file": file,
            }
            autoProcAttachmentContainerList.append(attachmentContainer)

        autoProcContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList

        return autoProcContainer

    def getIsa(self):
        """extract the last estimated Isa value from _SCALE.log file."""
        Isa = None
        for scaleLogFile in self.resultsDirectory.glob("*SCALE.log"):
            with open(scaleLogFile, "r") as fp:
                for line in fp:
                    if "estimated I/sigma asymptotic limit" in line:
                        Isa = line.split()[-1]
        return Isa


class Xia2DialsExecTask(AbstractTask):
    def run(self, inData):
        logger.info("xia2DIALS Execution started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        outData = {}
        logger.debug(f"working directory is {self.getWorkingDirectory()}")
        self.onlineAutoProcessing = inData["onlineAutoProcessing"]
        self.masterFilePath = inData["masterFilePath"]
        self.imageNoStart = inData["imageNoStart"]
        self.imageNoEnd = inData["imageNoEnd"]
        self.proteinAcronym = inData["proteinAcronym"]
        self.sampleName = inData["sampleName"]

        self.spaceGroupNumber = inData.get("spaceGroupNumber", 0)
        self.spaceGroupString = inData.get("spaceGroupString", "")

        self.unitCell = inData.get("unitCell", None)
        self.lowResLimit = inData.get("lowResolutionLimit", None)
        self.highResLimit = inData.get("highResolutionLimit", None)
        self.anomalous = inData.get("anomalous", False)

        outData["workingDirectory"] = str(self.getWorkingDirectory())

        xia2DialsSetup = UtilsConfig.get("Xia2DialsTask", "xia2DialsSetup", None)
        xia2DialsExecutable = UtilsConfig.get(
            "Xia2DialsTask", "xia2DialsExecutable", "xia2"
        )
        maxNoProcessors = UtilsConfig.get(
            "Xia2DialsTask", "maxNoProcessors", os.cpu_count()
        )
        xia2DialsFastMode = distutils.util.strtobool(
            UtilsConfig.get("Xia2DialsTask", "xia2DialsFastMode", "false")
        )

        # prepare nproc, njobs for dials.integrate
        dialsIntegratePhil = self.getWorkingDirectory() / "dials_integrate.phil"
        cpusPerJob = 8
        numJobs = 4
        with open(dialsIntegratePhil, "w") as fp:
            fp.write(
                f"""integration {{
    block {{
        size = Auto
        units = *degrees radians frames
    }}
    mp {{
        nproc={cpusPerJob}
        njobs={numJobs}
    }}
}}"""
            )

        # set up command line
        if xia2DialsSetup:
            commandLine = f". {xia2DialsSetup} \n"
        else:
            commandLine = ""
        commandLine += f" {xia2DialsExecutable}"
        # add flags, if present
        commandLine += " pipeline=dials"
        if self.anomalous:
            commandLine += " atom=X"
        commandLine += (
            f" image={self.masterFilePath}:{self.imageNoStart}:{self.imageNoEnd}"
        )
        if maxNoProcessors:
            commandLine += f" nproc={int(maxNoProcessors)}"

        commandLine += f" project={self.proteinAcronym} crystal={self.sampleName}"

        if xia2DialsFastMode:
            commandLine += " dials.fast_mode=True"
        if self.spaceGroupNumber != 0:
            commandLine += f" space_group={self.spaceGroupString}"
            commandLine += " unit_cell={cell_a},{cell_b},{cell_c},{cell_alpha},{cell_beta},{cell_gamma}".format(
                **self.unitCell
            )

        if self.lowResLimit is not None or self.highResLimit is not None:
            low = self.lowResLimit if self.lowResLimit else 1000.0
            high = self.highResLimit if self.highResLimit else 0.1
            commandLine += f" resolution.d_min={low}"
            commandLine += f" resolution.d_max={high}"
        commandLine += (
            f" integrate.mosaic=new dials.integrate.phil_file={dialsIntegratePhil}"
        )
        logger.info("xia2Dials command is {}".format(commandLine))

        if self.onlineAutoProcessing:
            returncode = self.submitCommandLine(
                commandLine,
                jobName="EDNA2_xia2",
                timeout=UtilsConfig.get("Xia2DialsTask","slurm_timeout"),
                partition=UtilsConfig.get("Xia2DialsTask","slurm_partition"),
                ignoreErrors=False,
            )
            if returncode != 0:
                self.setFailure()
                return
        else:
            try:
                self.runCommandLine(commandLine, listCommand=[])
            except RuntimeError:
                self.setFailure()
                return
        return outData


class Xia2JsonIspybTask(AbstractTask):
    def run(self, inData):
        xia2DialsExecDir = inData["xia2DialsExecDir"]
        xia2DialsSetup = UtilsConfig.get("Xia2DialsTask", "xia2DialsSetup", None)

        if xia2DialsSetup:
            commandLine = f". {xia2DialsSetup} \n"
        else:
            commandLine = ""

        commandLine += f" cd {xia2DialsExecDir}; \n"
        commandLine += f" xia2.ispyb_json"

        logger.info("xia2.ispyb_json command is {}".format(commandLine))
        outData = {"ispyb_json": None}
        try:
            self.runCommandLine(commandLine, listCommand=[])
        except RuntimeError:
            self.setFailure()
            return outData
        out_file = xia2DialsExecDir + "/ispyb.json"
        outData["ispyb_json"] = str(out_file)
        return outData
