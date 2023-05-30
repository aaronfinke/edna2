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
from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, UploadGPhLResultsToISPyB
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
                timeEnd=datetime.now().isoformat(timespec="seconds")
            )
    
    def run(self, inData):
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.processingPrograms="autoproc2"
        self.processingProgramStaraniso = "autoproc_staraniso2"
        self.processingCommandLine = ""

        self.setLogFileName(f"autoPROC_{self.startDateTimeFormatted}.log")
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
                unitCellList = [float(x) for x in self.unitCell.split(",")]
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
        """ 
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
        """
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

        if inData.get("dataCollectionId") is not None:
            #set ISPyB to running
            self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine = self.processingCommandLine,
                processingPrograms = self.processingPrograms,
                isAnom = self.doAnom,
                timeStart = self.timeStart)
            self.integrationIdStaraniso, self.programIdStaraniso = ISPyBStoreAutoProcResults.setIspybToRunning(
                dataCollectionId=self.dataCollectionId,
                processingCommandLine = self.processingCommandLine,
                processingPrograms = self.processingProgramStaraniso,
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
        self.autoPROCExecDir = self.getWorkingDirectory() / "AutoPROCExec_0"
        inc_x = 1
        while self.autoPROCExecDir.is_dir():
            self.autoPROCExecDir = self.getWorkingDirectory() / "AutoPROC_exec_{0}".format(inc_x)
            inc_x += 1

        autoPROCSetup = UtilsConfig.get(self,"autoPROCSetup", None)
        autoPROCExecutable = UtilsConfig.get(self,"autoPROCExecutable", "process")
        maxNoProcessors = UtilsConfig.get(self, "maxNoProcessors", None)
        pathToNeggiaPlugin = UtilsConfig.get(self, "pathToNeggiaPlugin", None)
        autoPROCmacro = UtilsConfig.get(self, "macro", None)

        if autoPROCSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + autoPROCSetup + "\n"
        commandLine += " {0}".format(autoPROCExecutable)
        #add flags, if present
        commandLine += " -B"
        commandLine += " -d {0}".format(str(self.autoPROCExecDir))
        commandLine += " -nthreads {0}".format(maxNoProcessors)
        
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
        self.endDateTime = datetime.now().isoformat(timespec="seconds")

        ispybXml = self.autoPROCExecDir / "autoPROC.xml"
        if ispybXml.is_file():
            self.outData["ispybXml"] = str(ispybXml)
            autoProcContainer = self.autoPROCXMLtoISPyBdict(ispybXml, data_collection_id=self.dataCollectionId, 
                                                            program_id=self.programId, 
                                                            integration_id=self.integrationId)

        ispybXmlStaraniso = self.autoPROCExecDir / "autoPROC_staraniso.xml"
        if ispybXmlStaraniso.is_file():
            self.outData["ispybXml_staraniso"] = str(ispybXmlStaraniso)
            autoProcContainerStaraniso = self.autoPROCXMLtoISPyBdict(ispybXmlStaraniso, data_collection_id=self.dataCollectionId, 
                                                                     program_id=self.programIdStaraniso, 
                                                                     integration_id=self.integrationIdStaraniso)

        #get CIF Files and gzip them
        autoPROCStaranisoAllCif = self.autoPROCExecDir / "Data_1_autoPROC_STARANISO_all.cif"
        autoPROCStaranisoAllCifGz = self.resultsDirectory / f"{self.pyarchPrefix}_autoPROC_STARANISO_all.cif.gz"
        autoPROCTruncateAllCif = self.autoPROCExecDir / "Data_2_autoPROC_TRUNCATE_all.cif"
        autoPROCTruncateAllCifGz = self.resultsDirectory / f"{self.pyarchPrefix}_autoPROC_TRUNCATE_all.cif.gz"
        autoPROCXdsAsciiHkl = self.autoPROCExecDir / "XDS_ASCII.HKL"
        autoPROCXdsAsciiHklGz = self.resultsDirectory / f"{self.pyarchPrefix}_XDS_ASCII.HKL.gz"

        try:
            logger.debug(f"gzip'ing {autoPROCStaranisoAllCif}")
            with open(autoPROCStaranisoAllCif,"rb") as fp_in:
                with gzip.open(autoPROCStaranisoAllCifGz, "wb") as fp_out:
                    shutil.copyfileobj(fp_in, fp_out)
        except:
            logger.error(f"gzip'ing {autoPROCStaranisoAllCif} failed.")
        try:
            logger.debug(f"gzip'ing {autoPROCTruncateAllCif}")
            with open(autoPROCTruncateAllCif,"rb") as fp_in:
                with gzip.open(autoPROCTruncateAllCifGz, "wb") as fp_out:
                    shutil.copyfileobj(fp_in, fp_out)
        except:
            logger.error(f"gzip'ing {autoPROCTruncateAllCif} failed.")
        try:
            logger.debug(f"gzip'ing {autoPROCXdsAsciiHkl}")
            with open(autoPROCXdsAsciiHkl,"rb") as fp_in:
                with gzip.open(autoPROCXdsAsciiHklGz, "wb") as fp_out:
                    shutil.copyfileobj(fp_in, fp_out)
        except:
            logger.error(f"gzip'ing {autoPROCXdsAsciiHkl} failed.")
        
        

        #copy files to results directory
        autoPROCLogFile = self.getLogPath()
        autoPROCReportPdf = self.autoPROCExecDir / "report.pdf"
        autoPROCStaranisoReportPdf = self.autoPROCExecDir / "report_staraniso.pdf"
        autoPROCStaranisoAllDataUniqueMtz = self.autoPROCExecDir / "staraniso_alldata-unique.mtz"
        autoPROCStaranisoAllDataUniqueStats = self.autoPROCExecDir / "staraniso_alldata-unique.stats"
        autoPROCStaranisoAllDataUniqueTable1 = self.autoPROCExecDir / "staraniso_alldata-unique.table1"
        autoPROCSummaryInlinedHtml = self.autoPROCExecDir / "summary_inlined.html"
        autoPROCSummaryTarGz = self.autoPROCExecDir / "summary.tar.gz"
        autoPROCTruncateUniqueMtz = self.autoPROCExecDir / "truncate-unique.mtz"
        autoPROCTruncateUniqueStats = self.autoPROCExecDir / "truncate-unique.stats"
        autoPROCTruncateUniqueTable1 = self.autoPROCExecDir / "truncate-unique.table1"

        autoPROCLogFile_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_autoPROC.log"
        autoPROCReportPdf_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_report.pdf"
        autoPROCStaranisoReportPdf_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_report_staraniso.pdf"
        autoPROCStaranisoAllDataUniqueMtz_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_staraniso_alldata-unique.mtz"
        autoPROCStaranisoAllDataUniqueStats_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_staraniso_alldata-unique.stats"
        autoPROCStaranisoAllDataUniqueTable1_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_staraniso_alldata-unique.table1"
        autoPROCSummaryInlinedHtml_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_summary_inlined.html"
        autoPROCSummaryTarGz_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_summary.tar.gz"
        autoPROCTruncateUniqueMtz_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_truncate-unique.mtz"
        autoPROCTruncateUniqueStats_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_truncate-unique.stats"
        autoPROCTruncateUniqueTable1_resultsDir = self.resultsDirectory / f"{self.pyarchPrefix}_truncate-unique.table1"

        autoProcAttachmentContainerList = []
        autoProcAttachmentContainerStaranisoList = []

        for files in [(autoPROCLogFile, autoPROCLogFile_resultsDir),
                      (autoPROCReportPdf, autoPROCReportPdf_resultsDir),
                      (autoPROCStaranisoReportPdf,autoPROCStaranisoReportPdf_resultsDir),
                      (autoPROCStaranisoAllDataUniqueMtz,autoPROCStaranisoAllDataUniqueMtz_resultsDir),
                      (autoPROCStaranisoAllDataUniqueStats,autoPROCStaranisoAllDataUniqueStats_resultsDir),
                      (autoPROCStaranisoAllDataUniqueTable1, autoPROCStaranisoAllDataUniqueTable1_resultsDir),
                      (autoPROCSummaryInlinedHtml,autoPROCSummaryInlinedHtml_resultsDir),
                      (autoPROCSummaryTarGz,autoPROCSummaryTarGz_resultsDir),
                      (autoPROCTruncateUniqueMtz,autoPROCTruncateUniqueMtz_resultsDir),
                      (autoPROCTruncateUniqueStats,autoPROCTruncateUniqueStats_resultsDir),
                      (autoPROCTruncateUniqueTable1,autoPROCTruncateUniqueTable1_resultsDir)]:
            UtilsPath.systemCopyFile(files[0],files[1])
            # pyarchFile = UtilsPath.createPyarchFilePath(files[1])
            pyarchFile = files[1]
            attachmentContainer = {
                "file" : pyarchFile,
            }
            if "staraniso" in str(files[1]):
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
            elif "truncate-unique" in str(files[1]):
                autoProcAttachmentContainerList.append(attachmentContainer)
            else:
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
                autoProcAttachmentContainerList.append(attachmentContainer)
        
        for file in [autoPROCStaranisoAllCifGz, autoPROCTruncateAllCifGz, autoPROCXdsAsciiHklGz]:
            # pyarchFile = UtilsPath.createPyarchFilePath(files[1])
            pyarchFile = file
            attachmentContainer = {
                "file" : pyarchFile,
            }
            if "STARANISO" in str(file):
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
            elif "TRUNCATE" in str(file):
                autoProcAttachmentContainerList.append(attachmentContainer)
            else:
                autoProcAttachmentContainerStaranisoList.append(attachmentContainer)
                autoProcAttachmentContainerList.append(attachmentContainer)

        
        autoProcContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList
        autoProcContainerStaraniso["autoProcProgramAttachment"] = autoProcAttachmentContainerStaranisoList
        
        #save as json
        autoProcContainerJson = self.resultsDirectory / "autoPROC.json"
        autoProcContainerStaranisoJson = self.resultsDirectory / "autoPROC_staraniso.json"

        with open(autoProcContainerJson,'w') as fp:
            json.dump(autoProcContainer,fp, indent=2, default=lambda o:str(o))
        
        with open(autoProcContainerStaranisoJson,'w') as fp:
            json.dump(autoProcContainerStaraniso,fp, indent=2, default=lambda o:str(o))

        ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=autoProcContainer, workingDirectorySuffix="uploadFinal")
        ispybStoreAutoProcResults.execute()
        ispybStoreAutoProcResultsStaraniso = ISPyBStoreAutoProcResults(inData=autoProcContainerStaraniso, workingDirectorySuffix="uploadFinal_staraniso")
        ispybStoreAutoProcResultsStaraniso.execute()

        if self.tmpdir is not None:
            self.tmpdir.cleanup()
        return 

    @staticmethod
    def autoPROCXMLtoISPyBdict(xml_path, data_collection_id=None, program_id=None, integration_id=None, trunc_len=256):
        dict_data = UtilsXML.dictfromXML(xml_path)
        autoProcXMLContainer = dict_data["AutoProcContainer"]
        autoProcContainer = {
            "dataCollectionId": data_collection_id
        }
        autoProcContainer["autoProcProgram"] = {}
        for k,v in autoProcXMLContainer["AutoProcProgramContainer"]["AutoProcProgram"].items():
                if isinstance(v,str) and len(v) > trunc_len:
                    autoProcContainer["autoProcProgram"][k] = v[:trunc_len-1]
                    logger.warning(f"string {k} truncated for loading to ISPyB to {trunc_len} characters: \"{v}\"")
                else:
                    autoProcContainer["autoProcProgram"][k] = v
        #fix some entries in autoProcProgram
        autoProcContainer["autoProcProgram"]["processingStatus"] = "SUCCESS"

        timeStart = autoProcContainer["autoProcProgram"]["processingStartTime"]
        timeEnd = autoProcContainer["autoProcProgram"]["processingEndTime"]
        timeStart_structtime = time.strptime(timeStart, "%a %b %d %H:%M:%S %Z %Y")
        timeEnd_structtime = time.strptime(timeEnd, "%a %b %d %H:%M:%S %Z %Y")
        autoProcContainer["autoProcProgram"]["processingStartTime"] = time.strftime('%Y-%m-%dT%H:%M:%SZ',timeStart_structtime)
        autoProcContainer["autoProcProgram"]["processingEndTime"] = time.strftime('%Y-%m-%dT%H:%M:%SZ',timeEnd_structtime)


        autoProcContainer["autoProcProgram"]["autoProcProgramId"] = program_id
        autoProcContainer["autoProc"] = autoProcXMLContainer["AutoProc"]
        #fix some entries in autoProc
        autoProcContainer["autoProc"]["refinedCellA"] = autoProcContainer["autoProc"].pop("refinedCell_a")
        autoProcContainer["autoProc"]["refinedCellB"] = autoProcContainer["autoProc"].pop("refinedCell_b")
        autoProcContainer["autoProc"]["refinedCellC"] = autoProcContainer["autoProc"].pop("refinedCell_c")
        autoProcContainer["autoProc"]["refinedCellAlpha"] = autoProcContainer["autoProc"].pop("refinedCell_alpha")
        autoProcContainer["autoProc"]["refinedCellBeta"] = autoProcContainer["autoProc"].pop("refinedCell_beta")
        autoProcContainer["autoProc"]["refinedCellGamma"] = autoProcContainer["autoProc"].pop("refinedCell_gamma")

        autoProcContainer["autoProc"]["autoProcProgramId"] = program_id
        autoProcContainer["autoProcScaling"] = autoProcXMLContainer["AutoProcScalingContainer"]["AutoProcScaling"]
        autoProcContainer["autoProcScalingStatistics"] = autoProcXMLContainer["AutoProcScalingContainer"]["AutoProcScalingStatistics"]
        # fix some entries in autoProcScalingStatistics
        for shell in autoProcContainer["autoProcScalingStatistics"]:
            shell["ccAno"] = shell.pop("ccAnomalous")
            shell["sigAno"] = shell.pop("DanoOverSigDano")

        autoProcContainer["AutoProcScalingIntegration"] = autoProcXMLContainer["AutoProcScalingContainer"]["AutoProcIntegrationContainer"]["AutoProcIntegration"]
        autoProcContainer["AutoProcScalingIntegration"]["autoProcProgramId"] = program_id
        autoProcContainer["AutoProcScalingIntegration"]["autoProcIntegrationId"] = integration_id
        #fix some entries in AutoProcScalingIntegration
        autoProcContainer["AutoProcScalingIntegration"]["cellA"] = autoProcContainer["AutoProcScalingIntegration"].pop("cell_a")
        autoProcContainer["AutoProcScalingIntegration"]["cellB"] = autoProcContainer["AutoProcScalingIntegration"].pop("cell_b")
        autoProcContainer["AutoProcScalingIntegration"]["cellC"] = autoProcContainer["AutoProcScalingIntegration"].pop("cell_c")
        autoProcContainer["AutoProcScalingIntegration"]["cellAlpha"] = autoProcContainer["AutoProcScalingIntegration"].pop("cell_alpha")
        autoProcContainer["AutoProcScalingIntegration"]["cellBeta"] = autoProcContainer["AutoProcScalingIntegration"].pop("cell_beta")
        autoProcContainer["AutoProcScalingIntegration"]["cellGamma"] = autoProcContainer["AutoProcScalingIntegration"].pop("cell_gamma")
        autoProcContainer["autoProcIntegration"] = autoProcContainer.pop("AutoProcScalingIntegration")
        return autoProcContainer


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


