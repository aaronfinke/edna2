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

from edna2.tasks.XDSTasks import XDSIndexing, XDSIntegration, XDSRerunCorrect
from edna2.tasks.SubWedgeAssembly import SubWedgeAssembly
from edna2.tasks.CCP4Tasks import PointlessTask, AimlessTask, TruncateTask, UniqueifyTask
from edna2.tasks.XSCALETasks import XSCALETask
from edna2.tasks.PhenixTasks import PhenixXTriageTask
from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, ISPyBStoreAutoProcStatus
from edna2.tasks.WaitFileTask import WaitFileTask


class Edna2ProcTask(AbstractTask):
    """
    Runs XDS both with and without anomalous
    treament of reflections, given the path
    to a master file.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "required": ["dataCollectionId"],
            "properties":
            {
                "dataCollectionId": {"type":["integer","null"]},
                "masterFilePath": {"type":["string","null"]},
                "spaceGroup": {"type":["integer","string"]},
                "unitCell": {"type":["string","null"]},
                "residues": {"type":["integer","null"]},
                "anomalous": {"type":["boolean","null"]},
                "imageNoStart": {"type":["integer","null"]},
                "imageNoEnd": {"type":["integer","null"]},
                "workingDirectory": {"type":["string","null"]},
            }
        }
    
    #sets ISPyB to FAILED if it's already logged
    def setFailure(self):
        self._dictInOut["isFailure"] = True
        if self.integrationId is not None and self.programId is not None:
            ISPyBStoreAutoProcResults.setIspybToFailed(
                dataCollectionId=self.dataCollectionId,
                autoProcProgramId=self.programId,
                autoProcIntegrationId=self.integrationId,
                processingCommandLine=self.processingCommandLine, 
                processingPrograms=self.processingPrograms, 
                isAnom=self.anomalous, 
                timeStart=self.startDateTime, 
                timeEnd=datetime.now().isoformat(timespec='seconds')
            )



    def run(self, inData):
        UtilsLogging.addLocalFileHandler(logger, self.getWorkingDirectory()/"EDNA2Proc.log")
        self.tmpdir = None
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec='seconds')
        self.processingPrograms="EDNA2_proc"
        self.processingCommandLine = ""
        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.anomalous = inData.get("anomalous",False)
        self.anomalousFlag = False
        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        self.imageNoStart = inData.get("imageNoStart",None)
        self.imageNoEnd = inData.get("imageNoEnd",None)
        self.masterFilePath = inData.get("masterFilePath",None)
        outData = {}
        self.resultFilePaths = []

        logger.info("EDNA2 Auto Processing started")
        logger.info(f"Running on {socket.gethostname()}")
        try:
            logger.info(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        #set up SG and unit cell
        self.spaceGroupNumber, self.spaceGroupString = UtilsCCTBX.parseSpaceGroup(self.spaceGroup)

        # set up unit cell
        if self.unitCell is not None:
            self.unitCell = UtilsCCTBX.parseUnitCell(self.unitCell)
        else:
            logger.info("No unit cell supplied")

        if self.onlineAutoProcessing:
            if self.dataCollectionId is None:
                logger.error("No dataCollectionId, exiting.")
                self.setFailure()
                return 
            try:
                dataCollectionWS3VO = UtilsIspyb.findDataCollection(self.dataCollectionId)
                ispybDataCollection = dict(dataCollectionWS3VO)
                logger.debug("ispybDataCollection: {}".format(ispybDataCollection))
                self.imageDirectory = ispybDataCollection.get("imageDirectory")
                if UtilsConfig.isEMBL():
                    self.fileTemplate = ispybDataCollection["fileTemplate"].replace("%05d", "#" * 5)
                elif UtilsConfig.isMAXIV():
                    self.fileTemplate = ispybDataCollection["fileTemplate"]
                else:
                    self.fileTemplate = ispybDataCollection["fileTemplate"].replace("%04d", "####")
                self.imageNoStart = ispybDataCollection["startImageNumber"]
                numImages = ispybDataCollection["numberOfImages"]
                self.imageNoEnd = numImages - self.imageNoStart + 1
                pathToStartImage = os.path.join(self.imageDirectory,
                                                UtilsImage.eiger_template_to_image(self.fileTemplate, self.imageNoStart))
                pathToEndImage = os.path.join(self.imageDirectory,
                                                UtilsImage.eiger_template_to_image(self.fileTemplate, self.imageNoEnd))
            except:
                logger.warning("Retrieval of data from ISPyB Failed, trying manually...")
                if self.masterFilePath is None:
                    logger.error("No dataCollectionId or masterFilePath found, exiting")
                    self.setFailure()
                    return
                self.masterFilePath = Path(self.masterFilePath)
                self.imageDirectory = self.masterFilePath.parent
                masterFileName = self.masterFilePath.name
                try:
                    self.fileTemplate = re.search(r"[A-Za-z0-9-_]+(?=_master\.h5)", masterFileName)[0]
                except:
                    logger.error(f"File template not found: {masterFileName}")
                    self.setFailure()
                    return
                numImages = UtilsImage.getNumberOfImages(masterFilePath=self.masterFilePath)
                self.imageNoStart = inData.get("imageNoStart",1)
                self.imageNoEnd = inData.get("imageNoEnd",numImages)
                logger.debug(f"imageNoStart: {self.imageNoStart}, imageNoEnd: {self.imageNoEnd}")
                pathToStartImage = self.imageDirectory / Path(self.fileTemplate + f"_data_{self.imageNoStart:06d}.h5")
                lastFileNumber = int(math.ceil(self.imageNoEnd / 100.0))
                pathToEndImage = self.imageDirectory / Path(self.fileTemplate + f"_data_{lastFileNumber:06d}.h5")
        else:
            if self.masterFilePath is None:
                logger.error("No dataCollectionId or masterFilePath found, exiting")
                self.setFailure()
                return
            self.masterFilePath = Path(self.masterFilePath)
            self.imageDirectory = self.masterFilePath.parent
            masterFileName = self.masterFilePath.name
            try:
                self.fileTemplate = re.search(r"[A-Za-z0-9-_]+(?=_master\.h5)", masterFileName)[0]
            except:
                logger.error(f"File template not found: {masterFileName}")
                self.setFailure()
                return
            numImages = UtilsImage.getNumberOfImages(masterFilePath=self.masterFilePath)
            self.imageNoStart = inData.get("imageNoStart",1)
            self.imageNoEnd = inData.get("imageNoEnd",numImages)
            logger.debug(f"imageNoStart: {self.imageNoStart}, imageNoEnd: {self.imageNoEnd}")
            pathToStartImage = self.imageDirectory / Path(self.fileTemplate + f"_data_{self.imageNoStart:06d}.h5")
            lastFileNumber = int(math.ceil(self.imageNoEnd / 100.0))
            pathToEndImage = self.imageDirectory / Path(self.fileTemplate + f"_data_{lastFileNumber:06d}.h5")

        if self.imageNoEnd - self.imageNoStart < 8:
            #if self.imageNoEnd - self.imageNoStart < -1:
            logger.error("There are fewer than 8 images, aborting")
            self.setFailure()
            return
        logger.info(f"dataCollectionId: {self.dataCollectionId}")

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
            isH5 = True
        elif UtilsConfig.isMAXIV():
            minSizeFirst = 100000
            minSizeLast = 100000
            isH5 = True
        else:
            minSizeFirst = 1000000
            minSizeLast = 1000000        
        
        listPrefix = self.fileTemplate.split("_")
        if UtilsConfig.isALBA():
            self.pyarchPrefix = "ap_{0}_{1}".format("_".join(listPrefix[:-2]),
                                                       listPrefix[-2])
        elif UtilsConfig.isMAXIV():
            self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])
        else:
            if len(listPrefix) > 2:
                self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[-3], listPrefix[-2])
            elif len(listPrefix) > 1:
                self.pyarchPrefix = "ap_{0}_run{1}".format(listPrefix[:-2], listPrefix[-2])
            else:
                self.pyarchPrefix = "ap_{0}_run".format(listPrefix[0])

        if self.masterFilePath is None:
            self.masterFilePath = os.path.join(self.imageDirectory,
                    UtilsImage.eiger_template_to_master(self.fileTemplate))



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

        self.integrationId = None
        self.programId = None
        
        #get integrationIDs and programIDs, set them to running
        try: 
            self.integrationId,self.programId = self.createIntegrationId(
                                    "Creating integration ID", 
                                    isAnom=self.anomalous)
            logger.info(f"integrationID: {self.integrationId}, programId: {self.programId}")
        except Exception as e:
            logger.error("Could not get integration ID: \n{0}".format(traceback.format_exc(e)))
        
        #set working directory, results directory, log file
        workingDirectory = self.getWorkingDirectory()
        self.resultsDirectory = Path(workingDirectory / "results")
        self.resultsDirectory.mkdir(parents=True, exist_ok=True)
        self.setLogFileName(f'edna2_proc.log')

        #XXX todo: I hate this
        if not inData.get("subWedge"):  
            logger.info("Generating subwedge...")
            #num_of_images, image_list = Edna2ProcTask.generateImageListFromH5Master(inData=inData)
            self.imgNumLow, self.imgNumHigh, self.imageList = self.generateImageListFromH5Master_fast(masterFilePath=self.masterFilePath)
            self.subWedgeAssembly = SubWedgeAssembly(inData=self.imageList)
            self.subWedgeAssembly.execute()
            self.xdsIndexingInData = self.subWedgeAssembly.outData
        else:
            self.xdsIndexingInData = inData
        self.numImages = self.imgNumHigh - self.imgNumLow + 1
        logger.info(f"Number of images: {self.numImages}")
        if self.numImages < 8:
            logger.error("There are fewer than 8 images. Aborting.")
            self.setFailure()
            return
        
        self.xdsIndexingInData["unitCell"] = self.unitCell
        self.xdsIndexingInData["spaceGroupNumber"] = self.spaceGroupNumber
        self.xdsIndexingInData["isAnom"] = self.anomalous
        self.xdsIndexingInData["onlineAutoProcessing"] = self.onlineAutoProcessing

        self.indexing = XDSIndexing(inData=self.xdsIndexingInData, workingDirectorySuffix="init")

        logger.info('XDS Indexing started')

        self.logToIspyb(self.integrationId,
                            'Indexing', 'Launched', 'XDS started')

        self.indexing.execute()
        
        time1 = time.perf_counter()
        self.timeXdsIndexing = time1 - self.timeStart

        if self.indexing.isFailure():
            logger.error("Error at indexing step. Stopping.")
            self.logToIspyb(self.integrationId,
                'Indexing',
                'Failed',
                'XDS failed after {0:.1f}s'.format(self.timeXdsIndexing))
            self.setFailure()
            return 
        else: 
            self.logToIspyb(self.integrationId,
                'Indexing',
                'Successful',
                'XDS finished after {0:.1f}s'.format(self.timeXdsIndexing))
            logger.info(f"XDS indexing time took {self.timeXdsIndexing:0.1f} seconds")
            logger.info("Indexing successful. a= {cell_a}, b= {cell_b}, c= {cell_c}, al = {cell_alpha}, be = {cell_beta}, ga = {cell_gamma}" \
                        .format(**self.indexing.outData["idxref"]["unitCell"]))

        #Now set up integration
        self.integrationInData = self.indexing.outData
        del self.integrationInData["workingDirectory"]
        self.integrationInData["subWedge"] = self.indexing.inData["subWedge"]
        self.integrationInData["onlineAutoProcessing"] = self.onlineAutoProcessing
        self.integration = XDSIntegration(inData=self.integrationInData, workingDirectorySuffix="init")
        logger.info('Integration started')

        self.logToIspyb(self.integrationId,
                            'Integration', 'Launched', 'XDS started')

        self.integration.execute()

        time2 = time.perf_counter()
        self.timeXdsIntegration = time2-time1

        if self.integration.isFailure():
            logger.error("Error at integration step. Stopping.")
            self.logToIspyb(self.integrationId,
                'Integration',
                'Failed',
                'XDS failed after {0:.1f}s'.format(self.timeXdsIntegration))
            self.setFailure()
            return 
        else: 
            self.logToIspyb(self.integrationId,
                'Integration',
                'Successful',
                'XDS finished after {0:.1f}s'.format(self.timeXdsIntegration))
            logger.info(f"XDS integration time took {self.timeXdsIntegration:0.1f} seconds")
            logger.info("Integration Successful.")

        #copy the XDS.INP file from the successful run into the results directory.
        xds_INP_result_path = self.resultsDirectory / f"{self.pyarchPrefix}_successful_XDS.INP"
        integrateLp_path = self.resultsDirectory / f"{self.pyarchPrefix}_INTEGRATE.LP"
        correctLp_path = self.resultsDirectory / f"{self.pyarchPrefix}_CORRECT.LP"
        integrateHkl_path = self.resultsDirectory / f"{self.pyarchPrefix}_INTEGRATE.HKL"
        xdsAsciiHkl_path = self.resultsDirectory / f"{self.pyarchPrefix}_XDS_ASCII.HKL"
        try:
            shutil.copy(Path(self.integration.outData["integrateHkl"]), integrateHkl_path)
            shutil.copy(Path(self.integration.outData["xdsInp"]), xds_INP_result_path)
            shutil.copy(Path(self.integration.outData["integrateLp"]), integrateLp_path)
            shutil.copy(Path(self.integration.outData["correctLp"]), correctLp_path)
            shutil.copy(Path(self.integration.outData["xdsAsciiHkl"]), xdsAsciiHkl_path)
        except Exception as e:
            logger.warning("Couldn't copy file to result directory")
            logger.warning(e)

        self.resultFilePaths.extend([xds_INP_result_path,integrateLp_path,correctLp_path,integrateHkl_path,xdsAsciiHkl_path])

        #calculate resolution cutoff

        logger.info("Starting first resolution cutoff...")
        self.completenessEntries = self.integration.outData['completenessEntries']
        self.firstResCutoff = self.getResCutoff(self.completenessEntries)
        if self.firstResCutoff is None:
            logger.error('No bins with CC1/2 greater than 30%')
            logger.error("Something could be wrong, or the completeness could be too low!")
            logger.error("bravais lattice/SG could be incorrect or something more insidious like")
            logger.error("incorrect parameters in XDS.INP like distance, X beam, Y beam, etc.")
            logger.error("Stopping")
            self.setFailure()
            return
        logger.info(f"Resolution cutoff is {self.firstResCutoff}")

        # run pointless 
        pointlessTaskinData = {
            'input_file' :  self.integration.outData["xdsAsciiHkl"],
            'output_file' : "ep_pointless_unmerged.mtz"
        }
        logger.info("Starting pointless task...")
        self.pointlessTask = PointlessTask(inData=pointlessTaskinData, workingDirectorySuffix="init")
        self.pointlessTask.execute()

        logger.debug(f"Pointless output: {self.pointlessTask.outData}")

        #now rerun CORRECT with corrected parameters
        rerunCor_data = {
            "xdsInp" : self.integration.outData["xdsInp"],
            "spotXds": self.integration.inData["spotXds"],
            "gxParmXds" : self.integration.outData["gxParmXds"],
            "gainCbf": self.integration.inData["gainCbf"],
            "blankCbf":self.integration.inData["blankCbf"],
            "xCorrectionsCbf": self.integration.inData["xCorrectionsCbf"],
            "yCorrectionsCbf": self.integration.inData["yCorrectionsCbf"],
            "bkginitCbf": self.integration.inData["bkginitCbf"],
            "integrateLp": self.integration.outData["integrateLp"],
            "integrateHkl": self.integration.outData["integrateHkl"],
            "sg_nr_from_pointless": self.pointlessTask.outData["sgnumber"],
            "cell_from_pointless" : self.pointlessTask.outData["cell"],
            "subWedge" : self.integration.inData["subWedge"],
            "resCutoff": self.firstResCutoff,
            "onlineAutoProcessing": self.onlineAutoProcessing,
            "isAnom": self.anomalous
        }


        logger.info("Rerunning CORRECT with the unit cell/SG from POINTLESS...")
        self.xdsRerun = XDSRerunCorrect(inData=rerunCor_data, workingDirectorySuffix='0')

        self.xdsRerun.execute()

        time3 = time.perf_counter()
        self.timeRerunCorrect = time3 - time2

        if self.xdsRerun.isFailure():
            logger.error('Rerun of CORRECT failed')
            self.setFailure()
            self.logToIspyb(self.integrationId,
                         'Scaling',
                         'Failed',
                         'Scaling failed after {0:.1}s'.format(self.timeRerunCorrect))
            return
        else:
            logger.info(f"Rerun of CORRECT finished.")
            self.logToIspyb(self.integrationId,
                         'Scaling',
                         'Successful',
                         'Scaling finished in {0:.1f}s'.format(self.timeRerunCorrect))

        logger.info("Starting second resolution cutoff...")
        self.completenessEntries = self.xdsRerun.outData['completenessEntries']
        
        self.resCutoff = self.getResCutoff(self.completenessEntries)
        if self.resCutoff is None:
            logger.error("Error in determining resolution after CORRECT rerun.")
            logger.error('No bins with CC1/2 greater than 30%')
            logger.error("Something could be wrong, or the completeness could be too low!")
            logger.error("bravais lattice/SG could be incorrect or something more insidious like")
            logger.error("incorrect parameters in XDS.INP like distance, X beam, Y beam, etc.")
            logger.error("Stopping")
            self.setFailure()
            self.logToIspyb(self.integrationId,
                            'Scaling',
                            'Failed',
                            'resolution cutoffs failed')
            self.setFailure()
            return
        logger.info(f"Resolution cutoff is {self.firstResCutoff}")
        self.logToIspyb(self.integrationId,
                        'Scaling',
                        'Successful',
                        'Resolution cutoffs finished')

        self.bins = [x["res"] for x in self.xdsRerun.outData["completenessEntries"] if x["include_res_based_on_cc"] is True]

        self.pointlessTaskReruninData = {
            'input_file' :  self.xdsRerun.outData["xdsAsciiHkl"],
            'output_file' : "ep__pointless_unmerged.mtz"
        }
        logger.info("Starting pointless tasks...")
        self.pointlessTaskRerun = PointlessTask(inData=self.pointlessTaskReruninData, workingDirectorySuffix="rerun")

        self.pointlessTaskRerun.execute()
        
        if self.pointlessTaskRerun.isFailure():
            logger.error("Pointless task failed.")
            self.setFailure()
            return

        self.aimlessTaskinData = {
            "input_file" : self.pointlessTaskRerun.outData["pointlessUnmergedMtz"],
            "output_file" : "ep__aimless.mtz",
            "start_image" : self.indexing.outData["start_image"],
            "end_image" : self.indexing.outData["end_image"],
            "dataCollectionId" : self.dataCollectionId,
            "res" : self.resCutoff,
            "anomalous" : True,
        }
        logger.info("Starting aimless...")
        self.aimlessTask = AimlessTask(inData=self.aimlessTaskinData, workingDirectorySuffix="0")
        self.aimlessTask.execute()

        logger.info(f"Aimless finished.")

        highAnomSignalFound = self.if_anomalous_signal(self.aimlessTask.outData["aimlessLog"],threshold=1.0)
        
        if highAnomSignalFound:
            logger.info("Significant anomalous signal for this dataset.")
            outData["HighAnomSignal"] = True
        else:
            logger.info("Insufficient anomalous signal for this dataset.")
            outData["HighAnomSignal"] = False

        if highAnomSignalFound and not self.anomalous:
            self.anomalousFlag = True
            logger.info("Rerunning CORRECT/AIMLESS with anomalous flags due to significant anomalous signal found.")
            self.anomalous = True
            rerunCor_Anomdata = {
                "xdsInp" : self.integration.outData["xdsInp"],
                "spotXds": self.integration.inData["spotXds"],
                "gxParmXds" : self.integration.outData["gxParmXds"],
                "gainCbf": self.integration.inData["gainCbf"],
                "blankCbf":self.integration.inData["blankCbf"],
                "xCorrectionsCbf": self.integration.inData["xCorrectionsCbf"],
                "yCorrectionsCbf": self.integration.inData["yCorrectionsCbf"],
                "bkginitCbf": self.integration.inData["bkginitCbf"],
                "integrateLp": self.integration.outData["integrateLp"],
                "integrateHkl": self.integration.outData["integrateHkl"],
                "sg_nr_from_pointless": self.pointlessTask.outData["sgnumber"],
                "cell_from_pointless" : self.pointlessTask.outData["cell"],
                "subWedge" : self.integration.inData["subWedge"],
                "resCutoff": self.firstResCutoff,
                "onlineAutoProcessing": self.onlineAutoProcessing,
                "isAnom": True
            }


            logger.info("Rerunning CORRECT with the unit cell/SG from POINTLESS and with anomalous flag on...")
            self.xdsRerunAnom = XDSRerunCorrect(inData=rerunCor_Anomdata, workingDirectorySuffix='anom')

            self.xdsRerunAnom.execute()

            if self.xdsRerun.isFailure():
                logger.error('Rerun of CORRECT failed')
                self.setFailure()
                self.logToIspyb(self.integrationId,
                            'Scaling',
                            'Failed',
                            'Scaling failed after {0:.1}s'.format(self.timeRerunCorrect))
                return
            else:
                logger.info(f"Rerun of CORRECT finished.")
                self.logToIspyb(self.integrationId,
                            'Scaling',
                            'Successful',
                            'Scaling finished in {0:.1f}s'.format(self.timeRerunCorrect))

            logger.info("Starting third resolution cutoff...")
            self.completenessEntries = self.xdsRerun.outData['completenessEntries']
            
            self.resCutoff = self.getResCutoff(self.completenessEntries)
            if self.resCutoff is None:
                logger.error("Error in determining resolution after CORRECT rerun.")
                logger.error('No bins with CC1/2 greater than 30%')
                logger.error("Something could be wrong, or the completeness could be too low!")
                logger.error("bravais lattice/SG could be incorrect or something more insidious like")
                logger.error("incorrect parameters in XDS.INP like distance, X beam, Y beam, etc.")
                logger.error("Stopping")
                self.setFailure()
                self.logToIspyb(self.integrationId,
                                'Scaling',
                                'Failed',
                                'resolution cutoffs failed')
                self.setFailure()
                return
            logger.info(f"Resolution cutoff is {self.firstResCutoff}")
            self.logToIspyb(self.integrationId,
                            'Scaling',
                            'Successful',
                            'Resolution cutoffs finished')

            self.bins = [x["res"] for x in self.xdsRerunAnom.outData["completenessEntries"] if x["include_res_based_on_cc"] is True]

            self.pointlessTaskRerunAnominData = {
                'input_file' :  self.xdsRerun.outData["xdsAsciiHkl"],
                'output_file' : "ep__pointless_unmerged.mtz"
            }
            logger.info("Starting pointless tasks...")
            self.pointlessTaskRerunAnom = PointlessTask(inData=self.pointlessTaskRerunAnominData, workingDirectorySuffix="rerunAnom")

            self.pointlessTaskRerunAnom.execute()
            
            self.aimlessTaskInDataAnom = {
                "input_file" : self.pointlessTaskRerunAnom.outData["pointlessUnmergedMtz"],
                "output_file" : "ep__aimless.mtz",
                "start_image" : self.indexing.outData["start_image"],
                "end_image" : self.indexing.outData["end_image"],
                "dataCollectionId" : self.dataCollectionId,
                "res" : self.resCutoff,
                "anomalous" : True,
            }
            logger.info("Starting aimless with anomalous flag on...")
            self.aimlessTaskAnom = AimlessTask(inData=self.aimlessTaskInDataAnom, workingDirectorySuffix="anom")
            self.aimlessTaskAnom.execute()

            logger.info(f"Aimless with anomalous flag finished.")

            self.xdsRerun = self.xdsRerunAnom
            self.pointlessTaskRerun = self.pointlessTaskRerunAnom
            self.aimlessTask = self.aimlessTaskAnom


        pointlessUnmergedMtzPath = self.resultsDirectory / (f"{self.pyarchPrefix}_ep__pointless_unmerged.mtz")
        self.resultFilePaths.append(pointlessUnmergedMtzPath)

        try:
            shutil.copy(Path(self.pointlessTaskRerun.outData["pointlessUnmergedMtz"]), pointlessUnmergedMtzPath)
        except Exception as e:
            logger.warning("Couldn't copy file to result directory")
            logger.warning(e)


        aimlessMergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless.mtz"
        aimlessUnmergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless_unmerged.mtz.gz"
        aimlessLogPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless.log"

        self.resultFilePaths.extend([aimlessMergedMtzPath,
                                     aimlessUnmergedMtzPath,aimlessLogPath])

        try:
            shutil.copy(Path(self.aimlessTask.outData["aimlessMergedMtz"]), aimlessMergedMtzPath)
            shutil.copy(Path(self.aimlessTask.outData["aimlessUnmergedMtz"]), aimlessUnmergedMtzPath)
            shutil.copy(Path(self.aimlessTask.outData["aimlessLog"]), aimlessLogPath)
        except Exception as e:
            logger.warning("Couldn't copy file to result directory")
            logger.warning(e)


        # self.timeXscaleStart = time.perf_counter()
        # self.xscaleTaskData = {
        #     "xdsAsciiPath": self.xdsRerun.outData["xdsAsciiHkl"],
        #     "bins" : self.bins,
        #     "sgNumber": self.pointlessTask.outData["sgnumber"],
        #     "cell" : self.pointlessTask.outData["cell"],
        # }

        # self.xscaleTaskData_merge = self.xscaleTaskData
        # self.xscaleTaskData_merge['isAnom'] = self.anomalous
        # self.xscaleTaskData_merge['merge'] = True
        # self.xscaleTaskData_merge['res'] = self.resCutoff
        # self.xscaleTaskData_merge['onlineAutoProcessing'] = self.onlineAutoProcessing


        # self.xscaleTask_merge = XSCALETask(inData=self.xscaleTaskData_merge, workingDirectorySuffix="merged")

        # self.xscaleTaskData_unmerge = self.xscaleTaskData
        # self.xscaleTaskData_unmerge['isAnom'] = self.anomalous
        # self.xscaleTaskData_unmerge['merge'] = False
        # self.xscaleTaskData_unmerge['res'] = self.resCutoff
        # self.xscaleTaskData_unmerge['onlineAutoProcessing'] = self.onlineAutoProcessing

        # self.xscaleTask_unmerge = XSCALETask(inData=self.xscaleTaskData_unmerge, workingDirectorySuffix="unmerged")

        # logger.info("Starting XSCALE merging...")
        # self.logToIspyb(self.integrationId,
        #              'Scaling', 'Launched', 'Start of XSCALE')

        # self.xscaleTask_merge.start()
        # self.xscaleTask_unmerge.start()

        # self.xscaleTask_merge.join()
        # self.xscaleTask_unmerge.join()

        # time4 = time.perf_counter()
        # self.timeXscale = time4-time3

        # for task in [self.xscaleTask_merge,self.xscaleTask_unmerge]:
        #     if task.isFailure():
        #         logger.error("XSCALE generation failed")
        #         self.logToIspyb(self.integrationId,
        #                  'Scaling',
        #                  'Failed',
        #                  'XSCALE failed after {0:.1f}s'.format(self.timeXscale))
        #         self.setFailure()
        #         return
        
        # logger.info("XSCALE generation finished.")
        # self.logToIspyb(self.integrationId,
        #         'Scaling',
        #         'Successful',
        #         'XSCALE finished in {0:.1f}s'.format(self.timeXscale))

        # xscaleTask_mergeLPFile = self.resultsDirectory / "ep__merged_XSCALE.LP"
        # xscaleTask_unmergeLPFile = self.resultsDirectory / "ep__unmerged_XSCALE.LP"
        # self.resultFilePaths.extend([xscaleTask_mergeLPFile,
        #                              xscaleTask_unmergeLPFile])
        # try:
        #     shutil.copy(self.xscaleTask_merge.outData["xscaleLp"], xscaleTask_mergeLPFile)
        #     shutil.copy(self.xscaleTask_unmerge.outData["xscaleLp"], xscaleTask_unmergeLPFile)
        # except Exception as e:
        #     logger.warning("Couldn't copy file to result directory")
        #     logger.warning(e)


        # logger.debug(f"XSCALE output: {self.xscaleTask_merge.outData}")

        
        logger.info("Start phenix.xtriage run...")
        self.phenixXTriageTaskData = {
            "input_file" : aimlessUnmergedMtzPath
        }
        self.phenixXTriageTask = PhenixXTriageTask(inData=self.phenixXTriageTaskData)
        self.phenixXTriageTask.start()

        #now run truncate/unique
        truncateData = {}
        truncateData["inputFile"] = aimlessMergedMtzPath
        tempFile = tempfile.NamedTemporaryFile(suffix='.mtz',
                                                prefix='tmp2-',
                                                dir=self.aimlessTask.getWorkingDirectory(),
                                                delete=False)
        truncateOut = tempFile.name
        tempFile.close()
        truncateData["outputFile"] = truncateOut
        truncateData["nres"] = self.inData.get("nres")
        truncateData["isAnom"] = self.anomalous
        truncateData["res"] = self.resCutoff
        truncateData["outputFile"] = truncateOut

        os.chmod(truncateOut, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

        logger.info("Start ccp4/truncate...")
        self.truncate = TruncateTask(inData=truncateData)
        self.truncate.execute()

        uniqueifyData = {"inputFile":self.truncate.outData["truncateOutputMtz"],
                            "outputFile": "truncate_unique.mtz"}
        self.uniqueify = UniqueifyTask(inData=uniqueifyData)
        
        self.uniqueify.execute()
        self.phenixXTriageTask.join()

        truncateLog = self.resultsDirectory / f"{self.pyarchPrefix}_truncate.log"
        uniqueMtz = self.resultsDirectory / f"{self.pyarchPrefix}_truncate.mtz"
        phenixXTriageTaskLog = self.resultsDirectory / f"{self.pyarchPrefix}_phenix_xtriage_anom.mtz"


        try: 
            shutil.copy(Path(self.truncate.outData["truncateLogPath"]),truncateLog)
            shutil.copy(Path(self.uniqueify.outData["uniqueifyOutputMtz"]),uniqueMtz)
            shutil.copy(Path(self.phenixXTriageTask.outData["logPath"]), phenixXTriageTaskLog)
        except:
            logger.warning("Couldn't copy file to result directory")
            logger.warning(e)

        logger.info("Phenix.xtriage finished.")

        if self.phenixXTriageTask.outData["hasTwinning"] and self.phenixXTriageTask.outData["hasPseudotranslation"]:
            logger.warning("Pseudotranslation and twinning detected by phenix.xtriage!")
        elif self.phenixXTriageTask.outData["hasTwinning"]:
            logger.warning("Twinning detected by phenix.xtriage!")
        elif self.phenixXTriageTask.outData["hasPseudotranslation"]:
            logger.warning("Pseudotranslation detected by phenix.xtriage!")
        else:
            logger.info("No twinning or pseudotranslation detected by phenix.xtriage.")

        self.endDateTime = datetime.now().isoformat(timespec='seconds')

        if inData.get("test",False):
            self.tmpdir = tempfile.TemporaryDirectory() 
            self.pyarchDirectory = Path(self.tmpdir.name)
        else:
            self.pyarchDirectory = self.storeDataOnPyarch(self.resultFilePaths)

        # if inData.get("test",False):
        #     self.tmpdir = tempfile.TemporaryDirectory() 
        #     self.pyarchDirectory = Path(self.tmpdir.name)
        # else:
        #     reg = re.compile(r"(?:/gpfs/offline1/visitors/biomax/|/data/visitors/biomax/)")
        #     pyarchDirectory = re.sub(reg, "/data/staff/ispybstorage/visitors/biomax/", str(self.resultsDirectory))
        #     self.pyarchDirectory = Path(pyarchDirectory)
        #     try:
        #         self.pyarchDirectory.mkdir(exist_ok=True,parents=True, mode=0o755)
        #         logger.info(f"Created pyarch directory: {self.pyarchDirectory}")
        #         for file in self.resultsDirectory.iterdir():
        #             pyarchFile = UtilsPath.createPyarchFilePath(file)
        #             shutil.copy(file,pyarchFile)
        #     except OSError as e:
        #         logger.error(f"Error when creating pyarch_dir: {e}")
        #         self.tmpdir = tempfile.TemporaryDirectory() 
        #         self.pyarchDirectory = Path(self.tmpdir.name)



        # Let's get results into a container for ispyb
        self.autoProcResultsContainer = self.generateAutoProcScalingResultsContainer(
                programId=self.programId, integrationId=self.integrationId, isAnom=self.anomalous)

        with open(self.resultsDirectory / "ednaPROC.json","w") as fp:
            json.dump(self.autoProcResultsContainer,fp, indent=2, default=lambda o:str(o))

        #now send it to ISPyB
        logger.info("Sending data to ISPyB...")
        self.ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=self.autoProcResultsContainer, workingDirectorySuffix='final')
        self.ispybStoreAutoProcResults.execute()
        if self.ispybStoreAutoProcResults.isFailure():
            logger.error("ISPyB Store autoproc results failed.")
            # self.setFailure()
            # return

        comments = ""
        if self.anomalousFlag:
            comments += "Strong anomalous signal detected. "
        if self.phenixXTriageTask.outData["hasTwinning"] and self.phenixXTriageTask.outData["hasPseudotranslation"]:
            comments += "Pseudotranslation and twinning detected by phenix.xtriage! "
        elif self.phenixXTriageTask.outData["hasTwinning"]:
            comments += "Twinning detected by phenix.xtriage! "
        elif self.phenixXTriageTask.outData["hasPseudotranslation"]:
            comments += "Pseudotranslation detected by phenix.xtriage! "
        if comments:
            UtilsIspyb.updateDataCollectionGroupComments(self.dataCollectionId,comments)

        outData = self.autoProcResultsContainer

        self.timeEnd = time.perf_counter()
        logger.info(f"Time to process was {self.timeEnd-self.timeStart:0.4f} seconds")
        if self.tmpdir is not None:
            self.tmpdir.cleanup()

        return outData
        
    @classmethod
    def storeDataOnPyarch(cls,resultFilePaths, pyarchDirectory=None):
        #create paths on Pyarch
        if pyarchDirectory is None:
            pyarchDirectory = UtilsPath.createPyarchFilePath(resultFilePaths[0]).parent
            if not pyarchDirectory.exists():
                pyarchDirectory.mkdir(parents=True, exist_ok=True, mode=0o755)
        for resultFile in [f for f in resultFilePaths if f.exists()]:
            resultFilePyarchPath = UtilsPath.createPyarchFilePath(resultFile)
            try:
                logger.info(f"Copying {resultFile} to pyarch directory")
                shutil.copy(resultFile,resultFilePyarchPath)
            except Exception as e:
                logger.warning(f"Couldn't copy file {resultFile} to results directory {pyarchDirectory}")
                logger.warning(e)
        return pyarchDirectory
            

    def generateAutoProcScalingResultsContainer(self, programId, integrationId, isAnom):

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

        pointlessTask =  self.pointlessTaskRerun
        autoProcContainer = {
            "autoProcProgramId" : programId,
            "spaceGroup" : pointlessTask.outData["sgstr"],
            "refinedCellA" : pointlessTask.outData["cell"]["length_a"],
            "refinedCellB" : pointlessTask.outData["cell"]["length_b"],
            "refinedCellC" : pointlessTask.outData["cell"]["length_c"],
            "refinedCellAlpha" : pointlessTask.outData["cell"]["angle_alpha"],
            "refinedCellBeta" : pointlessTask.outData["cell"]["angle_beta"],
            "refinedCellGamma" : pointlessTask.outData["cell"]["angle_gamma"],
        }
        autoProcResultsContainer["autoProc"] = autoProcContainer

        autoProcAttachmentContainerList = []
        for file in self.pyarchDirectory.iterdir():
            attachmentContainer = {
                "file" : file,
            }
            autoProcAttachmentContainerList.append(attachmentContainer)
            
        autoProcResultsContainer["autoProcProgramAttachment"] = autoProcAttachmentContainerList
        xdsRerun = self.xdsRerun.outData


        autoProcIntegrationContainer = {
            "autoProcIntegrationId" : integrationId,
            "autoProcProgramId" : programId,
            "startImageNumber" : self.imgNumLow,
            "endImageNumber" : self.imgNumHigh,
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
        
        autoProcScalingHasIntContainer = {
            "autoProcIntegrationId" : integrationId,
        }
        autoProcResultsContainer["autoProcScalingHasInt"] = autoProcScalingHasIntContainer

        # add the scaling statistics to a list of containers... 
        autoProcScalingStatisticsContainerList = []
        aimlessResults = self.aimlessTask.outData.get("aimlessResults")
        
        for shell,result in aimlessResults.items():
            autoProcScalingStatisticsContainer = {}
            autoProcScalingStatisticsContainer["scalingStatisticsType"] = shell
            for k,v in result.items():
                autoProcScalingStatisticsContainer[k] = v
            if shell == "overall":
                autoProcScalingStatisticsContainer["isa"] = xdsRerun.get("ISa", 0.0)
            # need to make a few adjustments for ISPyB...
            autoProcScalingStatisticsContainer["rmerge"] *= 100
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
        image_path = Path(inData['imagePath'])
        m = re.search(r"\S+_\d{1,2}(?=_master.h5)",image_path.name)
        image_list_stem = m.group(0)

        image_list = []
        master_file = h5py.File(inData['imagePath'],'r')
        for data_file in master_file['/entry/data'].keys():
            image_nr_high = int(master_file['/entry/data'][data_file].attrs['image_nr_high'])
            image_nr_low = int(master_file['/entry/data'][data_file].attrs['image_nr_low'])
            for i in range(image_nr_low,image_nr_high+1):
                image_list.append(f"{str(image_path.parent)}/{image_list_stem}_{i:06}.h5")
        master_file.close()
        return len(image_list), {"imagePath": image_list}
    
    @staticmethod
    def generateImageListFromH5Master_fast(masterFilePath):
        """Given an h5 master file, generate an image list for SubWedgeAssembly."""
        masterFilePath = Path(masterFilePath)
        m = re.search(r"\S+_\d{1,2}(?=_master.h5)",masterFilePath.name)
        image_list_stem = m.group(0)

        image_list = []
        with h5py.File(masterFilePath,'r') as master_file:
            data_file_low = list(master_file['/entry/data'].keys())[0]
            data_file_high = list(master_file['/entry/data'].keys())[-1]        
            image_nr_high = int(master_file['/entry/data'][data_file_high].attrs['image_nr_high'])
            image_nr_low = int(master_file['/entry/data'][data_file_low].attrs['image_nr_low'])
            image_list.append(f"{str(masterFilePath.parent)}/{image_list_stem}_{image_nr_low:06}.h5")
        return image_nr_low, image_nr_high, {"imagePath": image_list}

    def getResCutoff(self,completeness_entries):
        """
        get resolution cutoff based on CORRECT.LP
        suggestion.
        """
        return min([x['res'] for x in completeness_entries if x['include_res_based_on_cc']], default=None)

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
    
    def createIntegrationId(self, comments, isAnom=False):
        """
        gets integrationID and programID, 
        sets processing status to RUNNING.
        """
        statusInput = {
            "dataCollectionId": self.dataCollectionId,
            "autoProcIntegration" : {
                "anomalous": isAnom,
            },
            "autoProcProgram": {
                "processingCommandLine": self.processingCommandLine,
                "processingPrograms": self.processingPrograms,
                "processingStatus": "RUNNING",
                "processingStartTime": self.startDateTime,
            },
            "autoProcStatus": {
                "step":  "Indexing",
                "status": "Launched",
                "comments": comments,
                "bltimeStamp": datetime.now().isoformat(timespec='seconds'),
            }
        }
        autoprocStatus = ISPyBStoreAutoProcStatus(inData=statusInput,workingDirectorySuffix="createIntegrationId")

        # get our EDNAproc status id
        autoprocStatus.execute()
        return (autoprocStatus.outData["autoProcIntegrationId"],
                autoprocStatus.outData["autoProcProgramId"])
    
    def if_anomalous_signal(self, aimless_log, threshold = 1.0):
        """Grab the anomalous CC RCR value and see if it is 
        sufficiently large to run fast_ep. Generally, a value 
        greater than 1 indicates a significant anomalous signal."""
        cc_rcr = 0.0
        try:
            with open(aimless_log,'r') as fp:
                for line in fp:
                    if "$TABLE:  Correlations CC(1/2) within dataset" in line:
                        while "Overall" not in line: 
                            line = next(fp)
                        cc_rcr = float(line.split()[3])
        except:
            pass
        if cc_rcr >= threshold:
            return True
        else:
            return False
    