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

from cctbx import sgtbx

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsDetector
from edna2.utils import UtilsSymmetry
from edna2.utils import UtilsIspyb


logger = UtilsLogging.getLogger()

from edna2.tasks.XDSTasks import XDSIndexing, XDSIntegration, XDSRerunCorrect
from edna2.tasks.SubWedgeAssembly import SubWedgeAssembly
from edna2.tasks.CCP4Tasks import PointlessTask, AimlessTask, TruncateTask, UniqueifyTask
# from edna2.tasks.XSCALETasks import XSCALETask
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
                "dataCollectionId": {"type":"integer"},
                "masterFilePath": {"type":"string"},
                "spaceGroup": {"type":["integer","string"]},
                "unitCell": {"type":"string"},
                "residues": {"type":"integer"},
                "doAnomandNoAnom": {"type":"boolean"},
                "workingDirectory": {"type":"string"},
            }
        }
    
    #sets ISPyB to FAILED if it's already logged
    def setFailure(self):
        self._dictInOut["isFailure"] = True
        if self.integrationIdAnom is not None and self.programIdAnom is not None:
            ISPyBStoreAutoProcResults.setIspybToFailed(
                dataCollectionId=self.dataCollectionId,
                autoProcProgramId=self.programIdAnom, 
                autoProcIntegrationId=self.integrationIdAnom, 
                processingCommandLine=self.processingCommandLine, 
                processingPrograms=self.processingPrograms, 
                isAnom=True, 
                timeStart=self.startDateTime, 
                timeEnd=datetime.now().isoformat(timespec='seconds')
            )
        if self.integrationIdNoAnom is not None and self.programIdNoAnom is not None:
            ISPyBStoreAutoProcResults.setIspybToFailed(
                dataCollectionId=self.dataCollectionId,
                autoProcProgramId=self.programIdNoAnom, 
                autoProcIntegrationId=self.integrationIdNoAnom, 
                processingCommandLine=self.processingCommandLine, 
                processingPrograms=self.processingPrograms, 
                isAnom=False, 
                timeStart=self.startDateTime, 
                timeEnd=datetime.now().isoformat(timespec='seconds')
            )



    def run(self, inData):
        self.tmpdir = None
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec='seconds')
        self.processingPrograms="EDNA2_proc"
        self.processingCommandLine = ""

        if inData.get("doAnomandNoAnom", True):
            self.doAnom = True
            self.doNoAnom = True
        else: 
            self.doAnom = True
            self.doNoAnom = False



        logger.info("EDNA2 Auto Processing started")
        logger.info(f"Running on {socket.gethostname()}")
        try:
            logger.info(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        spaceGroupNumber = None
        spaceGroup = inData.get("spaceGroup",0)
        unitCell = inData.get("unitCell",None)

        if spaceGroup != 0:
            try:
                spaceGroupInfo = sgtbx.space_group_info(spaceGroup).symbol_and_number()
                spaceGroupString = spaceGroupInfo.split("No. ")[0][:-2]
                spaceGroupNumber = int(spaceGroupInfo.split("No. ")[1][:-1])
                logger.info("Supplied space group is {}, number {}".format(spaceGroupString, spaceGroupNumber))
            except:
                logger.debug("Could not parse space group")
                spaceGroupNumber = 0
        else:
            logger.info("No space group supplied")

        # need both SG and unit cell
        if spaceGroup != 0 and unitCell is not None:
            try:
                unitCellList = [float(x) for x in unitCell.split(',')]
                #if there are zeroes parsed in, need to deal with it
                if 0.0 in unitCellList:
                    raise Exception
                unitCell = {"cell_a": unitCellList[0],
                            "cell_b": unitCellList[1],
                            "cell_c": unitCellList[2],
                            "cell_alpha": unitCellList[3],
                            "cell_beta": unitCellList[4],
                            "cell_gamma": unitCellList[5]}
                logger.info("Supplied unit cell is {cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**unitCell))
            except:
                logger.debug("could not parse unit cell")
                unitCell = None
        else:
            logger.info("No unit cell supplied")

        self.dataCollectionId = inData.get('dataCollectionId')

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

            # Determine pyarch prefix
            if UtilsConfig.isALBA():
                listPrefix = template.split("_")
                self.pyarchPrefix = "ep_{0}_{1}".format("_".join(listPrefix[:-2]),
                                                        listPrefix[-2])
            else:
                listPrefix = template.split("_")
                self.pyarchPrefix = "ep_{0}_{1}".format(listPrefix[-3], listPrefix[-2])


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
        
        if isH5 and inData.get("masterFilePath") is None:
            self.masterFilePath = os.path.join(directory,
                                self.eiger_template_to_master(template))
        elif isH5:
            self.masterFilePath = inData.get("masterFilePath")
        else:
            logger.error("Only supporing HDF5 data at this time. Stopping.")
            self.setFailure()
            return


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

        self.integrationIdAnom = None
        self.programIdAnom = None
        self.integrationIdNoAnom = None
        self.programIdNoAnom = None
        
        #get integrationIDs and programIDs, set them to running
        if self.doAnom:
            try: 
                self.integrationIdAnom,self.programIdAnom = self.createIntegrationId(
                                        "Creating anomalous integration ID", 
                                        isAnom=True)
            except Exception as e:
                logger.error("Could not get anom integration ID: \n{0}".format(traceback.format_exc(e)))
            logger.info(f"anom integrationID: {self.integrationIdAnom}, programIdAnom: {self.integrationIdAnom}")

        if self.doNoAnom:
            try: 
                self.integrationIdNoAnom,self.programIdNoAnom = self.createIntegrationId(
                                        "Creating non-anomalous integration ID", 
                                        isAnom=False)
                logger.info(f"noanom integrationID: {self.integrationIdNoAnom}, programIdNoAnom: {self.programIdNoAnom}")
            except Exception as e:
                logger.error("Could not get noanom integration ID: \n{0}".format(traceback.format_exc(e)))
        
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
        
        self.xdsIndexingInData["unitCell"] = unitCell
        self.xdsIndexingInData["spaceGroupNumber"] = spaceGroupNumber
        self.xdsIndexingInData["isAnom"] = self.doAnom

        self.indexing = XDSIndexing(inData=self.xdsIndexingInData, workingDirectorySuffix="init")

        if self.integrationIdAnom is not None:
            self.integrationId = self.integrationIdAnom
        elif self.integrationIdNoAnom is not None:
            self.integrationId = self.integrationIdNoAnom
        else:
            self.integrationId = None

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
        shutil.copy(Path(self.integration.outData["xdsInp"]), xds_INP_result_path)

        integrateLp_path = self.resultsDirectory / f"{self.pyarchPrefix}_INTEGRATE.LP"
        shutil.copy(Path(self.integration.outData["integrateLp"]), integrateLp_path)

        correctLp_path = self.resultsDirectory / f"{self.pyarchPrefix}_CORRECT.LP"
        shutil.copy(Path(self.integration.outData["correctLp"]), correctLp_path)

        integrateHkl_path = self.resultsDirectory / f"{self.pyarchPrefix}_INTEGRATE.HKL"
        shutil.copy(Path(self.integration.outData["integrateHkl"]), integrateHkl_path)

        xdsAsciiHkl_path = self.resultsDirectory / f"{self.pyarchPrefix}_XDS_ASCII.HKL"
        shutil.copy(Path(self.integration.outData["xdsAsciiHkl"]), xdsAsciiHkl_path)


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
        }

        xdsRerun_anomData = rerunCor_data
        xdsRerun_anomData['isAnom'] = True
        
        xdsRerun_noAnomData = rerunCor_data
        xdsRerun_noAnomData['isAnom'] = False

        logger.info("Rerunning CORRECT with the unit cell/SG from POINTLESS...")
        self.xdsRerun_anom = XDSRerunCorrect(inData=xdsRerun_anomData, workingDirectorySuffix='anom')
        self.xdsRerun_noAnom = XDSRerunCorrect(inData=xdsRerun_noAnomData, workingDirectorySuffix='noAnom')

        self.xdsRerun_anom.start()
        self.xdsRerun_noAnom.start()
        self.xdsRerun_anom.join()
        self.xdsRerun_noAnom.join()

        time3 = time.perf_counter()
        self.timeRerunCorrect = time3 - time2

        if self.xdsRerun_anom.isFailure() or self.xdsRerun_noAnom.isFailure():
            logger.error('Anom, noAnom rerun of CORRECT failed')
            self.setFailure()
            self.logToIspyb(self.integrationId,
                         'Scaling',
                         'Failed',
                         'Scaling failed after {0:.1}s'.format(self.timeRerunCorrect))
            return
        else:
            logger.info(f"Anom, noAnom rerun of CORRECT finished.")
            self.logToIspyb(self.integrationId,
                         'Scaling',
                         'Successful',
                         'Scaling finished in {0:.1f}s'.format(self.timeRerunCorrect))

        logger.info("Starting second resolution cutoff...")
        self.anomCompletenessEntries = self.xdsRerun_anom.outData['completenessEntries']
        self.noAnomCompletenessEntries = self.xdsRerun_noAnom.outData['completenessEntries']
        
        self.anomResCutoff = self.getResCutoff(self.anomCompletenessEntries)
        self.noAnomResCutoff = self.getResCutoff(self.noAnomCompletenessEntries)
        if self.anomResCutoff is None or self.noAnomResCutoff is None:
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
                            'Anomalous/nonAnom resolution cutoffs failed')
            self.setFailure()
            return
        logger.info(f"Resolution cutoff is {self.firstResCutoff}")
        self.logToIspyb(self.integrationId,
                        'Scaling',
                        'Successful',
                        'Anomalous resolution cutoffs finished')

        self.bins = [x["res"] for x in self.xdsRerun_anom.outData["completenessEntries"] if x["include_res_based_on_cc"] is True]
        # self.timeXscaleStart = time.perf_counter()
        # self.xscaleTaskData = {
        #     "xdsAsciiPath_anom" : self.xdsRerun_anom.outData["xdsAsciiHkl"],
        #     "xdsAsciiPath_noAnom": self.xdsRerun_noAnom.outData["xdsAsciiHkl"],
        #     "bins" : self.bins,
        #     "sgNumber": self.pointlessTask.outData["sgnumber"],
        #     "cell" : self.pointlessTask.outData["cell"],
        # }

        # self.xscaleTaskData_mergeAnom = self.xscaleTaskData
        # self.xscaleTaskData_mergeAnom['isAnom'] = True
        # self.xscaleTaskData_mergeAnom['merge'] = True
        # self.xscaleTaskData_mergeAnom['res'] = self.anomResCutoff

        # self.xscaleTask_mergeAnom = XSCALETask(inData=self.xscaleTaskData_mergeAnom, workingDirectorySuffix="mergeAnom")
        
        # self.xscaleTaskData_mergenoAnom = self.xscaleTaskData
        # self.xscaleTaskData_mergenoAnom['isAnom'] = False
        # self.xscaleTaskData_mergenoAnom['merge'] = True
        # self.xscaleTaskData_mergenoAnom['res'] = self.noAnomResCutoff

        # self.xscaleTask_mergenoAnom = XSCALETask(inData=self.xscaleTaskData_mergenoAnom, workingDirectorySuffix="mergenoAnom")

        # self.xscaleTaskData_unmergeAnom = self.xscaleTaskData
        # self.xscaleTaskData_unmergeAnom['isAnom'] = True
        # self.xscaleTaskData_unmergeAnom['merge'] = False
        # self.xscaleTaskData_unmergeAnom['res'] = self.anomResCutoff

        # self.xscaleTask_unmergeAnom = XSCALETask(inData=self.xscaleTaskData_unmergeAnom, workingDirectorySuffix="unmergeAnom")

        # self.xscaleTaskData_unmergenoAnom = self.xscaleTaskData
        # self.xscaleTaskData_unmergenoAnom['isAnom'] = False
        # self.xscaleTaskData_unmergenoAnom['merge'] = False
        # self.xscaleTaskData_unmergenoAnom['res'] = self.noAnomResCutoff

        # self.xscaleTask_unmergenoAnom = XSCALETask(inData=self.xscaleTaskData_unmergenoAnom, workingDirectorySuffix="unmergenoAnom")

        # logger.info("Starting XSCALE merging...")
        # self.logToIspyb(self.integrationId,
        #              'Scaling', 'Launched', 'Start of XSCALE')

        # self.xscaleTask_mergeAnom.start()
        # self.xscaleTask_mergenoAnom.start()
        # self.xscaleTask_unmergeAnom.start()
        # self.xscaleTask_unmergenoAnom.start()

        # self.xscaleTask_mergeAnom.join()
        # self.xscaleTask_mergenoAnom.join()
        # self.xscaleTask_unmergeAnom.join()
        # self.xscaleTask_unmergenoAnom.join()

        # time4 = time.perf_counter()
        # self.timeXscale = time4-time3

        # for task in [self.xscaleTask_mergeAnom,self.xscaleTask_mergenoAnom,
        #              self.xscaleTask_unmergeAnom,self.xscaleTask_unmergenoAnom]:
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

        # self.xscaleTask_mergeAnomLPFile = self.resultsDirectory / "ep__merged_anom_XSCALE.LP"
        # shutil.copy(self.xscaleTask_mergeAnom.outData["xscaleLp"], self.xscaleTask_mergeAnomLPFile)

        # self.xscaleTask_mergenoAnomLPFile = self.resultsDirectory / "ep__merged_noanom_XSCALE.LP"
        # shutil.copy(self.xscaleTask_mergenoAnom.outData["xscaleLp"], self.xscaleTask_mergenoAnomLPFile)

        # self.xscaleTask_unmergeAnomLPFile = self.resultsDirectory / "ep__unmerged_anom_XSCALE.LP"
        # shutil.copy(self.xscaleTask_unmergeAnom.outData["xscaleLp"], self.xscaleTask_unmergeAnomLPFile)

        # self.xscaleTask_unmergenoAnomLPFile = self.resultsDirectory / "ep__unmerged_noanom_XSCALE.LP"
        # shutil.copy(self.xscaleTask_unmergenoAnom.outData["xscaleLp"], self.xscaleTask_unmergenoAnomLPFile)

        #logger.debug(f"XSCALE output: {self.xscaleTask_mergeAnom.outData}")

        self.pointlessTaskAnominData = {
            'input_file' :  self.xdsRerun_anom.outData["xdsAsciiHkl"],
            'output_file' : "ep__anom_pointless_unmerged.mtz"
        }
        logger.info("Starting pointless anom task...")
        self.pointlessTaskAnom = PointlessTask(inData=self.pointlessTaskAnominData, workingDirectorySuffix="anom")

        self.pointlessTaskNoAnominData = {
            'input_file' :  self.xdsRerun_anom.outData["xdsAsciiHkl"],
            'output_file' : "ep__noanom_pointless_unmerged.mtz"
        }
        logger.info("Starting pointless tasks...")
        self.pointlessTaskNoAnom = PointlessTask(inData=self.pointlessTaskNoAnominData, workingDirectorySuffix="noAnom")

        self.pointlessTaskNoAnom.start()
        self.pointlessTaskAnom.start()

        self.pointlessTaskNoAnom.join()
        self.pointlessTaskAnom.join()
        
        try:
            pointlessAnomUnmergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_" + self.pointlessTaskAnominData["output_file"]
            shutil.copy(Path(self.pointlessTaskAnom.outData["pointlessUnmergedMtz"]), pointlessAnomUnmergedMtzPath)
            pointlessNoAnomUnmergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_" + self.pointlessTaskNoAnominData["output_file"]
            shutil.copy(Path(self.pointlessTaskNoAnom.outData["pointlessUnmergedMtz"]), pointlessNoAnomUnmergedMtzPath)
        except:
            logger.error("error copying pointless results to directory")

        self.aimlessTask_anomData = {
            "input_file" : self.pointlessTaskAnom.outData["pointlessUnmergedMtz"],
            "output_file" : "ep__anom_aimless.mtz",
            "start_image" : self.indexing.outData["start_image"],
            "end_image" : self.indexing.outData["end_image"],
            "dataCollectionId" : self.dataCollectionId,
            "res" : self.anomResCutoff,
            "anom" : True
        }
        self.aimlessTask_anom = AimlessTask(inData=self.aimlessTask_anomData, workingDirectorySuffix="anom")
        self.aimlessTask_noAnomData = {
            "input_file" : self.pointlessTaskNoAnom.outData["pointlessUnmergedMtz"],
            "output_file" : "ep__noanom_aimless.mtz",
            "start_image" : self.indexing.outData["start_image"],
            "end_image" : self.indexing.outData["end_image"],
            "dataCollectionId" : self.dataCollectionId,
            "res" : self.noAnomResCutoff,
            "anom" : False,
        }
        logger.info("Starting aimless...")
        self.aimlessTask_noanom = AimlessTask(inData=self.aimlessTask_noAnomData, workingDirectorySuffix="noanom")
        self.aimlessTask_anom.start()
        self.aimlessTask_noanom.start()

        self.aimlessTask_anom.join()
        self.aimlessTask_noanom.join()


        logger.info(f"Aimless finished.")
        try:
            aimlessAnomMergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_anom_aimless.mtz"
            shutil.copy(Path(self.aimlessTask_anom.outData["aimlessMergedMtz"]), aimlessAnomMergedMtzPath)
            aimlessAnomUnmergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_anom_aimless_unmerged.mtz.gz"
            shutil.copy(Path(self.aimlessTask_anom.outData["aimlessUnmergedMtz"]), aimlessAnomUnmergedMtzPath)
            aimlessAnomLogPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless_anom.log"
            shutil.copy(Path(self.aimlessTask_anom.outData["aimlessLog"]), aimlessAnomLogPath)

            aimlessNoanomMergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_noanom_aimless.mtz"
            shutil.copy(Path(self.aimlessTask_noanom.outData["aimlessMergedMtz"]), aimlessNoanomMergedMtzPath)
            aimlessNoanomUnmergedMtzPath = self.resultsDirectory / f"{self.pyarchPrefix}_noanom_aimless_unmerged.mtz.gz"
            shutil.copy(Path(self.aimlessTask_noanom.outData["aimlessUnmergedMtz"]), aimlessNoanomUnmergedMtzPath)
            aimlessNoAnomLogPath = self.resultsDirectory / f"{self.pyarchPrefix}_aimless_noanom.log"
            shutil.copy(Path(self.aimlessTask_noanom.outData["aimlessLog"]), aimlessNoAnomLogPath)

        except Exception as e:
            logger.error("Error copying aimless results to directory: {0}".format(e))
            self.setFailure()
            return 
        
        logger.info("Start phenix.xtriage run...")
        self.phenixXTriageTaskData = {
            "input_file" : aimlessAnomUnmergedMtzPath
        }
        self.phenixXTriageTask = PhenixXTriageTask(inData=self.phenixXTriageTaskData)
        self.phenixXTriageTask.start()

        #now run truncate/unique
        truncateAnomData = {}
        truncateAnomData["inputFile"] = aimlessAnomMergedMtzPath
        tempFileAnom = tempfile.NamedTemporaryFile(suffix='.mtz',
                                                prefix='tmp2-',
                                                dir=self.aimlessTask_anom.getWorkingDirectory(),
                                                delete=False)
        truncateAnomOut = tempFileAnom.name
        tempFileAnom.close()
        truncateAnomData["outputFile"] = truncateAnomOut
        truncateAnomData["nres"] = self.inData.get("nres")
        truncateAnomData["isAnom"] = True
        truncateAnomData["res"] = self.anomResCutoff
        truncateAnomData["outputFile"] = truncateAnomOut
        os.chmod(truncateAnomOut, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

        truncateNoAnomData = {}
        truncateNoAnomData["inputFile"] = aimlessAnomMergedMtzPath
        tempFileNoAnom = tempfile.NamedTemporaryFile(suffix='.mtz',
                                                prefix='tmp2-',
                                                dir=self.aimlessTask_noanom.getWorkingDirectory(),
                                                delete=False)
        truncateOutNoAnom = tempFileNoAnom.name
        tempFileAnom.close()
        truncateNoAnomData["outputFile"] = truncateOutNoAnom
        truncateNoAnomData["nres"] = self.inData.get("nres")
        truncateNoAnomData["isAnom"] = False
        truncateNoAnomData["res"] = self.noAnomResCutoff
        truncateNoAnomData["outputFile"] = truncateOutNoAnom

        os.chmod(truncateOutNoAnom, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

        logger.info("Start ccp4/truncate...")
        self.truncateAnom = TruncateTask(inData=truncateAnomData)
        self.truncateNoAnom = TruncateTask(inData=truncateNoAnomData)
        self.truncateAnom.start()
        self.truncateNoAnom.start()
        self.truncateAnom.join()
        self.truncateNoAnom.join()

        uniqueifyAnomData = {"inputFile":self.truncateAnom.outData["truncateOutputMtz"],
                            "outputFile": "truncate_unique.mtz"}
        uniqueifyNoAnomData = {"inputFile":self.truncateNoAnom.outData["truncateOutputMtz"],
                            "outputFile": "truncate_unique.mtz"}
        self.uniqueifyAnom = UniqueifyTask(inData=uniqueifyAnomData)
        self.uniqueifyNoAnom = UniqueifyTask(inData=uniqueifyNoAnomData)
        
        self.uniqueifyAnom.start()
        self.uniqueifyNoAnom.start()

        self.uniqueifyAnom.join()
        self.uniqueifyNoAnom.join()
        self.phenixXTriageTask.join()

        try: 
            truncateAnomLog = self.resultsDirectory / f"{self.pyarchPrefix}_truncate_anom.log"
            truncateNoAnomLog = self.resultsDirectory / f"{self.pyarchPrefix}_truncate_noanom.log"
            shutil.copy(Path(self.truncateAnom.outData["truncateLogPath"]),truncateAnomLog)
            shutil.copy(Path(self.truncateNoAnom.outData["truncateLogPath"]),truncateNoAnomLog)
        except Exception as e:
            logger.error("Error copying truncate log to results directory: {0}".format(e))

        try:
            uniqueAnomMtz = self.resultsDirectory / f"{self.pyarchPrefix}_anom_truncate.mtz"
            uniqueNoAnomMtz = self.resultsDirectory / f"{self.pyarchPrefix}_noanom_truncate.mtz"
            shutil.copy(Path(self.uniqueifyAnom.outData["uniqueifyOutputMtz"]),uniqueAnomMtz)
            shutil.copy(Path(self.uniqueifyNoAnom.outData["uniqueifyOutputMtz"]),uniqueNoAnomMtz)
        except Exception as e:
            logger.error("Error copying uniqueify mtz to results directory: {0}".format(e))

        try:
            phenixXTriageTaskLog = self.resultsDirectory / f"{self.pyarchPrefix}_phenix_xtriage_anom.mtz"
            shutil.copy(Path(self.phenixXTriageTask.outData["logPath"]), phenixXTriageTaskLog)
        except:
            logger.error("Error copying phenix.xtriage results to directory")
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
            reg = re.compile(r"(?:/gpfs/offline1/visitors/biomax/|/data/visitors/biomax/)")
            pyarchDirectory = re.sub(reg, "/data/staff/ispybstorage/visitors/biomax/", str(self.resultsDirectory))
            self.pyarchDirectory = Path(pyarchDirectory)
            try:
                self.pyarchDirectory.mkdir(exist_ok=True,parents=True, mode=0o755)
                logger.info(f"Created pyarch directory: {self.pyarchDirectory}")
                for file in self.resultsDirectory.iterdir():
                    pyarchFile = UtilsPath.createPyarchFilePath(file)
                    shutil.copy(file,pyarchFile)
            except OSError as e:
                logger.error(f"Error when creating pyarch_dir: {e}")
                self.tmpdir = tempfile.TemporaryDirectory() 
                self.pyarchDirectory = Path(self.tmpdir.name)


        # Let's get results into a container for ispyb
        self.autoProcResultsContainerAnom = self.generateAutoProcScalingResultsContainer(
                programId=self.programIdAnom, integrationId=self.integrationIdAnom, isAnom=True)
        self.autoProcResultsContainerNoAnom = self.generateAutoProcScalingResultsContainer(
                programId=self.programIdNoAnom, integrationId=self.integrationIdNoAnom, isAnom=False)

        with open(self.resultsDirectory / "ednaPROC.json","w") as fp:
            json.dump(self.autoProcResultsContainerNoAnom,fp, indent=2, default=lambda o:str(o))

        #now send it to ISPyB
        logger.info("Sending anom data to ISPyB...")
        self.ispybStoreAutoProcResultsAnom = ISPyBStoreAutoProcResults(inData=self.autoProcResultsContainerAnom, workingDirectorySuffix='anom')
        self.ispybStoreAutoProcResultsAnom.execute()
        if self.ispybStoreAutoProcResultsAnom.isFailure():
            logger.error("ISPyB Store anom autoproc results failed.")
            self.setFailure()
            return
        
        logger.info("Sending noanom data to ISPyB...")
        self.ispybStoreAutoProcResultsNoAnom = ISPyBStoreAutoProcResults(inData=self.autoProcResultsContainerNoAnom, workingDirectorySuffix='noanom')
        self.ispybStoreAutoProcResultsNoAnom.execute()
        if self.ispybStoreAutoProcResultsNoAnom.isFailure():
            logger.error("ISPyB Store noanom autoproc results failed.")
            self.setFailure()
            return



        self.timeEnd = time.perf_counter()
        logger.info(f"Time to process was {self.timeEnd-self.timeStart:0.4f} seconds")
        if inData.get("test",False):
            self.tmpdir.cleanup()

        return 
        

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

        pointlessTask = self.pointlessTaskAnom if isAnom else self.pointlessTaskNoAnom
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
        xdsRerun = self.xdsRerun_anom.outData if isAnom else self.xdsRerun_noAnom.outData


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
        aimlessResults = self.aimlessTask_anom.outData.get("aimlessResults") if isAnom else self.aimlessTask_noanom.outData.get("aimlessResults")
        
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
    
    def createIntegrationId(self, comments, isAnom=True):
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
