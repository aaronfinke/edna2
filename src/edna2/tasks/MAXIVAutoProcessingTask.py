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
import numpy as np
from pathlib import Path
import re
import socket
import time
from datetime import datetime
STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

# for the os.chmod
from cctbx import sgtbx

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsCCTBX


logger = UtilsLogging.getLogger()

from edna2.tasks.Edna2ProcTask import Edna2ProcTask
from edna2.tasks.FastdpTask import FastdpTask
from edna2.tasks.AutoPROCTask import AutoPROCTask
from edna2.tasks.Xia2DIALSTasks import Xia2DialsTask
from edna2.tasks.ControlPyDozor import ControlPyDozor

class MAXIVAutoProcessingTask(AbstractTask):
    """
    Runs four autoprocessing pipelines.
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
                "workingDirectory": {"type":["string","null"]},
                "pdb": {"type":["string","null"]},
                "test": {"type":["boolean","null"]},
                "waitForFiles": {"type":["boolean","null"]},
                "doUploadIspyb": {"type":["boolean","null"]}
            }
        }
    
    def run(self, inData):
        UtilsLogging.addLocalFileHandler(logger, self.getWorkingDirectory()/"MAXIVAutoProcessing.log")
        logger.info("AutoPROC processing started")
        if os.environ.get('SLURM_JOB_ID'):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")

        self.anomFlag = False
        outData = {}
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.tmpdir = None
        self.imageNoStart = inData.get("imageNoStart", None)
        self.imageNoEnd = inData.get("imageNoEnd", None)

        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.masterFilePath = inData.get("masterFilePath", None)
        self.anomalous = inData.get("anomalous",False)
        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.residues = inData.get("residues",None)
        self.workingDirectory = inData.get("workingDirectory", self.getWorkingDirectory())
        self.doUploadIspyb = inData.get("doUploadIspyb",True)
        self.waitForFiles = inData.get("waitForFiles",True)

        self.pdb = inData.get("pdb",None)
        self.test = inData.get("test",False)

        self.proteinAcronym = "AUTOMATIC"
        self.sampleName = "DEFAULT"

        try:
            logger.debug(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass

        #set up SG and unit cell
        self.spaceGroupNumber, self.spaceGroupString = UtilsCCTBX.parseSpaceGroup(self.spaceGroup)

        # set up unit cell
        if self.unitCell is not None:
            self.unitCell = UtilsCCTBX.parseUnitCell_str(self.unitCell)
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
            #if self.imageNoEnd - self.imageNoStart < -1:
            logger.error("There are fewer than 8 images, aborting")
            self.setFailure()
            return
        
        imgQualityDozor = ControlPyDozor(inData={
            "dataCollectionId" : self.dataCollectionId,
            "masterFile" : self.masterFilePath,
            "startNo" : self.imageNoStart,
            "batchSize": numImages,
            "directory": str(self.getWorkingDirectory()/"ControlPyDozor_0"),
            "doISPyBUpload": True,
            "doSubmit": True,
            "returnSpotList": False
        }, workingDirectorySuffix="0"
        )

        edna2ProcTask = Edna2ProcTask(inData={
            "onlineAutoProcessing": False,
            "dataCollectionId": self.dataCollectionId,
            "masterFilePath": self.masterFilePath,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd": self.imageNoEnd,
            "anomalous": self.anomalous,
            "waitForFiles": False,
            "doUploadIspyb": True,
            "test": self.test
            }, workingDirectorySuffix="0")

        fastDpTask = FastdpTask(inData={
            "onlineAutoProcessing": True,
            "dataCollectionId": self.dataCollectionId,
            "masterFilePath": self.masterFilePath,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath":self.masterFilePath,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd": self.imageNoEnd,
            "waitForFiles": False,
            "doUploadIspyb": True,
            "anomalous": False,
            "test": self.test
            }, workingDirectorySuffix="0")
        
        imgQualityDozor.start()
        edna2ProcTask.start()
        fastDpTask.start()

        imgQualityDozor.join()
        fastDpTask.join()
        edna2ProcTask.join()
        
        if edna2ProcTask.isSuccess() and fastDpTask.isSuccess():
            outData = {
                "edna2Proc":edna2ProcTask.outData,
                "fastDp":fastDpTask.outData,
            }
            if edna2ProcTask.outData.get("HighAnomSignal",False) or fastDpTask.outData.get("HighAnomSignal",False):
                self.anomalous = True
            else:
                self.anomalous = False
        elif edna2ProcTask.isSuccess():
            outData = {
                "edna2Proc":edna2ProcTask.outData,
            }
            self.anomalous = edna2ProcTask.outData.get("HighAnomSignal",False)
        elif fastDpTask.isSuccess():
            outData = {
                "fastDp":fastDpTask.outData,
            }
            self.anomalous = fastDpTask.outData.get("HighAnomSignal",False)
        else:
            self.anomalous = False
         
        autoPROCTask = AutoPROCTask(inData={
            "onlineAutoProcessing": True,
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath":self.masterFilePath,
            "anomalous": self.anomalous,
            "test":self.test,
            "doUploadIspyb":self.doUploadIspyb,
            "waitForFiles":self.waitForFiles
            }, workingDirectorySuffix="0")
        

        xia2DialsTask = Xia2DialsTask(inData={
            "onlineAutoProcessing": True,
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath":self.masterFilePath,
            "anomalous": self.anomalous,
            "test":self.test,
            "doUploadIspyb":self.doUploadIspyb,
            "waitForFiles":self.waitForFiles
            }, workingDirectorySuffix="0")

        autoPROCTask.start()
        xia2DialsTask.start()

        autoPROCTask.join()
        xia2DialsTask.join()

        if autoPROCTask.isSuccess():
            outData["autoPROCTask"] = autoPROCTask.outData
        
        if xia2DialsTask.isSuccess():
            outData["xia2DialsTask"] = xia2DialsTask.outData

        return outData

