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
import json
STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

# for the os.chmod
from stat import *

from cctbx import sgtbx

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb


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
                "test": {"type":["boolean","null"]}
            }
        }
    
    def run(self, inData):
        self.anomFlag = False

        logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        outData = {}
        self.timeStart = time.perf_counter()
        self.startDateTime =  datetime.now().isoformat(timespec="seconds")
        self.startDateTimeFormatted = datetime.now().strftime("%y%m%d-%H%M%S")
        self.tmpdir = None

        self.dataCollectionId = inData.get("dataCollectionId", None)
        self.masterFilePath = inData.get("masterFilePath", None)
        self.anomalous = inData.get("anomalous",False)
        self.spaceGroup = inData.get("spaceGroup",0)
        self.unitCell = inData.get("unitCell",None)
        self.residues = inData.get("residues",None)
        self.workingDirectory = inData.get("workingDirectory", self.getWorkingDirectory())
        self.pdb = inData.get("pdb",None)
        self.test = inData.get("test",False)

        self.proteinAcronym = "AUTOMATIC"
        self.sampleName = "DEFAULT"

        logger.info("MAX IV Autoprocessing started")
        logger.info(f"Running on {socket.gethostname()}")
        try:
            logger.info(f"System load avg: {os.getloadavg()}")
        except OSError:
            pass


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
            self.spaceGroupNumber = 0
            self.spaceGroupString = ""            
            logger.info("No space group supplied")

        # need both SG and unit cell
        if self.spaceGroup != 0 and self.unitCell is not None:
            try:
                unitCellList = [float(x) for x in self.unitCell.split(",")]
                #if there are zeroes parsed in, need to deal with it
                if 0.0 in unitCellList:
                    raise Exception
                self.unitCell = {
                            "cell_a": unitCellList[0],
                            "cell_b": unitCellList[1],
                            "cell_c": unitCellList[2],
                            "cell_alpha": unitCellList[3],
                            "cell_beta": unitCellList[4],
                            "cell_gamma": unitCellList[5]
                            }
                logger.info("Supplied unit cell is {cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**self.unitCell))
            except:
                logger.debug("could not parse unit cell")
                self.unitCell = None
        else:
            logger.info("No unit cell supplied")

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
            if self.masterFilePath is None:
                self.masterFilePath = UtilsIspyb.getXDSMasterFilePath(self.dataCollectionId)

        except Exception as e:
            logger.warning("Retrieval of data from ISPyB Failed, trying manually...")
            logger.warning(e)
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
            "onlineAutoProcessing": True,
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd": self.imageNoEnd,
            "anomalous": False,
            "test": self.test
            }, workingDirectorySuffix="0")

        fastDpTask = FastdpTask(inData={
            "onlineAutoProcessing": True,
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath":self.masterFilePath,
            "imageNoStart": self.imageNoStart,
            "imageNoEnd": self.imageNoEnd,
            "anomalous": False,
            "test": self.test
            }, workingDirectorySuffix="0")
        
        imgQualityDozor.start()
        edna2ProcTask.start()
        fastDpTask.start()

        imgQualityDozor.join()
        edna2ProcTask.join()
        fastDpTask.join()
        return
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
            "test":self.test
            }, workingDirectorySuffix="0")
        

        xia2DialsTask = Xia2DialsTask(inData={
            "onlineAutoProcessing": True,
            "dataCollectionId": self.dataCollectionId,
            "unitCell": self.unitCell,
            "spaceGroup": self.spaceGroup,
            "masterFilePath":self.masterFilePath,
            "anomalous": self.anomalous,
            "test":self.test
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

