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

__authors__ = ["O. Svensson"]
__license__ = "MIT"
__date__ = "10/05/2019"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxv1/plugins/EDPluginControlImageQualityIndicators-v1.4/plugins/
#      EDPluginControlImageQualityIndicatorsv1_4.py

import os
import re
import time
import pathlib

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.WaitFileTask import WaitFileTask
from edna2.tasks.ControlDozor import ControlDozor
from edna2.tasks.PhenixTasks import DistlSignalStrengthTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()

DEFAULT_MIN_IMAGE_SIZE = 1000000
DEFAULT_WAIT_FILE_TIMEOUT = 300


class ImageQualityIndicators(AbstractTask):
    """
    This task controls the plugins that generate image quality indicators.
    """

    def __init__(self, inData, workingDirectorySuffix=None):
        AbstractTask.__init__(self, inData, workingDirectorySuffix)
        self.beamline = None
        self.directory = None
        self.template = None
        self.doSubmit = None
        self.doDozorM = None
        self.doDistlSignalStrength = None
        self.isFastMesh = None
        self.listImage = None
        self.batchSize = None
        self.minImageSize = None
        self.waitFileTimeout = None
        self.doIspybUpload = None
        self.dataCollectionId = None

    def run(self, inData):
        listImageQualityIndicators = []
        listControlDozorAllFile = []
        # Initialize parameters
        self.init(inData)
        # Set up batch list
        listOfBatches = self.createBatchList(inData)
        outData = dict()
        listDozorTask, listDistlTask = self.runDozorAndDistl(listOfBatches)
        if not self.isFailure():
            listDisltResult = self.synchronizeDislt(listDistlTask)
            (
                listImageQualityIndicators,
                listControlDozorAllFile,
            ) = self.synchronizeDozor(listDozorTask, listDisltResult)
            # Assemble all controlDozorAllFiles into one
            if self.doDozorM:
                imageQualityIndicatorsDozorAllFile = self.createDozorAllFile(
                    listControlDozorAllFile
                )
                outData["dozorAllFile"] = imageQualityIndicatorsDozorAllFile
        outData["imageQualityIndicators"] = listImageQualityIndicators
        return outData

    def init(self, inData):
        self.beamline = inData.get("beamline", None)
        self.doSubmit = inData.get("doSubmit", False)
        self.doDozorM = inData.get("doDozorM", False)
        self.doDistlSignalStrength = inData.get("doDistlSignalStrength", False)
        self.isFastMesh = inData.get("fastMesh", False)
        self.listImage = inData.get("image", [])
        self.batchSize = inData.get("batchSize", 1)
        self.doIspybUpload = inData.get("doIspybUpload", False)
        self.dataCollectionId = inData.get("dataCollectionId", None)
        # Configurations
        self.minImageSize = UtilsConfig.get(
            self, "minImageSize", defaultValue=DEFAULT_MIN_IMAGE_SIZE
        )
        self.waitFileTimeOut = UtilsConfig.get(
            self, "waitFileTimeOut", defaultValue=DEFAULT_WAIT_FILE_TIMEOUT
        )

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "beamline": {"type": "string"},
                "doDozorM": {"type": "boolean"},
                "doDistlSignalStrength": {"type": "boolean"},
                "doIndexing": {"type": "boolean"},
                "doIspybUpload": {"type": "boolean"},
                "processDirectory": {"type": "string"},
                "image": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "batchSize": {"type": "integer"},
                "fastMesh": {"type": "boolean"},
                "wedgeNumber": {"type": "integer"},
                "directory": {"type": "string"},
                "template": {"type": "string"},
                "startNo": {"type": "integer"},
                "endNo": {"type": "integer"},
                "dataCollectionId": {
                    "anyOf": [
                      {"type": "integer"},
                      {"type": "null"}
                    ]
                }
            },
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "required": ["imageQualityIndicators"],
            "properties": {
                "imageQualityIndicators": {
                    "type": "array",
                    "items": {"$ref": self.getSchemaUrl("imageQualityIndicators.json")},
                },
                "inputDozor": {"type": "number"},
                "dozorAllFile": {"type": "string"},
            },
        }

    def createBatchList(self, inData):
        listOfBatches = []
        listOfImagesInBatch = []
        if len(self.listImage) == 0:
            self.directory = pathlib.Path(inData["directory"])
            self.template = inData["template"]
            startNo = inData["startNo"]
            endNo = inData["endNo"]
            for index in range(startNo, endNo + 1):
                listOfImagesInBatch.append(index)
                if len(listOfImagesInBatch) == self.batchSize or index == 9999:
                    listOfBatches.append(listOfImagesInBatch)
                    listOfImagesInBatch = []
        else:
            firstImage = pathlib.Path(self.listImage[0])
            self.directory = firstImage.parent
            self.template = UtilsImage.getTemplate(firstImage)
            for image in self.listImage:
                imageNo = UtilsImage.getImageNumber(image)
                listOfImagesInBatch.append(imageNo)
                if len(listOfImagesInBatch) == self.batchSize:
                    listOfBatches.append(listOfImagesInBatch)
                    listOfImagesInBatch = []
        if len(listOfImagesInBatch) > 0:
            listOfBatches.append(listOfImagesInBatch)
        return listOfBatches

    def runDozorAndDistl(self, listOfBatches):
        #
        # Loop over batches:
        # - Wait for all files in batch
        # - Run Dozor and DistlSignalStrength (if required) in parallel
        #
        distl_tasks = []
        dozor_tasks = []
        template4d = re.sub("#+", "{0:04d}", self.template)
        for index, images_in_batch in enumerate(listOfBatches):
            listOfH5FilesInBatch = []
            image_no = images_in_batch[-1]
            # Wait for last image
            image_path = self.directory / template4d.format(image_no)
            logger.debug("Waiting for path: {0}".format(image_path))
            self.waitForImagePath(
                imagePath=image_path,
                batchSize=self.batchSize,
                isFastMesh=self.isFastMesh,
                minImageSize=self.minImageSize,
                waitFileTimeOut=self.waitFileTimeOut,
                listofH5FilesInBatch=listOfH5FilesInBatch,
            )
            logger.debug("Done waiting for path: {0}".format(image_path))
            if not self.isFailure():
                # Determine start and end image no
                batchStartNo = images_in_batch[0]
                batchEndNo = images_in_batch[-1]
                dozorTemplate = self.template
                # Run Control Dozor
                inDataControlDozor = {
                    "template": dozorTemplate,
                    "directory": self.directory,
                    "startNo": batchStartNo,
                    "endNo": batchEndNo,
                    "batchSize": self.batchSize,
                    "doSubmit": self.doSubmit,
                    "doDozorM": self.doDozorM,
                    "doIspybUpload": self.doIspybUpload,
                    "dataCollectionId": self.dataCollectionId
                }
                if self.beamline is not None:
                    inDataControlDozor["beamline"] = self.beamline
                controlDozor = ControlDozor(
                    inDataControlDozor,
                    workingDirectorySuffix="{0:04d}_{1:04d}".format(
                        batchStartNo, batchEndNo
                    ),
                )
                controlDozor.start()
                dozor_tasks.append(
                    (controlDozor, inDataControlDozor, list(images_in_batch))
                )
                # Check if we should run distl.signalStrength
                if self.doDistlSignalStrength:
                    for image_no in images_in_batch:
                        image_path = self.directory / template4d.format(image_no)
                        inDataDistl = {"referenceImage": str(image_path)}
                        distlTask = DistlSignalStrengthTask(
                            inData=inDataDistl,
                            workingDirectorySuffix=image_no,
                        )
                        distlTask.start()
                        distl_tasks.append((image_path, distlTask))
        return dozor_tasks, distl_tasks

    def synchronizeDislt(self, listDistlTask):
        listDistlResult = []
        # Synchronize all image quality indicator plugins and upload to ISPyB
        for (image, distlTask) in listDistlTask:
            imageQualityIndicators = {}
            if distlTask is not None:
                distlTask.join()
                if distlTask.isSuccess():
                    outDataDistl = distlTask.outData
                    if outDataDistl is not None:
                        imageQualityIndicators = outDataDistl["imageQualityIndicators"]
            imageQualityIndicators["image"] = str(image)
            listDistlResult.append(imageQualityIndicators)
        return listDistlResult

    def synchronizeDozor(self, listDozorTask, listDistlResult):
        listImageQualityIndicators = []
        listControlDozorAllFile = []
        for (
            controlDozor,
            inDataControlDozor,
            listBatch,
        ) in listDozorTask:
            controlDozor.join()
            # Check that we got at least one result
            if len(controlDozor.outData["imageQualityIndicators"]) == 0:
                # Run the dozor plugin again, this time synchronously
                firstImage = listBatch[0]
                lastImage = listBatch[-1]
                logger.warning(
                    "No dozor results! Re-executing Dozor for"
                    + " images {0} to {1}".format(firstImage, lastImage)
                )
                controlDozor = ControlDozor(
                    inDataControlDozor,
                    workingDirectorySuffix = "{0:04d}_{1:04d}_redo".format(
                        firstImage, lastImage
                    ),
                )
                controlDozor.execute()
            listOutDataControlDozor = list(
                controlDozor.outData["imageQualityIndicators"]
            )
            if self.doDistlSignalStrength:
                for outDataControlDozor in listOutDataControlDozor:
                    for distlResult in listDistlResult:
                        if outDataControlDozor["image"] == distlResult["image"]:
                            imageQualityIndicators = dict(outDataControlDozor)
                            imageQualityIndicators.update(distlResult)
                            listImageQualityIndicators.append(imageQualityIndicators)
            else:
                listImageQualityIndicators += listOutDataControlDozor
            # Check if dozorm
            if self.doDozorM:
                listControlDozorAllFile.append(controlDozor.outData["dozorAllFile"])
        return listImageQualityIndicators, listControlDozorAllFile

    def createDozorAllFile(self, listControlDozorAllFile):
        imageQualityIndicatorsDozorAllFile = str(
            self.getWorkingDirectory() / "dozor_all"
        )
        os.system("touch {0}".format(imageQualityIndicatorsDozorAllFile))
        for controlDozorAllFile in listControlDozorAllFile:
            command = (
                "cat "
                + controlDozorAllFile
                + " >> "
                + imageQualityIndicatorsDozorAllFile
            )
            os.system(command)
        return imageQualityIndicatorsDozorAllFile

    @classmethod
    def getH5FilePath(cls, filePath, batchSize=1, isFastMesh=False):
        imageNumber = UtilsImage.getImageNumber(filePath)
        prefix = UtilsImage.getPrefix(filePath)
        if isFastMesh:
            h5ImageNumber = int((imageNumber - 1) / 100) + 1
            h5FileNumber = 1
        else:
            h5ImageNumber = 1
            h5FileNumber = int((imageNumber - 1) / batchSize) * batchSize + 1
        h5MasterFileName = "{prefix}_{h5FileNumber}_master.h5".format(
            prefix=prefix, h5FileNumber=h5FileNumber
        )
        h5MasterFilePath = filePath.parent / h5MasterFileName
        h5DataFileName = "{prefix}_{h5FileNumber}_data_{h5ImageNumber:06d}.h5".format(
            prefix=prefix, h5FileNumber=h5FileNumber, h5ImageNumber=h5ImageNumber
        )
        h5DataFilePath = filePath.parent / h5DataFileName
        return h5MasterFilePath, h5DataFilePath, h5FileNumber

    def waitForImagePath(
        self,
        imagePath,
        batchSize,
        isFastMesh,
        minImageSize,
        waitFileTimeOut,
        listofH5FilesInBatch,
    ):
        # Force an 'ls' in parent directory - this sometimes helps to 'unblock'
        # the file system
        os.system("ls {0} > /dev/null".format(os.path.dirname(imagePath)))
        # If Eiger, just wait for the h5 file
        if imagePath.suffix == ".h5":
            h5MasterFilePath, h5DataFilePath, hdf5ImageNumber = self.getH5FilePath(
                imagePath, batchSize=batchSize, isFastMesh=isFastMesh
            )
            if h5DataFilePath not in listofH5FilesInBatch:
                listofH5FilesInBatch.append(h5DataFilePath)
                logger.info("Eiger data, waiting for master" + " and data files...")
                inDataWaitFileTask = {
                    "file": str(h5DataFilePath),
                    "size": minImageSize,
                    "timeOut": waitFileTimeOut,
                }
                workingDirectorySuffix = h5DataFilePath.name.split(".h5")[0]
                waitFileTask = WaitFileTask(
                    inData=inDataWaitFileTask,
                    workingDirectorySuffix=workingDirectorySuffix,
                )
                logger.info("Waiting for file {0}".format(h5DataFilePath))
                logger.debug("Wait file timeOut set to %f" % waitFileTimeOut)
                waitFileTask.execute()
                time.sleep(0.1)
            if not os.path.exists(h5DataFilePath):
                errorMessage = "Time-out while waiting for image %s" % h5DataFilePath
                logger.error(errorMessage)
                self.setFailure()
        else:
            if not imagePath.exists():
                logger.info("Waiting for file {0}".format(imagePath))
                inDataWaitFileTask = {
                    "file": str(imagePath),
                    "size": minImageSize,
                    "timeOut": waitFileTimeOut,
                }
                workingDirectorySuffix = imagePath.name.split(imagePath.suffix)[0]
                waitFileTask = WaitFileTask(
                    inData=inDataWaitFileTask,
                    workingDirectorySuffix=workingDirectorySuffix,
                )
                logger.debug("Wait file timeOut set to %.0f s" % waitFileTimeOut)
                waitFileTask.execute()
            if not imagePath.exists():
                errorMessage = "Time-out while waiting for image " + str(imagePath)
                logger.error(errorMessage)
                self.setFailure()
