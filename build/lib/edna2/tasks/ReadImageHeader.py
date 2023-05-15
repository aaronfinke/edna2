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
__date__ = "21/04/2019"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxv1/plugins/EDPluginGroupReadImageHeader-v1.0/plugins/
#      EDPluginControlReadImageHeaderv10.py

import os
import h5py
import numpy
import time
import pathlib

from edna2.utils import UtilsLogging

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsImage

logger = UtilsLogging.getLogger()

# Constants

# Default time out for wait file
DEFAULT_TIME_OUT = 30  # s

# Map between image suffix and image type
SUFFIX_ADSC = "img"
SUFFIX_MARCCD1 = "mccd"
SUFFIX_MARCCD2 = "marccd"
SUFFIX_Pilatus2M = "cbf"
SUFFIX_Pilatus6M = "cbf"
SUFFIX_Eiger4M = "cbf"
SUFFIX_Eiger9M = "cbf"
SUFFIX_Eiger16M = "cbf"

DICT_SUFFIX_TO_IMAGE = {
    SUFFIX_ADSC: "ADSC",
    SUFFIX_MARCCD1: "MARCCD",
    SUFFIX_MARCCD2: "MARCCD",
    SUFFIX_Pilatus2M: "Pilatus2M",
    SUFFIX_Pilatus6M: "Pilatus6M",
    SUFFIX_Eiger4M: "Eiger4M",
    SUFFIX_Eiger9M: "Eiger9M",
    SUFFIX_Eiger16M: "Eiger16M",
}


class ReadImageHeader(AbstractTask):

    def run(self, inData):
        listImagePath = inData["imagePath"]
        isFastMesh = inData.get("isFastMesh", False)
        hasOverlap = inData.get("hasOverlap", False)
        listSubWedge = []
        for imagePath in listImagePath:
            imageSuffix = os.path.splitext(imagePath)[1][1:]
            if imageSuffix == "cbf":
                subWedge = self.createCBFHeaderData(imagePath)
            elif imageSuffix == "h5":
                skipNumberOfImages = inData.get("skipNumberOfImages", False)
                subWedge = self.createHdf5HeaderData(
                    imagePath,
                    skipNumberOfImages,
                    isFastMesh=isFastMesh,
                    hasOverlap=hasOverlap,
                )
            else:
                raise RuntimeError(
                    "{0} cannot read image header from images with extension {1}".format(
                        self.__class__.__name__, imageSuffix
                    )
                )
            listSubWedge.append(subWedge)
        outData = {"subWedge": listSubWedge}
        return outData

    @classmethod
    def readCBFHeader(cls, filePath):
        """
        Returns an dictionary with the contents of a CBF image header.
        """
        dictHeader = None
        with open(filePath, "rb") as f:
            logger.info("Reading header from image " + filePath)
            f.seek(0, 0)
            doContinue = True
            iMax = 60
            index = 0
            while doContinue:
                line = f.readline().decode("utf-8")
                index += 1
                if "_array_data.header_contents" in line:
                    dictHeader = {}
                if "_array_data.data" in line or index > iMax:
                    doContinue = False
                if dictHeader is not None and line[0] == "#":
                    # Check for date
                    strTmp = line[2:].replace("\r\n", "")
                    if line[6] == "/" and line[10] == "/":
                        dictHeader["DateTime"] = strTmp
                    else:
                        strKey = strTmp.split(" ")[0]
                        strValue = strTmp.replace(strKey, "")[1:]
                        dictHeader[strKey] = strValue
        return dictHeader

    @classmethod
    def createCBFHeaderData(cls, imagePath):
        # Waiting for file
        timedOut, finalSize = UtilsPath.waitForFile(
            imagePath, expectedSize=100000, timeOut=DEFAULT_TIME_OUT
        )
        if timedOut:
            errorMessage = "Timeout when waiting for image %s" % imagePath
            logger.error(errorMessage)
            raise BaseException(errorMessage)
        dictHeader = cls.readCBFHeader(imagePath)
        detector = dictHeader["Detector:"]
        if (
            "PILATUS 3M" in detector
            or "PILATUS3 2M" in detector
            or "PILATUS 2M" in detector
            or "PILATUS2 3M" in detector
        ):
            detectorName = "PILATUS2 3M"
            detectorType = "pilatus2m"
            numberPixelX = 1475
            numberPixelY = 1679
        elif "PILATUS 6M" in detector or "PILATUS3 6M" in detector:
            detectorName = "PILATUS2 6M"
            detectorType = "pilatus6m"
            numberPixelX = 2463
            numberPixelY = 2527
        elif "eiger" in detector.lower() and "4m" in detector.lower():
            detectorName = "EIGER 4M"
            detectorType = "eiger4m"
            numberPixelX = 2070
            numberPixelY = 2167
        elif "eiger2" in detector.lower() and "9m" in detector.lower():
            detectorName = "EIGER2 9M"
            detectorType = "eiger9m"
            numberPixelX = 3108
            numberPixelY = 3262
        elif "eiger2" in detector.lower() and "16m" in detector.lower():
            detectorName = "EIGER2 16M"
            detectorType = "eiger16m"
            numberPixelX = 4148
            numberPixelY = 4362
        else:
            raise RuntimeError(
                "{0} cannot read image header from images with dector type {1}".format(
                    cls.__class__.__name__, detector
                )
            )
        experimentalCondition = {}
        detector = {"numberPixelX": numberPixelX, "numberPixelY": numberPixelY}
        # Pixel size
        listPixelSizeXY = dictHeader["Pixel_size"].split(" ")
        detector["pixelSizeX"] = float(listPixelSizeXY[0]) * 1000
        detector["pixelSizeY"] = float(listPixelSizeXY[3]) * 1000
        # Beam position
        listBeamPosition = (
            dictHeader["Beam_xy"]
            .replace("(", " ")
            .replace(")", " ")
            .replace(",", " ")
            .split()
        )
        detector["beamPositionX"] = float(listBeamPosition[0]) * detector["pixelSizeX"]
        detector["beamPositionY"] = float(listBeamPosition[1]) * detector["pixelSizeY"]
        distance = float(dictHeader["Detector_distance"].split(" ")[0]) * 1000
        detector["distance"] = distance
        detector["serialNumber"] = dictHeader["Detector:"]
        detector["name"] = detectorName
        detector["type"] = detectorType
        experimentalCondition["detector"] = detector
        # Beam object
        beam = {
            "wavelength": float(dictHeader["Wavelength"].split(" ")[0]),
            "exposureTime": float(dictHeader["Exposure_time"].split(" ")[0]),
        }
        experimentalCondition["beam"] = beam
        # Goniostat object
        goniostat = {}
        rotationAxisStart = float(dictHeader["Start_angle"].split(" ")[0])
        oscillationWidth = float(dictHeader["Angle_increment"].split(" ")[0])
        goniostat["rotationAxisStart"] = rotationAxisStart
        goniostat["rotationAxisEnd"] = rotationAxisStart + oscillationWidth
        goniostat["oscillationWidth"] = oscillationWidth
        experimentalCondition["goniostat"] = goniostat
        # Create the image object
        image = {"path": imagePath}
        if "DateTime" in dictHeader:
            image["date"] = dictHeader["DateTime"]
        imageNumber = UtilsImage.getImageNumber(imagePath)
        image["number"] = imageNumber
        subWedge = {"experimentalCondition": experimentalCondition, "image": [image]}
        return subWedge

    @classmethod
    def readHdf5Header(cls, filePath):
        """
        Returns an dictionary with the contents of an Eiger Hdf5 image header.
        """
        logger.info("Reading header from image " + str(filePath))
        f = h5py.File(filePath, "r")
        dictHeader = {
            "wavelength": f["entry"]["instrument"]["beam"]["incident_wavelength"][()],
            "beam_center_x": f["entry"]["instrument"]["detector"]["beam_center_x"][()],
            "beam_center_y": f["entry"]["instrument"]["detector"]["beam_center_y"][()],
            "count_time": f["entry"]["instrument"]["detector"]["count_time"][()],
            "detector_distance": f["entry"]["instrument"]["detector"][
                "detector_distance"
            ][()],
            "translation": list(
                f["entry"]["instrument"]["detector"]["geometry"]["translation"][
                    "distances"
                ]
            ),
            "x_pixel_size": f["entry"]["instrument"]["detector"]["x_pixel_size"][()],
            "y_pixel_size": f["entry"]["instrument"]["detector"]["y_pixel_size"][()],
            "omega_range_average": f["entry"]["sample"]["goniometer"][
                "omega_range_average"
            ][()],
            "detector_number": f["entry"]["instrument"]["detector"]["detector_number"][
                ()
            ].decode("utf-8"),
            "description": f["entry"]["instrument"]["detector"]["description"][
                ()
            ].decode("utf-8"),
            "data_collection_date": f["entry"]["instrument"]["detector"][
                "detectorSpecific"
            ]["data_collection_date"][()].decode("utf-8"),
            "data": list(f["entry"]["data"]),
        }
        # 'Old' Eiger files have just one entry for 'omega'
        omega = f["entry"]["sample"]["goniometer"]["omega"][()]
        if type(omega) == numpy.float32:
            dictHeader["omega_start"] = float(omega)
        else:
            dictHeader["omega_start"] = float(omega[0])
        f.close()
        return dictHeader

    @classmethod
    def createHdf5HeaderData(
        cls, imagePath, skipNumberOfImages=False, hasOverlap=False, isFastMesh=True
    ):
        h5MasterFilePath, h5DataFilePath, h5FileNumber = UtilsImage.getH5FilePath(
            pathlib.Path(imagePath), isFastMesh=isFastMesh, hasOverlap=hasOverlap
        )
        # Waiting for file
        timedOut, finalSize = UtilsPath.waitForFile(
            h5MasterFilePath, expectedSize=2000000, timeOut=DEFAULT_TIME_OUT
        )
        if timedOut:
            errorMessage = "Timeout when waiting for image %s" % imagePath
            logger.error(errorMessage)
            raise BaseException(errorMessage)
        logger.info("Final size for {0}: {1}".format(h5MasterFilePath, finalSize))
        noTrialsLeft = 5
        dictHeader = None
        while noTrialsLeft > 0:
            try:
                dictHeader = cls.readHdf5Header(h5MasterFilePath)
                noTrialsLeft = 0
            except Exception as e:
                logger.warning(
                    "Cannot read header from {0}, no trials left: {1}".format(
                        h5MasterFilePath, noTrialsLeft
                    )
                )
                time.sleep(5)
                noTrialsLeft -= 1
        if dictHeader is None:
            raise RuntimeError("Cannot read header from {0}!".format(h5MasterFilePath))
        description = dictHeader["description"]
        if "Eiger 4M" in description:
            detectorName = "EIGER 4M"
            detectorType = "eiger4m"
            numberPixelX = 2070
            numberPixelY = 2167
        elif "eiger2" in description.lower() and "9m" in description.lower():
            detectorName = "EIGER2 9M"
            detectorType = "eiger9m"
            numberPixelX = 3108
            numberPixelY = 3262
        elif "eiger" in description.lower() and "16M" in description:
            detectorName = "EIGER 16M"
            detectorType = "eiger16m"
            numberPixelX = 4148
            numberPixelY = 4362
        else:
            raise RuntimeError(
                "{0} cannot read image header from images with detector type {1}".format(
                    cls.__class__.__name__, description
                )
            )
        # Image number
        image_number = UtilsImage.getImageNumber(imagePath)
        # Find out size of data set
        prefix = str(h5MasterFilePath).split("master.h5")[0]
        listDataImage = []
        noImages = 0
        if not skipNumberOfImages:
            for data in dictHeader["data"]:
                dataFilePath = prefix + data + ".h5"
                timedOut, finalSize = UtilsPath.waitForFile(
                    dataFilePath, expectedSize=100000, timeOut=DEFAULT_TIME_OUT
                )
                if timedOut:
                    raise RuntimeError(
                        "Timeout waiting for file {0}".format(dataFilePath)
                    )
                # listDataImage.append({
                #     'path': dataFilePath
                # })
                f = h5py.File(dataFilePath, "r")
                dataShape = f["entry"]["data"]["data"].shape
                noImages += dataShape[0]
                f.close()
        experimentalCondition = {}
        # Pixel size and beam position
        detector = {
            "numberPixelX": int(numberPixelX),
            "numberPixelY": int(numberPixelY),
            "pixelSizeX": round(dictHeader["x_pixel_size"] * 1000, 3),
            "pixelSizeY": round(dictHeader["y_pixel_size"] * 1000, 3),
            "beamPositionX": round(
                float(dictHeader["beam_center_x"] * dictHeader["x_pixel_size"] * 1000),
                3,
            ),
            "beamPositionY": round(
                float(dictHeader["beam_center_y"] * dictHeader["y_pixel_size"] * 1000),
                3,
            ),
            "distance": round(float(dictHeader["detector_distance"]) * 1000, 3),
            "serialNumber": dictHeader["detector_number"],
            "name": detectorName,
            "type": detectorType,
        }
        experimentalCondition["detector"] = detector
        # Beam object
        beam = {
            "wavelength": round(float(dictHeader["wavelength"]), 6),
            "exposureTime": round(float(dictHeader["count_time"]), 6),
        }
        experimentalCondition["beam"] = beam
        # Goniostat object
        goniostat = {}
        rotationAxisStart = round(float(dictHeader["omega_start"]), 4)
        oscillationWidth = round(float(dictHeader["omega_range_average"]), 4)
        # Offset for the image number
        rotationAxisStart += (image_number - 1)*oscillationWidth
        goniostat["rotationAxisStart"] = rotationAxisStart
        goniostat["rotationAxisEnd"] = rotationAxisStart + oscillationWidth
        goniostat["oscillationWidth"] = oscillationWidth
        experimentalCondition["goniostat"] = goniostat
        # Create the image object
        image_dict = {
            "path": imagePath,
            "date": dictHeader["data_collection_date"],
            "number": 1,
        }
        # imageNumber = UtilsImage.getImageNumber(imagePath)
        # image['number'] = imageNumber
        subWedge = {
            "experimentalCondition": experimentalCondition,
            "image": [image_dict],
        }
        return subWedge
