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

__author__ = "Aaron Finke"
__contact__ = "aaron.finke@maxiv.lu.se"
__copyright__ = "MAX IV Laboratory"
__updated__ = "2023-11-03"

import os
import numpy
import time
import textwrap
from pathlib import Path
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt

from edna2.tasks.AbstractTask import AbstractTask
from edna2.utils import UtilsConfig
from edna2.utils import UtilsImage


from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb


logger = UtilsLogging.getLogger()

class Dozor2Task(AbstractTask):  # pylint: disable=too-many-instance-attributes
    """
    The Dozor2Task is responsible for executing the 'dozor2' program.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "nx" : {"type": "integer"},
                "ny" : {"type": "integer"},
                "pixel" : {"type": "number"},
                "exposure" : {"type": "number"},
                "detector_distance" : {"type": "number"},
                "X-ray_wavelength":  {"type": "number"},
                "fraction_polarization" : {"type": "number"},
                "orgx" :  {"type": "number"},
                "orgy" :  {"type": "number"},
                "oscillation_range" :  {"type": "number"},
                "starting_angle" :  {"type": "number"},
                "number_images":  {"type": "integer"},
                "first_image_number":  {"type": "integer"},
                "name_template_image":  {"type": "string"},
                "library":  {"type": "string"},
                "pixel_min" :  {"type": "number"},
                "pixel_max": {"type": "number"},
                "ix_max": {"type": "number"},
                "iy_max" : {"type": "number"},
                "ix_min": {"type": "number"},
                "iy_min": {"type": "number"},
                "spot_size":{"type": "integer"},
                "bad_zona":{"type": "string"},
                "wedge_number": {"type": "integer"},
                "spot_level": {"type": "integer"},
                "dataCollectionId":  {"type": "integer"},
                "masterFilePath" :  {"type": "string"}
            }
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "logPath": {"type": "string"},
                "workingDirectory": {"type": "string"},
            },
        }

    def run(self, inData):
        time_start = time.perf_counter()
        outData = {}
        self.masterFilePath = inData.get("masterFilePath")
        self.dataCollectionId = inData.get("dataCollectionId")

        dozorInData = self.getInDataSchema()["properties"]
        for k,v in dozorInData.items():
            dozorInData[k] = inData.get(k)
        del dozorInData["masterFilePath"]
        del dozorInData["dataCollectionId"]
        if self.dataCollectionId and not self.masterFilePath:
            self.masterFilePath = UtilsIspyb.getXDSMasterFilePath(self.dataCollectionId)
        if self.masterFilePath and Path(self.masterFilePath).exists():
            dozorMasterFileInData = self.parseDozorInDataFromMasterFile(self.masterFilePath)
            for k,v in dozorMasterFileInData.items():
                dozorInData[k] = v if dozorInData[k] is None else dozorInData[k]
            dozorInData["name_template_image"] = str(self.masterFilePath).replace('master','??????')
        dozorInData["library"] = UtilsConfig.get(self,"library")
        logger.debug(f'dozorInData: {dozorInData}')
        dozorInputFileData = self.generateDozorInputFileData(dozorInData)
        dozorDat = self.writeDozorInputFileData(dozorInputFileData,self.getWorkingDirectory())

        dozorSetup = UtilsConfig.get(self,"dozorSetup")
        cmd = ""
        cmd += f"source {dozorSetup}\n" if dozorSetup is not None else ""
        cmd += "dozor dozor.dat"

        self.runCommandLine(cmd)    
        time_end = time.perf_counter()
        logger.info(f"Process Time: {time_end - time_start} s")
        return outData
    
    @staticmethod
    def parseDozorInDataFromMasterFile(masterFile):
        imageHeader = UtilsImage.getImageMetadataFromH5MasterFile(masterFile)
        if not imageHeader:
            logger.error("Could not read image header")
            return
        
        dozorInData = {
                    "nx" : imageHeader['nx'],
                    "ny" :imageHeader['ny'],
                    "pixel" : imageHeader['x_pixel_size'],
                    "exposure" : imageHeader['count_time'],
                    "detector_distance" : imageHeader['detector_distance'],
                    "X-ray_wavelength": imageHeader['wavelength'],
                    "fraction_polarization" : UtilsConfig.get("Dozor2Task","polarization_fraction",1.0),
                    "orgx" : imageHeader['beam_center_x'],
                    "orgy" :  imageHeader['beam_center_y'],
                    "oscillation_range" : imageHeader['omega_range_average'],
                    "starting_angle" :  imageHeader['starting_angle'],
                    "number_images":  imageHeader['num_images'],
                    "first_image_number":  1,
                    "pixel_min" : 0,
                    "pixel_max": imageHeader['saturation_value'],
                    "ix_max": 0,
                    "iy_max" : 0,
                    "ix_min": 0,
                    "iy_min": 0
                }
        return dozorInData
    
    @staticmethod
    def generateDozorInputFileData(dozorInData):
        inputData = ""
        for k,v in dozorInData.items():
            inputData += f"{k} {v}\n" if v is not None else ""
        return inputData
    
    @staticmethod
    def writeDozorInputFileData(dozorInputFileData,dirf):
        dirf = Path(dirf)
        with open(dirf / "dozor.dat",'w') as fp:
            fp.write(dozorInputFileData)
        return dirf / "dozor.dat"
    

    @staticmethod
    def parseOutput(workingDir, logPath):
        outData = {}
        pathDozorm2Map = workingDir / "dozorm_001.map"
        if pathDozorm2Map.exists():
            dictMap = DozorM2.parseMap(pathDozorm2Map)
            nx = dictMap["nx"]
            ny = dictMap["ny"]
            dictCoord = DozorM2.parseDozorm2LogFile(logPath)
            # Fix for problem with 1D scans
            # listPositions = DozorM2.check1Dpositions(listPositions, nx, ny)
            crystalMapPath = DozorM2.makeCrystalPlot(dictMap["crystal"], workingDir)
            if nx != 1 and ny != 1:
                imageNumberMapPath = DozorM2.makeImageNumberMap(dictMap["imageNumber"], workingDir)
            else:
                imageNumberMapPath = None
            outData = {
                "dozorMap": str(pathDozorm2Map),
                "dictCoord": dictCoord,
                "nx": nx,
                "ny": ny,
                "score": dictMap["score"],
                "crystal": dictMap["crystal"],
                "imageNumber": dictMap["imageNumber"],
                "crystalMapPath": crystalMapPath,
                "imageNumberMapPath": imageNumberMapPath
            }
        return outData

    @staticmethod
    def parseDozorm2LogFile(logPath):
        #                    SCAN 1
        #                    ======
        #
        #    Total N.of crystals in Loop =  3
        # Cryst Aperture Central  Coordinate  Int/Sig  N.of Images CRsize Score   Dmin Helic   Start     Finish     Int/Sig
        # number size     image      X    Y          All  dX  dY   X   Y   sum                 x     y     x     y   helical
        # ------------------------------------------------------------------------------------------------------------------> X
        #     1   100.0     125   28.0    4.2  172.1  47  13   5  12   4  3846.2  3.06   NO
        #     2    20.0     133   20.0    4.0   47.5   5   3   2   3   2   147.8  3.46  YES    18     4    20     4   198.0
        #     3    20.0     198   31.0    6.0   37.5   2   2   1   2   1   112.2  3.68   NO
        with open(str(logPath)) as fd:
            listLogLines = fd.readlines()
        doParseLine = False
        do3dCoordinates = False
        scan1 = None
        scan2 = None
        listCoord = None
        for line in listLogLines:
            # print([line])
            if "SCAN 1" in line:
                listPositions = []
            elif "SCAN 2" in line:
                scan1 = listPositions
                listPositions = []
            elif "3D COORDINATES" in line:
                scan2 = listPositions
                do3dCoordinates = True
                listCoord = []
            if line.startswith("------"):
                doParseLine = True
            elif len(line) == 1:
                doParseLine = False
            elif doParseLine:
                listValues = line.split()
                if not do3dCoordinates:
                    try:
                        iOverSigma = float(listValues[5])
                    except:
                        iOverSigma = listValues[5]
                    position = {
                        "number": int(listValues[0]),
                        "apertureSize": str(int(float(listValues[1]))),
                        "imageNumber": int(listValues[2]),
                        "xPosition": float(listValues[3]),
                        "yPosition": float(listValues[4]),
                        "iOverSigma": iOverSigma,
                        "numberOfImagesTotal": int(listValues[6]),
                        "numberOfImagesTotalX": int(listValues[7]),
                        "numberOfImagesTotalY": int(listValues[8]),
                        "crSizeX": int(listValues[9]),
                        "crSizeY": int(listValues[10]),
                        "score": float(listValues[11]),
                        "dmin": float(listValues[12]),
                        "helical": listValues[13] == 'YES'
                    }
                    if position["helical"]:
                        position["helicalStartX"] = listValues[14]
                        position["helicalStartY"] = listValues[15]
                        position["helicalStopX"] = listValues[16]
                        position["helicalStopY"] = listValues[17]
                        position["helicalIoverSigma"] = listValues[18]
                    listPositions.append(position)
                else:
                    coord = {
                        "number": int(listValues[0]),
                        "averageScore": float(listValues[1]),
                        "dmin": float(listValues[2]),
                        "sc1": int(listValues[3]),
                        "sc2": int(listValues[4]),
                        "size": float(listValues[5]),
                        "scanX":  float(listValues[6]),
                        "scanY1": float(listValues[7]),
                        "scanY2": float(listValues[8]),
                        "dx": float(listValues[9]),
                        "dy1": float(listValues[10]),
                        "dy2": float(listValues[11]),
                        "sampx": float(listValues[12]),
                        "sampy": float(listValues[13]),
                        "phiy": float(listValues[14]),
                        # "alfa": float(listValues[15]),
                        # "sampx": float(listValues[16]),
                        # "sampy": float(listValues[17]),
                        # "phiy": float(listValues[18])
                    }
                    listCoord.append(coord)
        if scan1 is None:
            scan1 = listPositions
        dictCoord = {
            "scan1": scan1,
            "scan2": scan2,
            "coord": listCoord
        }
        return dictCoord

    @staticmethod
    def makeCrystalPlot(arrayCrystal, workingDirectory, debug=False):
        npArrayCrystal = numpy.array(arrayCrystal)
        ySize, xSize = npArrayCrystal.shape
        if xSize == 1:
            # Vertical line scan - transpose the matrix
            npArrayCrystal = numpy.transpose(npArrayCrystal)
            ySize, xSize = npArrayCrystal.shape
        npArrayCrystalAbs = numpy.abs(npArrayCrystal)
        # Make '999' be the max crystal number + 1
        maxNumber = numpy.amax(numpy.where(npArrayCrystalAbs < 999, npArrayCrystalAbs, 0))
        npArrayCrystalAbs = numpy.where(npArrayCrystalAbs == 999, maxNumber + 1, npArrayCrystalAbs)
        # minValue = numpy.amin(npArrayCrystal)
        # newZeroValue = minValue - 1
        # npArrayCrystal = numpy.where(npArrayCrystal == 0.0, newZeroValue, npArrayCrystal)

        maxSize = max(xSize, ySize)
        if maxSize < 10:
            fontSize = 12
            dpi = 75
        elif maxSize < 50:
            fontSize = 8
            dpi = 100
        else:
            fontSize =5
            dpi = 150

        font = {'family': 'normal',
                'weight': 'normal',
                'size': fontSize}

        matplotlib.rc('font', **font)

        fig, ax = plt.subplots()

        im = ax.imshow(
            npArrayCrystalAbs,
            cmap=matplotlib.cm.Spectral
        )

        ax.set_xticks(numpy.arange(len(range(xSize))))
        ax.set_yticks(numpy.arange(len(range(ySize))))

        ax.set_xticklabels(list(range(1, xSize+1)))
        ax.set_yticklabels(list(range(1, ySize+1)))

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(ySize):
            for j in range(xSize):
                if abs(npArrayCrystal[i, j]) > 0.001:
                    text = ax.text(j, i, npArrayCrystal[i, j],
                                   ha="center", va="center", color="b")

        ax.set_title("Crystal map")
        fig.tight_layout(pad=2)
        w, h = fig.get_size_inches()
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        x1 = float(x1)
        x2 = float(x2)
        y1 = float(y1)
        y2 = float(y2)
        fig.set_size_inches(w + 2, abs(y2 - y1) / (x2 - x1) * w + 2)
        if debug:
            plt.show()
        crystalMapPath = os.path.join(workingDirectory, "crystalMap.png")
        plt.savefig(crystalMapPath, dpi=dpi)

        return crystalMapPath

    @staticmethod
    def makeImageNumberMap(arrayImageNumber, workingDirectory, debug=False):
        npImageNumber = numpy.array(arrayImageNumber)
        npArrayImageNumber = numpy.zeros(npImageNumber.shape)
        ySize, xSize = npImageNumber.shape

        maxSize = max(xSize, ySize)
        if maxSize < 10:
            fontSize = 12
            dpi = 75
        elif maxSize < 50:
            fontSize = 8
            dpi = 100
        else:
            fontSize =5
            dpi = 150

        font = {'family': 'normal',
                'weight': 'normal',
                'size': fontSize}

        matplotlib.rc('font', **font)

        fig, ax = plt.subplots()
        im = ax.imshow(
            npArrayImageNumber,
            cmap=matplotlib.cm.Greys
        )

        ax.set_xticks(numpy.arange(len(range(xSize))))
        ax.set_yticks(numpy.arange(len(range(ySize))))

        ax.set_xticklabels(list(range(1, xSize+1)))
        ax.set_yticklabels(list(range(1, ySize+1)))

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        for i in range(ySize):
            for j in range(xSize):
                text = ax.text(j, i, arrayImageNumber[i][j],
                               ha="center", va="center", color="b")

        ax.set_title("Image numbers")
        fig.tight_layout(pad=2)
        w, h = fig.get_size_inches()
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        x1 = float(x1)
        x2 = float(x2)
        y1 = float(y1)
        y2 = float(y2)
        fig.set_size_inches(w + 2, abs(y2 - y1) / (x2 - x1) * w + 2)
        if debug:
            plt.show()
        imageNumberPath = os.path.join(workingDirectory, "imageNumber.png")
        plt.savefig(imageNumberPath, dpi=dpi)

        return imageNumberPath

    @staticmethod
    def parseMatrix(index, listLines, spacing, isFloat=True):
        arrayValues = []
        # Parse matrix - starts and ends with "---" line
        while not listLines[index].startswith("-------"):
            index += 1
        index += 1
        while not listLines[index].startswith("-------"):
            # listScores = textwrap.wrap(listLines[index][5:], spacing)
            listScores = []
            for line_pos in range(0, len(listLines[index]), spacing):
                sub_string = listLines[index][line_pos + 5:line_pos + spacing + 5].strip()
                if sub_string != "":
                    if isFloat:
                        listScores.append(float(sub_string))
                    else:
                        listScores.append(int(sub_string))
            arrayValues.append(listScores)
            index += 1
        index += 1
        return index, arrayValues



    @staticmethod
    def parseMap(mapPath):
        with open(str(mapPath)) as fd:
            listLines = fd.readlines()
        # Parse map dimensions
        index = 1
        nx, ny = map(int, listLines[index].split())
        # Parse scores
        index, arrayScore = DozorM2.parseMatrix(index, listLines, spacing=6, isFloat=True)
        # Parse rel. contamination
        index, relContamination = DozorM2.parseMatrix(index, listLines, spacing=6, isFloat=True)
        # Parse crystals
        index, arrayCrystal = DozorM2.parseMatrix(index, listLines, spacing=4, isFloat=False)
        # Parse image number
        index, arrayImageNumber = DozorM2.parseMatrix(index, listLines, spacing=5, isFloat=False)
        dictMap = {
            "nx": nx,
            "ny": ny,
            "score": arrayScore,
            "relContamination": relContamination,
            "crystal": arrayCrystal,
            "imageNumber": arrayImageNumber
        }
        return dictMap

    @staticmethod
    def updateMeshPositions(meshPositions, arrayScore):
        newMeshPositions = []
        for position in meshPositions:
            # pprint.pprint(position)
            indexY = position["indexY"]
            indexZ = position["indexZ"]
            # print(indexY, indexZ)
            dozormScore = arrayScore[indexZ][indexY]
            dozorScore = position["dozor_score"]
            # print(dozorScore, dozormScore)
            newPosition = dict(position)
            newPosition["dozor_score_orig"] = dozorScore
            newPosition["dozor_score"] = dozormScore
            newMeshPositions.append(newPosition)
        return newMeshPositions

    @staticmethod
    def check1Dpositions(listPositions, nx, ny):
        newListPositions = []
        for position in listPositions:
            newPosition = dict(position)
            if nx == 1:
                newPosition["xPosition"] = 1.0
            if ny == 1:
                newPosition["yPosition"] = 1.0
            newListPositions.append(newPosition)
        return newListPositions
