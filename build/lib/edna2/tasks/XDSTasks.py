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
__date__ = "20/04/2020"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx

# mxPluginExec/plugins/EDPluginGroupXDS-v1.0/plugins/EDPluginXDSv1_0.py
# mxPluginExec/plugins/EDPluginGroupXDS-v1.0/plugins/EDPluginXDSIndexingv1_0.py

import os
import math
import shutil
import numpy as np

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsDetector
from edna2.utils import UtilsSymmetry


logger = UtilsLogging.getLogger()


R2D = 180 / math.pi


class XDSTask(AbstractTask):
    """
    Common base class for all XDS tasks
    """

    def run(self, inData):
        commandLine = "xds_par"
        listXDS_INP = self.generateXDS_INP(inData)
        self.writeXDS_INP(listXDS_INP, self.getWorkingDirectory())
        self.setLogFileName("xds.log")
        self.runCommandLine(commandLine, listCommand=[])
        # Work in progress!
        outData = self.parseXDSOutput(self.getWorkingDirectory())
        return outData

    @staticmethod
    def generateImageLinks(in_data, working_directory=None):
        first_sub_wedge = in_data["subWedge"][0]
        first_image_path = first_sub_wedge["image"][0]["path"]
        prefix = UtilsImage.getPrefix(first_image_path)
        suffix = UtilsImage.getSuffix(first_image_path)
        if suffix == "h5":
            lowest_xds_image_number = 1
            highest_xds_image_number = 1
            h5MasterFilePath, h5DataFilePath, h5FileNumber = UtilsImage.getH5FilePath(
                first_image_path, hasOverlap=True, isFastMesh=False
            )
            h5MasterFile = os.path.basename((str(h5MasterFilePath)))
            h5DataFile = os.path.basename((str(h5DataFilePath)))
            list_image_link = [
                [str(h5MasterFilePath), h5MasterFile],
                [str(h5DataFilePath), h5DataFile],
            ]
            if working_directory is not None:
                os.symlink(str(h5MasterFilePath), h5MasterFile)
                os.symlink(str(h5DataFilePath), h5DataFile)
            template = h5MasterFile.replace("master", "??????")
            lowest_xds_image_number = None
            highest_xds_image_number = None
            for subwedge in in_data["subWedge"]:
                image_list = subwedge["image"]
                for image_dict in image_list:
                    image_path = image_dict["path"]
                    image_number = UtilsImage.getImageNumber(image_path)
                    if lowest_xds_image_number is None or lowest_xds_image_number > image_number:
                        lowest_xds_image_number = image_number
                    if highest_xds_image_number is None or highest_xds_image_number < image_number:
                        highest_xds_image_number = image_number
        else:
            template = "%s_xdslink_?????.%s" % (prefix, suffix)
            xds_lowest_image_number_global = 1
            # First we have to find the smallest goniostat rotation axis start:
            oscillation_start_min = 0

            # Loop through the list of sub wedges
            list_of_list = []
            lowest_xds_image_number = None
            highest_xds_image_number = None
            list_spot_range = []
            for sub_wedge in in_data["subWedge"]:
                list_image_link = []
                image_list = sub_wedge["image"]
                goniostat = sub_wedge["experimentalCondition"]["goniostat"]
                oscillation_start = goniostat["rotationAxisStart"]
                oscillation_range = goniostat["oscillationWidth"]

                # First find the lowest and highest image numbers
                lowest_image_number = None
                for image in image_list:
                    image_number = image["number"]
                    if lowest_image_number is None or image_number < lowest_image_number:
                        lowest_image_number = image_number

                # Loop through the list of images
                spot_range_min = None
                spot_range_max = None
                for image in image_list:
                    image_number = image["number"]
                    image_oscillation_start = (
                        oscillation_start
                        + (image_number - lowest_image_number) * oscillation_range
                    )
                    # if xdsLowestImageNumberGlobal is None:
                    #     xdsLowestImageNumberGlobal = 1 + int((imageOscillationStart - oscillationStartMin) / oscillationRange)
                    xds_image_number = xds_lowest_image_number_global + int(
                        (image_oscillation_start - oscillation_start_min) / oscillation_range + 0.5
                    )
                    print(
                        xds_image_number,
                        image_oscillation_start,
                        oscillation_start_min,
                        oscillation_range,
                    )
                    source_path = image["path"]
                    target = "%s_xdslink_%05d.%s" % (prefix, xds_image_number, suffix)
                    print([source_path, target])
                    list_image_link.append([source_path, target])
                    if working_directory is not None and not os.path.exists(target):
                        os.symlink(source_path, target)
                    if (
                        lowest_xds_image_number is None
                        or lowest_xds_image_number > xds_image_number
                    ):
                        lowest_xds_image_number = xds_image_number
                    if (
                        highest_xds_image_number is None
                        or highest_xds_image_number < xds_image_number
                    ):
                        highest_xds_image_number = xds_image_number
                    if spot_range_min is None or spot_range_min > xds_image_number:
                        spot_range_min = xds_image_number
                    if spot_range_max is None or spot_range_max < xds_image_number:
                        spot_range_max = xds_image_number
                list_spot_range.append([spot_range_min, spot_range_max])
                list_of_list.append(list_image_link)
        previous_exclude_data_range_max = 1
        list_exclude_data_range = []
        for spot_range_min, spot_range_max in list_spot_range:
            if spot_range_min > previous_exclude_data_range_max + 1:
                list_exclude_data_range.append([previous_exclude_data_range_max, spot_range_min - 1])
            previous_exclude_data_range_max = spot_range_max + 1
        dictImageLinks = {
            "imageLink": list_of_list,
            "spotRange": list_spot_range,
            "dataRange": [lowest_xds_image_number, highest_xds_image_number],
            "excludeDataRange": list_exclude_data_range,
            "template": template,
        }
        return dictImageLinks

    @staticmethod
    def generateXDS_INP(inData):
        """
        This method creates a list of XDS.INP commands
        """
        # Take the first sub webge in input as reference
        firstSubwedge = inData["subWedge"][0]
        listImage = firstSubwedge["image"]
        image = listImage[0]
        experimentalCondition = firstSubwedge["experimentalCondition"]
        detector = experimentalCondition["detector"]
        dictXDSDetector = XDSTask.getXDSDetector(detector)
        beam = experimentalCondition["beam"]
        goniostat = experimentalCondition["goniostat"]
        distance = round(detector["distance"], 3)
        wavelength = round(beam["wavelength"], 3)
        oscRange = goniostat["oscillationWidth"]
        startAngle = goniostat["rotationAxisStart"] - int(
            goniostat["rotationAxisStart"]
        )
        listXDS_INP = [
            "OVERLOAD=10048500",
            "DIRECTION_OF_DETECTOR_X-AXIS={0}".format(
                UtilsConfig.get("XDSTask", "DIRECTION_OF_DETECTOR_X-AXIS")
            ),
            "DIRECTION_OF_DETECTOR_Y-AXIS={0}".format(
                UtilsConfig.get("XDSTask", "DIRECTION_OF_DETECTOR_Y-AXIS")
            ),
            "ROTATION_AXIS={0}".format(UtilsConfig.get("XDSTask", "ROTATION_AXIS")),
            "INCIDENT_BEAM_DIRECTION={0}".format(
                UtilsConfig.get("XDSTask", "INCIDENT_BEAM_DIRECTION")
            ),
            "NX={0} NY={1} QX={2} QY={2}".format(
                dictXDSDetector["nx"], dictXDSDetector["ny"], dictXDSDetector["pixel"]
            ),
            "ORGX={0} ORGY={1}".format(
                dictXDSDetector["orgX"], dictXDSDetector["orgY"]
            ),
            "DETECTOR={0}  MINIMUM_VALID_PIXEL_VALUE={1}  OVERLOAD={2}".format(
                dictXDSDetector["name"],
                dictXDSDetector["minimumValidPixelValue"],
                dictXDSDetector["overload"],
            ),
            "SENSOR_THICKNESS={0}".format(dictXDSDetector["sensorThickness"]),
            "TRUSTED_REGION={0} {1}".format(
                dictXDSDetector["trustedRegion"][0], dictXDSDetector["trustedRegion"][1]
            ),
        ]
        # for trustedRegion in dictXDSDetector["untrustedRectangle"]:
        #     listXDS_INP.append('UNTRUSTED_RECTANGLE={0} {1} {2} {3}'.format(
        #         trustedRegion[0], trustedRegion[1],
        #         trustedRegion[2],trustedRegion[3]
        #     ))
        listXDS_INP += [
            "DETECTOR_DISTANCE={0}".format(distance),
            "X-RAY_WAVELENGTH={0}".format(wavelength),
            "OSCILLATION_RANGE={0}".format(oscRange),
            "STARTING_ANGLE={0}".format(startAngle),
            "INDEX_QUALITY= 0.25",
        ]
        if "spaceGroupNumber" in inData:
            spaceGroupNumber = inData["spaceGroupNumber"]
            cell = inData["cell"]
            unitCellConstants = "{a} {b} {c} {alpha} {beta} {gamma}".format(**cell)
            listXDS_INP += [
                "SPACE_GROUP_NUMBER={0}".format(spaceGroupNumber),
                "UNIT_CELL_CONSTANTS={0}".format(unitCellConstants),
            ]
        if (
            image["path"].endswith("h5")
            and UtilsConfig.get("XDSTask", "LIB") is not None
        ):
            listXDS_INP += ["LIB= {0}".format(UtilsConfig.get("XDSTask", "LIB"))]
        return listXDS_INP

    @staticmethod
    def createSPOT_XDS(listDozorSpotFile, oscRange):
        """
              implicit none
              integer nmax
              parameter(nmax=10000000)
              real*4 x(3),j
              integer n,i,k
              real*4 xa(nmax,3),ja(nmax)
              logical new
        c
              n=0
              do while(.true.)
                 read(*,*,err=1,end=1)x,j
                 new = .true.
                 do i = n,1,-1
                    if (abs(xa(i,3)-x(3)) .gt. 20.0 ) goto 3
                    do k = 1,2
                       if (abs(x(k)-xa(i,k)) .gt. 6.0) goto 2
                    enddo
                    new = .false.
                    xa(i,:)=(xa(i,:)*ja(i)+x*j)/(ja(i)+j)
                    ja(i)=ja(i)+j
          2         continue
                 enddo
          3       if (new) then
                    n=n+1
                    xa(n,:)=x
                    ja(n)=j
                 endif
              enddo
          1   continue
              do i=1,n
                 write(*,*)xa(i,:), ja(i)
              enddo
              end
        """
        listSpotXds = []
        n = 0
        for dozorSpotFile in listDozorSpotFile:
            # Read the file
            with open(str(dozorSpotFile)) as f:
                dozorLines = f.readlines()
            omega = float(dozorLines[2].split()[1]) % 360
            frame = int((omega - oscRange / 2) / oscRange) + 1
            print(omega, frame)
            for dozorLine in dozorLines[3:]:
                new = True
                listValues = dozorLine.split()
                n, xPos, yPos, intensity, sigma = list(map(float, listValues))
                # Subtracting 1 from X and Y: this is because for dozor the upper left pixel in the image is (1,1),
                # whereas for the rest of the world it is (0,0)
                xPos = xPos - 1
                yPos = yPos - 1
                index = 0
                for spotXds in listSpotXds:
                    frameOld = spotXds[2]
                    if abs(frameOld - frame) > 20:
                        break
                    xPosOld = spotXds[0]
                    yPosOld = spotXds[1]
                    intensityOld = spotXds[3]
                    if abs(xPosOld - xPos) <= 6 and abs(yPosOld - yPos) <= 6:
                        new = False
                        intensityNew = intensity + intensityOld
                        xPosNew = (
                            xPosOld * intensityOld + xPos * intensity
                        ) / intensityNew
                        yPosNew = (
                            yPosOld * intensityOld + yPos * intensity
                        ) / intensityNew
                        listSpotXds[index] = [xPosNew, yPosNew, frameOld, intensityNew]
                    index += 1

                if new:
                    spotXds = [xPos, yPos, frame, intensity]
                    listSpotXds.append(spotXds)

        strSpotXds = ""
        for spotXds in listSpotXds:
            strSpotXds += "{0:13.6f}{1:17.6f}{2:17.8f}{3:17.6f}    \n".format(*spotXds)
        return strSpotXds

    @staticmethod
    def writeSPOT_XDS(listDozorSpotFile, oscRange, workingDirectory):
        spotXds = XDSTask.createSPOT_XDS(listDozorSpotFile, oscRange)
        filePath = workingDirectory / "SPOT.XDS"
        with open(str(filePath), "w") as f:
            f.write(spotXds)

    def writeXDS_INP(self, listXDS_INP, workingDirectory):
        fileName = "XDS.INP"
        filePath = workingDirectory / fileName
        with open(str(filePath), "w") as f:
            for line in listXDS_INP:
                f.write(line + "\n")

    @staticmethod
    def getXDSDetector(dictDetector):
        untrustedRectangle = []
        dictXDSDetector = None
        detectorType = dictDetector["type"]
        nx = UtilsDetector.getNx(detectorType)
        ny = UtilsDetector.getNy(detectorType)
        pixel = UtilsDetector.getPixelsize(detectorType)
        orgX = round(dictDetector["beamPositionX"] / pixel, 3)
        orgY = round(dictDetector["beamPositionY"] / pixel, 3)
        if detectorType == "pilatus2m":
            untrustedRectangle = [
                [487, 495, 0, 1680],
                [981, 989, 0, 1680],
                [0, 1476, 195, 213],
                [0, 1476, 407, 425],
                [0, 1476, 619, 637],
                [0, 1476, 831, 849],
                [0, 1476, 1043, 1061],
                [0, 1476, 1255, 1273],
                [0, 1476, 1467, 1485],
            ]
            sensorThickness = 0.32
        elif detectorType == "pilatus6m":
            untrustedRectangle = [
                [487, 495, 0, 2528],
                [981, 989, 0, 2528],
                [1475, 1483, 0, 2528],
                [1969, 1977, 0, 2528],
                [0, 2464, 195, 213],
                [0, 2464, 407, 425],
                [0, 2464, 619, 637],
                [0, 2464, 831, 849],
                [0, 2464, 1043, 1061],
                [0, 2464, 1255, 1273],
                [0, 2464, 1467, 1485],
                [0, 2464, 1679, 1697],
                [0, 2464, 1891, 1909],
                [0, 2464, 2103, 2121],
                [0, 2464, 2315, 2333],
            ]
            sensorThickness = 0.32
        elif detectorType == "eiger4m":
            untrustedRectangle = [
                [1029, 1040, 0, 2167],
                [0, 2070, 512, 550],
                [0, 2070, 1063, 1103],
                [0, 2070, 1614, 1654],
            ]
            sensorThickness = 0.32
        elif detectorType == "eiger9m":
            untrustedRectangle = [
                [1029, 1040, 0, 3269],
                [2069, 2082, 0, 3269],
                [0, 3110, 513, 553],
                [0, 3110, 1064, 1104],
                [0, 3110, 1615, 1655],
                [0, 3110, 2166, 2206],
                [0, 3110, 2717, 2757],
            ]
        elif detectorType == "eiger16m":
            untrustedRectangle = [
                [0, 4149, 512, 549],
                [0, 4149, 1062, 1099],
                [0, 4149, 1612, 1649],
                [0, 4149, 2162, 2199],
                [0, 4149, 2712, 2749],
                [0, 4149, 3262, 3299],
                [0, 4149, 3812, 3849],
                [513, 514, 0, 4362],
                [1028, 1039, 0, 4362],
                [1553, 1554, 0, 4362],
                [2068, 2079, 0, 4362],
                [2593, 2594, 0, 4362],
                [3108, 3119, 0, 4362],
                [3633, 3634, 0, 4362],
            ]
            sensorThickness = 0.75
        else:
            raise RuntimeError("Unknown detector: {0}".format(detectorType))
        dictXDSDetector = {
            "name": "PILATUS",
            "nx": nx,
            "ny": ny,
            "orgX": orgX,
            "orgY": orgY,
            "pixel": pixel,
            "untrustedRectangle": untrustedRectangle,
            "trustedRegion": [0.0, 1.41],
            "trustedpixel": [7000, 30000],
            "minimumValidPixelValue": 0,
            "overload": 1048500,
            "sensorThickness": sensorThickness,
        }
        return dictXDSDetector


class XDSIndexing(XDSTask):
    def generateXDS_INP(self, inData):
        first_sub_wedge = inData["subWedge"][0]
        # listDozorSpotFile = inData["dozorSpotFile"]
        experimental_condition = first_sub_wedge["experimentalCondition"]
        goniostat = experimental_condition["goniostat"]
        oscRange = goniostat["oscillationWidth"]
        # XDSTask.writeSPOT_XDS(
        #     listDozorSpotFile,
        #     oscRange=oscRange,
        #     workingDirectory=self.getWorkingDirectory(),
        # )
        list_xds_inp = XDSTask.generateXDS_INP(inData)
        list_xds_inp.insert(0, "JOB= XYCORR INIT COLSPOT IDXREF")
        dict_image_links = self.generateImageLinks(inData, self.getWorkingDirectory())
        list_xds_inp.append(
            "NAME_TEMPLATE_OF_DATA_FRAMES= {0}".format(dict_image_links["template"])
        )
        list_spot_range = dict_image_links["spotRange"]
        for spot_range_min, spot_range_max in list_spot_range:
            list_xds_inp.append(
                "SPOT_RANGE= {0} {1}".format(spot_range_min, spot_range_max)
            )
        list_xds_inp.append(
            "DATA_RANGE= {0} {1}".format(
                dict_image_links["dataRange"][0], dict_image_links["dataRange"][1]
            )
        )
        list_spot_range = dict_image_links["excludeDataRange"]
        for exclude_range_min, exclude_range_max in list_spot_range:
            list_xds_inp.append(
                "EXCLUDE_DATA_RANGE= {0} {1}".format(exclude_range_min, exclude_range_max)
            )
        return list_xds_inp

    @staticmethod
    def parseXDSOutput(workingDirectory):
        idxref_path = workingDirectory / "IDXREF.LP"
        xparm_path = workingDirectory / "XPARM.XDS"
        spot_path = workingDirectory / "SPOT.XDS"
        out_data = {
            "idxref": XDSIndexing.readIdxrefLp(idxref_path),
            "xparm": XDSIndexing.parseXparm(xparm_path),
            "xparmXdsPath": xparm_path,
            "spotXdsPath": spot_path
        }
        return out_data

    @staticmethod
    def parseParameters(indexLine, listLines, resultXDSIndexing):
        if "MOSAICITY" in listLines[indexLine]:
            resultXDSIndexing["mosaicity"] = float(listLines[indexLine].split()[-1])
        elif "DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM" in listLines[indexLine]:
            resultXDSIndexing["xBeam"] = float(listLines[indexLine].split()[-1])
            resultXDSIndexing["yBeam"] = float(listLines[indexLine].split()[-2])
        elif "CRYSTAL TO DETECTOR DISTANCE" in listLines[indexLine]:
            resultXDSIndexing["distance"] = float(listLines[indexLine].split()[-1])

    @staticmethod
    def parseLattice(indexLine, listLines, resultXDSIndexing):
        if listLines[indexLine].startswith(" * ") and not listLines[
            indexLine + 1
        ].startswith(" * "):
            listLine = listLines[indexLine].split()
            latticeCharacter = int(listLine[1])
            bravaisLattice = listLine[2]
            spaceGroup = UtilsSymmetry.getMinimumSymmetrySpaceGroupFromBravaisLattice(
                bravaisLattice
            )
            spaceGroupNumber = UtilsSymmetry.getITNumberFromSpaceGroupName(spaceGroup)
            qualityOfFit = float(listLine[3])
            resultXDSIndexing.update(
                {
                    "latticeCharacter": latticeCharacter,
                    "spaceGroupNumber": spaceGroupNumber,
                    "qualityOfFit": qualityOfFit,
                    "a": float(listLine[4]),
                    "b": float(listLine[5]),
                    "c": float(listLine[6]),
                    "alpha": float(listLine[7]),
                    "beta": float(listLine[8]),
                    "gamma": float(listLine[9]),
                }
            )

    @staticmethod
    def readIdxrefLp(pathToIdxrefLp, resultXDSIndexing=None):
        if resultXDSIndexing is None:
            resultXDSIndexing = {}
        if pathToIdxrefLp.exists():
            with open(str(pathToIdxrefLp)) as f:
                listLines = f.readlines()
            indexLine = 0
            doParseParameters = False
            doParseLattice = False
            while indexLine < len(listLines):
                if (
                    "DIFFRACTION PARAMETERS USED AT START OF INTEGRATION"
                    in listLines[indexLine]
                ):
                    doParseParameters = True
                    doParseLattice = False
                elif (
                    "DETERMINATION OF LATTICE CHARACTER AND BRAVAIS LATTICE"
                    in listLines[indexLine]
                ):
                    doParseParameters = False
                    doParseLattice = True
                if doParseParameters:
                    XDSIndexing.parseParameters(indexLine, listLines, resultXDSIndexing)
                elif doParseLattice:
                    XDSIndexing.parseLattice(indexLine, listLines, resultXDSIndexing)
                indexLine += 1
        return resultXDSIndexing

    @staticmethod
    def volum(cell):
        """
        Calculate the cell volum from either:
         - the 6 standard cell parameters (a, b, c, alpha, beta, gamma)
         - or the 3 vectors A, B, C
        Inspired from XOconv written by Pierre Legrand:
        https://github.com/legrandp/xdsme/blob/67001a75f3c363bfe19b8bd7cae999f4fb9ad49d/XOconv/XOconv.py#L758"""
        if len(cell) == 6 and isinstance(cell[0], float):
            # expect a, b, c, alpha, beta, gamma (angles in degree).
            ca, cb, cg = map(XDSIndexing.cosd, cell[3:6])
            return (
                cell[0]
                * cell[1]
                * cell[2]
                * (1 - ca ** 2 - cb ** 2 - cg ** 2 + 2 * ca * cb * cg) ** 0.5
            )
        elif len(cell) == 3 and isinstance(cell[0], np.array):
            # expect vectors of the 3 cell parameters
            A, B, C = cell
            return A * B.cross(C)
        else:
            return "Can't parse input arguments."

    @staticmethod
    def cosd(a):
        return math.cos(a / R2D)

    @staticmethod
    def sind(a):
        return math.sin(a / R2D)

    @staticmethod
    def reciprocal(cell):
        """
        Calculate the 6 reciprocal cell parameters: a*, b*, c*, alpha*, beta*...
        Inspired from XOconv written by Pierre Legrand:
        https://github.com/legrandp/xdsme/blob/67001a75f3c363bfe19b8bd7cae999f4fb9ad49d/XOconv/XOconv.py#L776
        """
        sa, sb, sg = map(XDSIndexing.sind, cell[3:6])
        ca, cb, cg = map(XDSIndexing.cosd, cell[3:6])
        v = XDSIndexing.volum(cell)
        rc = (
            cell[1] * cell[2] * sa / v,
            cell[2] * cell[0] * sb / v,
            cell[0] * cell[1] * sg / v,
            math.acos((cb * cg - ca) / (sb * sg)) * R2D,
            math.acos((ca * cg - cb) / (sa * sg)) * R2D,
            math.acos((ca * cb - cg) / (sa * sb)) * R2D,
        )
        return rc

    @staticmethod
    def BusingLevy(rcell):
        """
        Inspired from XOconv written by Pierre Legrand:
        https://github.com/legrandp/xdsme/blob/67001a75f3c363bfe19b8bd7cae999f4fb9ad49d/XOconv/XOconv.py#L816
        """
        ex = np.array([1, 0, 0])
        ey = np.array([0, 1, 0])
        cosr = list(map(XDSIndexing.cosd, rcell[3:6]))
        sinr = list(map(XDSIndexing.sind, rcell[3:6]))
        Vr = XDSIndexing.volum(rcell)
        BX = ex * rcell[0]
        BY = rcell[1] * (ex * cosr[2] + ey * sinr[2])
        c = rcell[0] * rcell[1] * sinr[2] / Vr
        cosAlpha = (cosr[1] * cosr[2] - cosr[0]) / (sinr[1] * sinr[2])
        BZ = np.array([rcell[2] * cosr[1], -1 * rcell[2] * sinr[1] * cosAlpha, 1 / c])
        return np.array([BX, BY, BZ]).transpose()

    @staticmethod
    def parseXparm(pathToXparmXds):
        """
        Inspired from parse_xparm written by Pierre Legrand:
        https://github.com/legrandp/xdsme/blob/67001a75f3c363bfe19b8bd7cae999f4fb9ad49d/XOconv/XOconv.py#L372
        """
        if pathToXparmXds.exists():
            with open(str(pathToXparmXds)) as f:
                xparm = f.readlines()
            xparamDict = {
                "rot": list(map(float, xparm[1].split()[3:])),
                "beam": list(map(float, xparm[2].split()[1:])),
                "distance": float(xparm[8].split()[2]),
                "originXDS": list(map(float, xparm[8].split()[:2])),
                "A": list(map(float, xparm[4].split())),
                "B": list(map(float, xparm[5].split())),
                "C": list(map(float, xparm[6].split())),
                "cell": list(map(float, xparm[3].split()[1:])),
                "pixel_size": list(map(float, xparm[7].split()[3:])),
                "pixel_numb": list(map(float, xparm[7].split()[1:])),
                "symmetry": int(xparm[3].split()[0]),
                "num_init": list(map(float, xparm[1].split()[:3]))[0],
                "phi_init": list(map(float, xparm[1].split()[:3]))[1],
                "delta_phi": list(map(float, xparm[1].split()[:3]))[2],
                "detector_X": list(map(float, xparm[9].split())),
                "detector_Y": list(map(float, xparm[10].split())),
            }
        else:
            xparamDict = {}
        return xparamDict


class XDSGenerateBackground(XDSTask):
    def generateXDS_INP(self, inData):
        listXDS_INP = XDSTask.generateXDS_INP(inData)
        listXDS_INP.insert(0, "JOB= XYCORR INIT COLSPOT")
        dictImageLinks = self.generateImageLinks(inData, self.getWorkingDirectory())
        listXDS_INP.append(
            "NAME_TEMPLATE_OF_DATA_FRAMES= {0}".format(dictImageLinks["template"])
        )
        listXDS_INP.append(
            "DATA_RANGE= {0} {1}".format(
                dictImageLinks["dataRange"][0], dictImageLinks["dataRange"][1]
            )
        )
        return listXDS_INP

    @staticmethod
    def parseXDSOutput(workingDirectory):
        if (workingDirectory / "BKGINIT.cbf").exists():
            outData = {
                "gainCbf": str(workingDirectory / "GAIN.cbf"),
                "blankCbf": str(workingDirectory / "BLANK.cbf"),
                "bkginitCbf": str(workingDirectory / "BKGINIT.cbf"),
                "xCorrectionsCbf": str(workingDirectory / "X-CORRECTIONS.cbf"),
                "yCorrectionsCbf": str(workingDirectory / "Y-CORRECTIONS.cbf"),
            }
        else:
            outData = {}
        return outData


class XDSIntegration(XDSTask):
    def generateXDS_INP(self, inData):
        # Copy XPARM.XDS, GAIN.CBF file
        shutil.copy(inData["xparmXds"], self.getWorkingDirectory())
        shutil.copy(inData["gainCbf"], self.getWorkingDirectory())
        shutil.copy(inData["xCorrectionsCbf"], self.getWorkingDirectory())
        shutil.copy(inData["yCorrectionsCbf"], self.getWorkingDirectory())
        shutil.copy(inData["blankCbf"], self.getWorkingDirectory())
        shutil.copy(inData["bkginitCbf"], self.getWorkingDirectory())
        listXDS_INP = XDSTask.generateXDS_INP(inData)
        listXDS_INP.insert(0, "JOB= DEFPIX INTEGRATE CORRECT")
        dictImageLinks = self.generateImageLinks(inData, self.getWorkingDirectory())
        listXDS_INP.append(
            "NAME_TEMPLATE_OF_DATA_FRAMES= {0}".format(dictImageLinks["template"])
        )
        listXDS_INP.append(
            "DATA_RANGE= {0} {1}".format(
                dictImageLinks["dataRange"][0], dictImageLinks["dataRange"][1]
            )
        )
        return listXDS_INP

    @staticmethod
    def parseXDSOutput(workingDirectory):
        outData = {}
        if (workingDirectory / "XDS_ASCII.HKL").exists():
            outData = {
                "xdsAsciiHkl": str(workingDirectory / "XDS_ASCII.HKL"),
                "correctLp": str(workingDirectory / "CORRECT.LP"),
                "bkgpixCbf": str(workingDirectory / "BKGPIX.cbf"),
            }
        return outData


class XDSIndexAndIntegration(XDSTask):
    def generateXDS_INP(self, inData):
        listXDS_INP = XDSTask.generateXDS_INP(inData)
        listXDS_INP.insert(0, "JOB= XYCORR INIT IDXREF COLSPOT DEFPIX INTEGRATE CORRECT")
        dictImageLinks = self.generateImageLinks(inData, self.getWorkingDirectory())
        listXDS_INP.append(
            "NAME_TEMPLATE_OF_DATA_FRAMES= {0}".format(dictImageLinks["template"])
        )
        no_background_images = min(
            (dictImageLinks["dataRange"][1] - dictImageLinks["dataRange"][0]),
            4
        )
        listXDS_INP.append("BACKGROUND_RANGE= {0} {1}".format(
            dictImageLinks["dataRange"][0],
            dictImageLinks["dataRange"][0]+no_background_images-1
        ))
        listXDS_INP.append(
            "SPOT_RANGE= {0} {1}".format(
                dictImageLinks["dataRange"][0], dictImageLinks["dataRange"][1]
            )
        )
        listXDS_INP.append(
            "DATA_RANGE= {0} {1}".format(
                dictImageLinks["dataRange"][0], dictImageLinks["dataRange"][1]
            )
        )
        return listXDS_INP

    @staticmethod
    def parseXDSOutput(workingDirectory):
        outData = {}
        if (workingDirectory / "XDS_ASCII.HKL").exists():
            outData = {
                "xdsAsciiHkl": str(workingDirectory / "XDS_ASCII.HKL"),
                "correctLp": str(workingDirectory / "CORRECT.LP"),
                "bkgpixCbf": str(workingDirectory / "BKGPIX.cbf"),
            }
        return outData
