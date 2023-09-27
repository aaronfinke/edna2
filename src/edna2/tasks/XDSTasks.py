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
from pathlib import Path
import sys
import re
import h5py

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
        xdsSetup = UtilsConfig.get(self, 'xdsSetup')
        if xdsSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + xdsSetup + '\n'
        xdsExecutable = UtilsConfig.get(self, 'xdsExecutable',"xds_par")
        commandLine += xdsExecutable
        listXDS_INP = self.generateXDS_INP(inData)
        self.writeXDS_INP(listXDS_INP, self.getWorkingDirectory())
        self.setLogFileName("xds.log")
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        partition = UtilsConfig.get("XDSTask","slurm_partition",None)
        if self.onlineAutoProcessing:
            self.submitCommandLine(commandLine, jobName="EDNA2_XDS", partition=partition, ignoreErrors=False)
        else:
            self.runCommandLine(commandLine, listCommand=[])
        #check for errors
        errorList = self.checkLogForWarningAndErrors()
        # Work in progress!
        outData = self.parseXDSOutput(self.getWorkingDirectory())
        return outData

    @staticmethod
    def generateImageLinks(in_data, working_directory=None):
        first_sub_wedge = in_data["subWedge"][0]
        first_image_path = first_sub_wedge["image"][0]["path"]
        prefix = UtilsImage.getPrefix(first_image_path)
        suffix = UtilsImage.getSuffix(first_image_path)
        if suffix == "h5" and UtilsConfig.isMAXIV():
            lowest_xds_image_number = 1
            highest_xds_image_number = 1
            h5MasterFilePath, h5DataFilePath, h5FileNumber = UtilsImage.getH5FilePath(
                first_image_path, hasOverlap=False, isFastMesh=False
            )
            h5MasterFile = os.path.basename((str(h5MasterFilePath)))
            h5DataFile = os.path.basename((str(h5DataFilePath)))
            list_image_link = [
                [str(h5MasterFilePath), h5MasterFile],
                [str(h5DataFilePath), h5DataFile],
            ]
            if working_directory is not None:
                Path(working_directory / h5MasterFile).symlink_to(h5MasterFilePath)
                find_prefix = re.search(r"(\S+)_data_\d{6}.h5",str(h5DataFile))
                prefix = find_prefix.groups()[0]
                for datafile in h5DataFilePath.parent.glob(prefix + "_data_*"):
                    Path(working_directory / datafile.name).symlink_to(datafile)
            template = h5MasterFile.replace("master", "??????")
            lowest_xds_image_number = None
            highest_xds_image_number = None
            spot_range_min = None
            spot_range_max = None
            list_of_list = []
            list_spot_range = []
            # for subwedge in in_data["subWedge"]:
            #     image_list = subwedge["image"]
            #     for image_dict in image_list:
            #         image_path = image_dict["path"]
            #         lowest_xds_image_number = UtilsImage.getImageNumber(image_path)

            #grab last image number from master file
            with h5py.File(h5MasterFile,'r') as master_file:
                data_file = list(master_file['/entry/data'].keys())[0]
                lowest_xds_image_number = int(master_file['/entry/data'][data_file].attrs['image_nr_low'])
                last_data = list(master_file['/entry/data'].keys())[-1]
                highest_xds_image_number = int(master_file['/entry/data'][last_data].attrs['image_nr_high'])
            list_spot_range.append([lowest_xds_image_number,highest_xds_image_number])

        elif suffix == "h5":
            lowest_xds_image_number = 1
            highest_xds_image_number = 1
            h5MasterFilePath, h5DataFilePath, h5FileNumber = UtilsImage.getH5FilePath(
                first_image_path, hasOverlap=False, isFastMesh=False
            )
            h5MasterFile = os.path.basename((str(h5MasterFilePath)))
            h5DataFile = os.path.basename((str(h5DataFilePath)))
            list_image_link = [
                [str(h5MasterFilePath), h5MasterFile],
                [str(h5DataFilePath), h5DataFile],
            ]
            if working_directory is not None:
                Path(working_directory / h5MasterFile).symlink_to(h5MasterFilePath)
                find_prefix = re.search(r"(\S+)_data_\d{6}.h5",str(h5DataFile))
                prefix = find_prefix.groups()[0]
                for datafile in h5DataFilePath.parent.glob(prefix + "_data_*"):
                    Path(working_directory / datafile.name).symlink_to(datafile)

            
            template = h5MasterFile.replace("master", "??????")
            lowest_xds_image_number = None
            highest_xds_image_number = None
            spot_range_min = None
            spot_range_max = None
            list_spot_range = []
            list_of_list = []

            for subwedge in in_data["subWedge"]:
                image_list = subwedge["image"]
                for image_dict in image_list:
                    image_path = image_dict["path"]
                    image_number = UtilsImage.getImageNumber(image_path)
                    if lowest_xds_image_number is None or lowest_xds_image_number > image_number:
                        lowest_xds_image_number = image_number
                    if highest_xds_image_number is None or highest_xds_image_number < image_number:
                        highest_xds_image_number = image_number
                    if spot_range_min is None or spot_range_min > lowest_xds_image_number:
                        spot_range_min = lowest_xds_image_number
                    if spot_range_max is None or spot_range_max < highest_xds_image_number:
                        spot_range_max = highest_xds_image_number
                list_spot_range.append([spot_range_min, spot_range_max])

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
        startAngle = round(goniostat["rotationAxisStart"] - int(
            goniostat["rotationAxisStart"]
        ),4)
        listXDS_INP = [
            "DELPHI= {0}".format(
                UtilsConfig.get("XDSTask", "DELPHI")
            ),
            "NUMBER_OF_IMAGES_IN_CACHE= {0}".format(
                UtilsConfig.get("XDSTask", "NUMBER_OF_IMAGES_IN_CACHE")
            ),
            "MAXIMUM_NUMBER_OF_JOBS= {0}".format(
                UtilsConfig.get("XDSTask", "MAXIMUM_NUMBER_OF_JOBS")
            ),
            "MAXIMUM_NUMBER_OF_PROCESSORS= {0}".format(
                UtilsConfig.get("XDSTask", "MAXIMUM_NUMBER_OF_PROCESSORS")
            ),
            "INCLUDE_RESOLUTION_RANGE= 50.0 0.0",
            "OVERLOAD={0}".format(
                UtilsConfig.get("XDSTask", "OVERLOAD")
            ),
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
                UtilsConfig.get("XDSTask", "OVERLOAD"),
            ),
            "SENSOR_THICKNESS={0}".format(dictXDSDetector["sensorThickness"]),
            "TRUSTED_REGION={0} {1}".format(
                dictXDSDetector["trustedRegion"][0], dictXDSDetector["trustedRegion"][1]
            ),
            "FRACTION_OF_POLARIZATION= {0}".format(
                UtilsConfig.get("XDSTask", "FRACTION_OF_POLARIZATION")
            ),
            "POLARIZATION_PLANE_NORMAL= {0}".format(
                UtilsConfig.get("XDSTask", "POLARIZATION_PLANE_NORMAL")
            ),
            "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= {0}".format(
                UtilsConfig.get("XDSTask", "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS")
            ),
            "STRONG_PIXEL= {0}".format(
                UtilsConfig.get("XDSTask", "STRONG_PIXEL")
            ),
            "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= {0}".format(
                UtilsConfig.get("XDSTask", "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT")
            ),
            "SEPMIN= {0}".format(
                UtilsConfig.get("XDSTask", "SEPMIN")
            ),
            "CLUSTER_RADIUS= {0}".format(
                UtilsConfig.get("XDSTask", "CLUSTER_RADIUS")
            ),
            "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA= {0}".format(
                UtilsConfig.get("XDSTask", "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA_BETA")
            ),
            "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= {0}".format(
                UtilsConfig.get("XDSTask", "NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA")
            ),
            "REFINE(IDXREF)= {0}".format(
                UtilsConfig.get("XDSTask", "REFINE_IDXREF")
            ),
            "REFINE(INTEGRATE)= {0}".format(
                UtilsConfig.get("XDSTask", "REFINE_INTEGRATE")
            ),
            "REFINE(CORRECT)= {0}".format(
                UtilsConfig.get("XDSTask", "REFINE_CORRECT")
            ),
            ""

        ]
        for trustedRegion in dictXDSDetector["untrustedRectangle"]:
            listXDS_INP.append('UNTRUSTED_RECTANGLE={0} {1} {2} {3}'.format(
                trustedRegion[0], trustedRegion[1],
                trustedRegion[2],trustedRegion[3]
            ))
        listXDS_INP += [
            "DETECTOR_DISTANCE={0}".format(distance),
            "X-RAY_WAVELENGTH={0}".format(wavelength),
            "OSCILLATION_RANGE={0}".format(oscRange),
            "STARTING_ANGLE={0}".format(startAngle),
#            "INDEX_QUALITY= 0.3",
        ]
        if inData.get("idxref",None) is not None:
            if inData.get("idxref",None).get("spaceGroupNumber",0) != 0:
                spaceGroupNumber = inData.get("idxref").get("spaceGroupNumber")
                unitCell = inData.get("idxref").get("unitCell")
                unitCellConstants = "{cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**unitCell)
                listXDS_INP += [
                    "SPACE_GROUP_NUMBER={0}".format(spaceGroupNumber),
                    "UNIT_CELL_CONSTANTS={0}".format(unitCellConstants),
                ]

        elif inData.get("spaceGroupNumber",0) != 0:
            spaceGroupNumber = inData["spaceGroupNumber"]
            unitCell = inData["unitCell"]
            unitCellConstants = "{cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**unitCell)
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
                [0, 4150, 513, 553],
                [0, 4150, 1064, 1104],
                [0, 4150, 1615, 1655],
                [0, 4150, 2166, 2206],
                [0, 4150, 2717, 2757],
                [0, 4150, 3268, 3308],
                [0, 4150, 3819, 3859],
                [1029, 1042, 0, 4371],
                [2069, 2082, 0, 4371],
                [3109, 3122, 0, 4371],
            ]
            sensorThickness = UtilsDetector.getSensorThickness(detectorType)

        else:
            raise RuntimeError("Unknown detector: {0}".format(detectorType))
        dictXDSDetector = {
            "name": "EIGER",
            "nx": nx,
            "ny": ny,
            "orgX": orgX,
            "orgY": orgY,
            "pixel": pixel,
            "untrustedRectangle": untrustedRectangle,
            "trustedRegion": [0.0, 1.2],
            "trustedpixel": [7000, 30000],
            "minimumValidPixelValue": 0,
            "sensorThickness": sensorThickness,
        }
        return dictXDSDetector
    
    @staticmethod
    def parseCorrectLp(inData):
        outData = {}
        correctLPPath = Path(inData["correctLp"])
        try:
            with open(correctLPPath,'r') as fp:
                lines =  [l.strip('\n') for l in fp.readlines()]
        except IOError:
            logger.error("Could not open the specified XDS output file for reading: {0}".format(inData["correctLp"]))
            return None
        try:
            limits = ([i for i,s in enumerate(lines) if 'REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES' in s][0],
            [i for i,s in enumerate(lines) if 'MEAN INTENSITY AS FUNCTION OF SPINDLE POSITION WITHIN DATA IMAGE' in s][0])
        except IndexError:
            logger.error("Could not extract the data from the specified XDS output file: {0}".format(inData["correctLp"]))

        isaLine = [lines.index(x) for x in lines if '     a        b          ISa' in x][0]
        a,b,Isa = [float(x) for x in lines[isaLine+1].split()]


        refinedDiffractionParams = XDSTask._extractRefinedDiffractionParams(lines[limits[0]:limits[1]])
        completeness_entry_begin = [i for i,s in enumerate(lines) if 'LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected' in s][-1]
        completeness_entry_end = [i for i,s in enumerate(lines[completeness_entry_begin:]) if 'total' in s][0]
        completenessEntries = XDSTask._extractCompletenessEntries(lines[completeness_entry_begin+1:completeness_entry_begin+completeness_entry_end+1])

        # now for the last bit: check if we were given a path to the
        # gxparm file and if it exists get the space group and unit
        # cell constants from it
        gxparm_data = {}
        if Path(inData["gxParmXds"]).is_file():
            gxparm_data = XDSIndexing.parseXparm(Path(inData["gxParmXds"]))
        outData['ISa'] = Isa
        outData['wilsonA'] = a
        outData['wilsonB'] = b
        outData['refinedDiffractionParams'] = refinedDiffractionParams
        outData['completenessEntries'] = completenessEntries
        outData['gxparmData'] = gxparm_data

        return outData

    
    @staticmethod
    def _extractRefinedDiffractionParams(lines):
        """
        Get the refined diffraction parameters from CORRECT.LP.
        Returns None if it can't find all values.
        """
        outData = {}
        try:
            crystal_mosaicity =  list(filter(lambda element: 'CRYSTAL MOSAICITY (DEGREES)' in element, lines))[0].split()[-1]
            direct_beam_coordinates = list(filter(lambda element: 'DIRECT BEAM COORDINATES (REC. ANGSTROEM)' in element, lines))[0].split()[-3:]
            direct_beam_detector_coordinates = list(filter(lambda element: 'DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM' in element, lines))[0].split()[-2:]
            detector_origin = list(filter(lambda element: 'DETECTOR ORIGIN (PIXELS) AT' in element, lines))[0].split()[-2:]
            crystal_to_detector_distance = list(filter(lambda element: 'CRYSTAL TO DETECTOR DISTANCE (mm)' in element, lines))[0].split()[-1]
            coordinates_of_unit_cell_a_axis = list(filter(lambda element: 'COORDINATES OF UNIT CELL A-AXIS' in element, lines))[0].split()[-3:]
            coordinates_of_unit_cell_b_axis = list(filter(lambda element: 'COORDINATES OF UNIT CELL B-AXIS' in element, lines))[0].split()[-3:]
            coordinates_of_unit_cell_c_axis = list(filter(lambda element: 'COORDINATES OF UNIT CELL C-AXIS' in element, lines))[0].split()[-3:]
        except IndexError as idx:
            logger.error("Could not extract all refinement parameters from CORRECT.LP")

        outData['crystal_mosaicity'] = float(crystal_mosaicity)
        outData['direct_beam_coordinates'] = [float(x) for x in direct_beam_coordinates]
        outData['direct_beam_detector_coordinates'] = [float(x) for x in direct_beam_detector_coordinates]
        outData['detector_origin'] = [float(x) for x in detector_origin]
        outData['crystal_to_detector_distance'] = float(crystal_to_detector_distance)
        outData['coordinates_of_unit_cell_a_axis'] =  [float(x) for x in coordinates_of_unit_cell_a_axis]
        outData['coordinates_of_unit_cell_b_axis'] =  [float(x) for x in coordinates_of_unit_cell_b_axis]
        outData['coordinates_of_unit_cell_c_axis'] =  [float(x) for x in coordinates_of_unit_cell_c_axis]
        unit_cell_constants = list(filter(lambda element: 'UNIT_CELL_CONSTANTS=' in element, lines))[0].split()[-6:]
        outData['cell_a'], outData['cell_b'], outData['cell_c'],\
            outData['cell_alpha'], outData['cell_beta'], outData['cell_gamma'] = [float(x) for x in unit_cell_constants]

        return outData

    @staticmethod
    def _extractCompletenessEntries(lines):
        """
        Since the latest XDS version there's no guarantee the fields
        will be separated by whitespace. What's fixed is the size of
        each field. So we'll now use fixed offsets to extract the
        fields.

        The Fortran code uses this format statement:
        1130  FORMAT(F9.2,I12,I8,I10,F11.1,'%',F10.1,'%',F9.1,'%',I9,F8.2,      &
                F8.1,'%',F8.1,A1,I6,A1,F8.3,I8)
        """
        logger.debug(f'extracting completeness entries...')
        outData = {"completeness_entries": []}
        offsets = {
        'res': (0,9),
        'observed': (9, 21),
        'unique': (21, 29),
        'possible': (29, 39),
        'complete': (39, 50),
        'rfactor': (51, 61),
        'isig': (81, 89),
        'rmeas': (89, 98),
        'half_dataset_correlation': (99, 107),
        'include_res_based_on_cc': (107,109),
        }

        for line in lines:
            #deal with blank lines and total completeness
            if line.strip() == '':
                continue
            if 'total' in line:
                res_dict = {}
                for (name, (start, end)) in offsets.items():
                    if name == 'include_res_based_on_cc':
                        continue
                    value = float(line[start:end]) if not 'total' in line[start:end] else 'total'
                    res_dict[name] = value
                res_dict.pop('res', None)
                outData['total_completeness'] = res_dict
            else:
                res_dict = {}
                for (name, (start, end)) in offsets.items():
                    if name == 'include_res_based_on_cc':
                        res_dict[name] = True if '*' in line[start:end] else False
                    else:
                        value = float(line[start:end])
                        res_dict[name] = value
                outData['completeness_entries'].append(res_dict)
        return outData

    def checkLogForWarningAndErrors(self):
        """Checks the plugin/XDS log file for warning and error messages"""
        errorList = []
        if self.onlineAutoProcessing:
            strLog = self.getSlurmLog()
        else:
            strLog = self.getLog()
        listLog = strLog.split("\n")
        for strLogLine in listLog:
            # Check for Errors
            if "!!! ERROR " in strLogLine:
                errorList.append(strLogLine)
                logger.error(strLogLine)
        return errorList
                    

        

class XDSIndexing(XDSTask):
    def run(self, inData):
        xdsSetup = UtilsConfig.get('XDSTask', 'xdsSetup')
        if xdsSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + xdsSetup + '\n'
        xdsExecutable = UtilsConfig.get(self, 'xdsExecutable',"xds_par")
        commandLine += xdsExecutable
        listXDS_INP = self.generateXDS_INP(inData)
        self.writeXDS_INP(listXDS_INP, self.getWorkingDirectory())
        self.setLogFileName("xds.log")
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        partition = UtilsConfig.get("XDSTask","slurm_partition",None)
        if self.onlineAutoProcessing:
            self.submitCommandLine(commandLine, jobName="EDNA2_XDS", partition=partition, ignoreErrors=False)
            log = self.getSlurmLogFileName()
        else:
            self.runCommandLine(commandLine, listCommand=[])
            log = self.getLogFileName()
        #check for errors
        errorList = self.checkLogForWarningAndErrors()
        insufficientIndexing = False
        if errorList:
            # ignore this error
            insufficientIndexing = True if 'INSUFFICIENT PERCENTAGE (< 50%) OF INDEXED REFLECTIONS' in ' '.join(errorList) else False
            if not insufficientIndexing:
                self.setFailure()
        # Work in progress!
        outData = self.parseXDSOutput(self.getWorkingDirectory())
        outData['insufficientIndexing'] = insufficientIndexing
        return outData

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
                "SPOT_RANGE= {0} {1}".format(spot_range_min, spot_range_max/2)
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
        gainCbf_path = workingDirectory / "GAIN.cbf"
        blankCbf_path = workingDirectory / "BLANK.cbf"
        xCorrectionsCbf_path = workingDirectory / "X-CORRECTIONS.cbf"
        yCorrectionsCbf_path = workingDirectory / "Y-CORRECTIONS.cbf"
        bkginitCbf_path = workingDirectory / "BKGINIT.cbf"

        with open(workingDirectory / "XDS.INP", "r") as fp:
            for line in fp:
                if "DATA_RANGE=" in line:
                    first_image, last_image = [int(x) for x in line.split() if not "DATA_RANGE" in x]
        out_data = {
            "workingDirectory": workingDirectory,
            "xdsInp": str(workingDirectory / "XDS.INP"),
            "idxref": XDSIndexing.readIdxrefLp(idxref_path),
            "xparm": XDSIndexing.parseXparm(Path(xparm_path)),
            "xparmXds": xparm_path,
            "spotXds": spot_path,
            "gainCbf" : gainCbf_path,
            "blankCbf" : blankCbf_path,
            "xCorrectionsCbf" : xCorrectionsCbf_path,
            "yCorrectionsCbf" : yCorrectionsCbf_path,
            "bkginitCbf" : bkginitCbf_path,
            "start_image" : first_image,
            "end_image" : last_image,
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
                    "unitCell": {
                        "cell_a": float(listLine[4]),
                        "cell_b": float(listLine[5]),
                        "cell_c": float(listLine[6]),
                        "cell_alpha": float(listLine[7]),
                        "cell_beta": float(listLine[8]),
                        "cell_gamma": float(listLine[9]),
                        }
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
        if pathToXparmXds.is_file():
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

    def run(self, inData):
        xdsSetup = UtilsConfig.get('XDSTask', 'xdsSetup')
        if xdsSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + xdsSetup + '\n'
        xdsExecutable = UtilsConfig.get(self, 'xdsExecutable',"xds_par")
        commandLine += xdsExecutable
        listXDS_INP = self.generateXDS_INP(inData)
        self.writeXDS_INP(listXDS_INP, self.getWorkingDirectory())
        self.setLogFileName("xds.log")
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        partition = UtilsConfig.get("XDSTask","slurm_partition",None)
        if self.onlineAutoProcessing:
            self.submitCommandLine(commandLine, jobName="EDNA2_XDS", partition=partition, ignoreErrors=False)
        else:
            self.runCommandLine(commandLine, listCommand=[])
        #check for errors
        errorList = self.checkLogForWarningAndErrors()
        # xds succeeds if XDS_ASCII.HKL is generated
        if not (self.getWorkingDirectory() / "XDS_ASCII.HKL").exists():
            self.setFailure()
        outData = self.parseXDSOutput(self.getWorkingDirectory())
        return outData

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
                "workingDirectory": str(workingDirectory),
                "xdsInp": str(workingDirectory / "XDS.INP"),
                "xdsAsciiHkl": str(workingDirectory / "XDS_ASCII.HKL"),
                "integrateHkl": str(workingDirectory / "INTEGRATE.HKL"),
                "integrateLp": str(workingDirectory / "INTEGRATE.LP"),
                "correctLp": str(workingDirectory / "CORRECT.LP"),
                "bkgpixCbf": str(workingDirectory / "BKGPIX.cbf"),
                "gxParmXds":str(workingDirectory / "GXPARM.XDS"),
            }
            correctLpData = XDSTask.parseCorrectLp(outData)
            outData["ISa"] = correctLpData["ISa"]
            outData["refinedDiffractionParams"] = correctLpData["refinedDiffractionParams"]
            outData["completenessEntries"] = correctLpData["completenessEntries"]["completeness_entries"]
            outData["total_completeness"] = correctLpData["completenessEntries"]["total_completeness"]
            outData["gxparmData"] = correctLpData["gxparmData"]
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
                "workingDirectory": str(workingDirectory),
                "xdsAsciiHkl": str(workingDirectory / "XDS_ASCII.HKL"),
                "correctLp": str(workingDirectory / "CORRECT.LP"),
                "bkgpixCbf": str(workingDirectory / "BKGPIX.cbf"),
            }
        return outData

class XDSRerunCorrect(XDSTask):   
    def run(self, inData):
        xdsSetup = UtilsConfig.get('XDSTask', 'xdsSetup')
        if xdsSetup is None:
            commandLine = ""
        else:
            commandLine = ". " + xdsSetup + '\n'
        xdsExecutable = UtilsConfig.get(self, 'xdsExecutable',"xds_par")
        commandLine += xdsExecutable
        # Copy XPARM.XDS, GAIN.CBF file
        # recycle GXPARM.XDS to XPARM.XDS
        try:
            shutil.copy(inData["gxParmXds"], Path(self.getWorkingDirectory() / 'XPARM.XDS'))
            shutil.copy(inData["xdsInp"],self.getWorkingDirectory())
            shutil.copy(inData["gainCbf"], self.getWorkingDirectory())
            shutil.copy(inData["xCorrectionsCbf"], self.getWorkingDirectory())
            shutil.copy(inData["yCorrectionsCbf"], self.getWorkingDirectory())
            shutil.copy(inData["blankCbf"], self.getWorkingDirectory())
            shutil.copy(inData["bkginitCbf"], self.getWorkingDirectory())
            shutil.copy(inData["integrateHkl"], self.getWorkingDirectory())
        except:
            logger.error("Error copying files to rerun CORRECT.")
            self.setFailure()
            return

        listXDS_INP = self.generateXDS_INP(inData)
        self.writeXDS_INP(listXDS_INP, self.getWorkingDirectory())
        self.setLogFileName("xds.log")
        self.onlineAutoProcessing = inData.get("onlineAutoProcessing",False)
        partition = UtilsConfig.get("XDSTask","slurm_partition",None)
        if self.onlineAutoProcessing:
            self.submitCommandLine(commandLine, jobName="EDNA2_XDS", partition=partition, ignoreErrors=False)
        else:
            self.runCommandLine(commandLine, listCommand=[])
        #check for errors
        errorList = self.checkLogForWarningAndErrors()
        # xds succeeds if XDS_ASCII.HKL is generated
        if not (self.getWorkingDirectory() / "XDS_ASCII.HKL").exists():
            self.setFailure()
        outData = self.parseXDSOutput(self.getWorkingDirectory())
        return outData

     
    def generateXDS_INP(self, inData):
        
        anom = 'FALSE' if inData["isAnom"] else 'TRUE'

        dict_image_links = XDSTask.generateImageLinks(inData, self.getWorkingDirectory())

        with open(inData["xdsInp"],'r') as f:
            listXDS_INP = [x.strip('\n') for x in f.readlines()]

        # take the beam divergence and reflection range from INTEGRATE.LP
        # and put them in as starting parameters
        with open(inData["integrateLp"],'r') as f:
            for line in f:
                if "BEAM_DIVERGENCE=" in line:
                    beam_div = line[1:-1]
                    refl_range = next(f)[1:-1]
        listXDS_INP.append(beam_div)
        listXDS_INP.append(refl_range)

        # find the JOB= line and fix it
        jobLine = [listXDS_INP.index(x) for x in listXDS_INP if 'JOB=' in x][0]
        listXDS_INP[jobLine] = 'JOB=CORRECT'

        # now grab the UNIT_CELL and SPACE_GROUP_NUMBER lines if they exist
        cellLine = [listXDS_INP.index(x) for x in listXDS_INP if 'UNIT_CELL=' in x]
        sgLine = [listXDS_INP.index(x) for x in listXDS_INP if 'SPACE_GROUP_NUMBER=' in x]
        resRange = [listXDS_INP.index(x) for x in listXDS_INP if 'INCLUDE_RESOLUTION_RANGE=' in x]

        if resRange != []:
            listXDS_INP[resRange[0]] = f"INCLUDE_RESOLUTION_RANGE= 50.0 {inData['resCutoff']}"

        if cellLine != []:             
            listXDS_INP[cellLine[0]] = f'UNIT_CELL_CONSTANTS= {inData["cell_from_pointless"]["length_a"]:.2f} {inData["cell_from_pointless"]["length_b"]:.2f} {inData["cell_from_pointless"]["length_c"]:.2f}' \
                                       f'{inData["cell_from_pointless"]["angle_alpha"]:.2f} {inData["cell_from_pointless"]["angle_beta"]:.2f} {inData["cell_from_pointless"]["angle_gamma"]:.2f}'
        else:
            listXDS_INP.append(f'UNIT_CELL_CONSTANTS= {inData["cell_from_pointless"]["length_a"]:.2f} {inData["cell_from_pointless"]["length_b"]:.2f} {inData["cell_from_pointless"]["length_c"]:.2f} ' \
                               f'{inData["cell_from_pointless"]["angle_alpha"]:.2f} {inData["cell_from_pointless"]["angle_beta"]:.2f} {inData["cell_from_pointless"]["angle_gamma"]:.2f}')

        if sgLine != []:
            listXDS_INP[sgLine[0]] = f"SPACE_GROUP_NUMBER= {inData['sg_nr_from_pointless']}"
        else:
            listXDS_INP.append(f"SPACE_GROUP_NUMBER= {inData['sg_nr_from_pointless']}")
        
        listXDS_INP.append(f"FRIEDEL'S_LAW= {anom}")

        return listXDS_INP

    @staticmethod
    def parseXDSOutput(workingDirectory):
        outData = {}
        if (workingDirectory / "XDS_ASCII.HKL").exists():
            outData = {
                "workingDirectory": str(workingDirectory),
                "xdsInp": str(workingDirectory / "XDS.INP"),
                "xdsAsciiHkl": str(workingDirectory / "XDS_ASCII.HKL"),
                "integrateHkl": str(workingDirectory / "INTEGRATE.HKL"),
                "integrateLp": str(workingDirectory / "INTEGRATE.LP"),
                "correctLp": str(workingDirectory / "CORRECT.LP"),
                "bkgpixCbf": str(workingDirectory / "BKGPIX.cbf"),
                "gxParmXds":str(workingDirectory / "GXPARM.XDS"),
            }
            correctLpData = XDSTask.parseCorrectLp(outData)
            outData["ISa"] = correctLpData["ISa"]
            outData["refinedDiffractionParams"] = correctLpData["refinedDiffractionParams"]
            outData["completenessEntries"] = correctLpData["completenessEntries"]["completeness_entries"]
            outData["total_completeness"] = correctLpData["completenessEntries"]["total_completeness"]
            outData["gxparmData"] = correctLpData["gxparmData"]

        return outData
