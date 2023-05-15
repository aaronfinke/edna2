#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the 'Software'), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ['O. Svensson']
__license__ = 'MIT'
__date__ = '25/03/2021'

import os
import fabio
import pathlib

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.ReadImageHeader import ReadImageHeader

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


class H5ToBinnedCBFTask(AbstractTask):

    def getInDataSchema(self):
        return {
            "type": "object",
            "required": ["hdf5File"],
            "properties": {
                "imageNumber": {"type": "integer"},
                "startImageNumber": {"type": "integer"},
                "hdf5ImageNumber": {"type": "integer"},
                "hdf5File": {"type": "string"},
                "forcedOutputDirectory": {"type": "string"}
            }
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "outputCBFFile": {"type": "string"}
            }
        }

    def run(self, inData):
        outData = {}
        hdf5File = pathlib.Path(inData['hdf5File'])
        # Read the header
        inDataReadImageHeader = {
            "imagePath": [inData['hdf5File']],
            "isFastMesh": False
        }
        readImageHeader = ReadImageHeader(
            inData = inDataReadImageHeader
        )
        readImageHeader.execute()
        if readImageHeader.isSuccess():
            firstSubWedge = readImageHeader.outData["subWedge"][0]
            experimentalCondition = firstSubWedge["experimentalCondition"]
            detector = experimentalCondition["detector"]
            beam = experimentalCondition["beam"]
            goniostat = experimentalCondition["goniostat"]
            beamX = detector["beamPositionX"] / 2.0 / detector["pixelSizeX"]
            beamY = detector["beamPositionY"] / 2.0 / detector["pixelSizeY"]
            directory = hdf5File.parent
            prefix = UtilsImage.getPrefix(hdf5File)
            hdf5ImageNumber = 1
            if 'master.h5' in str(hdf5File):
                masterFile = hdf5File
            else:
                if UtilsConfig.isEMBL():
                    fileName = '{0}_master.h5'.format(prefix)
                else:
                    fileName = '{0}_{1}_master.h5'.format(prefix, hdf5ImageNumber)
                masterFile = directory / fileName
            image = fabio.open(str(masterFile))
            image.data = image.data.reshape(2181, 2, 2074, 2).sum(3).sum(1)
            cbfImage = image.convert("cbf")
            pilatus_headers = fabio.cbfimage.PilatusHeader("Silicon sensor, thickness 0.000750 m")
            pilatus_headers["Start_angle"] = goniostat["rotationAxisStart"]
            pilatus_headers["Angle_increment"] = goniostat["oscillationWidth"]
            pilatus_headers["Detector"] = "Eiger2 16M binned to 4M"
            pilatus_headers["Pixel_size"] = (detector["pixelSizeY"]*2, detector["pixelSizeX"]*2)
            pilatus_headers["Exposure_time"] = beam["exposureTime"]
            pilatus_headers["Wavelength"] = beam["wavelength"]
            pilatus_headers["Detector_distance"] = detector["distance"] / 1000.0
            pilatus_headers["Beam_xy"] = (beamY, beamX)
            pilatus_headers["Count_cutoff"] = 1009869
            cbfImage.pilatus_headers = pilatus_headers
            directory = inData.get("forcedOutputDirectory", str(self.getWorkingDirectory()))
            cbfImagePath = os.path.join(directory, os.path.basename(inData['hdf5File']).replace(".h5", ".cbf"))
            cbfImage.save(cbfImagePath)
            outData["outputCBFFile"] = cbfImagePath
        return outData
