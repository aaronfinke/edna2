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

__authors__ = ["O. Svensson"]
__license__ = "MIT"
__date__ = "26/03/2020"

DETECTOR_PARAMS = {
    "pilatus2m": {
        "nx": 1475,
        "ny": 1679,
        "pixel": 0.172,
        "sensorThickness": 0.32,
        "xdsUntrustedRectangle": [
            [487, 495, 0, 1680],
            [981, 989, 0, 1680],
            [0, 1476, 195, 213],
            [0, 1476, 407, 425],
            [0, 1476, 619, 637],
            [0, 1476, 831, 849],
            [0, 1476, 1043, 1061],
            [0, 1476, 1255, 1273],
            [0, 1476, 1467, 1485],
        ],
    },
    "pilatus6m": {
        "nx": 2463,
        "ny": 2527,
        "pixel": 0.172,
        "sensorThickness": 0.32,
        "xdsUntrustedRectangle": [
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
        ],
    },
    "eiger4m": {
        "nx": 2070,
        "ny": 2167,
        "pixel": 0.075,
        "sensorThickness": 0.32,
        "xdsUntrustedRectangle": [
            [1029, 1040, 0, 2167],
            [0, 2070, 512, 550],
            [0, 2070, 1063, 1103],
            [0, 2070, 1614, 1654],
        ],
    },
    "eiger9m": {
        "nx": 3108,
        "ny": 3262,
        "pixel": 0.075,
        "sensorThickness": 0.75,
        "xdsUntrustedRectangle": [],
    },
    "eiger16m": {
        "nx": 4148,
        "ny": 4362,
        "pixel": 0.075,
        "sensorThickness": 0.75,
        "xdsUntrustedRectangle": [],
    },
}


def __getDetectorValue(detectorType, key):
    if detectorType in DETECTOR_PARAMS:
        p = DETECTOR_PARAMS[detectorType]
        return p[key]
    else:
        raise RuntimeError(
            'Detector type "{0}" not defined in UtilsDetector.py!'.format(detectorType)
        )


def getNx(detectorType):
    return __getDetectorValue(detectorType, "nx")


def getNy(detectorType):
    return __getDetectorValue(detectorType, "ny")


def getPixelsize(detectorType):
    return __getDetectorValue(detectorType, "pixel")


def getXdsUntrustedRectangle(detectorType):
    return __getDetectorValue(detectorType, "xdsUntrustedRectangle")


def getSensorThickness(detectorType):
    return __getDetectorValue(detectorType, "sensorThickness")
