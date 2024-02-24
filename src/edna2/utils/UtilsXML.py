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
__date__ = "15/05/2023"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxv1/src/EDHandlerESRFPyarchv1_0.py

import os
import time
import pathlib
import tempfile
import xmltodict
import json

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()

def jsonFromXML(filePath) -> str:
    with open(filePath,"r") as fp:
        xmlFile = fp.read()
    orderedDict = xmltodict.parse(xmlFile)
    return json.dumps(orderedDict)

def dictfromXML(filePath) -> dict:
    with open(filePath,"r") as fp:
        xmlFile = fp.read()
    orderedDict = xmltodict.parse(xmlFile)
    return json.loads(json.dumps(orderedDict))
