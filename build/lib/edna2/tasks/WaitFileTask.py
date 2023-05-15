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
# mxPluginExec/plugins/EDPluginMXWaitFile-v1.1/plugins/EDPluginMXWaitFilev1_1.py

import pathlib

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()

DEFAULT_TIMEOUT = 120  # s


class WaitFileTask(AbstractTask):
    """
    This plugin waits for a file to appear on disk
    """

    def run(self, inData):
        # Wait for file if it's not already on disk'
        if "file" not in inData:
            raise BaseException("No expected file path in input!")
        filePath = pathlib.Path(self.inData["file"])
        expectedSize = inData.get("expectedSize", None)
        configTimeOut = UtilsConfig.get(self, "timeOut", DEFAULT_TIMEOUT)
        timeOut = inData.get("timeOut", configTimeOut)
        hasTimedOut, finalSize = UtilsPath.waitForFile(
            filePath, expectedSize=expectedSize, timeOut=timeOut
        )
        outData = {"timedOut": hasTimedOut}
        if finalSize is not None:
            outData["finalSize"] = finalSize
        return outData
