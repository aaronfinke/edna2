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

import unittest
import tempfile
import subprocess
import pathlib

from edna2.utils import UtilsTest
from edna2.utils import UtilsLogging

from edna2.tasks.WaitFileTask import WaitFileTask

logger = UtilsLogging.getLogger()


class MXWaitFileExecTest(unittest.TestCase):
    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    def test_WaitFile(self):
        referenceDataPath = self.dataPath / "WaitFile_reference.json"
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        waitFile = WaitFileTask(inData=inData)
        waitFile.execute()
        outData = waitFile.outData
        logger.info(outData)
        self.assertFalse(outData["timedOut"])
        self.assertEqual(outData["finalSize"], 2974985859)

    def test_WaitFileSizeTooBig(self):
        referenceDataPath = self.dataPath / "WaitFile_reference_sizetoobig.json"
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        waitFile = WaitFileTask(inData=inData)
        waitFile.execute()
        outData = waitFile.outData
        logger.info(outData)
        self.assertTrue(outData["timedOut"])

    def test_WaitFilesDownloading(self):
        tmpdir = tempfile.TemporaryDirectory()
        cmd = f"rsync -a /data/visitors/biomax/20170251/20231031/raw/Tau/Tau-natA3/Tau-natA3_1_*.h5 {tmpdir.name}/."
        proc = subprocess.Popen(cmd, shell=True)
        lastFile = pathlib.Path(tmpdir.name) / "Tau-natA3_1_data_000009.h5"
        inData = {"file": lastFile, "size": 400000000}
        waitFile = WaitFileTask(inData=inData)
        waitFile.execute()
        outData = waitFile.outData
        logger.info(outData)
        tmpdir.cleanup()
        self.assertFalse(outData["timedOut"])
        self.assertEqual(outData["finalSize"], 487765579)
