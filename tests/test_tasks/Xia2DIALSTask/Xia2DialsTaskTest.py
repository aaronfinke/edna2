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
__date__ = "23/01/2023"

import os
import unittest

from edna2.utils import UtilsTest
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

from edna2.tasks.Xia2DIALSTask import Xia2DialsTask

logger = UtilsLogging.getLogger()
import tracemalloc

class Xia2DialsExecTest(unittest.TestCase):

    def setUp(self):
        tracemalloc.start()
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    def test_execute_Edna2ProcTask(self):
        referenceDataPath = self.dataPath / 'inDataXia2Dials.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        inData["timeOut"] = 1800
        xia2DialsTask = Xia2DialsTask(inData=inData)
        xia2DialsTask.execute()
        self.assertTrue(xia2DialsTask.isSuccess())
        
    # def test_execute_Edna2ProcTask(self):
    #     referenceDataPath = self.dataPath / 'inDataXia2Dials_problems.json'
    #     problemIdFile = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
    #     problemIdList = problemIdFile["dataCollectionIds"]
    #     for problemId in problemIdList:
    #         xia2DialsTask = Xia2DialsTask(inData={
    #             "dataCollectionId": problemId,
    #             "test":True,
    #             "timeOut": 1800
    #         })
    #         xia2DialsTask.execute()
        # self.assertTrue(edna2proctask.isSuccess())

if __name__ == '__main__':
    unittest.main()