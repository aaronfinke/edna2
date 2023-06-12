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

__authors__ = ["O. Svensson"]
__license__ = "MIT"
__date__ = "12/05/2023"

import sys
import unittest
import inspect
from datetime import datetime
import time

from edna2.utils import UtilsTest
from edna2.utils import UtilsConfig
from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults

class ISPyBStoreAutoProcResultsTest(unittest.TestCase):
        
    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)
    
    @unittest.skipIf(UtilsConfig.getSite() == 'Default',
                     'Cannot run ispyb test with default config')

    def test_execute_storeAutoProcResultsAutoPROC(self):
        referenceDataPath = self.dataPath / 'ISPyBStoreAutoProcResults_autoPROC.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        self.dataCollectionId = inData["dataCollectionId"]
        self.processingCommandLine = ""
        self.processingPrograms = "autoPROC_test"
        self.timeStart = datetime.now().isoformat(timespec="seconds")

        self.integrationId, self.programId = ISPyBStoreAutoProcResults.setIspybToRunning(
            dataCollectionId=self.dataCollectionId,
            processingCommandLine = self.processingCommandLine,
            processingPrograms = self.processingPrograms,
            isAnom = False,
            timeStart = self.timeStart)

        self.assertFalse(self.integrationId is None)
        self.assertFalse(self.programId is None)
        t = 10
        while t:
            timer = '{:02d} s remaining'.format(t)
            print(timer, end="\r")
            time.sleep(1)
            t -= 1
        inData["autoProcProgram"]["autoProcProgramId"] = self.programId
        inData["autoProc"]["autoProcProgramId"] = self.programId
        inData["autoProcIntegration"]["autoProcProgramId"] = self.programId
        inData["autoProcIntegration"]["autoProcIntegrationId"] = self.integrationId

        ispybStoreAutoProcResults = ISPyBStoreAutoProcResults(inData=inData)
        ispybStoreAutoProcResults.execute()

        self.assertTrue(ispybStoreAutoProcResults)

if __name__ == '__main__':
    unittest.main()



