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

__authors__ = ["D. Fastus"]
__license__ = "MIT"
__date__ = "01/10/2023"

import unittest

# import utils to test tasks, config and logging
from edna2.utils import UtilsTest
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

# import task to test in this case the AlphaFoldTask
from edna2.tasks.AlphaFoldTask import AlphaFoldTask

# configure logger
logger = UtilsLogging.getLogger()

# set up execution class by inheriting from unittest.TestCase
class AlphaFoldTaskExecTest(unittest.TestCase):
    
    # set up test case with setting up data path from UtilsTest
    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    # decorator to skip test if site is default and need to be configured
    @unittest.skipIf(UtilsConfig.getSite() == 'Default',
                    'Cannot run AlphaFold test with default config')
    
    # method to test the execution of the AlphaFoldTask
    # takes the test data path in a json format 
    def test_execute_AlphaFoldPrediction(self):
        referenceDataPath = self.dataPath / 'inDataAlphaFoldTask.json'

        # create a temporary directory to store the output data
        tmpDir = UtilsTest.createTestTmpDirectory('AlphaFoldTask')
        # load the input data from the reference data path and define the temporary directory
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath,
                                                    tmpDir=tmpDir)
        # create an instance of the AlphaFoldTask and execute it within the test class
        alphaFoldTask = AlphaFoldTask(inData=inData)
        alphaFoldTask.execute() 

        # test run is ok when the task is successful and the output data is available
        # is Success item is defined though the parser structure in the AlphaFoldTask
        self.assertTrue(alphaFoldTask.isSuccess())
        outData = alphaFoldTask.outData
        self.assertTrue(outData['isSuccess'])
        
# instantiate the test case by calling unittest.main()
# an object of TestCase will be instantiated and the test methods will be run
if __name__ == '__main__':
    unittest.main()
