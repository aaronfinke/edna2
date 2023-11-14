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

# import task to test in this case the process predicted model task within the PhenixTasks
from edna2.tasks.PhenixTasks import PhenixProcessPredictedModelTask

# configure logger
logger = UtilsLogging.getLogger()

# set up execution class by inheriting from unittest.TestCase
class ProcPredModelExecTest(unittest.TestCase):
    
    # set up test case with setting up data path from UtilsTest
    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    # decorator to skip test if site is default and need to be configured
    @unittest.skipIf(UtilsConfig.getSite() == 'Default',
                    'Cannot run Phenix test with default config')
    
    # method to test the execution of the PhenixProcessPredictedModelTask
    # takes the test data path in a json format 
    def test_execute_ProcPredModel(self):
        referenceDataPath = self.dataPath / 'inDataProcPredTask.json'

        # load the input data from the reference data path and define the temporary directory
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)

        # create an instance of the PhenixProcessPredictedModelTask and execute it within the test class
        phenixProcessPredictedModelTask = PhenixProcessPredictedModelTask(inData=inData)
        phenixProcessPredictedModelTask.execute()

        # test run is ok when the task is successful and the output data is available
        # is Success item is defined though the parser structure in the PhenixProcessPredictedModelTask
        self.assertTrue(phenixProcessPredictedModelTask.isSuccess())
        outData = phenixProcessPredictedModelTask.outData
        self.assertTrue(outData['jobComplete'])
        
# instantiate the test case by calling unittest.main()
# an object of TestCase will be instantiated and the test methods will be run
if __name__ == '__main__':
    unittest.main()
