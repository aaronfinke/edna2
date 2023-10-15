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

from edna2.utils import UtilsTest
from edna2.utils import UtilsLogging

from edna2.tasks.HelloWorldTask import HelloWorldTask

logger = UtilsLogging.getLogger()


class HelloWorldTaskExecTest(unittest.TestCase):

    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    def test_execute(self):
        referenceDataPath = self.dataPath / 'inDataHelloWorldTask.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        helloWorldTask = HelloWorldTask(inData=inData)
        helloWorldTask.execute()
        self.assertTrue(helloWorldTask.isSuccess())
        outData = helloWorldTask.outData
        self.assertTrue('results' in outData)
        self.assertTrue('Hello world' in outData['results'])
        logger.info('Results from HelloWorldTask: {0}'.format(outData['results']))
if __name__ == "__main__":
    unittest.main()