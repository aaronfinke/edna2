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

__authors__ = ['O. Svensson']
__license__ = 'MIT'
__date__ = '21/04/2019'

import time, socket

from edna2.utils import UtilsLogging

from edna2.tasks.AbstractTask import AbstractTask

logger = UtilsLogging.getLogger()


class SlurmTestTask(AbstractTask):
    """
    Test task for testing Slurm execution
    """
    

    def run(self, inData):
        name = inData["name"]
        logger.info(f"name: {name}")
        cmd = "echo $HOSTNAME \n"
        cmd += "sleep 10"
        returnCode = self.submitCommandLine(cmd, jobName="EDNA2_test", partition="all", ignoreErrors=False)
        if returnCode != 0:
            self.setFailure()

        outData = {
            "hostName" : socket.gethostname(),
            "returnCode":returnCode
        }
        return outData
