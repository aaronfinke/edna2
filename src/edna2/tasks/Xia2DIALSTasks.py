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
__date__ = "25/04/2023"

import os
import shutil
import tempfile
import numpy as np
from pathlib import Path
import gzip
import time
import re
import json
from datetime import datetime

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsPath
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsXML


logger = UtilsLogging.getLogger()

from edna2.tasks.ISPyBTasks import ISPyBStoreAutoProcResults, ISPyBStoreAutoProcStatus
from edna2.tasks.WaitFileTask import WaitFileTask


STRF_TEMPLATE = "%a %b %d %H:%M:%S %Y"

class Xia2DialsTask(AbstractTask):
    def setFailure(self):
        self._dictInOut["isFailure"] = True
        if self.integrationId is not None and self.programId is not None:
            ISPyBStoreAutoProcResults.setIspybToFailed(
                dataCollectionId=self.dataCollectionId,
                autoProcProgramId=self.programId, 
                autoProcIntegrationId=self.integrationId, 
                processingCommandLine=self.processingCommandLine, 
                processingPrograms=self.processingPrograms, 
                isAnom=False, 
                timeStart=self.startDateTime, 
                timeEnd=datetime.now().isoformat(timespec="seconds")
            )

    pass

