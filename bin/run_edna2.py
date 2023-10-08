#!/usr/bin/env python3
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
__date__ = "13/05/2019"

import os
import sys
import json
import pathlib
import argparse
import importlib
from datetime import datetime

# Set up PYTHONPATH

filePath = pathlib.Path(__file__)
if filePath.parent == pathlib.Path("."):
    cwd = os.getcwd()
    projectHome = pathlib.Path(cwd).parent
else:
    projectHome = filePath.parents[2]
edna2TopLevelDir = projectHome / "src" / "edna2"
sys.path.insert(0, str(edna2TopLevelDir))
from edna2.utils import UtilsLogging

# Parse command line

parser = argparse.ArgumentParser()

parser.add_argument(action="store", dest="taskName", help="Name of EDNA2 task")

parser.add_argument(
    "--inData", action="store", dest="inData", help="Input data in JSON format"
)

parser.add_argument(
    "--inDataFile",
    action="store",
    dest="inDataFile",
    help="File with input data in JSON format",
)

parser.add_argument(
    "--outDataFile",
    action="store",
    dest="outDataFile",
    help="File where output data in JSON format will be stored",
)

parser.add_argument("--debug", action="store_true", help="Debug log level")

parser.add_argument("--warning", action="store_true", help="Warning log level")

parser.add_argument("--error", action="store_true", help="Error log level")

results = parser.parse_args()

taskName = results.taskName
inData = results.inData
inDataFile = results.inDataFile
outDataFile = results.outDataFile
debugLogLevel = results.debug
warningLogLevel = results.warning
errorLogLevel = results.error

if taskName is None:
    parser.print_help()
    sys.exit(1)
elif inData is None and inDataFile is not None:
    with open(inDataFile) as fd:
        inData = fd.read()
else:
    print("Error - no indata provided!")
    parser.print_help()
    sys.exit(1)

# Set up logging

if errorLogLevel:
    logger = UtilsLogging.getLogger("ERROR")
elif warningLogLevel:
    logger = UtilsLogging.getLogger("WARNING")
elif debugLogLevel:
    logger = UtilsLogging.getLogger("DEBUG")
else:
    logger = UtilsLogging.getLogger("INFO")

# Load and run EDNA2 task
timeStr = datetime.now().strftime('%y%m%d%H%M%S')

# edna2 = __import__("edna2.tasks.{0}".format(taskName))
edna2 = importlib.import_module("edna2.tasks.{0}".format(taskName))
# tasks = getattr(edna2, "tasks")
# tasksModule = getattr(edna2, taskName)
TaskClass = getattr(edna2, taskName)

task = TaskClass(inData=json.loads(inData),workingDirectorySuffix="0")
task.execute()
if task.isFailure():
    logger.error("Error when executing {0}!".format(taskName))
else:
    outData = json.dumps(task.outData, indent=4)
    if outDataFile is None:
        print(outData)
    else:
        with open(str(outDataFile), "w") as f:
            f.write(outData)
