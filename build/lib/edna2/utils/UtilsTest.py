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
# kernel/tests/src/EDTestCasePlugin.py (for loadTestImage)

import os
import re
import json
import pathlib
import datetime
import tempfile
import threading

from edna2.utils import UtilsLogging
from edna2.utils import UtilsImage

from urllib.request import urlopen, ProxyHandler, build_opener

logger = UtilsLogging.getLogger()

URL_EDNA_SITE = "http://www.edna-site.org/data/tests/images"
MAX_DOWNLOAD_TIME = 300


def __timeoutDuringDownload():
    """
    Function called after a timeout in the download part, raises RuntimeError.
    """
    raise RuntimeError(
        "Could not automatically download test images!\n"
        + "If you are behind a firewall,"
        + "please set the environment variable http_proxy.\n"
        + "Otherwise please try to download the images manually from\n"
        + "http://www.edna-site.org/data/tests/images"
    )


def getTestdataPath():
    pathFile = pathlib.Path(__file__)
    pathProjectBase = pathFile.parents[3]
    testdataPath = pathProjectBase / "tests" / "testdata"
    return testdataPath


def getTestRunPath():
    now = datetime.datetime.now()
    dateTime = now.strftime("%Y%m%d_%H%M%S")
    runPath = getTestdataPath() / "rundir" / dateTime
    if not runPath.exists():
        runPath.mkdir(0o755, parents=True)
    return runPath


def getTestImageDirPath():
    dirPath = getTestdataPath() / "images"
    if not dirPath.exists():
        dirPath.mkdir(0o755, parents=True)
    return dirPath


def loadTestImage(imageFileName):
    """
    This method tries to download images from
    http://www.edna-site.org/data/tests/images
    """
    imageDirPath = getTestImageDirPath()
    if not imageDirPath.exists():
        imageDirPath.mkdir(mode=0o777, parents=True)
    # Check if h5 data
    if imageFileName.endswith(".h5"):
        imagePath = imageDirPath / imageFileName
        hasOverlap = imagePath.name.startswith("ref-")
        h5MasterFilePath, h5DataFilePath, h5FileNumber = UtilsImage.getH5FilePath(
            imagePath, hasOverlap=hasOverlap
        )
        listImagePath = [h5MasterFilePath, h5DataFilePath]
    else:
        listImagePath = [imageDirPath / imageFileName]
    for imagePath in listImagePath:
        if not imagePath.exists():
            logger.info(
                "Trying to download image %s" % str(imagePath)
                + ", timeout set to %d s" % MAX_DOWNLOAD_TIME
            )
            if "http_proxy" in os.environ:
                dictProxies = {"http": os.environ["http_proxy"]}
                proxy_handler = ProxyHandler(dictProxies)
                opener = build_opener(proxy_handler).open
            else:
                opener = urlopen

            timer = threading.Timer(MAX_DOWNLOAD_TIME + 1, __timeoutDuringDownload)
            timer.start()
            data = opener(
                "%s/%s" % (URL_EDNA_SITE, imagePath.name),
                data=None,
                timeout=MAX_DOWNLOAD_TIME,
            ).read()
            timer.cancel()

            try:
                open(str(imagePath), "wb").write(data)
            except IOError:
                raise IOError(
                    "Unable to write downloaded data to disk at %s" % imagePath
                )

        if os.path.exists(str(imagePath)):
            logger.info("Image %s successfully downloaded." % imagePath)
        else:
            raise RuntimeError(
                "Could not automatically download test image %r!\n"
                + "If you are behind a firewall, "
                + "please set the environment variable http_proxy.\n"
                + "Otherwise please try to download the images manually from\n"
                + "http://www.edna-site.org/data/tests/images" % imageFileName
            )


def substitute(data, searchString, substituteString):
    dataString = json.dumps(data)
    newDataString = dataString.replace(searchString, substituteString)
    newData = json.loads(newDataString)
    return newData


def getSearchStringFileNames(searchString, data):
    listFileNames = []
    dataString = json.dumps(data)
    # With suffix
    expression = r'"\{0}/([\w-]+\.\w+)"'.format(searchString)
    for match in re.findall(expression, dataString):
        listFileNames.append(match)
    if len(listFileNames) == 0:
        # Try without suffix
        expression = r'"\{0}/([\w-]+)"'.format(searchString)
        for match in re.findall(expression, dataString):
            listFileNames.append(match)
    return listFileNames


def substitueTestData(inData, loadTestImages=True, taskDataPath=None, tmpDir=None):
    # $EDNA2_TESTDATA_IMAGES
    searchString = "$EDNA2_TESTDATA_IMAGES"
    substituteString = getTestImageDirPath().as_posix()
    listFileNames = getSearchStringFileNames(searchString, inData)
    if loadTestImages:
        for imageFileName in listFileNames:
            loadTestImage(imageFileName)
    newInData = substitute(inData, searchString, substituteString)
    # $EDNA2_TASK_DATA
    if taskDataPath is not None:
        searchString = "$EDNA2_TASK_DATA"
        substituteString = str(taskDataPath)
        newInData = substitute(newInData, searchString, substituteString)
    # $EDNA2_TMP_DATA
    if tmpDir is not None:
        searchString = "$EDNA2_TMP_DATA"
        substituteString = str(tmpDir)
        newInData = substitute(newInData, searchString, substituteString)
    # Any other environment variables...
    newInData = json.loads(os.path.expandvars(json.dumps(newInData)))
    return newInData


def loadAndSubstitueTestData(dataPath, loadTestImages=True, tmpDir=None):
    with open(str(dataPath)) as f:
        inData = json.loads(f.read())
    taskDataPath = dataPath.parent
    return substitueTestData(
        inData, loadTestImages=loadTestImages, taskDataPath=taskDataPath, tmpDir=tmpDir
    )


def createTestTmpDirectory(testName):
    return tempfile.mkdtemp(prefix=testName + "_")


def prepareTestDataPath(modulePath):
    dataPath = pathlib.Path(modulePath).parent / "data"
    os.chdir(str(getTestRunPath()))
    return dataPath
