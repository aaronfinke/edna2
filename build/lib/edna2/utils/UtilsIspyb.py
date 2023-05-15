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
__date__ = "05/09/2019"

import datetime
import os
import json
import time
import requests

from suds.client import Client
from suds.sax.date import DateTime
from suds.transport.https import HttpAuthenticated

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


def getDataFromURL(url):
    if "http_proxy" in os.environ:
        os.environ["http_proxy"] = ""
    response = requests.get(url)
    data = {"statusCode": response.status_code}
    if response.status_code == 200:
        data["data"] = json.loads(response.text)[0]
    else:
        data["text"] = response.text
    return data


def getRawDataFromURL(url):
    if "http_proxy" in os.environ:
        os.environ["http_proxy"] = ""
    response = requests.get(url)
    data = {"statusCode": response.status_code}
    if response.status_code == 200:
        data["content"] = response.content
    else:
        data["text"] = response.text
    return data


def getWdslRoot():
    dictConfig = UtilsConfig.getTaskConfig("ISPyB")
    wdslRoot = dictConfig["ispyb_ws_url"]
    return wdslRoot


def getTransport():
    transport = None
    logger = UtilsLogging.getLogger()
    if "ISPyB_user" not in os.environ:
        logger.error("No ISPyB user name defined as environment variable!")
    elif "ISPyB_pass" not in os.environ:
        logger.error("No ISPyB password defined as environment variable!")
    else:
        ispybUserName = os.environ["ISPyB_user"]
        ispybPassword = os.environ["ISPyB_pass"]
        transport = HttpAuthenticated(username=ispybUserName, password=ispybPassword)
    return transport


def getCollectionWebService():
    logger = UtilsLogging.getLogger()
    collectionWdsl = getToolsForCollectionWebService()
    transport = getTransport()
    if transport is None:
        logger.error(
            "No transport defined, ISPyB web service client cannot be instantiated."
        )
        collectionWSClient = None
    else:
        collectionWSClient = Client(collectionWdsl, transport=transport, cache=None)
    return collectionWSClient


def getToolsForCollectionWebService():
    return os.path.join(getWdslRoot(), "ispybWS", "ToolsForCollectionWebService?wsdl")


def getToolsForAutoprocessingWebService():
    return os.path.join(getWdslRoot(), "ispybWS", "ToolsForAutoprocessingWebService?wsdl")


def getAutoprocessingWebService():
    logger = UtilsLogging.getLogger()
    collectionWdsl = getToolsForAutoprocessingWebService()
    transport = getTransport()
    if transport is None:
        logger.error(
            "No transport defined, ISPyB web service client cannot be instantiated."
        )
        autoprocessingWSClient = None
    else:
        autoprocessingWSClient = Client(collectionWdsl, transport=transport, cache=None)
    return autoprocessingWSClient

def findDataCollection(dataCollectionId, client=None):
    e = None
    dataCollectionWS3VO = None
    noTrials = 5
    logger = UtilsLogging.getLogger()
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                "No web service client available, cannot contact findDataCollection web service."
            )
        elif dataCollectionId is None:
            logger.error(
                "No dataCollectionId given, cannot contact findDataCollection web service."
            )
        else:
            dataCollectionWS3VO = client.service.findDataCollection(dataCollectionId)
    except Exception as e:
        logger.error(
            "ISPyB error for findDataCollection: {0}, {1} trials left".format(
                e, noTrials
            )
        )
    return dataCollectionWS3VO


def findDataCollectionFromFileLocationAndFileName(imagePath, client=None):
    logger = UtilsLogging.getLogger()
    dataCollectionWS3VO = None
    noTrials = 10
    fileLocation = os.path.dirname(imagePath)
    fileName = os.path.basename(imagePath)
    if fileName.endswith(".h5"):
        prefix = UtilsImage.getPrefix(fileName)
        imageNumber = UtilsImage.getImageNumber(fileName)
        fileName = "{0}_{1:04d}.h5".format(prefix, imageNumber)
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                "No web service client available, cannot contact findDataCollectionFromFileLocationAndFileName web service."
            )
        elif fileLocation is None:
            logger.error(
                "No fileLocation given, cannot contact findDataCollectionFromFileLocationAndFileName web service."
            )
        elif fileName is None:
            logger.error(
                "No fileName given, cannot contact findDataCollectionFromFileLocationAndFileName web service."
            )
        else:
            dataCollectionWS3VO = (
                client.service.findDataCollectionFromFileLocationAndFileName(
                    fileLocation, fileName
                )
            )
    except Exception as e:
        logger.error(
            "ISPyB error for findDataCollectionFromFileLocationAndFileName: {0}, {1} trials left".format(
                e, noTrials
            )
        )
        raise e
    if dataCollectionWS3VO is None:
        time.sleep(1)
        if noTrials == 0:
            logger.error("No data collections found for path {0}".format(imagePath))
        else:
            logger.warning(
                "Cannot find {0} in ISPyB - retrying, {1} trials left".format(
                    imagePath, noTrials
                )
            )
    return dataCollectionWS3VO


def setImageQualityIndicatorsPlot(dataCollectionId, plotFile, csvFile):
    logger = UtilsLogging.getLogger()
    client = getCollectionWebService()
    if client is None:
        logger.error(
            "No web service client available, cannot contact setImageQualityIndicatorsPlot web service."
        )
    returnDataCollectionId = client.service.setImageQualityIndicatorsPlot(
        dataCollectionId, plotFile, csvFile
    )
    return returnDataCollectionId


def storeOrUpdateAutoProcProgram(
    programs,
    commandline,
    status,
    message=None,
    start_time=None,
    end_time=None,
    environment=None,
    record_time_stamp=None,
    autoProcProgramId=None,
):
    if start_time is None:
        start_time = datetime.datetime.now()
    if end_time is None:
        end_time = datetime.datetime.now()
    if record_time_stamp is None:
        record_time_stamp = datetime.datetime.now()
    client = getAutoprocessingWebService()
    autoProcProgramId = client.service.storeOrUpdateAutoProcProgram(
        arg0=autoProcProgramId,
        processingCommandLine=commandline,
        processingPrograms=programs,
        processingStatus=status,
        processingMessage=message,
        processingStartTime=DateTime(start_time),
        processingEndTime=DateTime(end_time),
        processingEnvironment=environment,
        recordTimeStamp=record_time_stamp
    )
    return autoProcProgramId

def storeOrUpdateAutoProcProgramAttachment(
    auto_proc_program_id,
    file_type,
    file_path,
    record_time_stamp=None,
    auto_proc_program_attachment_id=None,
):
    if record_time_stamp is None:
        record_time_stamp = DateTime(datetime.datetime.now())
    client = getAutoprocessingWebService()
    auto_proc_program_attachment_id = client.service.storeOrUpdateAutoProcProgramAttachment(
        arg0=auto_proc_program_attachment_id,
        fileType=file_type,
        fileName=os.path.basename(file_path),
        filePath=os.path.dirname(file_path),
        recordTimeStamp=record_time_stamp,
        autoProcProgramId=auto_proc_program_id
    )
    return auto_proc_program_attachment_id

def storeOrUpdateAutoProcIntegration(
    auto_proc_program_id,
    dataCollection_id,
    start_image_number=None,
    end_image_number=None,
    auto_proc_integration_id=None,
    refined_detector_distance=None,
    refined_x_beam=None,
    refined_y_beam=None,
    rotation_axis_x=None,
    rotation_axis_y=None,
    rotation_axis_z=None,
    beam_vector_x=None,
    beam_vector_y=None,
    beam_vector_z=None,
    cell_a=None,
    cell_b=None,
    cell_c=None,
    cell_alpha=None,
    cell_beta=None,
    cell_gamma=None,
    record_time_stamp=None,
    anomalous=None,
):
    if record_time_stamp is None:
        record_time_stamp = DateTime(datetime.datetime.now())
    client = getAutoprocessingWebService()
    auto_proc_integration_id = client.service.storeOrUpdateAutoProcIntegration(
        arg0=auto_proc_integration_id,
        autoProcProgramId=auto_proc_program_id,
        startImageNumber=start_image_number,
        endImageNumber=end_image_number,
        refinedDetectorDistance=refined_detector_distance,
        refinedXbeam=refined_x_beam,
        refinedYbeam=refined_y_beam,
        rotationAxisX=rotation_axis_x,
        rotationAxisY=rotation_axis_y,
        rotationAxisZ=rotation_axis_z,
        beamVectorX=beam_vector_x,
        beamVectorY=beam_vector_y,
        beamVectorZ=beam_vector_z,
        cellA=cell_a,
        cellB=cell_b,
        cellC=cell_c,
        cellAlpha=cell_alpha,
        cellBeta=cell_beta,
        cellGamma=cell_gamma,
        recordTimeStamp=record_time_stamp,
        anomalous=anomalous,
        dataCollectionId=dataCollection_id
    )
    return auto_proc_integration_id
