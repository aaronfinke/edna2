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

import os
import json
import time
import requests
from datetime import datetime
from pathlib import Path

from suds.client import Client
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
        collectionWSClient = Client(collectionWdsl, transport=transport, cache=None, location=collectionWdsl)
    return collectionWSClient


def getToolsForCollectionWebService():
    return os.path.join(getWdslRoot(), "ispybWS", "ToolsForCollectionWebService?wsdl")

def getAutoprocessingWebService():
    logger = UtilsLogging.getLogger()
    collectionWdsl = getToolsForAutoprocessingWebService()
    transport = getTransport()
    if transport is None:
        logger.error(
            "No transport defined, ISPyB web service client cannot be instantiated."
        )
        collectionWSClient = None
    else:
        collectionWSClient = Client(collectionWdsl, transport=transport, cache=None, location=collectionWdsl)
    return collectionWSClient


def getToolsForAutoprocessingWebService():
    return os.path.join(getWdslRoot(), "ispybWS", "ToolsForAutoprocessingWebService?wsdl")

def getBLSampleWebService():
    logger = UtilsLogging.getLogger()
    collectionWdsl = getToolsForBLSampleWebService()
    transport = getTransport()
    if transport is None:
        logger.error(
            "No transport defined, ISPyB web service client cannot be instantiated."
        )
        collectionWSClient = None
    else:
        collectionWSClient = Client(collectionWdsl, transport=transport, cache=None, location=collectionWdsl)
    return collectionWSClient

def getToolsForBLSampleWebService():
        return os.path.join(getWdslRoot(), "ispybWS", "ToolsForBLSampleWebService?wsdl")

def getProteinAcronymAndSampleNameFromDataCollectionId(dataCollectionId,client=None):
    """get the protein acronym and sample name from a given DataCollectionId."""
    proteinAcronym, sampleName = None, None
    if client is None:
        client = getBLSampleWebService()
    if client is None:
        logger.error(
                "No web service client available, cannot contact getBLSampleWebService web service."
            )
        return proteinAcronym, sampleName
    try:
        sampleId = getSampleIdFromDataCollectionId(dataCollectionId)
        if sampleId is None:
            return proteinAcronym, sampleName
        sampleInfo = client.service.getSampleInformation(sampleId)
        proteinAcronym = sampleInfo.proteinAcronym
        sampleName = sampleInfo.sampleName
    except Exception as e:
        logger.error(f"Could not retrieve protein acronym/sample name for {dataCollectionId}: {e}")

    return proteinAcronym, sampleName
        

def getSampleIdFromDataCollectionId(dataCollectionId,client=None):
    """ get the Sample Id from a given DataCollectionId."""
    blSampleId = None
    if client is None:
        client = getCollectionWebService()
    if client is None:
        logger.error(
                "No web service client available, cannot contact getBLSampleWebService web service."
            )
        return blSampleId
    try:
        blSampleVO = client.service.getDataCollectionInfo(dataCollectionId).blSampleVO
        blSampleId = blSampleVO.blSampleId
    except Exception as e:
        logger.error(f"Could not retrieve SampleId from datacollectionId {dataCollectionId}: {e}")
    return blSampleId

def getDataCollectionGroupId(dataCollectionId, client=None):
    """get dataCollectionGroupId from dataCollectionId"""
    dataCollectionGroupId = None
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
            return dataCollectionGroupId
        dataCollectionWS3VO = client.service.findDataCollection(dataCollectionId)
        dataCollectionGroupId = dataCollectionWS3VO.dataCollectionGroupId
    except Exception as e:
        logger.error(
            "ISPyB error for getDataCollectionGroupId: {0} trials left".format(
                e
            )
        )
    return dataCollectionGroupId

def findDataCollectionWS3VO(dataCollectionId,client=None):
    dataCollectionWS3VO = None
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
            return dataCollectionWS3VO
        dataCollectionWS3VO = client.service.findDataCollection(dataCollectionId)
    except Exception as e:
        logger.error(
            "ISPyB error for getDataCollectionGroupId: {0} trials left".format(
                e
            )
        )
    return dataCollectionWS3VO

def findDataCollectionGroupWS3VO(dataCollectionGroupId,client=None):
    dataCollectionGroupWS3VO = None
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
            return dataCollectionGroupWS3VO
        dataCollectionGroupWS3VO = client.service.findDataCollectionGroup(dataCollectionGroupId)
    except Exception as e:
        logger.error(
            "ISPyB error for getDataCollectionGroupId: {0} trials left".format(
                e
            )
        )
    return dataCollectionGroupWS3VO

def storeOrUpdateDataCollection(dataCollectionWS3VO, client=None):
    dataCollectionId = None
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
            return dataCollectionId
        dataCollectionId = client.service.storeOrUpdateDataCollection(dataCollectionWS3VO)
    except Exception as e:
        logger.error(
            "ISPyB error for getDataCollectionGroupId: {0} trials left".format(
                e
            )
        )
    return dataCollectionId

def storeOrUpdateDataCollectionGroup(dataCollectionGroupWS3VO, client=None):
    dataCollectionGroupId = None
    try:
        if client is None:
            client = getCollectionWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
            return dataCollectionGroupId
        dataCollectionGroupId = client.service.storeOrUpdateDataCollectionGroup(dataCollectionGroupWS3VO)
    except Exception as e:
        logger.error(
            "ISPyB error for getDataCollectionGroupId: {0} trials left".format(
                e
            )
        )
    return dataCollectionGroupId


def updateDataCollectionGroupComments(dataCollectionId, comments):
    client = getCollectionWebService()
    iDataCollectionGroupId = None
    dataCollectionGroupId = getDataCollectionGroupId(dataCollectionId, client=client)
    dataCollectionGroupWS3VO = findDataCollectionGroupWS3VO(dataCollectionGroupId, client=client)
    if not comments in str(dataCollectionGroupWS3VO.comments):
        dataCollectionGroupWS3VO.comments = comments
        iDataCollectionGroupId = storeOrUpdateDataCollectionGroup(dataCollectionGroupWS3VO, client=client)
        dataCollectionWS3VO = findDataCollectionWS3VO(dataCollectionId, client=client)
        if hasattr(dataCollectionWS3VO, "comments"):
            if not comments in dataCollectionWS3VO.comments:
                dataCollectionWS3VO.comments += " " + comments
                dataCollectionId = storeOrUpdateDataCollection(dataCollectionWS3VO, client=client)
        else:
            dataCollectionWS3VO.comments = comments
            dataCollectionId = storeOrUpdateDataCollection(dataCollectionWS3VO, client=client)
    logger.debug(f"dataCollectionGroupId: {iDataCollectionGroupId}")





def storeOrUpdateAutoProcProgram(
        autoProcProgramId=None,  
        processingCommandLine=None,
        processingPrograms=None,
        processingStatus=None,
        processingMessage=None,
        processingStartTime=None,
        processingEndTime=None,
        processingEnvironment=None,
        client=None):
    """
    Store Autoproc status: RUNNING,SUCCESS,FAILED,TIMEOUT
    """
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcProgramId = client.service.storeOrUpdateAutoProcProgram(   
            arg0= autoProcProgramId,
            processingCommandLine = processingCommandLine,
            processingPrograms=processingPrograms,
            processingStatus=processingStatus,
            processingMessage=processingMessage,
            processingStartTime=processingStartTime,
            processingEndTime=processingEndTime,
            processingEnvironment=processingEnvironment,
            recordTimeStamp=datetime.now()
        )
    except Exception as e:
        logger.error(
            "ISPyB error for storeOrUpdateAutoProcProgram: {0} trials left".format(
                e
            )
        )
    logger.debug(f"autoProcProgramId: {autoProcProgramId}")
    return autoProcProgramId

def storeOrUpdateAutoProcIntegration(
        autoProcIntegrationId=None,
        autoProcProgramId=None,
        startImageNumber=None,
        endImageNumber=None,
        refinedDetectorDistance= None,
        refinedXbeam=None,
        refinedYbeam=None,
        rotationAxisX=None,
        rotationAxisY=None,
        rotationAxisZ=None,
        beamVectorX=None,
        beamVectorY=None,
        beamVectorZ=None,
        cellA=None,
        cellB=None,
        cellC=None,
        cellAlpha=None,
        cellBeta=None,
        cellGamma=None,
        anomalous=None,
        dataCollectionId=None,
        client=None):
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcIntegrationId = client.service.storeOrUpdateAutoProcIntegration(
        arg0=autoProcIntegrationId,
        autoProcProgramId=autoProcProgramId,
        startImageNumber=startImageNumber,
        endImageNumber=endImageNumber,
        refinedDetectorDistance= refinedDetectorDistance,
        refinedXbeam=refinedXbeam,
        refinedYbeam=refinedYbeam,
        rotationAxisX=rotationAxisX,
        rotationAxisY=rotationAxisY,
        rotationAxisZ=rotationAxisZ,
        beamVectorX=beamVectorX,
        beamVectorY=beamVectorY,
        beamVectorZ=beamVectorZ,
        cellA=cellA,
        cellB=cellB,
        cellC=cellC,
        cellAlpha=cellAlpha,
        cellBeta=cellBeta,
        cellGamma=cellGamma,
        recordTimeStamp=datetime.now(),
        anomalous=anomalous,
        dataCollectionId=dataCollectionId
        )
    except Exception as e:
        logger.error(
            "ISPyB error for storeOrUpdateAutoProcProgram: {0} trials left".format(
                e
            )
    )
    logger.debug(f"autoProcIntegrationId: {autoProcIntegrationId}")
    return autoProcIntegrationId

    
def storeOrUpdateAutoProcProgramAttachment(
        file=None,
        autoProcProgramAttachmentId=None,
        autoProcProgramId=None,
        client=None
        ):
    if not isinstance(file,Path):
        file=Path(file)
    if not file.is_file():
        logger.error(f"File {file} does not exist, skipping storeOrUpdateAutoProcProgramAttachment")
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )

        strFileType = "Result" if file.suffix == ".mtz" or file.suffix == ".gz" else "Log"
        strFileName = file.name
        strFilePath = str(file.parent)
        timeStamp = datetime.now()
        autoProcProgramAttachmentId = client.service.storeOrUpdateAutoProcProgramAttachment(
            arg0=autoProcProgramAttachmentId, 
            fileType=strFileType, 
            fileName=strFileName, 
            filePath=strFilePath, 
            recordTimeStamp=timeStamp, 
            autoProcProgramId=autoProcProgramId
            )    
    except Exception as e:
        logger.error(
            "ISPyB error for autoProcProgramAttachmentId: {0}".format(
                e
            )
    )
    logger.debug(f"autoProcProgramAttachmentId: {autoProcProgramAttachmentId}")
    return autoProcProgramAttachmentId

def storeOrUpdateAutoProc(
        client=None,
        autoProcId=None,
        autoProcProgramId=None,
        spaceGroup=None,
        refinedCellA=None,
        refinedCellB=None,
        refinedCellC=None,
        refinedCellAlpha=None,
        refinedCellBeta=None,
        refinedCellGamma=None,
    ):
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcId = client.service.storeOrUpdateAutoProc(
        arg0=autoProcId,
        autoProcProgramId=autoProcProgramId,
        spaceGroup=spaceGroup,
        refinedCellA=refinedCellA,
        refinedCellB=refinedCellB,
        refinedCellC=refinedCellC,
        refinedCellAlpha=refinedCellAlpha,
        refinedCellBeta=refinedCellBeta,
        refinedCellGamma=refinedCellGamma,
        recordTimeStamp=datetime.now()
        )
    except Exception as e:
        logger.error(
            "ISPyB error for storeOrUpdateAutoProc: {0}".format(
                e
            )
        )

    logger.debug(f"autoProcId: {autoProcId}")
    return autoProcId

def storeOrUpdateAutoProcScalingHasInt(
    autoProcScalingHasIntId=None,
    autoProcIntegrationId=None,
    autoProcScalingId=None,
    client=None,
    ):
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcScalingHasIntId = client.service.storeOrUpdateAutoProcScalingHasInt(
        arg0=autoProcScalingHasIntId,
        autoProcIntegrationId=autoProcIntegrationId,
        autoProcScalingId=autoProcScalingId,
        recordTimeStamp=datetime.now().isoformat(timespec='seconds')
        )
    except Exception as e:
        logger.error(
            "ISPyB error for autoProcScalingHasIntId: {0}".format(
                e
            )
        )
    logger.debug(f"autoProcScalingHasIntId: {autoProcScalingHasIntId}")
    return autoProcScalingHasIntId

def storeOrUpdateAutoProcScaling(
    client=None,
    autoProcScalingId=None,
    autoProcId=None,
    resolutionEllipsoidAxis11=None,
    resolutionEllipsoidAxis12=None,
    resolutionEllipsoidAxis13=None,
    resolutionEllipsoidAxis21=None,
    resolutionEllipsoidAxis22=None,
    resolutionEllipsoidAxis23=None,
    resolutionEllipsoidAxis31=None,
    resolutionEllipsoidAxis32=None,
    resolutionEllipsoidAxis33=None,
    resolutionEllipsoidValue1=None,
    resolutionEllipsoidValue2=None,
    resolutionEllipsoidValue3=None,
    ):

    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcScalingId = client.service.storeOrUpdateAutoProcScaling(
                arg0=autoProcScalingId,
                autoProcId=autoProcId,
                recordTimeStamp=datetime.now().isoformat(timespec='seconds'),
                resolutionEllipsoidAxis11=resolutionEllipsoidAxis11,
                resolutionEllipsoidAxis12=resolutionEllipsoidAxis12,
                resolutionEllipsoidAxis13=resolutionEllipsoidAxis13,
                resolutionEllipsoidAxis21=resolutionEllipsoidAxis21,
                resolutionEllipsoidAxis22=resolutionEllipsoidAxis22,
                resolutionEllipsoidAxis23=resolutionEllipsoidAxis23,
                resolutionEllipsoidAxis31=resolutionEllipsoidAxis31,
                resolutionEllipsoidAxis32=resolutionEllipsoidAxis32,
                resolutionEllipsoidAxis33=resolutionEllipsoidAxis33,
                resolutionEllipsoidValue1=resolutionEllipsoidValue1,
                resolutionEllipsoidValue2=resolutionEllipsoidValue2,
                resolutionEllipsoidValue3=resolutionEllipsoidValue3,
                )
    except Exception as e:
        logger.error(
            "ISPyB error for autoProcScalingId: {0}".format(
                e
            )
        )
    logger.debug(f"autoProcScalingId: {autoProcScalingId}")
    return autoProcScalingId

def storeOrUpdateAutoProcScalingStatistics(
        client=None,
        autoProcScalingStatisticsId=None,
        scalingStatisticsType=None,
        comments=None,
        resolutionLimitLow=None,
        resolutionLimitHigh=None,
        rmerge=None,
        rmeasWithinIplusIminus=None,
        rmeasAllIplusIminus=None,
        rpimWithinIplusIminus=None,
        rpimAllIplusIminus=None,
        fractionalPartialBias=None,
        nTotalObservations=None,
        nTotalUniqueObservations=None,
        meanIoverSigI=None,
        completeness=None,
        multiplicity=None,
        anomalousCompleteness=None,
        anomalousMultiplicity=None,
        anomalous=None,
        autoProcScalingId=None,
        ccHalf=None,
        ccAno=None,
        sigAno=None,
        isa=None,
        completenessSpherical=None,
        anomalousCompletenessSpherical=None,
        completenessEllipsoidal=None,
        anomalousCompletenessEllipsoidal=None,
        recordTimeStamp=datetime.now().isoformat(timespec='seconds'),
):
    
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcScalingStatisticsId = client.service.storeOrUpdateAutoProcScalingStatistics(
                arg0=autoProcScalingStatisticsId,
                scalingStatisticsType=scalingStatisticsType,
                comments=comments,
                resolutionLimitLow=resolutionLimitLow,
                resolutionLimitHigh=resolutionLimitHigh,
                rmerge=rmerge,
                rmeasWithinIplusIminus=rmeasWithinIplusIminus,
                rmeasAllIplusIminus=rmeasAllIplusIminus,
                rpimWithinIplusIminus=rpimWithinIplusIminus,
                rpimAllIplusIminus=rpimAllIplusIminus,
                fractionalPartialBias=fractionalPartialBias,
                nTotalObservations=nTotalObservations,
                nTotalUniqueObservations=nTotalUniqueObservations,
                meanIoverSigI=meanIoverSigI,
                completeness=completeness,
                multiplicity=multiplicity,
                anomalousCompleteness=anomalousCompleteness,
                anomalousMultiplicity=anomalousMultiplicity,
                recordTimeStamp=recordTimeStamp,
                anomalous=anomalous,
                autoProcScalingId=autoProcScalingId,
                ccHalf=ccHalf,
                ccAno=ccAno,
                sigAno=sigAno,
                isa=isa,
                completenessSpherical=completenessSpherical,
                anomalousCompletenessSpherical=anomalousCompletenessSpherical,
                completenessEllipsoidal=completenessEllipsoidal,
                anomalousCompletenessEllipsoidal=anomalousCompletenessEllipsoidal,
                )
    except Exception as e:
        logger.error(
            "ISPyB error for autoProcScalingStatisticsId: {0}".format(
                e
            )
        )
    logger.debug(f"autoProcScalingStatisticsId: {autoProcScalingStatisticsId}")
    return autoProcScalingStatisticsId

def storeOrUpdateAutoProcStatus(
        client=None,
        autoProcStatusId=None,
        autoProcIntegrationId=None,
        step=None,
        status=None,
        comments=None,
        bltimeStamp=None):
    
    try:
        if client is None:
            client = getAutoprocessingWebService()
        if client is None:
            logger.error(
                    "No web service client available, cannot contact findDataAutoprocessing web service."
                )
        autoProcStatusId = client.service.storeOrUpdateAutoProcStatus(
                arg0=autoProcStatusId, 
                autoProcIntegrationId=autoProcIntegrationId, 
                step=step, 
                status=status, 
                comments=comments, 
                bltimeStamp=bltimeStamp, 
                )
    except Exception as e:
        logger.error(
            "ISPyB error for autoProcStatusId: {0}".format(
                e
            )
        )
    logger.debug(f"autoProcStatusId: {autoProcStatusId}")
    return autoProcStatusId




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

def getXDSInfo(dataCollectionId, client=None):
    """
    Collects data from getXDSInfo collection service.
    """
    if dataCollectionId is None:
        logger.error(
                "No dataCollectionId given, cannot contact findDataCollection web service."
            )
        return None
    
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
        dataCollectionWS3VO = client.service.getXDSInfo(dataCollectionId)
    except Exception as e:
        logger.error(
            "ISPyB error for findDataCollection: {0}, {1} trials left".format(
                e, noTrials
            )
        )
    return dataCollectionWS3VO

def getXDSMasterFilePath(dataCollectionId) -> Path:

    if dataCollectionId is None:
        logger.error(
                "No dataCollectionId given, cannot contact findDataCollection web service."
            )
        return None
    
    dataCollectionXDSWS3VO = getXDSInfo(dataCollectionId)
    if dataCollectionXDSWS3VO is None:
        logger.error(
            "No dataCollectionId given, cannot contact findDataCollection web service."
        )

    XDSInfoDict = Client.dict(dataCollectionXDSWS3VO)
    imageDirectory = XDSInfoDict.get("imageDirectory")
    fileTemplate = XDSInfoDict.get("fileTemplate")
    if imageDirectory is None or fileTemplate is None:
        logger.error("No master file found in XDSInfo!")
        return None
    
    imageDirectory = Path(imageDirectory)
    fileTemplate = fileTemplate.replace("%06d","master")
    return imageDirectory / fileTemplate
