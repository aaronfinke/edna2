#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the 'Software'), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ['O. Svensson']
__license__ = 'MIT'
__date__ = '21/04/2019'

from datetime import datetime
import pprint
import time

import xmltodict
# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxPluginExec/plugins/EDPluginGroupISPyB-v1.4/plugins/
#     EDPluginISPyBRetrieveDataCollectionv1_4.py


from suds.client import Client
from suds.transport.http import HttpAuthenticated

import os
import gzip
import pathlib

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb

from edna2.tasks.AbstractTask import AbstractTask

logger = UtilsLogging.getLogger()


class ISPyBRetrieveDataCollection(AbstractTask):

    def run(self, inData):
        dictConfig = UtilsConfig.getTaskConfig('ISPyB')
        username = dictConfig['username']
        password = dictConfig['password']
        httpAuthenticated = HttpAuthenticated(username=username,
                                              password=password)
        wdsl = dictConfig['ispyb_ws_url'] + '/ispybWS/ToolsForCollectionWebService?wsdl'
        client = Client(wdsl, transport=httpAuthenticated, cache=None)
        if 'image' in inData:
            path = pathlib.Path(inData['image'])
            indir = path.parent.as_posix()
            infilename = path.name
            dataCollection = client.service.findDataCollectionFromFileLocationAndFileName(
                indir,
                infilename
            )
        elif 'dataCollectionId' in inData:
            dataCollectionId = inData['dataCollectionId']
            dataCollection = client.service.findDataCollection(dataCollectionId)
        else:
            errorMessage = 'Neither image nor data collection id given as input!'
            logger.error(errorMessage)
            raise BaseException(errorMessage)
        if dataCollection is not None:
            outData = Client.dict(dataCollection)
        else:
            outData = {}
        return outData


class GetListAutoprocIntegration(AbstractTask):

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "token": {"type": "string"},
                "proposal": {"type": "string"},
                "dataCollectionId": {"type": "integer"}
            }
        }

    # def getOutDataSchema(self):
    #     return {
    #         "type": "array",
    #         "items": {
    #             "type": "object",
    #             "properties": {
    #                 "AutoProcIntegration_autoProcIntegrationId": {"type": "integer"}
    #             }
    #         }
    #     }

    def run(self, inData):
        # urlExtISPyB, token, proposal, dataCollectionId
        token = inData['token']
        proposal = inData['proposal']
        dataCollectionId = inData['dataCollectionId']
        dictConfig = UtilsConfig.getTaskConfig('ISPyB')
        restUrl = dictConfig['ispyb_ws_url'] + '/rest'
        ispybWebServiceURL = os.path.join(
            restUrl, token, 'proposal', str(proposal), 'mx',
            'autoprocintegration', 'datacollection', str(dataCollectionId),
            'view')
        dataFromUrl = UtilsIspyb.getDataFromURL(ispybWebServiceURL)
        outData = {}
        if dataFromUrl['statusCode'] == 200:
            outData['autoprocIntegration'] = dataFromUrl['data']
        else:
            outData['error'] = dataFromUrl
        return outData


class GetListAutoprocAttachment(AbstractTask):

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "token": {"type": "string"},
                "proposal": {"type": "string"},
                "autoProcProgramId": {"type": "integer"}
            }
        }

    # def getOutDataSchema(self):
    #     return {
    #         "type": "array",
    #         "items": {
    #             "type": "object",
    #             "properties": {
    #                 "AutoProcIntegration_autoProcIntegrationId": {"type": "integer"}
    #             }
    #         }
    #     }

    def run(self, inData):
        # urlExtISPyB, token, proposal, autoProcProgramId
        token = inData['token']
        proposal = inData['proposal']
        autoProcProgramId = inData['autoProcProgramId']
        dictConfig = UtilsConfig.getTaskConfig('ISPyB')
        restUrl = dictConfig['ispyb_ws_url'] + '/rest'
        ispybWebServiceURL = os.path.join(
            restUrl, token, 'proposal', str(proposal), 'mx',
            'autoprocintegration', 'attachment', 'autoprocprogramid',
            str(autoProcProgramId), 'list')
        dataFromUrl = UtilsIspyb.getDataFromURL(ispybWebServiceURL)
        outData = {}
        if dataFromUrl['statusCode'] == 200:
            outData['autoprocAttachment'] = dataFromUrl['data']
        else:
            outData['error'] = dataFromUrl
        return outData


class GetListAutoprocessingResults(AbstractTask):
    """
    This task receives a list of data collection IDs and returns a list
    of dictionaries with all the auto-processing results and file attachments
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "token": {"type": "string"},
                "proposal": {"type": "string"},
                "dataCollectionId": {
                    "type": "array",
                    "items": {
                        "type": "integer",
                    }
                }
            }
        }

    # def getOutDataSchema(self):
    #     return {
    #         "type": "object",
    #         "required": ["dataForMerge"],
    #         "properties": {
    #             "dataForMerge": {
    #                 "type": "object",
    #                 "items": {
    #                     "type": "object",
    #                     "properties": {
    #                         "spaceGroup": {"type": "string"}
    #                     }
    #                 }
    #             }
    #         }
    #     }

    def run(self, inData):
        urlError = None
        token = inData['token']
        proposal = inData['proposal']
        listDataCollectionId = inData['dataCollectionId']
        dictForMerge = {}
        dictForMerge['dataCollection'] = []
        for dataCollectionId in listDataCollectionId:
            dictDataCollection = {
                'dataCollectionId': dataCollectionId
            }
            inDataGetListIntegration = {
                'token': token,
                'proposal': proposal,
                'dataCollectionId': dataCollectionId
            }
            getListAutoprocIntegration = GetListAutoprocIntegration(
                inData=inDataGetListIntegration
            )
            getListAutoprocIntegration.setPersistInOutData(False)
            getListAutoprocIntegration.execute()
            resultAutoprocIntegration = getListAutoprocIntegration.outData
            if 'error' in resultAutoprocIntegration:
                urlError = resultAutoprocIntegration['error']
                break
            else:
                listAutoprocIntegration = resultAutoprocIntegration['autoprocIntegration']
                # Get v_datacollection_summary_phasing_autoProcProgramId
                for autoprocIntegration in listAutoprocIntegration:
                    if 'v_datacollection_summary_phasing_autoProcProgramId' in autoprocIntegration:
                        autoProcProgramId = autoprocIntegration[
                            'v_datacollection_summary_phasing_autoProcProgramId'
                        ]
                        inDataGetListAttachment = {
                            'token': token,
                            'proposal': proposal,
                            'autoProcProgramId': autoProcProgramId
                        }
                        getListAutoprocAttachment = GetListAutoprocAttachment(
                            inData=inDataGetListAttachment
                        )
                        getListAutoprocAttachment.setPersistInOutData(False)
                        getListAutoprocAttachment.execute()
                        resultAutoprocAttachment = getListAutoprocAttachment.outData
                        if 'error' in resultAutoprocAttachment:
                            urlError = resultAutoprocAttachment['error']
                        else:
                            autoprocIntegration['autoprocAttachment'] = resultAutoprocAttachment['autoprocAttachment']
                    dictDataCollection['autoprocIntegration'] = listAutoprocIntegration
            dictForMerge['dataCollection'].append(dictDataCollection)
            # dictForMerge[dataCollectionId] = dictDataCollection
        if urlError is None:
            outData = dictForMerge
        else:
            outData = {
                'error': urlError
            }
        return outData


class RetrieveAttachmentFiles(AbstractTask):
    """
    This task receives a list of data collection IDs and returns a list
    of dictionaries with all the auto-processing results and file attachments
    """

    # def getInDataSchema(self):
    #     return {
    #         "type": "object",
    #         "properties": {
    #             "token": {"type": "string"},
    #             "proposal": {"type": "string"},
    #             "dataCollectionId": {
    #                 "type": "array",
    #                 "items": {
    #                     "type": "integer",
    #                 }
    #             }
    #         }
    #     }

    # def getOutDataSchema(self):
    #     return {
    #         "type": "object",
    #         "required": ["dataForMerge"],
    #         "properties": {
    #             "dataForMerge": {
    #                 "type": "object",
    #                 "items": {
    #                     "type": "object",
    #                     "properties": {
    #                         "spaceGroup": {"type": "string"}
    #                     }
    #                 }
    #             }
    #         }
    #     }

    def run(self, inData):
        urlError = None
        listPath = []
        token = inData['token']
        proposal = inData['proposal']
        listAttachment = inData['attachment']
        dictConfig = UtilsConfig.getTaskConfig('ISPyB')
        restUrl = dictConfig['ispyb_ws_url'] + '/rest'
        # proposal/MX2112/mx/autoprocintegration/autoprocattachmentid/21494689/get
        for dictAttachment in listAttachment:
            attachmentId = dictAttachment['id']
            fileName = dictAttachment['fileName']
            ispybWebServiceURL = os.path.join(
                restUrl, token, 'proposal', str(proposal), 'mx',
                'autoprocintegration', 'autoprocattachmentid', str(attachmentId),
                'get')
            rawDataFromUrl = UtilsIspyb.getRawDataFromURL(ispybWebServiceURL)
            if rawDataFromUrl['statusCode'] == 200:
                rawData = rawDataFromUrl['content']
                if fileName.endswith('.gz'):
                    rawData = gzip.decompress(rawData)
                    fileName = fileName.split('.gz')[0]
                with open(fileName, "wb") as f:
                    f.write(rawData)
                listPath.append(str(self.getWorkingDirectory() / fileName))
            else:
                urlError = rawDataFromUrl
        if urlError is None:
            outData = {
                'filePath': listPath
            }
        else:
            outData = {
                'error': urlError
            }
        return outData


class ISPyBFindDetectorByParam(AbstractTask):

    def run(self, inData):
        dictConfig = UtilsConfig.getTaskConfig('ISPyB')
        username = dictConfig['username']
        password = dictConfig['password']
        httpAuthenticated = HttpAuthenticated(username=username,
                                              password=password)
        wdsl = dictConfig['ispyb_ws_url'] + '/ispybWS/ToolsForCollectionWebService?wsdl'
        client = Client(wdsl, transport=httpAuthenticated, cache=None)
        manufacturer = inData['manufacturer']
        model = inData['model']
        mode = inData['mode']
        detector = client.service.findDetectorByParam(
            "",
            manufacturer,
            model,
            mode
        )
        if detector is not None:
            outData = Client.dict(detector)
        else:
            outData = {}
        return outData


class UploadGPhLResultsToISPyB(AbstractTask):

    def run(self, in_data):
        # Load XML file
        xml_path = in_data["autoPROCXML"]
        programId = in_data.get("programId",None)
        with open(xml_path) as f:
            xml_string = f.read()
        auto_proc_dict = xmltodict.parse(xml_string)
        # pprint.pprint(auto_proc_dict)
        auto_proc_container = auto_proc_dict["AutoProcContainer"]
        auto_proc_program_container = auto_proc_container["AutoProcProgramContainer"]
        # 1. Create AutoProcProgram entry
        auto_proc_program = auto_proc_program_container["AutoProcProgram"]
        auto_proc_program_id = UtilsIspyb.storeOrUpdateAutoProcProgram(
            autoProcProgramId = programId,
            processingPrograms=self.check_length(auto_proc_program["processingPrograms"]),
            processingCommandLine=self.check_length(auto_proc_program["processingCommandLine"]),
            processingStartTime=self.get_time(auto_proc_program["processingStartTime"]),
            processingEndTime=self.get_time(auto_proc_program["processingEndTime"]),
            processingEnvironment=self.check_length(auto_proc_program["processingEnvironment"]),
            processingMessage=self.check_length(auto_proc_program["processingMessage"]),
            processingStatus=auto_proc_program["processingStatus"]
        )
        pprint.pprint(auto_proc_program_id)
        out_data = {}
        return out_data

    def check_length(self, parameter, max_string_length=255):
        if type(parameter) == str and len(parameter) > max_string_length:
            old_parameter = parameter
            parameter = parameter[0:max_string_length - 3] + "..."
            logger.warning(
                "String truncated to %d characters for ISPyB! Original string: %s" % (max_string_length, old_parameter))
            logger.warning("Truncated string: %s" % parameter)
        return parameter

    def get_time(self, time_value):
        # Fri May 12 08:31:54 CEST 2023
        return datetime.datetime.strptime(time_value, "%a %b %d %H:%M:%S %Z %Y")
    
class ISPyBStoreAutoProcResults(AbstractTask):
    """
    Stores the contents of autoprocessing in ISPyB.
    """
    
    @staticmethod
    def getAutoProcProgramContainer():
        return  {
                "autoProcProgramId": None,
                "processingCommandLine": None,
                "processingPrograms": None,
                "processingStatus": None,
                "processingStartTime": None,
                "processingEndTime": None,
                "processingEnvironment": None,
                }

    @staticmethod
    def getAutoProcContainer():
        return {
                "autoProcId": None,
                "autoProcProgramId": None,
                "spaceGroup": None,
                "refinedCellA": None,
                "refinedCellB": None,
                "refinedCellC": None,
                "refinedCellAlpha": None,
                "refinedCellBeta": None,
                "refinedCellGamma": None,
                }
    
    @staticmethod
    def getAutoProcProgramAttachmentContainer():
        return {
                "file": None,
                "autoProcProgramAttachmentId": None,
                "autoProcProgramId": None,
                }
    
    @staticmethod
    def getAutoProcIntegrationContainer():
        return {
                "autoProcIntegrationId": None,
                "autoProcProgramId": None,
                "startImageNumber": None,
                "endImageNumber": None,
                "refinedDetectorDistance": None,
                "refinedXbeam": None,
                "refinedYbeam": None,
                "rotationAxisX": None,
                "rotationAxisY": None,
                "rotationAxisZ": None,
                "beamVectorX": None,
                "beamVectorY": None,
                "beamVectorZ": None,
                "cellA": None,
                "cellB": None,
                "cellC": None,
                "cellAlpha": None,
                "cellBeta": None,
                "cellGamma": None,
                "anomalous": None,
                "dataCollectionId": None,
        }
    
    @staticmethod
    def getAutoProcScalingHasIntContainer():
        return  {                 
                "autoProcScalingHasIntId": None,
                "autoProcIntegrationId": None,
                "autoProcScalingId": None,
                }
    
    @staticmethod
    def getAutoProcScalingContainer():
        return  {
                "autoProcScalingId": None,
                "autoProcId": None,
                "resolutionEllipsoidAxis11": None,
                "resolutionEllipsoidAxis12": None,
                "resolutionEllipsoidAxis13": None,
                "resolutionEllipsoidAxis21": None,
                "resolutionEllipsoidAxis22": None,
                "resolutionEllipsoidAxis23": None,
                "resolutionEllipsoidAxis31": None,
                "resolutionEllipsoidAxis32": None,
                "resolutionEllipsoidAxis33": None,
                "resolutionEllipsoidValue1": None,
                "resolutionEllipsoidValue2": None,
                "resolutionEllipsoidValue3": None,
                }
    def getAutoProcStatusContainer():
        return {
            "autoProcStatusId": None,
            "autoProcIntegrationId": None,
            "step": None,
            "status": None,
            "comments": None,
            "bltimeStamp": None,
        }
    
    @staticmethod
    def getAutoProcScalingStatisticsContainer():
        return {                            
                "autoProcScalingStatisticsId": None,
                "scalingStatisticsType": None,
                "resolutionLimitLow": None,
                "resolutionLimitHigh": None,
                "rmerge": None,
                "rmeasWithinIplusIminus": None,
                "rmeasAllIplusIminus": None,
                "rpimWithinIplusIminus": None,
                "rpimAllIplusIminus": None,
                "fractionalPartialBias": None,
                "nTotalObservations": None,
                "nTotalUniqueObservations": None,
                "meanIoverSigI": None,
                "completeness": None,
                "multiplicity": None,
                "anomalousCompleteness": None,
                "anomalousMultiplicity": None,
                "anomalous": None,
                "autoProcScalingId": None,
                "ccHalf": None,
                "ccAno": None,
                "sigAno": None,
                "isa": None,
                "completenessSpherical": None,
                "anomalousCompletenessSpherical": None,
                "completenessEllipsoidal": None,
                "anomalousCompletenessEllipsoidal": None,
                }
    
    # def getInDataSchema(self):
    #     return {
    #          "$ref": self.getSchemaUrl("ispybAutoprocIntegration.json")
    #     }
    
    # def getOutDataSchema(self):
    #     return {
    #         "type":"object",
    #         "properties": {
    #             "autoProcId": {"type": ["integer","null"] },
    #             "autoProcIntegrationId": {"type": ["integer","null"] },
    #             "autoProcScalingId": {"type": ["integer","null"] },
    #             "autoProcProgramId": {"type":["integer","null"] },
    #         }
    #     }
    
    def run(self,inData):
        outData = {
            "autoProcId": None,
            "autoProcIntegrationId":None,
            "autoProcScalingId": None,
            "autoProcProgramId": None,
        }
        autoProcId = None
        autoProcIntegrationId = None
        autoProcScalingId = None
        autoProcProgramId = None

        client = UtilsIspyb.getAutoprocessingWebService()
        if client is None:
            logger.error("Cannot connect to ISPyB web service")
            self.setFailure()
            return outData
        
        dataCollectionId = inData.get("dataCollectionId", None)
        autoProcProgramData = inData.get("autoProcProgram", self.getAutoProcProgramContainer())
        autoProcProgramId = UtilsIspyb.storeOrUpdateAutoProcProgram(
            autoProcProgramId=autoProcProgramData.get("autoProcProgramId"),  
            processingCommandLine=autoProcProgramData.get("processingCommandLine"),
            processingPrograms=autoProcProgramData.get("processingPrograms"),
            processingStatus=autoProcProgramData.get("processingStatus"),
            processingStartTime=autoProcProgramData.get("processingStartTime"),
            client=client)
        if autoProcProgramId is None:
            logger.error("Couldn't create entry for AutoProcProgram in ISPyB!")
            self.setFailure()
            return outData
        outData["autoProcProgramId"] = autoProcProgramId
        listAutoProcProgramAttachment = inData.get("autoProcProgramAttachment", None)
        if listAutoProcProgramAttachment is not None:
            for program in listAutoProcProgramAttachment:
                autoProcProgramAttachmentId = UtilsIspyb.storeOrUpdateAutoProcProgramAttachment(
                    file=program.get("file"),
                    autoProcProgramAttachmentId=program.get("autoProcProgramAttachmentId"),
                    autoProcProgramId=autoProcProgramId,
                    client=client
                )
                if autoProcProgramAttachmentId is None:
                    logger.error("Error creating attachment point in ISPyB!")
                program["autoProcProgramAttachmentId"] = autoProcProgramAttachmentId
                program["autoProcProgramId"] = autoProcProgramId

        autoProcIntegrationData = inData.get("autoProcIntegration", self.getAutoProcIntegrationContainer())
        autoProcIntegrationId = UtilsIspyb.storeOrUpdateAutoProcIntegration(
                autoProcIntegrationId=autoProcIntegrationData.get("autoProcIntegrationId"),
                autoProcProgramId=autoProcProgramId,
                startImageNumber=autoProcIntegrationData.get("startImageNumber"),
                endImageNumber=autoProcIntegrationData.get("endImageNumber"),
                refinedDetectorDistance=autoProcIntegrationData.get("refinedDetectorDistance"),
                refinedXbeam=autoProcIntegrationData.get("refinedXbeam"),
                refinedYbeam=autoProcIntegrationData.get("refinedYbeam"),
                rotationAxisX=autoProcIntegrationData.get("rotationAxisX"),
                rotationAxisY=autoProcIntegrationData.get("rotationAxisY"),
                rotationAxisZ=autoProcIntegrationData.get("rotationAxisZ"),
                beamVectorX=autoProcIntegrationData.get("beamVectorX"),
                beamVectorY=autoProcIntegrationData.get("beamVectorY"),
                beamVectorZ=autoProcIntegrationData.get("beamVectorZ"),
                cellA=autoProcIntegrationData.get("cellA"),
                cellB=autoProcIntegrationData.get("cellB"),
                cellC=autoProcIntegrationData.get("cellC"),
                cellAlpha=autoProcIntegrationData.get("cellAlpha"),
                cellBeta=autoProcIntegrationData.get("cellBeta"),
                cellGamma=autoProcIntegrationData.get("cellGamma"),
                anomalous=autoProcIntegrationData.get("anomalous", False),
                dataCollectionId=dataCollectionId,
                client=client)
        if autoProcIntegrationId is None:
            logger.warning("Couldn't create entry for AutoProcIntegration in ISPyB!")
        outData["autoProcIntegrationId"] = autoProcIntegrationId
        if autoProcProgramData.get("processingStatus", None) == "FAILED":
            return outData
        autoProcData = inData.get("autoProc", None)
        if autoProcData is not None:
            autoProcId = UtilsIspyb.storeOrUpdateAutoProc(
                    client=client,
                    autoProcId=autoProcData.get("autoProcId"),
                    autoProcProgramId=autoProcProgramId,
                    spaceGroup=autoProcData.get("spaceGroup"),
                    refinedCellA=autoProcData.get("refinedCellA"),
                    refinedCellB=autoProcData.get("refinedCellB"),
                    refinedCellC=autoProcData.get("refinedCellC"),
                    refinedCellAlpha=autoProcData.get("refinedCellAlpha"),
                    refinedCellBeta=autoProcData.get("refinedCellBeta"),
                    refinedCellGamma=autoProcData.get("refinedCellGamma"),
            )
        if autoProcId is None:
            logger.debug("Couldn't create entry for AutoProc in ISPyB. Stopping here.")
            return outData
        outData["autoProcId"] = autoProcId
        autoProcScalingData = inData.get("autoProcScaling", self.getAutoProcScalingContainer())
        autoProcScalingId = UtilsIspyb.storeOrUpdateAutoProcScaling(
            client=client,
            autoProcScalingId=autoProcScalingData.get("autoProcScalingId"),
            autoProcId=autoProcId,
            resolutionEllipsoidAxis11=autoProcScalingData.get("resolutionEllipsoidAxis11"),
            resolutionEllipsoidAxis12=autoProcScalingData.get("resolutionEllipsoidAxis12"),
            resolutionEllipsoidAxis13=autoProcScalingData.get("resolutionEllipsoidAxis13"),
            resolutionEllipsoidAxis21=autoProcScalingData.get("resolutionEllipsoidAxis21"),
            resolutionEllipsoidAxis22=autoProcScalingData.get("resolutionEllipsoidAxis22"),
            resolutionEllipsoidAxis23=autoProcScalingData.get("resolutionEllipsoidAxis23"),
            resolutionEllipsoidAxis31=autoProcScalingData.get("resolutionEllipsoidAxis31"),
            resolutionEllipsoidAxis32=autoProcScalingData.get("resolutionEllipsoidAxis32"),
            resolutionEllipsoidAxis33=autoProcScalingData.get("resolutionEllipsoidAxis33"),
            resolutionEllipsoidValue1=autoProcScalingData.get("resolutionEllipsoidValue1"),
            resolutionEllipsoidValue2=autoProcScalingData.get("resolutionEllipsoidValue2"),
            resolutionEllipsoidValue3=autoProcScalingData.get("resolutionEllipsoidValue3"),
            )
        if autoProcScalingId is None:
            logger.error("Couldn't create entry for AutoProcScaling in ISPyB!")
            self.setFailure()
            return outData
        outData["autoProcScalingId"] = autoProcScalingId
        #autoProcScalingHasIntData = inData.get("autoProcScalingHasInt", self.getAutoProcScalingHasIntContainer())
        autoProcScalingHasIntId = UtilsIspyb.storeOrUpdateAutoProcScalingHasInt(
            autoProcScalingHasIntId=None,
            autoProcIntegrationId=autoProcIntegrationId,
            autoProcScalingId=autoProcScalingId,
            client=client,
            )
        if autoProcScalingHasIntId is None:
            logger.error("Couldn't create entry for AutoProcScalingHasInt in ISPyB!")
            self.setFailure()
            return outData
        
        autoProcScalingStatisticsDataList = inData.get("autoProcScalingStatistics", None)
        if autoProcScalingStatisticsDataList is None:
            logger.debug("Couldn't create entry for autoProcScalingStatisticsData in ISPyB. Stopping here.")
            return outData

        for autoProcScalingStatisticsData in autoProcScalingStatisticsDataList:
            autoProcScalingStatisticsId = UtilsIspyb.storeOrUpdateAutoProcScalingStatistics(
                client=client,
                autoProcScalingStatisticsId=autoProcScalingStatisticsData.get("autoProcScalingStatisticsId"),
                autoProcScalingId=autoProcScalingId,
                scalingStatisticsType=autoProcScalingStatisticsData.get("scalingStatisticsType"),
                comments=autoProcScalingStatisticsData.get("comments"),
                resolutionLimitLow=autoProcScalingStatisticsData.get("resolutionLimitLow"),
                resolutionLimitHigh=autoProcScalingStatisticsData.get("resolutionLimitHigh"),
                rmerge=autoProcScalingStatisticsData.get("rmerge"),
                rmeasWithinIplusIminus=autoProcScalingStatisticsData.get("rmeasWithinIplusIminus"),
                rmeasAllIplusIminus=autoProcScalingStatisticsData.get("rmeasAllIplusIminus"),
                rpimWithinIplusIminus=autoProcScalingStatisticsData.get("rpimWithinIplusIminus"),
                rpimAllIplusIminus=autoProcScalingStatisticsData.get("rpimAllIplusIminus"),
                fractionalPartialBias=autoProcScalingStatisticsData.get("fractionalPartialBias"),
                nTotalObservations=autoProcScalingStatisticsData.get("nTotalObservations"),
                nTotalUniqueObservations=autoProcScalingStatisticsData.get("nTotalUniqueObservations"),
                meanIoverSigI=autoProcScalingStatisticsData.get("meanIoverSigI"),
                completeness=autoProcScalingStatisticsData.get("completeness"),
                multiplicity=autoProcScalingStatisticsData.get("multiplicity"),
                anomalousCompleteness=autoProcScalingStatisticsData.get("anomalousCompleteness"),
                anomalousMultiplicity=autoProcScalingStatisticsData.get("anomalousMultiplicity"),
                anomalous=autoProcScalingStatisticsData.get("anomalous"),
                ccHalf=autoProcScalingStatisticsData.get("ccHalf"),
                ccAno=autoProcScalingStatisticsData.get("ccAno"),
                sigAno=autoProcScalingStatisticsData.get("sigAno"),
                isa=autoProcScalingStatisticsData.get("isa"),
                completenessSpherical=autoProcScalingStatisticsData.get("completenessSpherical"),
                anomalousCompletenessSpherical=autoProcScalingStatisticsData.get("anomalousCompletenessSpherical"),
                completenessEllipsoidal=autoProcScalingStatisticsData.get("completenessEllipsoidal"),
                anomalousCompletenessEllipsoidal=autoProcScalingStatisticsData.get("anomalousCompletenessEllipsoidal"),
            )
            if autoProcScalingStatisticsId is None:
                logger.error("Couldn't create entry for autoProcScalingStatistics in ISPyB!")
                self.setFailure()
                return outData
        return outData
    
    @staticmethod
    def setIspybToRunning(dataCollectionId=None, processingCommandLine=None, processingPrograms=None, isAnom=False, timeStart=None):
        inputStoreAutoProcAnom = {
            "dataCollectionId": dataCollectionId,
            "autoProcProgram":  {
                "processingCommandLine": processingCommandLine,
                "processingPrograms": processingPrograms,
                "processingStatus": "RUNNING",
                "processingStartTime": timeStart,
                },
            "autoProcIntegration" : {
                "anomalous":isAnom,
            }
        }
        autoProcStoreIspybResults = ISPyBStoreAutoProcResults(inData=inputStoreAutoProcAnom, workingDirectorySuffix="setRunning")
        autoProcStoreIspybResults.execute()

        return autoProcStoreIspybResults.outData["autoProcIntegrationId"], autoProcStoreIspybResults.outData["autoProcProgramId"]

    @staticmethod
    def setIspybToFailed(dataCollectionId=None, autoProcProgramId=None, autoProcIntegrationId=None, processingCommandLine=None, processingPrograms=None, isAnom=False, timeStart=None, timeEnd=None):
        inputStoreAutoProcAnom = {
            "dataCollectionId": dataCollectionId,
            "autoProcProgram":  {
                "autoProcProgramId": autoProcProgramId,
                "processingCommandLine": processingCommandLine,
                "processingPrograms": processingPrograms,
                "processingStatus": "FAILED",
                "processingStartTime": timeStart,
                "processingEndTime": timeEnd,
                },
            "autoProcIntegration" : {
                "anomalous":isAnom,
                "autoProcIntegrationId" : autoProcIntegrationId,
            }
        }
        autoProcStoreIspybResults = ISPyBStoreAutoProcResults(inData=inputStoreAutoProcAnom, workingDirectorySuffix="setFailed")
        autoProcStoreIspybResults.execute()

        return autoProcStoreIspybResults.outData["autoProcIntegrationId"], autoProcStoreIspybResults.outData["autoProcProgramId"]

class ISPyBStoreAutoProcStatus(AbstractTask):
    def getOutDataSchema(self):
        return {
            "autoProcIntegrationId": {"type": ["integer","null"]},
            "autoProcProgramId": {"type": ["integer","null"]},
            "autoProcStatusId": {"type": ["integer","null"]},
            }

    def run(self,inData):
        outData = {
            "autoProcIntegrationId": None,
            "autoProcProgramId": None,
            "autoProcStatusId": None,
            }
        dataCollectionId = inData.get("dataCollectionId")

        autoProcIntegration = inData.get("autoProcIntegration")
        autoProcIntegrationId = autoProcIntegration.get("autoProcIntegrationId")
        autoProcProgram = inData.get("autoProcProgram")
        autoProcProgramId = autoProcProgram.get("autoProcProgramId")
        autoProcStatus = inData.get("autoProcStatus")
        autoProcStatusId = inData.get("autoProcStatusId")

        if (autoProcIntegrationId is None) and (dataCollectionId is None):
            logger.error("Either data collection id or auto proc integration id must be given as input!")
            self.setFailure()
            return outData

        client = UtilsIspyb.getAutoprocessingWebService()
        if client is None:
            logger.error("Cannot connect to ISPyB web service")
            self.setFailure()
            return outData
        
        if autoProcIntegrationId is None:
            if autoProcProgram is not None:
                autoProcProgramId = UtilsIspyb.storeOrUpdateAutoProcProgram(
                    autoProcProgramId=autoProcProgramId,
                    processingCommandLine=autoProcProgram.get("processingCommandLine"),
                    processingPrograms=autoProcProgram.get("processingPrograms"),
                    processingStatus=autoProcProgram.get("processingStatus"),
                    processingStartTime=autoProcProgram.get("processingStartTime"),
                    processingEndTime=autoProcProgram.get("processingEndTime"),
                    processingEnvironment=autoProcProgram.get("processingEnvironment"),
                    client=client
                    )
                
            else:
                autoProcProgramId=None
            logger.debug(f"autoProcProgramId: {autoProcProgramId}")
            # If no autoProcessingId is given create a dummy entry in the integration table
            autoProcIntegrationId = UtilsIspyb.storeOrUpdateAutoProcIntegration(
                client=client,
                dataCollectionId=dataCollectionId,
                autoProcProgramId=autoProcProgramId,
                anomalous=autoProcIntegration.get("anomalous", False)
            )
        autoProcStatusId = UtilsIspyb.storeOrUpdateAutoProcStatus(
            client=client,
            autoProcStatusId=autoProcStatusId,
            autoProcIntegrationId=autoProcIntegrationId,
            step=autoProcStatus.get("step"),
            status=autoProcStatus.get("status"),
            comments=autoProcStatus.get("comments"),
            bltimeStamp=datetime.now()
        )
        
        outData = {
            "autoProcIntegrationId": autoProcIntegrationId,
            "autoProcProgramId": autoProcProgramId,
            "autoProcStatusId": autoProcStatusId,
            }
        return outData

def createIntegrationId(task, comments, isAnom=False):
    """
    gets integrationID and programID, 
    sets processing status to RUNNING.
    """
    statusInput = {
        "dataCollectionId": task.dataCollectionId,
        "autoProcIntegration" : {
            "anomalous": isAnom,
        },
        "autoProcProgram": {
            "processingCommandLine": task.processingCommandLine,
            "processingPrograms": task.processingPrograms,
            "processingStatus": "RUNNING",
            "processingStartTime": task.startDateTime,
        },
        "autoProcStatus": {
            "step":  "Indexing",
            "status": "Launched",
            "comments": comments,
            "bltimeStamp": datetime.now().isoformat(timespec='seconds'),
        }
    }
    autoprocStatus = ISPyBStoreAutoProcStatus(inData=statusInput,workingDirectorySuffix="createIntegrationId")

    # get our EDNAproc status id
    autoprocStatus.execute()
    return (autoprocStatus.outData["autoProcIntegrationId"],
            autoprocStatus.outData["autoProcProgramId"])
