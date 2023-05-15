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

import datetime
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
        with open(xml_path) as f:
            xml_string = f.read()
        auto_proc_dict = xmltodict.parse(xml_string)
        # pprint.pprint(auto_proc_dict)
        auto_proc_container = auto_proc_dict["AutoProcContainer"]
        auto_proc_program_container = auto_proc_container["AutoProcProgramContainer"]
        # 1. Create AutoProcProgram entry
        auto_proc_program = auto_proc_program_container["AutoProcProgram"]
        auto_proc_program_id = UtilsIspyb.storeOrUpdateAutoProcProgram(
            programs=self.check_length(auto_proc_program["processingPrograms"]),
            commandline=self.check_length(auto_proc_program["processingCommandLine"]),
            start_time=self.get_time(auto_proc_program["processingStartTime"]),
            end_time=self.get_time(auto_proc_program["processingEndTime"]),
            environment=self.check_length(auto_proc_program["processingEnvironment"]),
            message=self.check_length(auto_proc_program["processingMessage"]),
            status=auto_proc_program["processingStatus"]
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