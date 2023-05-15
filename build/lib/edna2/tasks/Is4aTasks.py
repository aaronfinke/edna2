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
__date__ = "10/05/2019"

import json
import shutil

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.ISPyBTasks import GetListAutoprocessingResults


class FindHklAsciiForMerge(AbstractTask):
    """
    This task receives a list of data collection IDs and returns a
    json schema for EXI2
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
        inDataGetListAutoprocessingResults = {
            'token': token,
            'proposal': proposal,
            'dataCollectionId': listDataCollectionId
        }
        getListAutoprocessingResults = GetListAutoprocessingResults(
            inData=inDataGetListAutoprocessingResults
        )
        getListAutoprocessingResults.execute()
        outDataAutoprocessing = getListAutoprocessingResults.outData
        if 'error' in outDataAutoprocessing:
            urlError = outDataAutoprocessing['error']
        else:
            index = 1
            properties = {}
            listOrder = []
            for dataCollection in outDataAutoprocessing['dataCollection']:
                dataCollectionId = dataCollection['dataCollectionId']
                dictEntry = {}
                listEnumNames = []
                listEnumValues = []
                proteinAcronym = None
                blSampleName = None
                if 'error' in dataCollection['autoprocIntegration']:
                    urlError = dataCollection['autoprocIntegration']['error']
                else:
                    for autoProcResult in dataCollection['autoprocIntegration']:
                        if proteinAcronym is None:
                            proteinAcronym = autoProcResult['Protein_acronym']
                            blSampleName = autoProcResult['BLSample_name']
                        for autoProcAttachment in autoProcResult['autoprocAttachment']:
                            if 'XDS_ASCII' in autoProcAttachment['fileName']:
                                fileName = autoProcAttachment['fileName']
                                program = autoProcResult['v_datacollection_processingPrograms']
                                attachmentId = autoProcAttachment['autoProcProgramAttachmentId']
                                enumName = '{0:30s} {1}'.format(program, fileName)
                                listEnumNames.append(enumName)
                                enumValue = attachmentId
                                listEnumValues.append(enumValue)
                if urlError is None:
                    entryKey = 'hkl_' + str(dataCollectionId)
                    if entryKey not in properties:
                        dictEntry['title'] = 'Select HKL for data Collection #{0} {2} {1}-{2}'.format(
                            index,
                            proteinAcronym,
                            blSampleName
                        )
                        dictEntry['enum'] = listEnumValues
                        dictEntry['enumNames'] = listEnumNames
                        properties[entryKey] = dictEntry
                        listOrder.append(entryKey)
                    entryKey = 'minimum_I/SIGMA_' + str(dataCollectionId)
                    if entryKey not in properties:
                        # Minimum sigma
                        dictEntry = {
                            'integer': 'string',
                            'type': 'string',
                            'title': 'minimum_I/SIGMA for data Collection #{0} {2} {1}-{2}'.format(
                                index,
                                proteinAcronym,
                                blSampleName
                            )
                        }
                        properties[entryKey] = dictEntry
                        listOrder.append(entryKey)
                        index += 1
        if urlError is None:
            schema = {
                'properties': properties,
                'type': 'object',
                'title': 'User input needed'
            }
            uiSchema = {
                'ui:order': listOrder
            }
            outData = {
                "schema": schema,
                "uiSchema": uiSchema
            }
        else:
            outData = {
                'error': urlError
            }
        return outData


class FindPipelineForMerge(AbstractTask):
    """
    This task receives a list of data collection IDs and returns a
    json schema for EXI2
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
        inDataGetListAutoprocessingResults = {
            'token': token,
            'proposal': proposal,
            'dataCollectionId': listDataCollectionId
        }
        getListAutoprocessingResults = GetListAutoprocessingResults(
            inData=inDataGetListAutoprocessingResults
        )
        getListAutoprocessingResults.execute()
        outDataAutoprocessing = getListAutoprocessingResults.outData
        if 'error' in outDataAutoprocessing:
            urlError = outDataAutoprocessing['error']
        else:
            index = 1
            properties = {}
            listOrder = []
            dictEntry = {}
            for dataCollection in outDataAutoprocessing['dataCollection']:
                dataCollectionId = dataCollection['dataCollectionId']
                listEnumValues = []
                proteinAcronym = None
                blSampleName = None
                if 'error' in dataCollection['autoprocIntegration']:
                    urlError = dataCollection['autoprocIntegration']['error']
                else:
                    for autoProcResult in dataCollection['autoprocIntegration']:
                        if len(autoProcResult['autoprocAttachment']) > 0:
                            if proteinAcronym is None:
                                proteinAcronym = autoProcResult['Protein_acronym']
                                blSampleName = autoProcResult['BLSample_name']
                            if '1' in autoProcResult['anomalous']:
                                anom = True
                            else:
                                anom = False
                            for autoProcAttachment in autoProcResult['autoprocAttachment']:
                                if 'XDS_ASCII' in autoProcAttachment['fileName']:
                                    fileName = autoProcAttachment['fileName']
                                    program = autoProcResult['v_datacollection_processingPrograms']
                                    attachmentId = autoProcAttachment['autoProcProgramAttachmentId']
                                    if anom:
                                        entryKey = program + '_anom'
                                    else:
                                        entryKey = program + '_noanom'
                                    if entryKey not in dictEntry:
                                        dictEntry[entryKey] = []
                                    dictEntry[entryKey].append({'id': attachmentId, 'fileName': fileName})
        if urlError is None:
            listEnumNames = []
            dictInput = {}
            for entryKey, listAttachment in dictEntry.items():
                if len(listAttachment) == len(outDataAutoprocessing['dataCollection']):
                    listEnumNames.append(entryKey)
                    dictInput[entryKey] = listAttachment
                    index += 1
            if len(listEnumNames) > 0:
                dictSchema = {
                    'title': 'Select processing pipeline for data Collection {0}-{1}'.format(
                        proteinAcronym,
                        blSampleName
                    ),
                    'type': 'string',
                    'enum': listEnumNames,
                    'enumNames': listEnumNames
                }
                key = "pipeline"
                properties[key] = dictSchema
                listOrder.append(key)
                # Minimum sigma
                dictSchema = {
                    'integer': 'string',
                    'type': 'string',
                    'title': 'minimum_I/SIGMA for data Collection {0}-{1}'.format(
                        proteinAcronym,
                        blSampleName
                    )
                }
                key = 'minimum_I/SIGMA'
                properties[key] = dictSchema
                listOrder.append(key)
        if urlError is None:
            schema = {
                'properties': properties,
                'type': 'object',
                'title': 'User input needed'
            }
            uiSchema = {
                'ui:order': listOrder
            }
            outData = {
                'schema': {
                    "schema": schema,
                    "uiSchema": uiSchema
                },
                'input': dictInput
            }
        else:
            outData = {
                'error': urlError
            }
        return outData


class MergeUtls(AbstractTask):
    """
    This task will run the Merge_utls.py program written by Shibom Basu
    """

    def run(self, inData):
        listHklLp = inData['listHklLp']
        workingDir = self.getWorkingDirectory()
        index = 1
        for hklLp in listHklLp:
            dataDir = workingDir / "data{0}".format(index)
            dataDir.mkdir(exist_ok=True)
            shutil.copy(hklLp['hkl'], str(dataDir / 'XDS_ASCII.HKL'))
            index += 1
        commandLine = 'Merge_utls.py --root {0} --expt serial-xtal'.format(str(workingDir))
        self.runCommandLine(commandLine, logPath=None)
        # Find Mergeing_results.json
        resultPath = self.getWorkingDirectory() / 'adm_serial-xtal' / 'adm_3' / 'Mergeing_results.json'
        if resultPath.exists():
            with open(str(resultPath)) as f:
                mergeResult = json.loads(f.read())
        outData = {'mergeResult': mergeResult}
        return outData
