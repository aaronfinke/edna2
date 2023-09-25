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

import os
import json
from builtins import RuntimeError

import numpy
import shlex
# import distro
import shutil
import base64
import pathlib
import tempfile

import matplotlib
import matplotlib.pyplot as plt

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.H5ToCBFTask import H5ToCBFTask
from edna2.tasks.ReadImageHeader import ReadImageHeader
from edna2.tasks.ISPyBTasks import ISPyBRetrieveDataCollection

from edna2.utils import UtilsPath
from edna2.utils import UtilsImage
from edna2.utils import UtilsIspyb
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsDetector

# Simplification of ControlDozor.py that instead runs a middle Python layer that
# relies on Dozor python bindings.

logger = UtilsLogging.getLogger()

# Default parameters

MAX_BATCH_SIZE = 5000


class ExecPyDozor(AbstractTask):  # pylint: disable=too-many-instance-attributes
    """
    The ExecPyDozor is responsible for executing the python code relying on
    'dozor' wrappers. Much of the meta data collection is done at the python
    level since the use case relies upon available Eiger master hdf5 files.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "required": ["masterFile", "firstImageNumber", "lastImageNumber"],
            "properties": {
                "masterFile": {"type": "string"},
                "maskFile": {"type": "string"},                
                "firstImageNumber": {"type": "integer"},
                "lastImageNumber": {"type": "integer"},
                "outputDirectory": {"type": "string"},
                "numOfProcesses": {"type": "integer"},
                "doSubmit": {"type" : "boolean"},
                "dozorCutoff": {"type": "integer"}
            }
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "required": ["imageDozor"],
            "properties": {
                "imageDozor": {
                    "type": "array",
                    "items": {
                        "$ref": self.getSchemaUrl("imageDozor.json")
                    }
                },
                "halfDoseTime": {"type": "number"},
                "dozorPlot":  {"type": "string"},
                "plotmtvFile":  {"type": "string"},
                "pngPlots":  {
                    "type": "array",
                    "items": {"type": "string"}
                },
                "dozorAllFile": {"type": "string"},
            },
        }

    def run(self, inData):
        doSubmit = inData.get('doSubmit', False)
        # Create PyDozor command line
        path = UtilsConfig.get(self,'slurm_path','dozor')
        noProcesses = UtilsConfig.get(self, 'core', 10)
        module_imports = UtilsConfig.get(self,'module_import',None)        
        cutoff = inData.get("dozorCutoff",5)

        executable = module_imports + "\n\n"
        executable += UtilsConfig.get(self, 'slurm_executable', 'dozor')
        executable += " -s {} ".format(inData['firstImageNumber'])
        executable += " -e {} ".format(inData['lastImageNumber'])                                
        executable += " -m {} ".format(inData['masterFile'])                                
        executable += " -o {} ".format(inData['outputDirectory'])
        executable += " -n {} ".format(noProcesses)
        executable += " -c {} ".format(cutoff)        
        #Don't create individual spot lists to avoid spamming file system.
        executable += " --skip_spots" 

        if doSubmit:
            partition = UtilsConfig.get(self,'slurm_partition',None)
        else:
            partition = None
        commandLine = executable 
        if inData.get('maskFile',"") != "":
            commandLine += " -M {} ".format(inData['maskFile'])
        self.setLogFileName('pydozor.log')
        print("runCommandLine = \n{}\n".format(commandLine))
        self.runCommandLine(commandLine, doSubmit=doSubmit, partition=partition)
        log = self.getLog()
        with open(inData['outputDirectory']+"/dozor.log",'r') as dozorOut:
            outData = self.parseOutput(inData, dozorOut, workingDir=self.getWorkingDirectory())
        return outData



    def parseOutput(self, inData, output, workingDir=None):
        """
        This method parses the output of dozor
        """
        resultDozor = {
            'imageDozor': []  # list of dict. each dict contains spotFile and Image_path
        }
        masterDict = ReadImageHeader.readHdf5Header(inData['masterFile'])
        angle_start = masterDict['omega_start']
        angle_step = round(masterDict['omega_range_average'],3)
        counter = 0

        for line in output:
            # Remove '|'
            listLine = shlex.split(line)
            if len(listLine) > 0:
                imageDozor = {}
                try: 
                    imageDozor['angle'] = angle_start + counter*angle_step
                    imageDozor['image'] = str(listLine[0])
                    imageDozor['number'] = int(listLine[1])
                    imageDozor['spotsNumOf'] = int(listLine[2])
                    imageDozor['spotsIntAver'] = 0 # TODO                
                    imageDozor['spotScore'] = float(listLine[3])
                    imageDozor['mainScore'] = float(listLine[3]) # TODO Difference compared to dozorSpotScore?
                    imageDozor['visibleResolution'] = float(listLine[4])
                    imageDozor['spotsResolution'] = float(listLine[4]) # TODO Difference compared to visibleResolution?
                except Exception as e:
                    logger.warning('Exception caught when parsing Dozor log!')
                    logger.warning(e)
                """ 
                Spot files might be processed at this point, but for many images, a more intelligent storage strategy
                is needed, in particular when we move to SSX datasets of the 10k order of magnitude. The dozor_offline
                code should create a h5 file containing all spot data in a single file.
                 
                if workingDir is not None:
                    spotFile = os.path.join(str(workingDir),
                                            '%05d.spot' % imageDozor['number'])
                    if os.path.exists(spotFile):
                        imageDozor['spotFile'] = spotFile
                """
                resultDozor['imageDozor'].append(imageDozor)
                counter += 1
        return resultDozor


class ControlPyDozor(AbstractTask):

    def __init__(self, inData, workingDirectorySuffix=None):
        AbstractTask.__init__(self, inData, workingDirectorySuffix=workingDirectorySuffix)
        self.directory = inData.get('masterFile',"")
        self.hasOverlap = False
        self.overlap = 0.0
        self.dataCollection = None
    

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "dataCollectionId": {"type": "integer"},
                "processDirectory": {"type": "string"},
                "dozorOutputDirectory": {"type": "string"},                
                "masterFile": {"type": "string"},
                "maskFile": {"type": "string"},                
                "directory": {"type": "string"},
                "beamline": {"type": "string"},
                "template": {"type": "string"},
                "startNo": {"type": "integer"},
                "endNo": {"type": "integer"},
                "batchSize": {"type": "integer"},
                "doISPyBUpload": {"type": "boolean"},
                "doDozorm": {"type": "boolean"},
                "doSubmit": {"type": "boolean"},
                "dozorCutoff": {"type": "integer"},
                "returnSpotList": {"type": "boolean"}
            }
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "required": ["imageQualityIndicators", "detectorType"],
            "properties": {
                "imageQualityIndicators": {
                    "type": "array",
                    "items": {
                        "$ref": self.getSchemaUrl("imageQualityIndicators.json")
                    }
                },
                "detectorType": {"type": "string"},
                "halfDoseTime": {"type": "number"},
                "dozorPlot":  {"type": "string"},
                "pathToCbfDirectory":  {"type": "string"},
                "pngPlots":  {"type": "string"},
                "dozorAllFile": {"type": "string"},
            },
        }

    def run(self, inData):
        # workingDirectory = self.getWorkingDirectory()
        outData = {}
        controlDozorAllFile = None
        listDozorAllFile = []
        hasHdf5Prefix = False
        detectorType = "Eiger16M_test"
        returnSpotList = inData.get('returnSpotList', False)
        batchSize = None
        # Check doDozorm
        doDozorm = inData.get('doDozorm', False)
        # Collect metadata from ISPyB if possible
        self.collectMetaData(inData)
        batchSize = self.determineBatchsize(inData)
        if not 'masterFile' in inData:
            inData['masterFile'] = self.determineMasterfile(inData)
            self.directory=inData['masterFile']
        if not 'maskFile' in inData:
            inData['maskFile'] = self.determineMaskfile(inData)
        if not 'endNo' in inData:
            inData['endNo'] = -1 #Negative means full range.
        # Check overlap
        overlap = inData.get('overlap', self.overlap)
        if overlap != 0:
            self.hasOverlap = True
        logger.debug("ExecPyDozor batch size: {0}".format(batchSize))
        outData['imageQualityIndicators'] = []
        #for listBatch in listAllBatches:
        outDataDozor = self.runPyDozorTask(
            inData=inData,
            batchSize=batchSize,
            overlap=overlap,
            workingDirectory=str(self.getWorkingDirectory()),
            hasOverlap=self.hasOverlap
        )
        #############################################################
        #Output will be much smaller for the dozor_offline script.
        ###########################################################
        if outDataDozor is not None:
            for imageDozor in outDataDozor['imageDozor']:
                imageQualityIndicators = {
                    'angle': imageDozor['angle'],
                    'number': imageDozor['number'],
                    'image': imageDozor['image'],
                    'dozorScore': imageDozor['spotScore'],
                    'dozorSpotScore': imageDozor['mainScore'],
                    'dozorSpotsNumOf': imageDozor['spotsNumOf'],
                    'dozorVisibleResolution': imageDozor['visibleResolution'],
                    'dozorSpotsResolution': imageDozor['spotsResolution']
                }
                """ 
                Spotfiles might be needed later in some form. Keeping this codeblock for now.

                if 'spotFile' in imageDozor:
                    if os.path.exists(imageDozor['spotFile']):
                        spotFile = imageDozor['spotFile']
                        imageQualityIndicators['dozorSpotFile'] = spotFile
                        if returnSpotList:
                            numpyArray = numpy.loadtxt(spotFile, skiprows=3)
                            imageQualityIndicators['dozorSpotList'] = \
                                base64.b64encode(numpyArray.tostring()).decode('utf-8')
                            imageQualityIndicators['dozorSpotListShape'] = \
                                list(numpyArray.shape)
                """
                outData['imageQualityIndicators'].append(imageQualityIndicators)
            if doDozorm:
                listDozorAllFile.append(outDataDozor['dozorAllFile'])
        # Assemble all dozorAllFiles into one
        if doDozorm:
            controlDozorAllFile = str(self.getWorkingDirectory() / "dozor_all")
            os.system('touch {0}'.format(controlDozorAllFile))
            for dozorAllFile in listDozorAllFile:
                command = 'cat ' + dozorAllFile + ' >> ' + controlDozorAllFile
                os.system(command)
        # Make plot if we have a data collection id
        if 'dataCollectionId' in inData:
            if "processDirectory" in inData:
                processDirectory = pathlib.Path(inData["processDirectory"])
            else:
                processDirectory = self.getWorkingDirectory()
            dozorPlotPath, dozorCsvPath = self.makePlot(inData['dataCollectionId'], outData, self.getWorkingDirectory())
            doIspybUpload = inData.get("doISPyBUpload", False)
            if doIspybUpload:
                self.storeDataOnPyarch(inData["dataCollectionId"],
                                       dozorPlotPath=dozorPlotPath,
                                       dozorCsvPath=dozorCsvPath,
                                       workingDirectory=processDirectory)
        # Read the header from the first image in the batch
        outData['detectorType'] = detectorType
        if doDozorm:
            outData['dozorAllFile'] = controlDozorAllFile
        return outData

    def collectMetaData(self, inData):
        if 'dataCollectionId' in inData:
            if inData['dataCollectionId']!=0:
                ispybInData = {
                    'dataCollectionId': inData['dataCollectionId']
                }
                ispybTask = ISPyBRetrieveDataCollection(inData=ispybInData)
                self.dataCollection = ispybTask.run(ispybInData)

    def determineBatchsize(self, inData):
        if self.dataCollection != None:
            batchSize = self.dataCollection['numberOfImages']
            if 'overlap' in self.dataCollection and \
                    abs(self.dataCollection['overlap']) > 1:
                self.hasOverlap = True
                self.overlap = self.dataCollection['overlap']
            else:
                self.overlap = 0.0
        else:
            # No connection to ISPyB, take parameters from input
            if 'batchSize' in inData:
                batchSize = inData['batchSize']
            else:
                batchSize = UtilsConfig.get('ControlPyDozor', 'batchSize')
      
        if batchSize > MAX_BATCH_SIZE or batchSize == None:
            batchSize = MAX_BATCH_SIZE
        return batchSize

    def determineMasterfile(self, inData):
        masterFile=""
        if self.dataCollection != None:
            imgDir = self.dataCollection['imageDirectory']
            if imgDir[-1] != '/': 
                imgDir = imgDir + '/'
            dataCollNo = self.dataCollection['dataCollectionNumber']
            imgPrefix = self.dataCollection['imagePrefix']
            imgSuffix = self.dataCollection['imageSuffix']
            masterFile = imgDir + imgPrefix + "_" + str(dataCollNo) + "_" + "master." + imgSuffix
        return masterFile
 
    def determineMaskfile(self, inData):
        return UtilsConfig.get('ControlPyDozor','mask_file') 

    @classmethod
    def runPyDozorTask(cls, inData, batchSize, overlap, workingDirectory,
                     hasOverlap):
        doSubmit = inData.get('doSubmit', False)
        doDozorm = inData.get('doDozorm', False)
        cutoff = inData.get('dozorCutoff',5)
        outputDozorDir = inData.get('dozorOutputDirectory',"output")
        outDataDozor = None
        inDataPyDozor = {'masterFile': inData['masterFile'],
                'maskFile': inData['maskFile'],
                'firstImageNumber': inData['startNo'],
                'lastImageNumber': inData['endNo'],
                'outputDirectory': outputDozorDir,
                'numOfProcesses': int(UtilsConfig.get('ExecPyDozor','core',10)),
                'doSubmit': doSubmit,
                'doDozorm': doDozorm,
                'dozorCutoff': cutoff
                }
        #Make sure dozor has existing directories to write to.
        outDirectory = pathlib.Path(workingDirectory) / outputDozorDir
        print("DEBUG: outDirectory = {}".format(outDirectory)) #ALEK DEBUG
        try:
            if not outDirectory.exists():
                outDirectory.mkdir(parents=True)
        except Exception as e:
            logger.warning("Couldn't create dozor output dirs: {0}".format(outDirectory))

        dozor = ExecPyDozor(inData=inDataPyDozor)
        dozor.execute()
        if not dozor.isFailure():
            outDataDozor = dozor.outData
        return outDataDozor

    def makePlot(self, dataCollectionId, outDataImageDozor, workingDirectory):
        minImageNumber = None
        maxImageNumber = None
        minAngle = None
        maxAngle = None
        minDozorValue = None
        maxDozorValue = None
        minResolution = None
        maxResolution = None
        plotFileName = 'dozor_{0}.png'.format(dataCollectionId)
        csvFileName = 'dozor_{0}.csv'.format(dataCollectionId)
        with open(str(workingDirectory / csvFileName), 'w') as gnuplotFile:
            gnuplotFile.write(
                '# Data directory: {0}\n'.format(self.directory)
            )
            gnuplotFile.write(
                '# {0:>9s}{1:>16s}{2:>16s}{3:>16s}{4:>16s}{5:>16s}\n'.format(
                    "'Image no'",
                    "'Angle'",
                    "'No of spots'",
                    "'Main score (*10)'",
                    "'Spot score'",
                    "'Visible res.'",
                )
            )
            for imageQualityIndicators in outDataImageDozor['imageQualityIndicators']:
               gnuplotFile.write(
                   '{0:10d},{1:15.3f},{2:15d},{3:15.3f},{4:15.3f},{5:15.3f}\n'.format(
                        imageQualityIndicators['number'],
                        imageQualityIndicators['angle'],
                        imageQualityIndicators['dozorSpotsNumOf'],
                        10 * imageQualityIndicators['dozorScore'],
                        imageQualityIndicators['dozorSpotScore'],
                        imageQualityIndicators['dozorVisibleResolution'],
                    )
               )
               if minImageNumber is None or minImageNumber > imageQualityIndicators['number']:
                   minImageNumber = imageQualityIndicators['number']
                   minAngle = imageQualityIndicators['angle']
               if maxImageNumber is None or maxImageNumber < imageQualityIndicators['number']:
                   maxImageNumber = imageQualityIndicators['number']
                   maxAngle = imageQualityIndicators['angle']
               if minDozorValue is None or minDozorValue > imageQualityIndicators['dozorScore']:
                   minDozorValue = imageQualityIndicators['dozorScore']
               if maxDozorValue is None or maxDozorValue < imageQualityIndicators['dozorScore']:
                   maxDozorValue = imageQualityIndicators['dozorScore']

               # Min resolution: the higher the value the lower the resolution
               if minResolution is None or minResolution < imageQualityIndicators['dozorVisibleResolution']:
                   # Disregard resolution worse than 10.0
                   if imageQualityIndicators['dozorVisibleResolution'] < 10.0:
                       minResolution = imageQualityIndicators['dozorVisibleResolution']

               # Max resolution: the lower the number the better the resolution
               if maxResolution is None or maxResolution > imageQualityIndicators['dozorVisibleResolution']:
                   maxResolution = imageQualityIndicators['dozorVisibleResolution']


        xtics = ''
        if minImageNumber is not None and minImageNumber == maxImageNumber:
           minAngle -= 1.0
           maxAngle += 1.0
        noImages = maxImageNumber - minImageNumber + 1
        if noImages <= 4:
           minImageNumber -= 0.1
           maxImageNumber += 0.1
           deltaAngle = maxAngle - minAngle
           minAngle -= deltaAngle * 0.1 / noImages
           maxAngle += deltaAngle * 0.1 / noImages
           xtics = '1'

        if maxResolution is None or maxResolution > 0.8:
           maxResolution = 0.8
        else:
           maxResolution = int(maxResolution * 10.0) / 10.0

        if minResolution is None or minResolution < 4.5:
           minResolution = 4.5
        else:
           minResolution = int(minResolution * 10.0) / 10.0 + 1

        if maxDozorValue < 0.001 and minDozorValue < 0.001:
           yscale = 'set yrange [-0.5:0.5]\n    set ytics 1'
        else:
           yscale = 'set autoscale  y'

        gnuplotScript = \
"""#
set terminal png
set output '{dozorPlotFileName}'
set title '{title}' offset 0,1 font 'helvetica bold,8' noenhanced
set grid x2 y2
set xlabel "Image number"
set x2label 'Angle (degrees)'
set y2label 'Resolution (A)'
set ylabel 'Number of spots / ExecDozor score (*10)'
set xtics {xtics} nomirror
set x2tics
set ytics nomirror
set y2tics
set xrange [{minImageNumber}:{maxImageNumber}]
set x2range [{minAngle}:{maxAngle}]
{yscale}
set y2range [{minResolution}:{maxResolution}]
set key below font 'helvetica bold,10'
plot '{dozorCsvFileName}' using 1:3 title 'Number of spots' axes x1y1 with points linetype rgb 'goldenrod' pointtype 7 pointsize 0.5, \
    '{dozorCsvFileName}' using 1:4 title 'ExecDozor score' axes x1y1 with points linetype 3 pointtype 7 pointsize 0.5, \
    '{dozorCsvFileName}' using 1:6 title 'Visible resolution' axes x1y2 with points linetype 1 pointtype 7 pointsize 0.5
""".format(title=self.directory, #Change directory to masterfile or something better
          dozorPlotFileName=plotFileName,
          dozorCsvFileName=csvFileName,
          minImageNumber=minImageNumber,
          maxImageNumber=maxImageNumber,
          minAngle=minAngle,
          maxAngle=maxAngle,
          minResolution=minResolution,
          maxResolution=maxResolution,
          xtics=xtics,
          yscale=yscale,
          )
        pathGnuplotScript = str(workingDirectory / 'gnuplot.sh')
        with open(pathGnuplotScript, 'w') as f:
            f.write(gnuplotScript)
        oldCwd = os.getcwd()
        os.chdir(str(workingDirectory))
        gnuplot = UtilsConfig.get(self, 'gnuplot', 'gnuplot')
        os.system('{0} {1}'.format(gnuplot, pathGnuplotScript))
        os.chdir(oldCwd)
        dozorPlotPath = workingDirectory / plotFileName
        dozorCsvPath = workingDirectory / csvFileName
        return dozorPlotPath, dozorCsvPath

    @classmethod
    def storeDataOnPyarch(cls, dataCollectionId,  dozorPlotPath, dozorCsvPath, workingDirectory):
        resultsDirectory = pathlib.Path(workingDirectory) / 'results'
        try:
            if not resultsDirectory.exists():
                resultsDirectory.mkdir(parents=True, mode=0o755)
            dozorPlotResultPath = resultsDirectory / dozorPlotPath.name
            dozorCsvResultPath = resultsDirectory / dozorCsvPath.name
            shutil.copy(dozorPlotPath, dozorPlotResultPath)
            shutil.copy(dozorCsvPath, dozorCsvResultPath)
        except Exception as e:
            logger.warning(
                "Couldn't copy files to results directory: {0}".format(
                    resultsDirectory))
            logger.warning(e)
        try:
            # Create paths on pyarch
            dozorPlotPyarchPath = UtilsPath.createPyarchFilePath(dozorPlotResultPath)
            dozorCsvPyarchPath = UtilsPath.createPyarchFilePath(dozorCsvResultPath)
            if not os.path.exists(os.path.dirname(dozorPlotPyarchPath)):
                os.makedirs(os.path.dirname(dozorPlotPyarchPath), 0o755)
            shutil.copy(dozorPlotResultPath, dozorPlotPyarchPath)
            shutil.copy(dozorCsvResultPath, dozorCsvPyarchPath)
            # Upload to data collection
            dataCollectionId = UtilsIspyb.setImageQualityIndicatorsPlot(
                dataCollectionId, dozorPlotPyarchPath, dozorCsvPyarchPath)
        except Exception as e:
            logger.warning("Couldn't copy files to pyarch.")
            logger.warning(e)
            
