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

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxPluginExec/plugins/EDPluginGroupMOSFLM-v1.0/plugins/EDPluginMOSFLMv10.py
# mxPluginExec/plugins/EDPluginGroupMOSFLM-v1.0/plugins/EDPluginMOSFLMIndexingv10.py

import os
import json

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.ReadImageHeader import ReadImageHeader

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsDnaTables
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


class AbstractMosflmTask(AbstractTask):
    """
    Common base class for all MOSLFM tasks
    """

    def __init__(self, inData):
        AbstractTask.__init__(self, inData)
        self.matrixFileName = None

    def run(self, inData):
        commandLine = 'ipmosflm DNA dnaTables.xml'
        # Translate from generic data model to specific MOSFLM data model
        mosflmInData = self.generateMOSFLMInData(inData)
        # Save MOSFLM inData
        with open(str(self.getWorkingDirectory() / "mosflmInData.json"), "w") as f:
            f.write(json.dumps(mosflmInData, indent=4))
        listCommand = self.generateMOSFLMCommands(mosflmInData, self.getWorkingDirectory())
        self.setLogFileName('mosflm.log')
        self.runCommandLine(commandLine, listCommand=listCommand)
        # Work in progress!
        outData = self.parseMosflmOutput(self.getWorkingDirectory())
        return outData

    @classmethod
    def generateMOSFLMInData(cls, inData):
        if "imagePath" in inData:
            inDataReadImageHeader = {
                "imagePath": inData["imagePath"]
            }
            readImageHeader = ReadImageHeader(inData=inDataReadImageHeader)
            readImageHeader.execute()
            listSubWedges = readImageHeader.outData["subWedge"]
        elif "subWedge" in inData:
            listSubWedges = inData["subWedge"]
        else:
            # Here we assume that we have inData in MOSFLM format...
            return inData
        firstSubWedge = listSubWedges[0]
        experimentalCondition = firstSubWedge["experimentalCondition"]
        detector = experimentalCondition["detector"]
        beam = experimentalCondition["beam"]
        goniostat = experimentalCondition["goniostat"]
        detectorType = detector["type"]
        if detectorType.startswith("pilatus"):
            mosflmDetector = "PILATUS"
        elif detectorType.startswith("eiger"):
            mosflmDetector = "EIGER"
        else:
            raise RuntimeError("Unknown detecor type for MOSFLM: " + detectorType)
        listImage = firstSubWedge["image"]
        firstImage = listImage[0]
        firstPath = firstImage["path"]
        mosflmInData = {
            # The MOSFLM beam position convention is opposite to the XDS / CBF/ Dozor convention
            "beam": {
                "x": detector["beamPositionY"],
                "y": detector["beamPositionX"]
            },
            "detector": {
                "numberPixelX": detector["numberPixelX"],
                "numberPixelY": detector["numberPixelY"],
                "pixelSizeX": detector["pixelSizeX"],
                "pixelSizeY": detector["pixelSizeY"],
                "type": mosflmDetector
            },
            "directory": os.path.dirname(firstPath),
            "distance": detector["distance"],
            "template": UtilsImage.getTemplate(firstPath),
            "wavelength": beam["wavelength"]
        }
        listMosflmImage = []
        for subWedge in listSubWedges:
            goniostat = subWedge["experimentalCondition"]["goniostat"]
            rotationAxisStart = goniostat['rotationAxisStart']
            oscillationWidth = goniostat['oscillationWidth']
            imageIndex = 0
            for image in subWedge["image"]:
                number = UtilsImage.getImageNumber(image["path"])
                rotationAxisStart = rotationAxisStart + oscillationWidth * imageIndex
                rotationAxisEnd = rotationAxisStart + oscillationWidth
                mosflmImage = {
                    "number":               number,
                    "rotationAxisStart":    rotationAxisStart,
                    "rotationAxisEnd":      rotationAxisEnd
                }
                listMosflmImage.append(mosflmImage)
        mosflmInData["image"] = listMosflmImage
        return mosflmInData
        #
        #                                                          "pixelSizeY"],
        #                                                      "image": [
        #     {
        #         "number": 1,
        #         "rotationAxisStart": 54.5,
        #         "rotationAxisEnd": 55.5
        #     },
        #     {
        #         "number": 2,
        #         "rotationAxisStart": 144.5,
        #         "rotationAxisEnd": 145.5
        #     }
        # ],
        # "matrix": {
        #     "matrixA": [
        #         [
        #             -0.00520639,
        #             0.00409075,
        #             0.00598777
        #         ],
        #         [
        #             -0.00922857,
        #             -0.00696139,
        #             0.00233786
        #         ],
        #         [
        #             0.00580022,
        #             -0.00738658,
        #             0.00518602
        #         ]
        #     ],
        #     "missettingsAngles": [
        #         0.0,
        #         0.0,
        #         0.0
        #     ],
        #     "matrixU": [
        #         [
        #             -0.4310086,
        #             0.3741439,
        #             0.8211259
        #         ],
        #         [
        #             -0.7639828,
        #             -0.6355417,
        #             -0.1114313
        #         ],
        #         [
        #             0.4801684,
        #             -0.6753541,
        #             0.5597635
        #         ]
        #     ],
        #     "cell": {
        #         "a": 74.4515,
        #         "b": 84.8724,
        #         "c": 115.476,
        #         "alpha": 70.0227,
        #         "beta": 76.0723,
        #         "gamma": 85.3234
        #     }
        # },
        # "mosaicityEstimation": 0.8500000000000001,
        # "deviationAngular": 1.906892,
        # "refinedDistance": 305.222015,
        # "symmetry": "P1"
        # }
        # outData = {
        #     "xds": xdsIndexingOutData
        # }

    def generateMOSFLMCommands(self, inData, workingDirectory):
        """
        This method creates a list of MOSFLM commands
        """
        detector = inData['detector']
        detectorType = detector['type']
        detectorCommand = detectorType
        # Check if reversephi and omega are configured
        omega = UtilsConfig.get('MOSFLM', 'omega')
        if omega is not None:
            detectorCommand += ' OMEGA {0}'.format(omega)
        reversePhi = UtilsConfig.get('MOSFLM', 'reversePhi')
        if reversePhi is not None:
            detectorCommand += ' REVERSEPHI'
        listCommand = [
            'WAVELENGTH {0}'.format(inData['wavelength']),
            'DISTANCE {0}'.format(inData['distance']),
            'BEAM {0} {1}'.format(inData['beam']['x'], inData['beam']['y']),
            'DIRECTORY {0}'.format(inData['directory']),
            'TEMPLATE {0}'.format(inData['template'])
        ]
        if 'symmetry' in inData:
            listCommand.append('SYMMETRY {0}'.format(inData['symmetry']))
        if 'mosaicity' in inData:
            listCommand.append('MOSAIC {0}'.format(inData['mosaicity']))

        newmatFileName = self.getNewmatFileName()
        newmatMatrix = inData.get('matrix', None)
        if newmatMatrix is not None:
            matrixFileName = self.getMatrixFileName()
            newmatPath = workingDirectory / matrixFileName
            self.writeNewmat(newmatMatrix, newmatPath)
            listCommand.append('MATRIX ' + matrixFileName)

        # Add exclude regions if Pilatus
        if detectorType == 'PILATUS':
            listCommand.append('DETECTOR {0}'.format(detectorCommand))
            if detector['numberPixelX'] == 1475 and \
               detector['numberPixelY'] == 1679:
                # Pilatus 2M
                listCommand.append('LIMITS EXCLUDE    0.0  83.76  288.96   85.14')
                listCommand.append('LIMITS EXCLUDE    0.0 168.73  288.96  170.10')
                listCommand.append('LIMITS EXCLUDE   33.54   0.0   36.63  253.87')
                listCommand.append('LIMITS EXCLUDE   70.00   0.0   73.1   253.87')
                listCommand.append('LIMITS EXCLUDE  106.46   0.0  109.56  253.87')
                listCommand.append('LIMITS EXCLUDE  142.93   0.0  146.02  253.87')
                listCommand.append('LIMITS EXCLUDE  179.39   0.0  182.49  253.87')
                listCommand.append('LIMITS EXCLUDE  215.86   0.0  218.95  253.87')
                listCommand.append('LIMITS EXCLUDE  252.32   0.0  255.42  253.87')
            elif detector['numberPixelX'] == 2463 and \
               detector['numberPixelY'] == 2527:
                # Pilatus 6M
                listCommand.append('LIMITS EXCLUDE  0.0 338.77 434.6 340.24')
                listCommand.append('LIMITS EXCLUDE  0.0 253.80 434.6 255.28')
                listCommand.append('LIMITS EXCLUDE  0.0 168.83 434.6 170.21')
                listCommand.append('LIMITS EXCLUDE  0.0  83.86 434.6  85.24')
                listCommand.append('LIMITS EXCLUDE 398.18  0.0 401.28 423.6')
                listCommand.append('LIMITS EXCLUDE 361.72  0.0 364.81 423.6')
                listCommand.append('LIMITS EXCLUDE 325.25  0.0 328.35 423.6')
                listCommand.append('LIMITS EXCLUDE 288.79  0.0 291.88 423.6')
                listCommand.append('LIMITS EXCLUDE 252.32  0.0 255.42 423.6')
                listCommand.append('LIMITS EXCLUDE 215.86  0.0 218.96 423.6')
                listCommand.append('LIMITS EXCLUDE 179.40  0.0 182.49 423.6')
                listCommand.append('LIMITS EXCLUDE 142.93  0.0 145.86 423.6')
                listCommand.append('LIMITS EXCLUDE 106.47  0.0 109.56 423.6')
                listCommand.append('LIMITS EXCLUDE  70.00  0.0 73.10  423.6')
                listCommand.append('LIMITS EXCLUDE  33.54  0.0 36.64  423.6')

        # Check if raster is configured
        raster = UtilsConfig.get('MOSFLM', 'raster')
        if raster is not None:
            listCommand.append('RASTER {0}'.format(raster))
        # Check if polarization is configured
        polarization = UtilsConfig.get('MOSFLM', 'polarization')
        if polarization is not None:
            listCommand.append('POLARIZATION {0}'.format(polarization))
        return listCommand

    def getNewmatFileName(self):
        newmatFileName = self.__class__.__name__ + '_newmat.mat'
        return newmatFileName

    def setNewmatFileName(self, newmatFileName):
        self.newmatFileName = newmatFileName

    def getMatrixFileName(self):
        if self.matrixFileName is None:
            self.matrixFileName = self.__class__.__name__ + '_matrix.mat'
        return self.matrixFileName

    def setMatrixFileName(self, matrixFileName):
        self.matrixFileName = matrixFileName

    @classmethod
    def splitStringIntoListOfFloats(cls, inputString):
        listFloats = []
        listString = inputString.split()
        for element in listString:
            if element != '':
                listFloats.append(float(element))
        return listFloats

    def getDataMOSFLMMatrix(self, matrixFileName):
        if matrixFileName is None:
            strMatrixFileName = self.getMatrixFileName()
        else:
            strMatrixFileName = matrixFileName
        xsDataMOSFLMNewmatMatrix = self.getNewmat(strMatrixFileName)
        return xsDataMOSFLMNewmatMatrix

    @classmethod
    def getNewmat(cls, newMatFilePath):
        listOfListOfFloat = []
        with open(str(newMatFilePath)) as f:
            listLine = f.readlines()
        # Convert into list of lists of float
        for line in listLine:
            if not line.startswith('SYMM'):
                listOfListOfFloat.append(
                    AbstractMosflmTask.splitStringIntoListOfFloats(line))
        # Fill in the data
        newmat = {}

        matrixA = [
            listOfListOfFloat[0],
            listOfListOfFloat[1],
            listOfListOfFloat[2],
        ]
        newmat['matrixA'] = matrixA
        newmat['missettingsAngles'] = listOfListOfFloat[3]
        newmat['matrixU'] = [
            listOfListOfFloat[4],
            listOfListOfFloat[5],
            listOfListOfFloat[6],
        ]
        cell = {
            'a': listOfListOfFloat[7][0],
            'b': listOfListOfFloat[7][1],
            'c': listOfListOfFloat[7][2],
            'alpha': listOfListOfFloat[7][3],
            'beta': listOfListOfFloat[7][4],
            'gamma': listOfListOfFloat[7][5]
        }
        newmat['cell'] = cell
        return newmat

    @classmethod
    def writeNewmat(cls, newmat, newmatPath):
        matrixA = newmat['matrixA']
        strNewmat = ''
        strNewmat += ' {0:11.8f} {1:11.8f} {2:11.8f}\n'.format(*matrixA[0])
        strNewmat += ' {0:11.8f} {1:11.8f} {2:11.8f}\n'.format(*matrixA[1])
        strNewmat += ' {0:11.8f} {1:11.8f} {2:11.8f}\n'.format(*matrixA[2])
        missettingsAngles = newmat['missettingsAngles']
        strNewmat += ' {0:11.3f} {0:11.3f} {0:11.3f}\n'.format(*missettingsAngles)
        matrixU = newmat['matrixU']
        strNewmat += ' {0:11.7f} {1:11.7f} {2:11.7f}\n'.format(*matrixU[0])
        strNewmat += ' {0:11.7f} {1:11.7f} {2:11.7f}\n'.format(*matrixU[1])
        strNewmat += ' {0:11.7f} {1:11.7f} {2:11.7f}\n'.format(*matrixU[2])
        cell = newmat['cell']
        strNewmat += ' {a:11.4f} {b:11.4f} {c:11.4f}'.format(**cell)
        strNewmat += ' {alpha:11.4f} {beta:11.4f} {gamma:11.4f}\n'.format(**cell)
        strNewmat += ' {0:11.3f} {0:11.3f} {0:11.3f}\n'.format(*missettingsAngles)
        with open(str(newmatPath), 'w') as f:
            f.write(strNewmat)

    # def generateExecutiveSummary(self, _edPlugin):
    #     """
    #     Generates a summary of the execution of the plugin.
    #     This method is common to all MOSFLM plugins.
    #     """
    #     self.DEBUG("EDPluginMOSFLMv10.generateExecutiveSummary")
    #     if (self.getStringVersion() is not None):
    #         self.addExecutiveSummaryLine(self.getStringVersion())


class MosflmIndexingTask(AbstractMosflmTask):

    def generateMOSFLMCommands(self, inData, workingDirectory):
        """
        This method creates a list of MOSFLM indexing commands given a valid
        inData
        """
        listCommand = AbstractMosflmTask.generateMOSFLMCommands(self,
                            inData, workingDirectory)
        listCommand.append('NEWMAT ' + self.getNewmatFileName())
        listImage = inData['image']
        for image in listImage:
            imageNumber = image['number']
            rotationAxisStart = image['rotationAxisStart']
            rotationAxisEnd = image['rotationAxisEnd']
            listCommand.append('AUTOINDEX DPS REFINE IMAGE {0} PHI {1} {2}'.format(
                imageNumber, rotationAxisStart, rotationAxisEnd))
        listCommand.append('GO')
        for image in listImage:
            listCommand.append('MOSAIC ESTIMATE {0}'.format(image['number']))
            listCommand.append('GO')
        return listCommand

    def parseMosflmOutput(self, workingDirectory):
        newMatFilePath = workingDirectory / self.getNewmatFileName()
        dnaTablesFilePath = workingDirectory / 'dnaTables.xml'
        return self.parseIndexingMosflmOutput(newMatFilePath, dnaTablesFilePath)

    def parseIndexingMosflmOutput(self, newMatFilePath, dnaTablesFilePath):
        outData = {}
        # Read Newmat file
        newmat = self.getNewmat(newMatFilePath)
        if newmat is None:
            errorMessage = "MOSFLM indexing error : No solution was obtained!"
            logger.error(errorMessage)
            self.setFailure()
        else:
            outData['newmat'] = self.getNewmat(newMatFilePath)
            # Then read the XML file
            dictDnaTables = UtilsDnaTables.getDict(dnaTablesFilePath)
            # Mosaicity estimation
            listTablesMosaicity = UtilsDnaTables.getTables(
                dictDnaTables, 'mosaicity_estimation')
            mosaicityValueSum = 0.0
            nValues = 0
            for table in listTablesMosaicity:
                listParam = UtilsDnaTables.getListParam(table)
                mosaicityValue = UtilsDnaTables.getListValue(
                    listParam, 'mosaicity', 'value')
                mosaicityValueSum += mosaicityValue
                nValues += 1
            outData['mosaicityEstimation'] = mosaicityValueSum / nValues
            # Refinement
            listTablesRefinement = UtilsDnaTables.getTables(
                dictDnaTables, 'refinement')
            tableRefinement = listTablesRefinement[0]
            listParamRefinement = UtilsDnaTables.getListParam(tableRefinement)
            outData['deviationAngular'] = UtilsDnaTables.getListValue(
                listParamRefinement, 'deviations', 'angular')
            outData['refinedDistance'] = UtilsDnaTables.getListValue(
                listParamRefinement, 'results', 'detector_distance')
            outData['spotsUsed'] = UtilsDnaTables.getListValue(
                listParamRefinement, 'parameters', 'used')
            outData['spotsTotal'] = UtilsDnaTables.getListValue(
                listParamRefinement, 'parameters', 'out_of')
            # Solution refinement
            listTablesSolutionRefinement = UtilsDnaTables.getTables(
                dictDnaTables, 'solution_refinement')
            tableSolutionRefinement = listTablesSolutionRefinement[0]
            listSolutionRefinement = UtilsDnaTables.getListParam(tableSolutionRefinement)
            outData['selectedSolutionNumber'] = UtilsDnaTables.getListValue(
                listSolutionRefinement, 'selection', 'number')
            outData['selectedSolutionSpaceGroup'] = UtilsDnaTables.getListValue(
                listSolutionRefinement, 'selection', 'spacegroup')
            outData['selectedSolutionSpaceGroupNumber'] = UtilsDnaTables.getListValue(
                listSolutionRefinement, 'selection', 'spacegroup_number')
            # Autoindex solutions
            listTablesAutoIndexSolution = UtilsDnaTables.getTables(
                dictDnaTables, 'autoindex_solutions')
            listIndexingSolution = []
            tablesAutoIndexSolution = listTablesAutoIndexSolution[0]
            listAutoindexSolutions = UtilsDnaTables.getListParam(tablesAutoIndexSolution)
            for autoIndexSolution in listAutoindexSolutions:
                indexingSolution = {
                    'index': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'index'),
                    'penalty': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'penalty'),
                    'lattice': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'lattice'),
                    'a': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'a'),
                    'b': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'b'),
                    'c': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'c'),
                    'alpha': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'alpha'),
                    'beta': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'beta'),
                    'gamma': UtilsDnaTables.getItemValue(
                        autoIndexSolution, 'gamma')
                }
                listIndexingSolution.append(indexingSolution)
            outData['indexingSolution'] = listIndexingSolution
            #  Beam refinement
            listTablesBeamRefinement = UtilsDnaTables.getTables(
                dictDnaTables, 'beam_refinement')
            tableBeamRefinement = listTablesBeamRefinement[0]
            listRefinedBeam = UtilsDnaTables.getListParam(tableBeamRefinement)
            initialBeamPositionX = UtilsDnaTables.getListValue(
                listRefinedBeam, 'initial_beam', 'x')
            initialBeamPositionY = UtilsDnaTables.getListValue(
                listRefinedBeam, 'initial_beam', 'y')
            refinedBeamPositionX = UtilsDnaTables.getListValue(
                listRefinedBeam, 'refined_beam', 'x')
            refinedBeamPositionY = UtilsDnaTables.getListValue(
                listRefinedBeam, 'refined_beam', 'y')
            beamPositionRefined = {
                'x': refinedBeamPositionX,
                'y': refinedBeamPositionY
            }
            outData['refinedBeam'] = beamPositionRefined
            beamPositionShift = {
                'x': initialBeamPositionX - refinedBeamPositionX,
                'y': initialBeamPositionY - refinedBeamPositionY,
            }
            outData['beamShift'] = beamPositionShift
        return outData
    #
    #
    # def generateExecutiveSummary(self, _edPlugin):
    #     """
    #     Generates a summary of the execution of the plugin.
    #     """
    #     EDPluginMOSFLMv10.generateExecutiveSummary(self, _edPlugin)
    #     self.DEBUG("EDPluginMOSFLMIndexingv10.createDataMOSFLMOutputIndexing")
    #     xsDataMOSFLMInputIndexing = self.getDataInput()
    #     xsDataMOSFLMOutputIndexing = self.getDataOutput()
    #     if not self.isFailure():
    #         self.addExecutiveSummaryLine("Execution of MOSFLM indexing successful.")
    #     self.addExecutiveSummaryLine("Image directory         : %s" % xsDataMOSFLMInputIndexing.getDirectory().getValue())
    #     self.addExecutiveSummaryLine("Image template          : %s" % xsDataMOSFLMInputIndexing.getTemplate().getValue())
    #     strImagesUsed = "Images used in indexing : "
    #     for xsDataMOSFLMImage in xsDataMOSFLMInputIndexing.getImage():
    #         strImagesUsed += "%5d" % xsDataMOSFLMImage.getNumber().getValue()
    #     self.addExecutiveSummaryLine(strImagesUsed)
    #     if (xsDataMOSFLMInputIndexing.getSymmetry() is not None):
    #         self.addExecutiveSummaryLine("Target symmetry     : %s" % xsDataMOSFLMInputIndexing.getSymmetry().getValue())
    #         self.addExecutiveSummaryLine("")
    #     if (xsDataMOSFLMOutputIndexing is not None):
    #         self.addExecutiveSummaryLine("")
    #         self.addExecutiveSummaryLine("List of possible solutions (index, penalty, lattice and cell):")
    #         iSelectedSolutionNumber = xsDataMOSFLMOutputIndexing.getSelectedSolutionNumber().getValue()
    #         # Add all solutions with penalty < 100 + 1 solution with penalty > 100
    #         bAddSolution = False
    #         listMOSFLMOutputIndexing = xsDataMOSFLMOutputIndexing.getPossibleSolutions()
    #         iNumberOfSolutions = len(listMOSFLMOutputIndexing)
    #         for i in range(iNumberOfSolutions):
    #
    #             iPenalty = listMOSFLMOutputIndexing[i].getPenalty().getValue()
    #             if (i < (iNumberOfSolutions - 1)):
    #                 iPenaltyNext = listMOSFLMOutputIndexing[i + 1].getPenalty().getValue()
    #                 if ((iPenalty >= 100) and (iPenaltyNext <= 100)):
    #                     bAddSolution = True
    #             if (bAddSolution):
    #                 iIndex = listMOSFLMOutputIndexing[i].getIndex().getValue()
    #                 xsDataCell = listMOSFLMOutputIndexing[i].getCell()
    #                 strLattice = listMOSFLMOutputIndexing[i].getLattice().getValue()
    #                 strPossibleSolution = "%3d %4d %2s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" % \
    #                                             (iIndex, iPenalty, strLattice, \
    #                                               xsDataCell.getLength_a().getValue(),
    #                                               xsDataCell.getLength_b().getValue(),
    #                                               xsDataCell.getLength_c().getValue(),
    #                                               xsDataCell.getAngle_alpha().getValue(),
    #                                               xsDataCell.getAngle_beta().getValue(),
    #                                               xsDataCell.getAngle_gamma().getValue(),
    #                                              )
    #                 self.addExecutiveSummaryLine(strPossibleSolution)
    #             iPenaltyOld = iPenalty
    #         self.addExecutiveSummaryLine("")
    #         self.addExecutiveSummaryLine("Choosen solution number   : %14d" % iSelectedSolutionNumber)
    #         strSelectedSpaceGroup = xsDataMOSFLMOutputIndexing.getSelectedSolutionSpaceGroup().getValue()
    #         self.addExecutiveSummaryLine("Selected space group      : %14s" % (strSelectedSpaceGroup))
    #         xsDataCellRefined = xsDataMOSFLMOutputIndexing.getRefinedNewmat().getRefinedCell()
    #         self.addExecutiveSummaryLine("Refined cell              : %6.2f %7.2f %7.2f %5.1f %5.1f %5.1f" % (\
    #                                       xsDataCellRefined.getLength_a().getValue(),
    #                                       xsDataCellRefined.getLength_b().getValue(),
    #                                       xsDataCellRefined.getLength_c().getValue(),
    #                                       xsDataCellRefined.getAngle_alpha().getValue(),
    #                                       xsDataCellRefined.getAngle_beta().getValue(),
    #                                       xsDataCellRefined.getAngle_gamma().getValue()
    #                                       ))
    #         iSpotsTotal = xsDataMOSFLMOutputIndexing.getSpotsTotal().getValue()
    #         iSpotsUsed = xsDataMOSFLMOutputIndexing.getSpotsUsed().getValue()
    #         self.addExecutiveSummaryLine("Number of spots used      : %14d " % (iSpotsUsed))
    #         self.addExecutiveSummaryLine("Number of spots total     : %14d " % (iSpotsTotal))
    #         fDeviationPositional = xsDataMOSFLMOutputIndexing.getDeviationPositional().getValue()
    #         fDeviationAngular = xsDataMOSFLMOutputIndexing.getDeviationAngular().getValue()
    #         self.addExecutiveSummaryLine("Spot deviation positional : %14.2f [mm]" % (fDeviationPositional))
    #         self.addExecutiveSummaryLine("Spot deviation angular    : %14.2f [degrees]" % (fDeviationAngular))
    #         fBeamshiftX = xsDataMOSFLMOutputIndexing.getBeamShift().getX().getValue()
    #         fBeamshiftY = xsDataMOSFLMOutputIndexing.getBeamShift().getY().getValue()
    #         self.addExecutiveSummaryLine("Beam shift (X, Y)         : %6.3f, %6.3f [mm]" % \
    #                                       (fBeamshiftX, fBeamshiftY))
    #         fMosaicityEstimated = xsDataMOSFLMOutputIndexing.getMosaicityEstimation().getValue()
    #         self.addExecutiveSummaryLine("Estimated mosaicity       : %14.2f [degrees]" % fMosaicityEstimated)


class MosflmGeneratePredictionTask(AbstractMosflmTask):

    def __init__(self, inData):
        AbstractMosflmTask.__init__(self, inData)
        self.predictionFileName = None

    def generateMOSFLMCommands(self, inData, workingDirectory):
        """
        This method creates a list of MOSFLM generate prediction commands
        given a valid inData
        """
        listCommand = AbstractMosflmTask.generateMOSFLMCommands(self,
                                                                inData,
                                                                workingDirectory)
        template = inData['template']
        imageNumber = inData['image'][0]['number']
        self.predictionFileName = \
            MosflmGeneratePredictionTask.getImageFileNameFromTemplate(
                template, imageNumber)
        if self.predictionFileName is not None:
            self.predictionFileName += '_pred.jpg'
        rotationAxisStart = inData['image'][0]['rotationAxisStart']
        rotationAxisEnd = inData['image'][0]['rotationAxisEnd']
        listCommand.append('XGUI ON')
        listCommand.append('IMAGE %d PHI %f TO %f' % \
                           (imageNumber, rotationAxisStart, rotationAxisEnd))
        listCommand.append('GO')
        listCommand.append('PREDICT_SPOTS')
        listCommand.append(
            'CREATE_IMAGE PREDICTION ON BINARY TRUE FILENAME %s' % \
            self.predictionFileName)
        listCommand.append('RETURN')
        listCommand.append('EXIT')
        return listCommand

    @classmethod
    def getImageFileNameFromTemplate(cls, template, imageNumber):
        hashFound = False
        hasFinished = False
        firstHash = None
        noHashes = 0
        try:
            for index, character in enumerate(template):
                if (not hashFound) and (not hasFinished):
                    if character == '#':
                        firstHash = index
                        hashFound = True
                else:
                    if (character != '#') and (not hasFinished):
                        hasFinished = True
                if hashFound and (not hasFinished):
                    noHashes += 1
            imageFileName = template[0:firstHash] + \
                               str(imageNumber).rjust(noHashes, '0')
        except Exception:
            logger.warning("Couldn't transform template {}".format(template) +
                           " to file name")
            imageFileName = None
        return imageFileName

    def parseMosflmOutput(self, workingDirectory):
        predictionPath = None
        if self.predictionFileName is not None:
            predictionPath = workingDirectory / self.predictionFileName
        outData = {
            'predictionImage': predictionPath
        }
        return outData
