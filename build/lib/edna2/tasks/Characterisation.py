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
__date__ = "23/11/2020"

import math

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.Best import Best
from edna2.tasks.ControlIndexing import ControlIndexing
from edna2.tasks.XDSTasks import XDSGenerateBackground
from edna2.tasks.XDSTasks import XDSIntegration
from edna2.tasks.ReadImageHeader import ReadImageHeader
from edna2.tasks.Raddose import Raddose
from edna2.tasks.H5ToCBFTask import H5ToCBFTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsSymmetry


class Characterisation(AbstractTask):
    """
    This task receives a list of images or data collection ids and
    returns result of indexing
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "dataCollectionId": {"type": "integer"},
                "diffractionPlan": {
                    "$ref": self.getSchemaUrl("ispybDiffractionPlan.json")
                },
                "experimentalCondition": {
                    "$ref": self.getSchemaUrl("experimentalCondition.json")
                },
                "imagePath": {
                    "type": "array",
                    "items": {
                        "type": "string",
                    },
                },
                "sample": {"$ref": self.getSchemaUrl("sample.json")},
                "token": {"anyOf": [{"type": "string"}, {"type": "null"}]},
            },
        }

    def run(self, inData):
        listImagePath = inData["imagePath"]
        prefix = UtilsImage.getPrefix(listImagePath[0])
        listSubWedge = self.getListSubWedge(inData)
        diffractionPlan = inData.get("diffractionPlan", {})
        sample = inData.get("sample", {})
        # Check if flux is present
        flux = None
        absorbedDoseRate = None
        experimentalCondition = inData.get("experimentalCondition", None)
        if experimentalCondition is not None:
            beam = experimentalCondition.get("beam", None)
            if beam is not None:
                flux = beam.get("flux", None)
                # if not "size" in beam:
                #     beam["size"] = listSubWedge[0]["experimentalCondition"]["beam"]["size"]
                beam["exposureTime"] = listSubWedge[0]["experimentalCondition"]["beam"][
                    "exposureTime"
                ]
                beam["wavelength"] = listSubWedge[0]["experimentalCondition"]["beam"][
                    "wavelength"
                ]
        # Convert images to CBF
        firstImage = listSubWedge[0]["image"][0]["path"]
        import os
        import pprint

        listH5ToCBF = []
        suffix = os.path.splitext(firstImage)[1]
        if suffix == ".h5":
            for subWedge in listSubWedge:
                imageList = subWedge["image"]
                for image in imageList:
                    imagePath = image["path"]
                    hdf5ImageNumber = UtilsImage.getImageNumber(imagePath)
                    inDataH5ToCBF = {
                        "hdf5File": imagePath,
                        "hdf5ImageNumber": hdf5ImageNumber,
                        "imageNumber": 1,
                        "forcedOutputDirectory": str(self.getWorkingDirectory()),
                        "forcedOutputImageNumber": hdf5ImageNumber,
                    }
                    h5ToCBF = H5ToCBFTask(inData=inDataH5ToCBF)
                    h5ToCBF.start()
                    listH5ToCBF.append((image, h5ToCBF))
        # Join CBF creation
        for image, h5ToCBF in listH5ToCBF:
            h5ToCBF.join()
            image["path"] = h5ToCBF.outData["outputCBFFile"]
        pprint.pprint(listSubWedge)
        # Start indexing
        outDataIndexing, outDataGB, listSubWedge = self.indexing(
            prefix, listSubWedge, self.getWorkingDirectory()
        )
        if outDataIndexing is not None and outDataGB is not None:
            listXdsAsciiHkl, correctLp, bkgpixCbf = self.integration(
                prefix, listSubWedge, outDataIndexing, outDataGB
            )
            if listXdsAsciiHkl is not None:
                # Check if Raddose should be run
                estimateRadiationDamage = self.checkEstimateRadiationDamage(
                    inData, flux
                )
                if estimateRadiationDamage:
                    # Check if forced space group
                    forcedSpaceGroup = None
                    numOperators = None
                    cell = None
                    if (
                        "diffractionPlan" in inData
                        and "forcedSpaceGroup" in inData["diffractionPlan"]
                    ):
                        forcedSpaceGroup = inData["diffractionPlan"]["forcedSpaceGroup"]
                        if forcedSpaceGroup is not None:
                            if forcedSpaceGroup != "":
                                forcedSpaceGroup = forcedSpaceGroup.replace(" ", "")
                                numOperators = UtilsSymmetry.getNumberOfSymmetryOperatorsFromSpaceGroupName(
                                    forcedSpaceGroup
                                )
                            else:
                                forcedSpaceGroup = None
                    if forcedSpaceGroup is None:
                        # Get indexing space group IT number
                        if "resultIndexing" in outDataIndexing:
                            resultIndexing = outDataIndexing["resultIndexing"]
                            cell = resultIndexing["cell"]
                            if "spaceGroupNumber" in resultIndexing:
                                spaceGroupNumber = resultIndexing["spaceGroupNumber"]
                                numOperators = UtilsSymmetry.getNumberOfSymmetryOperatorsFromSpaceGroupITNumber(
                                    spaceGroupNumber
                                )
                    if numOperators is None:
                        raise RuntimeError(
                            "Error when trying to determine number of symmetry operators!"
                        )
                    chemicalComposition = self.getDefaultChemicalComposition(
                        cell, numOperators
                    )
                    numberOfImages = self.getNumberOfImages(listSubWedge)
                    sample = inData["sample"]
                    inDataRaddose = {
                        "experimentalCondition": experimentalCondition,
                        "chemicalComposition": chemicalComposition,
                        "sample": sample,
                        "cell": cell,
                        "numberOfImages": numberOfImages,
                        "numOperators": numOperators,
                    }
                    # import pprint
                    # pprint.pprint(inDataRaddose)
                    raddose = Raddose(
                        inData=inDataRaddose, workingDirectorySuffix=prefix
                    )
                    raddose.execute()
                    if raddose.isSuccess():
                        absorbedDoseRate = raddose.outData["absorbedDoseRate"]
                inDataBest = {
                    "diffractionPlan": diffractionPlan,
                    "sample": sample,
                    "subWedge": listSubWedge,
                    "xdsAsciiHkl": listXdsAsciiHkl,
                    "bkgpixCbf": bkgpixCbf,
                    "correctLp": correctLp,
                    "crystalAbsorbedDoseRate": absorbedDoseRate,
                }
                bestTask = Best(inData=inDataBest, workingDirectorySuffix=prefix)
                bestTask.execute()

    @staticmethod
    def indexing(prefix, listSubWedge, workingDirectory):
        outDataIndexing = None
        outDataGB = None
        # Start indexing
        inDataIndexing = {
            "subWedge": listSubWedge,
        }
        indexingTask = ControlIndexing(
            inData=inDataIndexing, workingDirectorySuffix=prefix
        )
        indexingTask.start()
        # Start background esitmation
        inDataGenerateBackground = {
            "subWedge": [listSubWedge[0]],
        }
        generateBackground = XDSGenerateBackground(
            inData=inDataGenerateBackground, workingDirectorySuffix=prefix
        )
        generateBackground.start()
        generateBackground.join()
        # Check indexing
        indexingTask.join()
        if indexingTask.isSuccess():
            outDataIndexing = indexingTask.outData
            outDataGB = generateBackground.outData
        return outDataIndexing, outDataGB, listSubWedge

    @staticmethod
    def integration(prefix, listSubWedge, outDataIndexing, outDataGB):
        listXdsAsciiHkl = None
        correctLp = None
        bkgpixCbf = None
        resultIndexing = outDataIndexing["resultIndexing"]
        xparmXdsPath = outDataIndexing["xparmXdsPath"]
        if xparmXdsPath is not None:
            listTasks = []
            listXdsAsciiHkl = []
            for subWedge in listSubWedge:
                inDataIntergation = {
                    "subWedge": [subWedge],
                    "xparmXds": xparmXdsPath,
                    "spaceGroupNumber": resultIndexing["spaceGroupNumber"],
                    "cell": resultIndexing["cell"],
                    "gainCbf": outDataGB["gainCbf"],
                    "blankCbf": outDataGB["blankCbf"],
                    "bkginitCbf": outDataGB["bkginitCbf"],
                    "xCorrectionsCbf": outDataGB["xCorrectionsCbf"],
                    "yCorrectionsCbf": outDataGB["yCorrectionsCbf"],
                }
                imageNo = subWedge["image"][0]["number"]
                integrationTask = XDSIntegration(
                    inData=inDataIntergation,
                    workingDirectorySuffix=prefix + "_{0:04d}".format(imageNo),
                )
                integrationTask.start()
                listTasks.append(integrationTask)
            for task in listTasks:
                task.join()
                if "xdsAsciiHkl" in task.outData:
                    listXdsAsciiHkl.append(task.outData["xdsAsciiHkl"])
                if correctLp is None:
                    correctLp = task.outData["correctLp"]
                    bkgpixCbf = task.outData["bkgpixCbf"]
        return listXdsAsciiHkl, correctLp, bkgpixCbf

    @staticmethod
    def getListSubWedge(inData):
        listSubWedge = None
        # First check if we have data collection ids or image list
        # if "dataCollectionId" in inData:
        #     # TODO: get list of data collections from ISPyB
        #         logger.warning("Not implemented!")
        # el
        if "imagePath" in inData:
            listSubWedge = Characterisation.readImageHeaders(inData["imagePath"])
        else:
            raise RuntimeError("No dataCollectionId or imagePath in inData")
        return listSubWedge

    @staticmethod
    def readImageHeaders(listImagePath):
        # Read the header(s)
        inDataReadImageHeader = {"imagePath": listImagePath}
        readImageHeader = ReadImageHeader(
            inData=inDataReadImageHeader,
            workingDirectorySuffix=UtilsImage.getPrefix(listImagePath[0]),
        )
        readImageHeader.execute()
        listSubWedge = readImageHeader.outData["subWedge"]
        return listSubWedge

    @staticmethod
    def getDefaultChemicalComposition(cell, numOperators):

        # For default chemical composition
        averageAminoAcidVolume = 135.49
        averageCrystalSolventContent = 0.47
        averageSulfurContentPerAminoacid = 0.05
        averageSulfurConcentration = 314

        a = cell["a"]
        b = cell["b"]
        c = cell["c"]
        alpha = math.radians(cell["alpha"])
        beta = math.radians(cell["beta"])
        gamma = math.radians(cell["gamma"])

        unitCellVolume = (
            a
            * b
            * c
            * (
                math.sqrt(
                    1
                    - math.cos(alpha) * math.cos(alpha)
                    - math.cos(beta) * math.cos(beta)
                    - math.cos(gamma) * math.cos(gamma)
                    + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
                )
            )
        )
        polymerVolume = unitCellVolume * (1 - averageCrystalSolventContent)
        numberOfMonomersPerUnitCell = round(polymerVolume / averageAminoAcidVolume)
        numberOfMonomersPerAsymmetricUnit = round(
            numberOfMonomersPerUnitCell / numOperators
        )
        numberOfSulfurAtom = int(
            round(numberOfMonomersPerAsymmetricUnit * averageSulfurContentPerAminoacid)
        )

        chemicalCompositionMM = {
            "solvent": {
                "atom": [
                    {
                        "symbol": "S",
                        "concentration": averageSulfurConcentration,
                    }
                ]
            },
            "structure": {
                "chain": [
                    {
                        "type": "protein",
                        "numberOfCopies": 1,
                        "numberOfMonomers": numberOfMonomersPerAsymmetricUnit,
                        "heavyAtoms": [{"symbol": "S", "numberOf": numberOfSulfurAtom}],
                    }
                ],
                "numberOfCopiesInAsymmetricUnit": 1,
            },
        }

        return chemicalCompositionMM

    @staticmethod
    def checkEstimateRadiationDamage(inData, flux=None):
        estimateRadiationDamage = None
        # Check if radiation damage estimation is required or not in the diffraction plan
        diffractionPlan = inData.get("diffractionPlan", None)
        if diffractionPlan is not None:
            if "estimateRadiationDamage" in diffractionPlan:
                # Estimate radiation damage is explicitly set in the diffraction plan
                estimateRadiationDamage = diffractionPlan["estimateRadiationDamage"]
            else:
                strategyOption = diffractionPlan.get("strategyOption", None)
                if strategyOption is not None and "-DamPar" in strategyOption:
                    # The "-DamPar" option requires estimation of radiation damage
                    estimateRadiationDamage = True

        # Check if not set by the diffraction plan
        if estimateRadiationDamage is None:
            experimentalCondition = inData.get("experimentalCondition", None)
            if experimentalCondition is not None:
                beam = experimentalCondition.get("beam", None)
                if beam is not None:
                    flux = beam.get("flux")
                    if flux is not None:
                        estimateRadiationDamage = True

        if estimateRadiationDamage is None:
            estimateRadiationDamage = False

        return estimateRadiationDamage

    @staticmethod
    def getNumberOfImages(listSubwedge):
        noImages = 0
        for subWedge in listSubwedge:
            goniostat = subWedge["experimentalCondition"]["goniostat"]
            oscStart = goniostat["rotationAxisStart"]
            oscEnd = goniostat["rotationAxisEnd"]
            oscWidth = goniostat["oscillationWidth"]
            noImages += int(round((oscEnd - oscStart) / oscWidth, 0))
        return noImages
