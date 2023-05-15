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
__date__ = "11/01/2021"



from edna2.tasks.AbstractTask import AbstractTask




class Raddose(AbstractTask):
    """
    This task runs the program raddose
    """
    HENDERSON_LIMIT = 2e7

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "exposureTime": {"type": "number"},
                "flux": {"type": "number"},
                "beamSizeX": {"type": "number"},
                "beamSizeY": {"type": "number"},
                "wavelength": {"type": "number"},
                "numberOfImages": {"type": "number"},
            }
        }

    def run(self, inData):
        commandLine, listCommand = self.createCommandLine(inData)
        self.runCommandLine(commandLine, listCommand=listCommand)
        dictResults = self.parseLogFile(self.getLogPath())
        outData = self.createOutData(inData, dictResults, pathToLogFile=self.getLogPath())
        return outData

    @staticmethod
    def createCommandLine(inData):
        commandLine = "raddose"
        experimentalCondition = inData["experimentalCondition"]
        sample = inData["sample"]
        beam = experimentalCondition["beam"]
        numOperators = inData["numOperators"]
        # Analyse chemical composition
        chemicalComposition = inData["chemicalComposition"]
        # Solvent
        totalSATM = []
        solvent = chemicalComposition.get("solvent")
        if solvent is not None:
            totalSATM = solvent["atom"]
        # Structure
        structure = chemicalComposition.get("structure", None)
        totalNRESInStructure = 0
        totalNDNAInStructure = 0
        totalNRNAInStructure = 0
        totalPATM = []
        if structure is not None:
            listChain = structure["chain"]

            for chain in listChain:
                # heavy atoms of each chain to be added in the PATM
                atomicCompositionHeavyAtoms = chain.get("heavyAtoms", None)
                if atomicCompositionHeavyAtoms is not None:
                    iterator = 1
                    while iterator <= chain["numberOfCopies"]:
                        totalPATM = Raddose.mergeAtomicComposition(totalPATM, atomicCompositionHeavyAtoms)
                        iterator = iterator + 1

                type = chain["type"]
                numberOfMonomers = chain["numberOfMonomers"] * chain["numberOfCopies"]

                if (type == "protein"):
                    totalNRESInStructure = totalNRESInStructure + numberOfMonomers
                elif(type == "dna"):
                    totalNDNAInStructure = totalNDNAInStructure + numberOfMonomers
                elif(type == "rna"):
                    totalNRNAInStructure = totalNRNAInStructure + numberOfMonomers

            numStructInAU = structure["numberOfCopiesInAsymmetricUnit"]
            numStructInUC = int(numStructInAU * numOperators)

            # Ligands - not implemented yet
            if "ligand" in structure:
                raise RuntimeError("Ligand present in structure but support for ligands not yet implemented!")
            # xsDataLigands = xsDataStructure.getLigand()
            # for ligand in xsDataLigands:
            #
            #     # Light atoms to be added to the NRES
            #     nres = ligand.getNumberOfLightAtoms().getValue() * ligand.getNumberOfCopies().getValue() / 7.85
            #     totalNRESInStructure = totalNRESInStructure + nres
            #
            #     # Heavy atoms to be added to the PATM
            #     if (ligand.getHeavyAtoms() is not None):
            #         iterator = 1
            #         while iterator <= ligand.getNumberOfCopies().getValue():
            #             totalPATM = self.mergeAtomicComposition(totalPATM, ligand.getHeavyAtoms())
            #             iterator = iterator + 1

        listCommand = [
            "BEAM {x} {y}".format(**beam["size"]),
            "PHOSEC {flux}".format(**beam),
            "WAVELENGTH {wavelength}".format(**beam),
            # "CRYSTAL {x} {y} {z}".format(**sample["size"]),
            "CELL {a} {b} {c} {alpha} {beta} {gamma}".format(**inData["cell"]),
            "EXPOSURE {exposureTime}".format(**beam),
            "IMAGES {numberOfImages}".format(**inData),
            "NMON {0}".format(numStructInUC)
        ]

        if totalNRESInStructure != 0:
            listCommand.append("NRES {0}".format(int(round(totalNRESInStructure))))
        if totalNDNAInStructure != 0:
            listCommand.append("NDNA {0}".format(int(round(totalNDNAInStructure))))
        if totalNRNAInStructure != 0:
            listCommand.append("NRA {0}".format(int(round(totalNRNAInStructure))))
        if len(totalPATM) != 0:
            patmLine = "PATM"
            for patm in totalPATM:
                patmLine += " {symbol} {numberOf}".format(**patm)
            listCommand.append(patmLine)
        if len(totalSATM) != 0:
            satmLine = "SATM"
            for satm in totalSATM:
                satmLine += " {symbol} {concentration}".format(**satm)
            listCommand.append(satmLine)

        #
        # patmLine = "PATM"
        # for patm in inData["crystalPATM"]:
        #     patmLine += " {symbol} {numberOf}".format(**patm)
        # listCommand.append(patmLine)
        listCommand.append("END")
        return commandLine, listCommand

    @staticmethod
    def parseLogFile(logFilePath):
        dictResults = {}
        with open(str(logFilePath), "rb") as fd:
            content = fd.read()
        # Cut off the last 150 bytes
        content = content[0:len(content) - 150]
        for line in content.decode('utf-8').split("\n"):
            if "Dose in Grays" in line:
                dictResults["doseInGrays"] = round(float(line.split()[-1]), 1)
            elif "Total absorbed dose (Gy)" in line:
                dictResults["totalAbsorbedDose"] = round(float(line.split()[-1]), 1)
            elif "Solvent Content (%)" in line:
                dictResults["solventContent"] = round(float(line.split()[-1]), 1)
        return dictResults

    @staticmethod
    def createOutData(inData, dictResults, pathToLogFile=None):
        outData = {}
        numberOfImages = inData["numberOfImages"]
        exposureTime = inData["experimentalCondition"]["beam"]["exposureTime"]
        totalExposureTime = numberOfImages * exposureTime
        if "doseInGrays" in dictResults:
            absorbedDose = dictResults["doseInGrays"]
        elif "totalAbsorbedDose" in dictResults:
            absorbedDose = dictResults["totalAbsorbedDose"]
        else:
            raise RuntimeError["Neither doseInGrays nor totalAbsorbedDose in results!"]
        absorbedDoseRate = round(absorbedDose / totalExposureTime, 1)
        timeToReachHendersonLimit = round(Raddose.HENDERSON_LIMIT / absorbedDoseRate, 1)
        outData = {
            "absorbedDose": absorbedDose,
            "absorbedDoseRate": absorbedDoseRate,
            "pathToLogFile": pathToLogFile,
            "timeToReachHendersonLimit": timeToReachHendersonLimit
        }
        return outData

    @staticmethod
    def mergeAtomicComposition(atomicComposition1, atomicComposition2):

        mergedAtomicComposition = []
        dictionary = {}

        for atom2 in atomicComposition2:
            dictionary[atom2["symbol"]] = atom2["numberOf"]

        for atom1 in atomicComposition1:
            symbol = atom1["symbol"]
            if Raddose.exists(symbol, atomicComposition2):
                mergedAtom = {
                    "symbol": symbol,
                    "numberOf": atom1["numberOf"] + dictionary[symbol]
                }
                mergedAtomicComposition.append(mergedAtom)
            else:
                mergedAtomicComposition.append(atom1)

        for atom2 in atomicComposition2:
            symbol = atom2["symbol"]
            if not Raddose.exists(symbol, atomicComposition1):
                mergedAtomicComposition.append(atom2)

        return mergedAtomicComposition

    @staticmethod
    def exists(symbol, atomicComposition):
        exists = False
        for atom in atomicComposition:
            if atom["symbol"] == symbol:
                exists = True
                break
        return exists
