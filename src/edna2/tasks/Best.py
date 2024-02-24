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


from edna2.tasks.AbstractTask import AbstractTask

# aimedCompleteness : XSDataDouble optional
# aimedIOverSigma : XSDataDouble optional
# aimedRedundancy : XSDataDouble optional
# aimedResolution : XSDataDouble optional
# anomalousData : XSDataBoolean optional
# beamExposureTime : XSDataTime
# beamMaxExposureTime : XSDataTime optional
# beamMinExposureTime : XSDataTime optional
# bestFileContentDat : XSDataString
# bestFileContentHKL : XSDataString []
# bestFileContentPar : XSDataString
# complexity : XSDataString optional
# crystalAbsorbedDoseRate : XSDataAbsorbedDoseRate optional
# crystalShape : XSDataDouble optional
# crystalSusceptibility : XSDataDouble optional
# detectorDistanceMax : XSDataLength optional
# detectorDistanceMin : XSDataLength optional
# detectorType : XSDataString
# doseLimit : XSDataDouble optional
# goniostatMaxRotationSpeed : XSDataAngularSpeed optional
# goniostatMinRotationWidth : XSDataAngle optional
# minTransmission : XSDataDouble optional
# numberOfCrystalPositions : XSDataInteger optional
# radiationDamageModelBeta : XSDataDouble optional
# radiationDamageModelGamma : XSDataDouble optional
# rFriedel : XSDataDouble optional
# strategyOption : XSDataString optional
# transmission : XSDataDouble optional
# userDefinedRotationRange : XSDataAngle optional
# userDefinedRotationStart : XSDataAngle optional
# xdsBackgroundImage : XSDataFile optional
# xdsCorrectLp : XSDataFile optional
# xdsBkgpixCbf : XSDataFile optional
# xdsAsciiHkl : XSDataFile [] optional


class Best(AbstractTask):
    """
    This task runs the program BEST
    http://www.embl-hamburg.de/BEST
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "diffractionPlan": {
                    "$ref": self.getSchemaUrl("ispybDiffractionPlan.json")
                },
                "subWedge": {
                    "type": "array",
                    "items": {"$ref": self.getSchemaUrl("subWedge.json")},
                },
                "xdsBackgroundImage": {"type": "string"},
                "correctLp": {"type": "string"},
                "xdsBkgpixCbf": {"type": "string"},
                "xdsAsciiHkl": {
                    "type": "array",
                    "items": {
                        "type": "string",
                    },
                },
            },
        }

    def run(self, inData):
        commandLine = self.createBestCommandLine(inData)
        self.runCommandLine(commandLine, listCommand=[])
        outData = {}
        return outData

    @staticmethod
    def createBestCommandLine(inData):
        diffractionPlan = inData.get("diffractionPlan", {})
        firstSubWedge = inData["subWedge"][0]
        experimentalCondition = firstSubWedge["experimentalCondition"]
        detector = experimentalCondition["detector"]
        beam = experimentalCondition["beam"]
        commandLine = "best__v4.1.1_test_20201220"
        commandLine += Best.addOption(detector, "type", "-f")
        commandLine += Best.addOption(beam, "exposureTime", "-t")
        commandLine += Best.addOption(inData, "crystalAbsorbedDoseRate", "-GpS")
        commandLine += Best.addOption(diffractionPlan, "aimedCompleteness", "-C")
        commandLine += Best.addOption(diffractionPlan, "aimedIOverSigma", "-i2s")
        commandLine += Best.addOption(diffractionPlan, "aimedRedundancy", "-R")
        commandLine += Best.addOption(diffractionPlan, "aimedResolution", "-r")
        if diffractionPlan.get("anomalousData", False):
            if "numberOfCrystalPositions" in diffractionPlan:
                commandLine += " -a -p 0 360"
            elif "absorbedDoseRate" in inData:
                commandLine += " -asad"
            else:
                commandLine += " -a"
        commandLine += Best.addOption(diffractionPlan, "complexity", "-e")
        commandLine += Best.addOption(inData, "crystalShape", "-sh")
        commandLine += Best.addOption(diffractionPlan, "crystalSusceptibility", "-su")
        commandLine += Best.addOption(
            diffractionPlan, "detectorDistanceMax", "-DIS_MAX"
        )
        commandLine += Best.addOption(
            diffractionPlan, "detectorDistanceMin", "-DIS_MIN"
        )
        commandLine += Best.addOption(diffractionPlan, "doseLimit", "-DMAX")
        commandLine += Best.addOption(
            diffractionPlan, "goniostatMaxRotationSpeed", "-S"
        )
        commandLine += Best.addOption(
            diffractionPlan, "goniostatMinRotationWidth", "-w"
        )
        commandLine += Best.addOption(
            diffractionPlan, "maxExposureTimePerDataCollection", "-T"
        )
        commandLine += Best.addOption(diffractionPlan, "minExposureTimePerImage", "-M")
        commandLine += Best.addOption(diffractionPlan, "minTransmission", "-TRmin")
        commandLine += Best.addOption(
            diffractionPlan, "numberOfCrystalPositions", "-Npos"
        )
        commandLine += Best.addOption(
            diffractionPlan, "radiationDamageModelBeta", "-beta"
        )
        commandLine += Best.addOption(
            diffractionPlan, "radiationDamageModelGamma", "-gama"
        )
        commandLine += Best.addOption(diffractionPlan, "rFriedel", "-Rf")
        commandLine += Best.addOption(diffractionPlan, "strategyOption", "")
        commandLine += Best.addOption(diffractionPlan, "transmission", "-Trans")
        if (
            "userDefinedRotationRange" in diffractionPlan
            and "userDefinedRotationStart" in diffractionPlan
        ):
            commandLine += (
                " -phi {userDefinedRotationStart} {userDefinedRotationRange}".format(
                    **diffractionPlan
                )
            )
        # Output of GLE files and plotmtv files
        commandLine += " -g -o plot.mtv"
        # Integration data
        commandLine += " -xds " + inData["bkgpixCbf"]
        for xdsAsciiHklPath in inData["xdsAsciiHkl"]:
            commandLine += " " + xdsAsciiHklPath
        return commandLine

    @staticmethod
    def addOption(inData, name, option):
        returnValue = ""
        value = inData.get(name)
        if value is not None:
            returnValue = " {0} {1}".format(option, value)
        return returnValue
