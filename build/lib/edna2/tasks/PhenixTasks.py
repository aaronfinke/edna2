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
__date__ = "26/07/2019"

import pathlib

from edna2.tasks.AbstractTask import AbstractTask


class DistlSignalStrengthTask(AbstractTask):
    """
    This task runs the labelit.distl command for pre-screening reference images.
    """

    def run(self, inData):
        commandLine = 'distl.signal_strength '
        commandLine += inData['referenceImage']
        logPath = self.getWorkingDirectory() / 'distl.log'
        self.runCommandLine(commandLine, logPath=logPath)
        with open(str(logPath)) as f:
            logText = f.read()
        imageQualityIndicators = self.parseLabelitDistlOutput(logText)
        outData = {'imageQualityIndicators': imageQualityIndicators}
        return outData

    def parseLabelitDistlOutput(self, logText):
        imageQualityIndicators = {}
        for line in logText.split("\n"):
            if line.find("Spot Total") != -1:
                spotTotal = int(line.split()[3])
                imageQualityIndicators['spotTotal'] = spotTotal
            elif line.find("In-Resolution Total") != -1:
                inResTotal = int(line.split()[3])
                imageQualityIndicators['inResTotal'] = inResTotal
            elif line.find("Good Bragg Candidates") != -1:
                goodBraggCandidates = int(line.split()[4])
                imageQualityIndicators['goodBraggCandidates'] = goodBraggCandidates
            elif line.find("Ice Rings") != -1:
                iceRings = int(line.split()[3])
                imageQualityIndicators['iceRings'] = iceRings
            elif line.find("Method 1 Resolution") != -1:
                method1Res = float(line.split()[4])
                imageQualityIndicators['method1Res'] = method1Res
            elif line.find("Method 2 Resolution") != -1:
                if line.split()[4] != "None":
                    fMethod2Res = float(line.split()[4])
                    imageQualityIndicators['method2Res'] = fMethod2Res
            elif line.find("Maximum unit cell") != -1:
                if line.split()[4] != "None":
                    fMaxUnitCell = float(line.split()[4])
                    imageQualityIndicators['maxUnitCell'] = fMaxUnitCell
            elif line.find("%Saturation, Top 50 Peaks") != -1:
                pctSaturationTop50Peaks = float(line.split()[5])
                imageQualityIndicators['pctSaturationTop50Peaks'] = pctSaturationTop50Peaks
            elif line.find("In-Resolution Ovrld Spots") != -1:
                inResolutionOvrlSpots = int(line.split()[4])
                imageQualityIndicators['inResolutionOvrlSpots'] = inResolutionOvrlSpots
            elif line.find("Bin population cutoff for method 2 resolution") != -1:
                binPopCutOffMethod2Res = float(line.split()[7][:-1])
                imageQualityIndicators['binPopCutOffMethod2Res'] = binPopCutOffMethod2Res
            elif line.find("Total integrated signal, pixel-ADC units above local background (just the good Bragg candidates)") != -1:
                totalIntegratedSignal = float(line.split()[-1])
                imageQualityIndicators['totalIntegratedSignal'] = totalIntegratedSignal
            elif line.find("signals range from") != -1:
                listStrLine = line.split()
                signalRangeMin = float(listStrLine[3])
                signalRangeMax = float(listStrLine[5])
                signalRangeAverage = float(listStrLine[-1])
                imageQualityIndicators['signalRangeMin'] = signalRangeMin
                imageQualityIndicators['signalRangeMax'] = signalRangeMax
                imageQualityIndicators['signalRangeAverage'] = signalRangeAverage
            elif line.find("Saturations range from") != -1:
                listStrLine = line.split()
                saturationRangeMin = float(listStrLine[3][:-1])
                saturationRangeMax = float(listStrLine[5][:-1])
                saturationRangeAverage = float(listStrLine[-1][:-1])
                imageQualityIndicators['saturationRangeMin'] = saturationRangeMin
                imageQualityIndicators['saturationRangeMax'] = saturationRangeMax
                imageQualityIndicators['saturationRangeAverage'] = saturationRangeAverage
        return imageQualityIndicators
