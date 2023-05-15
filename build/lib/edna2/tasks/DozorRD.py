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

__author__ = "Olof Svensson"
__contact__ = "svensson@esrf.eu"
__copyright__ = "ESRF"
__updated__ = "2022-02-13"

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsLogging


logger = UtilsLogging.getLogger()


class DozorRD(AbstractTask):  # pylint: disable=too-many-instance-attributes
    """
    The DozorM2 is responsible for executing the 'dozorm2' program.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "list_dozor_all": {
                    "type": "array",
                    "items": {"type": "string"},
                },
                "wavelength": {"type": "number"},
                "exposureTime": {"type": "number"},
                "numberOfImages": {"type": "number"},
            },
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "logPath": {"type": "string"},
                "noSpotDecrease": {"type": "number"},
                "mainScore": {"type": "number"},
                "spotIntensity": {"type": "number"},
                "sumIntensity": {"type": "number"},
                "average": {"type": "number"},
            },
        }

    def run(self, inData):
        outData = {}
        commandLine = "dozorrd dozorrd.dat"
        for index, dozor_all_file in enumerate(inData["list_dozor_all"]):
            os.symlink(
                dozor_all_file, str(self.getWorkingDirectory() / "d_{0:03d}".format(index + 1))
            )
        commands = self.generateCommands(inData)
        with open(str(self.getWorkingDirectory() / "dozorrd.dat"), "w") as f:
            f.write(commands)
        self.setLogFileName("dozorrd.log")
        self.runCommandLine(commandLine)
        # Send log and error log to the logger
        logger.info("\n\n" + self.getLog())
        logger.info("")
        errorLog = self.getErrorLog()
        if errorLog != "":
            logger.info("\n\n" + self.getErrorLog())
            logger.info("")
        outData = self.parseDozorRDLogFile(self.getLogPath())
        # Generate plots
        plotMtvFile = self.getWorkingDirectory() / "dozor_rd.mtv"
        if plotMtvFile.exists():
            listPlotFile = self.generatePngPlots(
                plotmtvFile=plotMtvFile,
                workingDir=self.getWorkingDirectory()
            )
            outData["listPlotFile"] = listPlotFile
        return outData

    @staticmethod
    def generateCommands(inData):
        """
        This method creates the input file for dozorm
        """
        nameTemplateScan = "d_0??"
        command = "X-ray_wavelength {0}\n".format(inData["wavelength"])
        command += "exposure {0}\n".format(inData["exposureTime"])
        command += "number_images {0}\n".format(inData["numberOfImages"])
        command += "first_scan_number 1\n"
        command += "number_scans {0}\n".format(len(inData["list_dozor_all"]))
        command += "name_template_scan {0}\n".format(nameTemplateScan)
        command += "end\n"
        return command

    @staticmethod
    def parseDozorRDLogFile(logPath):
        #  Program dozorrd /A.Popov,G.Bourenkov /
        #  Version 1.1 //  31.01.2022
        #  Copyright 2020 by Alexander Popov and Gleb Bourenkov
        #   T1/2 estimates
        #  N.of spot decrease= 1.751
        #       main score   =10.719
        #   spot intensity   =67.269
        #   sum  intensity   = 2.218
        #   average   =26.735
        outData = {}
        outData["logPath"] = str(logPath)
        with open(str(logPath)) as fd:
            listLogLines = fd.readlines()
        dictParse = {
            "N.of spot decrease": "noSpotDecrease",
            "main score": "mainScore",
            "spot intensity": "spotIntensity",
            "sum  intensity": "sumIntensity",
            "average": "average",
        }
        for line in listLogLines:
            if "=" in line:
                key, value = line.split("=")
                key = key.strip()
                if key in dictParse:
                    try:
                        outData[dictParse[key]] = float(value)
                    except ValueError:
                        outData[dictParse[key]] = -9999
                        errorMessage = "Error when parsing '{0}': '{1}'".format(key, value.strip())
                        if "errorMessage" in outData:
                            outData["errorMessage"] += " " + errorMessage
                        else:
                            outData["errorMessage"] = errorMessage
        return outData

    def generatePngPlots(self, plotmtvFile, workingDir):
        listPlotFile = []
        # Create plot dictionary
        with open(plotmtvFile) as f:
            listLines = f.readlines()
        dictPlot = None
        plotData = None
        dictPlotList = None
        listPlots = []
        index = 0
        while index < len(listLines):
            # print("0" + listLines[index])
            if listLines[index].startswith("$"):
                dictPlot = {}
                dictPlotList = []
                listPlots.append(dictPlot)
                dictPlot["plotList"] = dictPlotList
                index += 1
                dictPlot["name"] = listLines[index].split("'")[1]
                index += 1
                # print(listLines[index])
                while listLines[index].startswith("%"):
                    listLine = listLines[index].split("=")
                    # print(listLine)
                    label = listLine[0][1:].strip()
                    # print("label: " + str([label]))
                    if "'" in listLine[1]:
                        value = listLine[1].split("'")[1]
                    else:
                        value = listLine[1]
                    value = value.replace("\n", "").strip()
                    # print("value: " + str([value]))
                    dictPlot[label] = value
                    index += 1
                    # print(listLines[index])
            elif listLines[index].startswith("#"):
                dictSubPlot = {}
                dictPlotList.append(dictSubPlot)
                plotName = listLines[index].split("#")[1].replace("\n", "").strip()
                dictSubPlot["name"] = plotName
                index += 1
                # print("1" + listLines[index])
                while listLines[index].startswith("%"):
                    listLine = listLines[index].split("=")
                    # print(listLine)
                    label = listLine[0][1:].strip()
                    # print("label: " + str([label]))
                    if "'" in listLine[1]:
                        value = listLine[1].split("'")[1]
                    else:
                        value = listLine[1]
                    value = value.replace("\n", "").strip()
                    # print("value: " + str([value]))
                    dictSubPlot[label] = value
                    index += 1
                    # print(listLines[index])
                dictSubPlot["xValues"] = []
                dictSubPlot["yValues"] = []
            else:
                listData = listLines[index].replace("\n", "").split()
                dictSubPlot["xValues"].append(float(listData[0]))
                dictSubPlot["yValues"].append(float(listData[1]))
                index += 1
        # pprint.pprint(listPlots)
        # Generate the plots
        for mtvplot in listPlots:
            listLegend = []
            xmin = None
            xmax = None
            ymin = None
            ymax = None
            for subPlot in mtvplot["plotList"]:
                xmin1 = min(subPlot["xValues"])
                if xmin is None or xmin > xmin1:
                    xmin = xmin1
                xmax1 = max(subPlot["xValues"])
                if xmax is None or xmax < xmax1:
                    xmax = xmax1
                ymin1 = min(subPlot["yValues"])
                if ymin is None or ymin > ymin1:
                    ymin = ymin1
                ymax1 = max(subPlot["yValues"])
                if ymax is None or ymax < ymax1:
                    ymax = ymax1
            if "xmin" in mtvplot:
                xmin = float(mtvplot["xmin"])
            if "ymin" in mtvplot:
                ymin = float(mtvplot["ymin"])
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            plt.xlabel(mtvplot["xlabel"])
            plt.ylabel(mtvplot["ylabel"])
            plt.title(mtvplot["name"])
            for subPlot in mtvplot["plotList"]:
                if "markercolor" in subPlot:
                    style = "bs-."
                else:
                    style = "r"
                plt.plot(subPlot["xValues"], subPlot["yValues"], style, linewidth=2)
                listLegend.append(subPlot["linelabel"])
            plt.legend(listLegend, loc='lower right')
            plotPath = os.path.join(str(workingDir), mtvplot["name"].replace(" ", "").replace(".", "_") + ".png")
            plt.savefig(plotPath, bbox_inches='tight', dpi=75)
            plt.close()
            listPlotFile.append(plotPath)
        return listPlotFile
