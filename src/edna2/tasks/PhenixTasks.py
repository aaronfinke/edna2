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
import gzip
import os

from typing import Tuple
import io
from mmtbx.scaling import printed_output
import html

from edna2.tasks.AbstractTask import AbstractTask
from edna2.utils import UtilsLogging
from edna2.utils import UtilsConfig

logger = UtilsLogging.getLogger()


class MmtbxXtriageTask(AbstractTask):
    """
    This task runs mmtbx.xtriage_analyses
    No command line required!
    Input should be unmerged AIMLESS MTZ output
    """

    def run(self, inData):
        if os.environ.get("DIALS", None) is None:
            xia2DialsSetup = UtilsConfig.get("Xia2DialsTask", "xia2DialsSetup", None)

            commandLine = "source {} \n".format(xia2DialsSetup)
        else:
            commandLine = ""
            logger.info(f"DIALS version is {os.environ.get('DIALS_VERSION', None)}")

        input_file, was_unzipped = self.gunzipInputFile(input_file=pathlib.Path(inData["input_file"]))
        logPath = self.getWorkingDirectory() / "Xtriage.log"

        from iotbx.reflection_file_reader import any_reflection_file
        from mmtbx.scaling.xtriage import master_params as xtriage_master_params
        from mmtbx.scaling.xtriage import xtriage_analyses

        xtriage_params = xtriage_master_params.fetch(sources=[]).extract()
        xtriage_params.scaling.input.xray_data.skip_sanity_checks = True

        try:
            reader = any_reflection_file(str(input_file))
        except Exception as e:
            logger.error(f"file could not be read: {e}")
            self.setFailure()

        if reader.file_type() != "ccp4_mtz":
            logger.error("input file type is not ccp4_mtz")
            self.setFailure()
            return
        logger.debug(f"Reading in {input_file}")
        arrays = reader.as_miller_arrays(merge_equivalents=False)
        for ma in arrays:
            if ma.info().labels == ["BATCH"]:
                batches = ma
            elif ma.info().labels == ["I", "SIGI"]:
                intensities = ma
            elif ma.info().labels == ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]:
                intensities = ma
            elif ma.info().labels == ["SCALEUSED"]:
                scales = ma
        if intensities is None or batches is None:
            logger.error("No intensities or batches found")
            self.setFailure()
            return

        mtz_object = reader.file_content()
        indices = mtz_object.extract_original_index_miller_indices()
        intensities = intensities.customized_copy(indices=indices, info=intensities.info())
        merged_intensities = intensities.merge_equivalents().array()

        try:
            xanalysis = xtriage_analyses(
                miller_obs=merged_intensities,
                unmerged_obs=intensities,
                text_out=open(logPath, "w"),
                params=xtriage_params,
            )
        except Exception as e:
            logger.error(f"xtriage analysis failed: {e}")

        issues = xanalysis.summarize_issues()
        xtriage_success = []
        xtriage_warnings = []
        xtriage_danger = []
        xs = io.StringIO()
        xout = _xtriage_output(xs)
        xanalysis.show(out=xout)
        xout.flush()
        sub_header_to_out = xout._sub_header_to_out

        for level, text, sub_header in issues._issues:
            summary = sub_header_to_out.get(sub_header, io.StringIO()).getvalue()
            summary = html.escape(summary)
            d = {"level": level, "text": text, "summary": summary, "header": sub_header}
            if level == 0:
                xtriage_success.append(d)
            elif level == 1:
                xtriage_warnings.append(d)
            elif level == 2:
                xtriage_danger.append(d)

        if was_unzipped:
            input_file.unlink(missing_ok=True)

        outData = {
            "logFile": str(logPath),
            "xtriage_success": xtriage_success,
            "xtriage_warnings": xtriage_warnings,
            "xtriage_danger": xtriage_danger,
        }

        return outData

    def gunzipInputFile(self, input_file: pathlib.Path) -> Tuple[pathlib.Path, bool]:
        """unzips input file and whether it was unzipped"""
        if ".gz" in input_file.suffixes:
            fIn = gzip.open(input_file, "rb")
            fileContent = fIn.read()
            fIn.close()
            strMtzFileName = input_file.stem
            self.strPathToLocalMtz = pathlib.Path(self.getWorkingDirectory() / strMtzFileName)
            logger.debug("unzipped .mtz.gz file")
            fOut = open(self.strPathToLocalMtz, "wb")
            fOut.write(fileContent)
            fOut.close()
            return (self.strPathToLocalMtz, True)
        else:
            return (input_file, False)


class _xtriage_output(printed_output):
    def __init__(self, out):
        super().__init__(out)
        self.gui_output = True
        self._out_orig = self.out
        self.out = io.StringIO()
        self._sub_header_to_out = {}

    def show_big_header(self, text):
        pass

    def show_header(self, text):
        self._out_orig.write(self.out.getvalue())
        self.out = io.StringIO()
        super().show_header(text)

    def show_sub_header(self, title):
        self._out_orig.write(self.out.getvalue())
        self.out = io.StringIO()
        self._current_sub_header = title
        assert title not in self._sub_header_to_out
        self._sub_header_to_out[title] = self.out

    def flush(self):
        self._out_orig.write(self.out.getvalue())
        self.out.flush()
        self._out_orig.flush()


class PhenixXTriageTask(AbstractTask):
    """
    This task runs phenix.xtriage
    """

    def run(self, inData):
        if os.environ.get("PHENIX", None) is None:
            commandLine = "source /mxn/groups/sw/mxsw/env_setup/phenix_env.sh \n"
        else:
            commandLine = ""
            logger.info(f"PHENIX version is {os.environ.get('PHENIX_VERSION', None)}")

        input_file, was_unzipped = self.gunzipInputFile(input_file=pathlib.Path(inData["input_file"]))
        commandLine += "phenix.xtriage "
        commandLine += str(input_file)
        commandLine += " obs=I,SIGI,merged "
        logPath = self.getWorkingDirectory() / "PhenixXtriage.log"
        self.runCommandLine(commandLine, logPath=logPath)
        if was_unzipped:
            input_file.unlink(missing_ok=True)

        outData = self.parseXtriageLogFile(logPath)

        return outData

    def gunzipInputFile(self, input_file: pathlib.Path) -> Tuple[pathlib.Path, bool]:
        """unzips input file and whether it was unzipped"""
        if ".gz" in input_file.suffixes:
            fIn = gzip.open(input_file, "rb")
            fileContent = fIn.read()
            fIn.close()
            strMtzFileName = input_file.stem
            self.strPathToLocalMtz = pathlib.Path(self.getWorkingDirectory() / strMtzFileName)
            fOut = open(self.strPathToLocalMtz, "wb")
            fOut.write(fileContent)
            fOut.close()
            return (self.strPathToLocalMtz, True)
        else:
            return (input_file, False)

    def parseXtriageLogFile(self, pathToLogFile: pathlib.Path):
        outData = {"logPath": str(pathToLogFile)}
        outData["hasTwinning"] = False
        outData["hasPseudotranslation"] = False
        outData["TwinLawsStatistics"] = []

        if pathToLogFile.is_file():
            with open(pathToLogFile, "r") as f:
                strLog = f.readlines()
            iIndex = 0
            listLines = [x.strip("\n") for x in strLog]
            bContinue = True
            while bContinue:
                if listLines[iIndex].startswith("Statistics depending on twin laws"):
                    # ------------------------------------------------------
                    # | Operator | type | R obs. | Britton alpha | H alpha |
                    # ------------------------------------------------------
                    # | k,h,-l   |  PM  | 0.025  | 0.458         | 0.478   |
                    # | -h,k,-l  |  PM  | 0.017  | 0.459         | 0.487   |
                    # ------------------------------------------------------
                    iIndex += 4
                    while not listLines[iIndex].startswith("---------"):
                        listLine = listLines[iIndex].split("|")
                        xsDataTwinLawsStatistics = {}
                        xsDataTwinLawsStatistics["operator"] = listLine[1].replace(" ", "")
                        xsDataTwinLawsStatistics["twinType"] = listLine[2].replace(" ", "")
                        xsDataTwinLawsStatistics["rObs"] = float(listLine[3])
                        xsDataTwinLawsStatistics["brittonAlpha"] = float(listLine[4])
                        xsDataTwinLawsStatistics["hAlpha"] = float(listLine[5])
                        xsDataTwinLawsStatistics["mlAlpha"] = float(listLine[6])
                        outData["TwinLawsStatistics"].append(xsDataTwinLawsStatistics)

                        iIndex += 1

                elif listLines[iIndex].startswith("Patterson analyses"):
                    # - Largest peak height   : 6.089
                    iIndex += 1
                    pattersonLargestPeakHeight = float(listLines[iIndex].split(":")[1])
                    outData["pattersonLargestPeakHeight"] = pattersonLargestPeakHeight
                    # (corresponding p value : 6.921e-01)
                    iIndex += 1
                    pattersonPValue = float(listLines[iIndex].split(":")[1].split(")")[0])
                    outData["pattersonPValue"] = pattersonPValue

                elif "indicating pseudo translational symmetry" in listLines[iIndex]:
                    #    The analyses of the Patterson function reveals a significant off-origin
                    #    peak that is 66.43 % of the origin peak, indicating pseudo translational symmetry.
                    #    The chance of finding a peak of this or larger height by random in a
                    #    structure without pseudo translational symmetry is equal to the 6.0553e-06.
                    #    The detected translational NCS is most likely also responsible for the elevated intensity ratio.
                    #    See the relevant section of the logfile for more details.
                    outData["hasPseudotranslation"] = True
                elif "As there are twin laws possible given the crystal symmetry, twinning could" in listLines[iIndex]:
                    #    The results of the L-test indicate that the intensity statistics
                    #    are significantly different than is expected from good to reasonable,
                    #    untwinned data.
                    #    As there are twin laws possible given the crystal symmetry, twinning could
                    #    be the reason for the departure of the intensity statistics from normality.
                    outData["hasTwinning"] = True
                iIndex += 1
                if iIndex == len(listLines):
                    bContinue = False
        return outData


class DistlSignalStrengthTask(AbstractTask):
    """
    This task runs the labelit.distl command for pre-screening reference images.
    """

    def run(self, inData):
        commandLine = "distl.signal_strength "
        commandLine += inData["referenceImage"]
        logPath = self.getWorkingDirectory() / "distl.log"
        self.runCommandLine(commandLine, logPath=logPath)
        with open(str(logPath)) as f:
            logText = f.read()
        imageQualityIndicators = self.parseLabelitDistlOutput(logText)
        outData = {"imageQualityIndicators": imageQualityIndicators}
        return outData

    def parseLabelitDistlOutput(self, logText):
        imageQualityIndicators = {}
        for line in logText.split("\n"):
            if line.find("Spot Total") != -1:
                spotTotal = int(line.split()[3])
                imageQualityIndicators["spotTotal"] = spotTotal
            elif line.find("In-Resolution Total") != -1:
                inResTotal = int(line.split()[3])
                imageQualityIndicators["inResTotal"] = inResTotal
            elif line.find("Good Bragg Candidates") != -1:
                goodBraggCandidates = int(line.split()[4])
                imageQualityIndicators["goodBraggCandidates"] = goodBraggCandidates
            elif line.find("Ice Rings") != -1:
                iceRings = int(line.split()[3])
                imageQualityIndicators["iceRings"] = iceRings
            elif line.find("Method 1 Resolution") != -1:
                method1Res = float(line.split()[4])
                imageQualityIndicators["method1Res"] = method1Res
            elif line.find("Method 2 Resolution") != -1:
                if line.split()[4] != "None":
                    fMethod2Res = float(line.split()[4])
                    imageQualityIndicators["method2Res"] = fMethod2Res
            elif line.find("Maximum unit cell") != -1:
                if line.split()[4] != "None":
                    fMaxUnitCell = float(line.split()[4])
                    imageQualityIndicators["maxUnitCell"] = fMaxUnitCell
            elif line.find("%Saturation, Top 50 Peaks") != -1:
                pctSaturationTop50Peaks = float(line.split()[5])
                imageQualityIndicators["pctSaturationTop50Peaks"] = pctSaturationTop50Peaks
            elif line.find("In-Resolution Ovrld Spots") != -1:
                inResolutionOvrlSpots = int(line.split()[4])
                imageQualityIndicators["inResolutionOvrlSpots"] = inResolutionOvrlSpots
            elif line.find("Bin population cutoff for method 2 resolution") != -1:
                binPopCutOffMethod2Res = float(line.split()[7][:-1])
                imageQualityIndicators["binPopCutOffMethod2Res"] = binPopCutOffMethod2Res
            elif (
                line.find(
                    "Total integrated signal, pixel-ADC units above local background (just the good Bragg candidates)"
                )
                != -1
            ):
                totalIntegratedSignal = float(line.split()[-1])
                imageQualityIndicators["totalIntegratedSignal"] = totalIntegratedSignal
            elif line.find("signals range from") != -1:
                listStrLine = line.split()
                signalRangeMin = float(listStrLine[3])
                signalRangeMax = float(listStrLine[5])
                signalRangeAverage = float(listStrLine[-1])
                imageQualityIndicators["signalRangeMin"] = signalRangeMin
                imageQualityIndicators["signalRangeMax"] = signalRangeMax
                imageQualityIndicators["signalRangeAverage"] = signalRangeAverage
            elif line.find("Saturations range from") != -1:
                listStrLine = line.split()
                saturationRangeMin = float(listStrLine[3][:-1])
                saturationRangeMax = float(listStrLine[5][:-1])
                saturationRangeAverage = float(listStrLine[-1][:-1])
                imageQualityIndicators["saturationRangeMin"] = saturationRangeMin
                imageQualityIndicators["saturationRangeMax"] = saturationRangeMax
                imageQualityIndicators["saturationRangeAverage"] = saturationRangeAverage
        return imageQualityIndicators


class PhenixProcessPredictedModelTask(AbstractTask):
    """
    This task runs phenix.process_predicted_model to replace B-factors and (optinally) break the model into domains
    """

    def run(self, inData):
        
        if os.environ.get("PHENIX", None) is None:
            phenixSetup = UtilsConfig.get("Phenix", "phenix_env")
            logger.debug(f"Phenix Setup: {phenixSetup}")
            if phenixSetup is None:
                logger.warning("Phenix setup not found!")
                commandLine = ""
            else:
                logger.info(f"PHENIX version is {os.environ.get('PHENIX_VERSION', None)}")
                commandLine = "source " + phenixSetup + "\n"
        else:
            commandLine = ""
            logger.info(f"PHENIX version is {os.environ.get('PHENIX_VERSION', None)}")

        commandLine += "phenix.process_predicted_model "
        commandLine += inData["PDB_file"]
        # commandLine += ' '

        logPath = self.getWorkingDirectory() / "PhenixProcPM.log"
        self.runCommandLine(commandLine, logPath=logPath)
        if logPath.exists():
            with open(str(logPath)) as f:
                logText = f.read()
        else:
            logger.error(f"Log file {logPath} does not exist")

        logger.info("Command line: {0}".format(commandLine))

        outData = self.parseProcessPredictedModel(logText)

        return outData

    def parseProcessPredictedModel(self, logText):
        outData = {}

        for line in logText.split("\n"):
            if line.find("Working directory:") != -1:
                outData["workingDirectory"] = line.split(":")[1].strip()
            elif line.find("Found model,") != -1:
                outData["processedModel"] = line.split(",")[1].strip()
            elif line.find("B-value field") != -1:
                outData["B-value"] = line.split(" ")[2:]
            elif line.find("Maximum B-Value") != -1:
                outData["maximumB-value"] = line.split(":")[1]
            elif line.find("Maximum rmsd") != -1:
                outData["maximumRmsd"] = line.split(" ")[3:4]
            elif line.find("Clusters:") != -1:
                outData["Clusters"] = line.split(":")[1]
            elif line.find("Threshold:") != -1:
                outData["Threshold"] = line.split(":")[1]
            elif line.find("Job complete") != -1:
                outData["jobComplete"] = True

        return outData
