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

__authors__ = ["A. Finke"]
__license__ = "MIT"
__date__ = "20/01/2024"

from jinja2 import Environment, PackageLoader
from cctbx import sgtbx, uctbx
import xmltodict
import html
import pandas as pd
import numpy as np
from pathlib import Path
from edna2.utils import UtilsLogging
import time

from edna2.tasks.AbstractTask import AbstractTask
logger = UtilsLogging.getLogger()

pd.options.mode.copy_on_write = True
SQRT2 = 1.414213


class Edna2ReportTask(AbstractTask):
    """
    Generates an HTML report from XDS,AIMLESS,CTRUNCATE,xtriage
    output. Template stolen without remorse from xia2/dials.
    Requirements:
    AIMLESS XML output
    POINTLESS XML output
    CTRUNCATE log file
    INTEGRATE.LP from XDS
    xdsstat output
    xtriage output
    """

    def run(self, inData):
        t1 = time.perf_counter()
        loader = PackageLoader("edna2", "templates")
        env = Environment(loader=loader)
        template = env.get_template("edna2.html")

        styles = {}
        self.datasetName = inData.get("datasetName","DEFAULT")
        pointlessLog = inData["pointlessLog"]
        aimlessXml = inData["aimlessXml"]
        aimlessLog = inData["aimlessLog"]
        cTruncateLog = inData["cTruncateLog"]
        integrateLp = inData["integrateLp"]
        xtriageOutput = inData["xtriageOutput"]
        pointlessOutput = inData["pointlessOutput"]
        xdsstatLp = inData["xdsstatLp"]
        XdsIndexingLp = inData["XdsIndexingLp"]
        XdsCorrectLp = inData["XdsCorrectLp"]
        files = inData["files"]

        refs = inData.get("refsList")
        xdsIndexOutput = self.generateHtmlFromText(XdsIndexingLp)
        xdsIntegrateOutput = self.generateHtmlFromText(integrateLp)
        xdsCorrectOutput = self.generateHtmlFromText(XdsCorrectLp)
        pointlessHtml = self.generateHtmlFromText(pointlessLog)
        aimlessOutput = self.generateHtmlFromText(aimlessLog)
        cTruncateOutput = self.generateHtmlFromText(cTruncateLog)
        try:
            aimlessResults = xmltodict.parse(open(aimlessXml).read())
        except Exception as e:
            logger.error(f"aimless xml could not be parsed: {e}")
            self.setFailure()
            return
        try:
            truncateDfDict, wilsonBFactor = self.generateTruncateDf(cTruncateLog)
        except Exception as e:
            logger.error(f"aimless xml could not be parsed: {e}")
            self.setFailure()
            return

        aimlessDf = self.generateAimlessDf(aimlessResults)
        overallStatsTable, spaceGroupString, unitCellString = self.generateOverallStatsTable(aimlessResults)
        resShells = self.resShellsFromAimless(aimlessDf)
        overallResultsTable = self.generateOverallResultsTableFromAimless(aimlessResults)

        analysisByResolutionPlots = self.analysisByResolutionPlotsPlotly(aimlessDf, truncateDfDict, wilsonBFactor)
        batchPlotDfs = self.batchPlotDfs(integrateLp, xdsstatLp)
        batchResultsPlotly = self.batchResultsPlotly(batchPlotDfs, aimlessResults)
        miscellaneousResultsPlotly = self.miscellaneous_plots_plotly(truncateDfDict)

        refs = [
            {
                "acta": "Winn, M. D. et al. (2011) Acta Cryst. D67, 235-242.",
                "url": "https://doi.org/10.1107/S0907444910045749",
            },
            {
                "acta": "Evans, P. R. and Murshudov, G. N. (2013) Acta Cryst. D69, 1204-1214.",
                "url": "http://journals.iucr.org/d/issues/2013/07/00/ba5190/index.html",
            },
            {"acta": "Evans, P. (2006) Acta Cryst. D62, 72-82.", "url": "https://doi.org/10.1107/S0907444905036693"},
            {"acta": "Kabsch, W. (2010) Acta Cryst. D66, 125-132.", "url": "https://doi.org/10.1107/S0907444909047337"},
        ]

        references = {cdict["acta"]: cdict.get("url") for cdict in refs}
        if pointlessOutput["alternativeSpacegroups"]:
            alternative_space_groups = ", ".join(x for x in pointlessOutput["alternativeSpacegroups"])
        else:
            alternative_space_groups = None
        html_source = template.render(
            page_title="edna2proc processing report",
            wname=datasetName,
            xdsIndexOutput=xdsIndexOutput,
            xdsIntegrateOutput=xdsIntegrateOutput,
            xdsCorrectOutput=xdsCorrectOutput,
            pointlessOutput=pointlessHtml,
            aimlessOutput=aimlessOutput,
            cTruncateOutput=cTruncateOutput,
            space_group=spaceGroupString,
            alternative_space_groups=alternative_space_groups,
            unit_cell=unitCellString,
            overall_stats_table=overallStatsTable,
            detailed_stats_table=overallResultsTable,
            resShells=resShells,
            mtz_files=files["mtz_files"],
            unmerged_files=files["unmerged_files"],
            other_files=files["other_files"],
            log_files_table=files["log_files"],
            ccHalfPlot=analysisByResolutionPlots,
            batchplots=batchResultsPlotly,
            xtriage_results=xtriageOutput,
            miscellaneous=miscellaneousResultsPlotly,
            references=references,
            styles=styles,
        )
        try:
            with open(Path(self.getWorkingDirectory()/"edna2proc_summary.html"), "wb") as f:
                f.write(html_source.encode("utf-8", "xmlcharrefreplace"))
        except Exception as e:
            logger.error(f"Couldn't write html file: {e}")
            return
        
        logger.info(f"process took {time.perf_counter() - t1:0.1f} s")

        outData = {
            "htmlFile": str(Path(self.getWorkingDirectory()/"edna2proc_summary.html"))
        }

        return outData


    @staticmethod
    def generateHtmlFromText(textLog):
        try:
            with open(textLog, "r") as fp:
                output = html.escape(fp.read())
        except Exception as e:
            logger.error(f"could not parse log to html: {e}")
            return
        return output


    def generateOverallStatsTable(self, aimlessResults):
        wavelengthString = str(float(aimlessResults["AIMLESS"]["ReflectionData"]["Dataset"]["Wavelength"]))

        spaceGroupName = aimlessResults["AIMLESS"]["Result"]["Dataset"]["SpacegroupName"]
        spaceGroupString = sgtbx.space_group_info(symbol=spaceGroupName).symbol_and_number()

        unitCell = aimlessResults["AIMLESS"]["Result"]["Dataset"]["cell"]
        unitCellString = "a = {a}, b = {b}, c = {c} , \u03B1 = {alpha}, \u03B2 = {beta}, \u03B3 = {gamma}".format(**unitCell)

        resLow = aimlessResults["AIMLESS"]["Result"]["Dataset"]["ResolutionLow"]
        resHigh = aimlessResults["AIMLESS"]["Result"]["Dataset"]["ResolutionHigh"]
        resolutionString = "{0} - {1} ({2} - {1})".format(resLow["Overall"], resHigh["Overall"], resLow["Outer"])

        completeness = aimlessResults["AIMLESS"]["Result"]["Dataset"]["Completeness"]
        completenessString = "{0} ({1})".format(completeness["Overall"], completeness["Outer"])

        multiplicity = aimlessResults["AIMLESS"]["Result"]["Dataset"]["Multiplicity"]
        multiplicityString = "{0} ({1})".format(multiplicity["Overall"], multiplicity["Outer"])

        CChalf = aimlessResults["AIMLESS"]["Result"]["Dataset"]["CChalf"]
        CChalfString = "{0} ({1})".format(CChalf["Overall"], CChalf["Outer"])

        IsigI = aimlessResults["AIMLESS"]["Result"]["Dataset"]["MeanIoverSD"]
        IsigIString = "{0} ({1})".format(IsigI["Overall"], IsigI["Outer"])

        rMerge = aimlessResults["AIMLESS"]["Result"]["Dataset"]["RmergeOverall"]
        rMergeString = "{0} ({1})".format(rMerge["Overall"], rMerge["Outer"])

        anomCompleteness = aimlessResults["AIMLESS"]["Result"]["Dataset"]["AnomalousCompleteness"]
        anomCompletenessString = "{0} ({1})".format(anomCompleteness["Overall"], anomCompleteness["Outer"])

        anomMultiplicity = aimlessResults["AIMLESS"]["Result"]["Dataset"]["AnomalousMultiplicity"]
        anomMultiplicityString = "{0} ({1})".format(anomMultiplicity["Overall"], anomMultiplicity["Outer"])

        columns = []
        columns.append(
            [
                "",
                "Wavelength (Å)",
                "Resolution range (Å)",
                "Completeness (%)",
                "Multiplicity",
                "CC-half",
                "I/sigma",
                "R<sub>merge</sub>(I)",
                # anomalous statistics
                "Anomalous completeness (%)",
                "Anomalous multiplicity",
            ]
        )
        columns.append(
            [
                self.datasetName,
                wavelengthString,
                resolutionString,
                completenessString,
                multiplicityString,
                CChalfString,
                IsigIString,
                rMergeString,
                anomCompletenessString,
                anomMultiplicityString,
            ],
        )
        overall_stats_table = [[c[i] for c in columns] for i in range(len(columns[0]))]
        return overall_stats_table, spaceGroupString, unitCellString

    @staticmethod
    def generateAimlessDf(aimlessResults):
        """
        generates pandas dataframe for output
        from aimless dict results
        """
        GraphStatsVsResolution = next(
            item for item in aimlessResults["AIMLESS"]["CCP4Table"] if item.get("@id") == "Graph-StatsVsResolution"
        )
        labels = GraphStatsVsResolution["headers"]["#text"].split()

        StatsVsResolution = np.genfromtxt(GraphStatsVsResolution["data"].split("\n"))
        df = pd.DataFrame(StatsVsResolution, columns=labels)
        ComplenessVsResolution = next(
            item
            for item in aimlessResults["AIMLESS"]["CCP4Table"]
            if item.get("@id") == "Graph-CompletenessVsResolution"
        )
        labels = ComplenessVsResolution["headers"]["#text"]
        labels = ComplenessVsResolution["headers"]["#text"].split()
        complenessVsResolution = np.genfromtxt(ComplenessVsResolution["data"].split("\n"))
        df1 = pd.DataFrame(complenessVsResolution, columns=labels)
        df = df.drop("Nmeas", axis=1)
        df = df.join(df1[["Nmeas", "Nref", "%poss", "Mlplct", "AnoCmp", "AnoMlt"]])
        CChalf = next(item for item in aimlessResults["AIMLESS"]["CCP4Table"] if item.get("@id") == "Graph-CChalf")
        labels = CChalf["headers"]["#text"].split()
        CChalfData = np.genfromtxt(CChalf["data"].split("\n"))
        df2 = pd.DataFrame(CChalfData, columns=labels)
        df = df.join(df2[["CC1/2", "CCanom"]])
        resLow = aimlessResults["AIMLESS"]["Result"]["Dataset"]["ResolutionLow"]
        resLow = float(resLow["Overall"])
        resLow

        df["dlow"] = df["Dmid"]
        df["dhigh"] = df["Dmid"].shift(1)
        df.at[0, "dhigh"] = resLow

        df["Dres"] = pd.Series(
            list(map(lambda x, y: f"{float(y):.2f}" + " - " + f"{float(x):.2f}", df["dlow"], df["dhigh"]))
        )

        return df

    @staticmethod
    def generateTruncateDf(cTruncateLog):
        wilsonData = []
        completenessData = []
        cumIntensDist = []
        twinLtest = []
        with open(cTruncateLog, "r") as fp:
            for line in fp:
                if "$TABLE: Intensity Completeness analysis" in line:
                    for _ in range(3):
                        line = next(fp)
                    labels_comp = line.strip("$").strip("/n").split()
                    for _ in range(2):
                        line = next(fp)
                    while "$$" not in line:
                        completenessData.append(line.strip("\n"))
                        line = next(fp)
                if "$TABLE: Wilson" in line:
                    line = next(fp)
                    try:
                        wilsonBFactor = line.split("=")[1].split(":")[0].strip()
                    except:
                        wilsonBFactor = "N/A"
                    line = next(fp)
                    labels_wilson = line.replace("$", "").strip("/n").split()
                    line = next(fp)
                    line = next(fp)
                    while "$$" not in line:
                        wilsonData.append(line.strip("\n"))
                        line = next(fp)
                if "$TABLE: Cumulative intensity distribution" in line:
                    for _ in range(2):
                        line = next(fp)
                    labels_cumul = line.replace("$", "").strip("/n").split()
                    for _ in range(2):
                        line = next(fp)
                    while "$$" not in line:
                        cumIntensDist.append(line.strip("\n"))
                        line = next(fp)
                if "$TABLE: L test for twinning" in line:
                    for _ in range(2):
                        line = next(fp)
                    labels_lTest = line.replace("$", "").strip("/n").split()
                    for _ in range(2):
                        line = next(fp)
                    while "$$" not in line:
                        twinLtest.append(line.strip("\n"))
                        line = next(fp)

        completenessDf = pd.DataFrame(np.genfromtxt(completenessData), columns=labels_comp)
        cumIntensDistDf = pd.DataFrame(np.genfromtxt(cumIntensDist), columns=labels_cumul)

        twinLtestDf = pd.DataFrame(np.genfromtxt(twinLtest), columns=labels_lTest)
        wilsonDf = pd.DataFrame(np.genfromtxt(wilsonData), columns=labels_wilson)
        dfDict = {
            "wilson": wilsonDf,
            "completeness": completenessDf,
            "cumIntensDist": cumIntensDistDf,
            "twinLtest": twinLtestDf,
        }
        return dfDict, wilsonBFactor

    @staticmethod
    def resShellsFromAimless(aimlessDf):
        df2 = aimlessDf[
            ["Dres", "Nmeas", "Nref", "Mlplct", "%poss", "AvI", "Mn(I/sd)", "Rmrg", "Rmeas", "Rpim", "CC1/2", "CCanom"]
        ]
        df2["Nmeas"] = df2["Nmeas"].astype("int32")
        df2["Nref"] = df2["Nref"].astype("int32")

        df2 = df2.rename(
            columns={
                "Dres": "Resolution (Å)",
                "Nmeas": "N(obs)",
                "Nref": "N(unique)",
                "Mlplct": "Multiplicity",
                "%poss": "Completeness (%)",
                "AvI": "Mean I",
                "Mn(I/sd)": "Mean I/\u03C3",
                "Rmrg": "R<sub>merge</sub>",
                "Rmeas": "R<sub>meas</sub>",
                "Rpim": "R<sub>pim</sub>",
                "CCanom": "CC<sub>ano</sub>",
                "CC1/2": "CC<sub>1/2</sub>",
            },
        )
        labels = list(df2.columns)
        resShells = df2.values.tolist()
        resShells.insert(0, labels)

        return resShells

    @staticmethod
    def generateOverallResultsTableFromAimless(aimlessDf):
        ds = aimlessDf["AIMLESS"]["Result"]["Dataset"]

        overall = []
        inner = []
        outer = []
        for k in [
            "ResolutionLow",
            "ResolutionHigh",
            "NumberObservations",
            "NumberReflections",
            "Multiplicity",
            "Completeness",
            "MeanIoverSD",
            "RmergeOverall",
            "RmeasOverall",
            "RpimOverall",
            "CChalf",
        ]:
            overall.append(ds[k]["Overall"])
            inner.append(ds[k]["Inner"])
            outer.append(ds[k]["Outer"])
        overall.insert(0, f"{overall.pop(0)} - {overall.pop(0)}")
        overall.insert(0, "Overall")
        outer.insert(0, f"{outer.pop(0)} - {outer.pop(0)}")
        outer.insert(0, "High Resolution")
        inner.insert(0, f"{inner.pop(0)} - {inner.pop(0)}")
        inner.insert(0, "Low Resolution")

        labels = [
            "",
            "Resolution (Å)",
            "Observations",
            "Unique reflections",
            "Multiplicity",
            "Completeness (%)",
            "Mean I/\u03C3(I)",
            "R<sub>merge</sub>",
            "R<sub>meas</sub>",
            "R<sub>pim</sub>",
            "CC<sub>1/2</sub>",
        ]

        columns = [labels, overall, inner, outer]
        columns = list(map(list, zip(*columns)))
        return columns

    @staticmethod
    def analysisByResolutionPlotsPlotly(aimlessDf, truncateDfDict, wilsonBFactor):
        tickvals, ticktext = Edna2ReportTask.d_star_sq_to_d_ticks(list(aimlessDf["1/d^2"]), nticks=8)

        return {
            "cc_one_half": {
                "data": [
                    {
                        "x": list(aimlessDf["1/d^2"]),
                        "y": list(aimlessDf["CC1/2"]),
                        "type": "scatter",
                        "name": "CC<sub>\u00bd</sub>",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                ],
                "layout": {
                    "title": "CC<sub>\u00bd</sub> vs resolution<br><sup>data from AIMLESS</sup>",
                    "xaxis": {
                        "title": "Resolution (\u00c5)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"title": "CC<sub>\u00bd</sub>", "rangemode": "tozero"},
                },
                "help": "The correlation coefficients, CC<sub>\u00bd</sub>, between random half-datasets. A correlation\ncoefficient of +1 indicates good correlation, and 0 indicates no correlation.\nCC<sub>\u00bd</sub> is typically close to 1 at low resolution, falling off to close to zero at\nhigher resolution. A typical resolution cutoff based on CC<sub>\u00bd</sub> is around 0.3-0.5.\n\n[1] Karplus, P. A., & Diederichs, K. (2012). Science, 336(6084), 1030-1033.\n    https://doi.org/10.1126/science.1218231\n[2] Diederichs, K., & Karplus, P. A. (2013). Acta Cryst D, 69(7), 1215-1222.\n    https://doi.org/10.1107/S0907444913001121\n[3] Evans, P. R., & Murshudov, G. N. (2013). Acta Cryst D, 69(7), 1204-1214.\n    https://doi.org/10.1107/S0907444913000061\n",
            },
            "i_over_sigi": {
                "data": [
                    {
                        "x": list(aimlessDf["1/d^2"]),
                        "y": list(aimlessDf["Mn(I/sd)"]),
                        "type": "scatter",
                        "name": "Mean I/σ(I)",
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                ],
                "layout": {
                    "title": "Mean I/\u03C3(I) vs resolution<br><sup>data from AIMLESS</sup>",
                    "xaxis": {
                        "title": "Resolution (\u00c5)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"title": "Mean I/\u03C3(I)", "rangemode": "tozero"},
                },
            },
            "completeness_vs_resolution": {
                "data": [
                    {
                        "x": list(aimlessDf["1/d^2"]),
                        "y": list(aimlessDf["%poss"]),
                        "type": "scatter",
                        "name": "Completeness",
                        "mode": "lines",
                    },
                    {
                        "x": list(aimlessDf["1/d^2"]),
                        "y": list(aimlessDf["AnoCmp"]),
                        "type": "scatter",
                        "name": "Anomalous Completeness",
                        "mode": "lines",
                    },
                ],
                "layout": {
                    "title": "Completeness vs resolution<br><sup>data from AIMLESS</sup>",
                    "xaxis": {
                        "title": "Resolution (\u00c5)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"title": "Completeness (%)", "rangemode": "tozero"},
                },
            },
            "multiplicity_vs_resolution": {
                "data": [
                    {
                        "x": list(aimlessDf["1/d^2"]),
                        "y": list(aimlessDf["Mlplct"]),
                        "type": "scatter",
                        "name": "Multiplicity",
                        "mode": "lines",
                    },
                    {
                        "x": list(aimlessDf["1/d^2"]),
                        "y": list(aimlessDf["AnoMlt"]),
                        "type": "scatter",
                        "name": "Anomalous Multiplicity",
                        "mode": "lines",
                    },
                ],
                "layout": {
                    "title": "Multiplicity vs resolution<br><sup>data from AIMLESS</sup>",
                    "xaxis": {
                        "title": "Resolution (\u00c5)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"title": "Multiplicity", "rangemode": "tozero"},
                },
            },
            "wilson_plot": {
                "data": [
                    {
                        "x": list(truncateDfDict["wilson"]["1/resol^2"]),
                        "y": list(truncateDfDict["wilson"]["ln(I/I_th)"]),
                        "type": "scatter",
                        "name": "Observed",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["wilson"]["1/resol^2"]),
                        "y": list(truncateDfDict["wilson"]["Reference_prot"]),
                        "type": "scatter",
                        "name": "Reference",
                        "mode": "lines",
                    },
                ],
                "layout": {
                    "title": f"Wilson Plot - estimated B-factor = {wilsonBFactor}<br><sup>data from CTRUNATE</sup>",
                    "xaxis": {
                        "title": "Resolution (\u00c5)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"title": "ln(I/I_th)"},
                },
                "help": """
    The isotropic Wilson temperature estimate (B-value) is an approximation to the 
    fall-off of scattering with resolution. This should be correlated with the refined 
    atomic B-values.  This averages out any anisotropy in the experimental observations.
    This approximation will be misleading for strongly anisotropic data.

    Computed using Popov & Bourenkov, Acta D (2003) D59, 1145

    The Wilson plot shows the fall off of the mean intensity with resolution. This is then
    used to calculate an absolute scale and temperature factor for a set of observed 
    intensities, using the theory of AC Wilson.  The reference plot is based upon an analysis
    of high resolution datasets in the PDB (BEST), which takes into account the nonrandom 
    distribution of atoms within the crystal. Some deviation from the reference plot is to 
    be expected, however, significant deviation may indicate problems, such as ice rings, 
    detector issues, or misprocessed data.""",
            },
            "completeness_vs_iSigIThresholds": {
                "data": [
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["Completeness"]),
                        "type": "scatter",
                        "name": "All Data",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["(I/s>15)"]),
                        "type": "scatter",
                        "name": "(I/\u03C3(I)>15)",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["(I/s>10)"]),
                        "type": "scatter",
                        "name": "(I/\u03C3(I)>10)",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["(I/s>5)"]),
                        "type": "scatter",
                        "name": "(I/\u03C3(I)>5)",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["(I/s>3)"]),
                        "type": "scatter",
                        "name": "(I/\u03C3(I)>3)",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["(I/s>2)"]),
                        "type": "scatter",
                        "name": "(I/\u03C3(I)>2)",
                        "mode": "lines",
                    },
                    {
                        "x": list(truncateDfDict["completeness"]["1/resol^2"]),
                        "y": list(truncateDfDict["completeness"]["(I/s>1)"]),
                        "type": "scatter",
                        "name": "(I/\u03C3(I)>1)",
                        "mode": "lines",
                    },
                ],
                "layout": {
                    "title": f"Intensity-Completeness Analysis<br><sup>data from CTRUNCATE</sup>",
                    "xaxis": {
                        "title": "Resolution (\u00c5)",
                        "tickvals": tickvals,
                        "ticktext": ticktext,
                    },
                    "yaxis": {"title": "Completeness", "rangemode": "tozero"},
                },
                "help": """
    The completeness at various resolution limit plots gives the completeness after applying a I/sigI
    cutoff.  The profiles give an indication of the quality of the data. 
    """,
            },
        }

    @staticmethod
    def d_star_sq_to_d_ticks(d_star_sq, nticks):
        min_d_star_sq = min(d_star_sq)
        dstep = (max(d_star_sq) - min_d_star_sq) / nticks
        tickvals = [min_d_star_sq + (i * dstep) for i in range(nticks)]
        ticktext = [f"{uctbx.d_star_sq_as_d(dsq):.2f}" for dsq in tickvals]
        return tickvals, ticktext

    @staticmethod
    def aimless_rcpDf(aimlessResults):
        GraphRaddamAnalysis = next(
            item for item in aimlessResults["AIMLESS"]["CCP4Table"] if item.get("@id") == "Graph-RadiationDamageAnalysis"
        )
        labels = GraphRaddamAnalysis["headers"]["#text"].split()
        RaddamAnalysis = np.genfromtxt(GraphRaddamAnalysis["data"].split("\n"))
        df = pd.DataFrame(RaddamAnalysis, columns=labels)
        return df

    @staticmethod
    def integrateDf(integrate_lp):
        bloop = []
        with open(integrate_lp, "r") as fp:
            for line in fp:
                if "IMAGE IER  SCALE     NBKG NOVL NEWALD NSTRONG  NREJ   SIGMAB   SIGMAR" in line:
                    line = next(fp)
                    linelen = len(line)
                    while len(line) == linelen:
                        bloop.append(line.strip())
                        line = next(fp)
        array = np.array([x.split() for x in bloop], dtype=float)
        integrate_labels = ["IMAGE", "IER", "SCALE", "NBKG", "NOVL", "NEWALD", "NSTRONG", "NREJ", "SIGMAB", "SIGMAR"]
        integrate_df = pd.DataFrame(array, columns=integrate_labels)

        return integrate_df

    @staticmethod
    def batchPlotDfs(integrateLp, xdsstatLp):
        integrate_df = Edna2ReportTask.integrateDf(integrateLp)

        xdsstatdata = []
        xdsstat_differencedata = []
        with open(xdsstatLp, "r") as fp:
            for line in fp:
                if "L\n" in line:
                    xdsstatdata.append(line.strip())
                if "DIFFERENCE\n" in line:
                    xdsstat_differencedata.append(line.strip())
        xdsstatlabels = [
            "frame",
            "nref",
            "nmisfits",
            "I",
            "sigmaI",
            "I/sigmaI",
            "fracReflObs",
            "correl",
            "Rmeas",
            "NrefusedforRMeas",
            "nuniqrefl",
            "l",
        ]
        xdsstat_difflabels = [
            "Framediff",
            "#refs",
            "R_d",
            "n-notfriedel",
            "Rd-notfriedel",
            "n-friedel",
            "Rd-friedel",
            "dummy",
        ]
        xdsstatdiff = pd.DataFrame(np.genfromtxt(xdsstat_differencedata), columns=xdsstat_difflabels)

        xdsstat_df = pd.DataFrame(np.genfromtxt(xdsstatdata), columns=xdsstatlabels)
        xdsstat_df["frame"] = xdsstat_df["frame"].astype(int)
        integrate_df = integrate_df.join(xdsstat_df, how="left")
        integrate_df = integrate_df.join(xdsstatdiff[["Framediff", "R_d"]], how="left")
        return integrate_df
    
    @staticmethod
    def batchResultsPlotly(df, aimlessResults):
        df1 = Edna2ReportTask.aimless_rcpDf(aimlessResults)

        # linear regression for Rcp Analysis
        x = np.array(df["Framediff"])
        x = x[~np.isnan(x)]
        y = np.array(df["R_d"])
        y = y[~np.isnan(y)]
        A = np.vstack([x, np.ones(len(x))]).T
        m, c = np.linalg.lstsq(A, y, rcond=None)[0]
        df["y_pred"] = c + m * df["Framediff"]
        limit = df["R_d"].iloc[0] * SQRT2

        return {
            "scale_rmerge_batch": {
                "data": [
                    {
                        "x": list(df["frame"]),
                        "y": list(df["SCALE"]),
                        "type": "scatter",
                        "name": "Scale",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                        "text": [f"#refl:{x}" for x in list(df["nref"])],
                    },
                    {
                        "x": list(df["frame"]),
                        "y": list(df["Rmeas"]),
                        "type": "scatter",
                        "name": "Rmeas",
                        "opacity": 0.8,
                        "mode": "lines",
                        "yaxis": "y2",
                        "line": {
                            "color": "rgb(255, 127, 14)",
                        },
                        "text": [f"#refl:{x}" for x in list(df["NrefusedforRMeas"])],
                    },
                ],
                "layout": {
                    "title": f"Scale and R<sub>meas</sub> vs frame<br><sub>data from XDS</sub>",
                    "xaxis": {"title": "frame"},
                    "yaxis": {"title": "Scale"},
                    "yaxis2": {"title": "Rmeas", "overlaying": "y", "rangemode": "tozero", "side": "right"},
                },
            },
            "IsigI_vs_batch": {
                "data": [
                    {
                        "x": list(df["frame"]),
                        "y": list(df["I/sigmaI"]),
                        "type": "scatter",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                ],
                "layout": {
                    "title": f"\u2039I/\u03C3(I)\u203A by frame<br><sub>data from XDS</sub>",
                    "xaxis": {"title": "frame"},
                    "yaxis": {"title": "\u2039I/\u03C3(I)\u203A"},
                },
            },
            "Rd_vs_batchdiff": {
                "data": [
                    {
                        "x": list(df["Framediff"]),
                        "y": list(df["R_d"]),
                        "type": "scatter",
                        "name": "Rd",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                    {
                        "x": list(df["Framediff"]),
                        "y": [limit for _ in range(len(df["Framediff"]))],
                        "type": "scatter",
                        "name": "Start value * \u221A(2)",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                    },
                    {
                        "x": list(df["Framediff"]),
                        "y": list(df["y_pred"]),
                        "type": "scatter",
                        "name": "linear regression",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(26,255,26)", "dash": "dot"},
                    },
                ],
                "layout": {
                    "title": f"Rd vs frame difference<br><sub>data from XDS</sub>",
                    "xaxis": {"title": "frame difference"},
                    "yaxis": {"title": "Rd"},
                    "showlegend": True,
                },
            },
            "completeness_vs_dose": {
                "data": [
                    {
                        "x": list(df1["Batch"]),
                        "y": list(df1["CmPoss"]),
                        "type": "scatter",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                ],
                "layout": {
                    "title": f"Completeness vs dose<br><sub>data from AIMLESS</sub>",
                    "xaxis": {"title": "Dose"},
                    "yaxis": {"title": "Completeness"},
                },
            },
            "Rcp_vs_dose": {
                "data": [
                    {
                        "x": list(df1["Batch"]),
                        "y": list(df1["Rcp"]),
                        "type": "scatter",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                ],
                "layout": {
                    "title": f"Rcp vs dose<br><sub>data from AIMLESS</sub>",
                    "xaxis": {"title": "Dose"},
                    "yaxis": {"title": "Rcp"},
                },
            },
        }

    @staticmethod
    def miscellaneous_plots_plotly(truncateDfDict):
        return {
            "cumIntensDist": {
                "data": [
                    {
                        "x": list(truncateDfDict["cumIntensDist"]["Z"]),
                        "y": list(truncateDfDict["cumIntensDist"]["Acent_theor"]),
                        "type": "scatter",
                        "name": "Acentric theory",
                        "opacity": 0.8,
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                    },
                    {
                        "x": list(truncateDfDict["cumIntensDist"]["Z"]),
                        "y": list(truncateDfDict["cumIntensDist"]["Acent_obser"]),
                        "type": "scatter",
                        "name": "Acentric observed",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)"},
                    },
                    {
                        "x": list(truncateDfDict["cumIntensDist"]["Z"]),
                        "y": list(truncateDfDict["cumIntensDist"]["Cent_theor"]),
                        "type": "scatter",
                        "name": "Centric theory",
                        "opacity": 0.8,
                        "line": {"color": "rgb(255, 127, 14)", "dash": "dot"},
                    },
                    {
                        "x": list(truncateDfDict["cumIntensDist"]["Z"]),
                        "y": list(truncateDfDict["cumIntensDist"]["Cent_obser"]),
                        "type": "scatter",
                        "name": "Centric observed",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)"},
                    },
                ],
                "layout": {
                    "title": f"Cumulative Intensity Distribution<br><sup>data from CTRUNCATE</sup>",
                    "xaxis": {
                        "title": "z",
                    },
                    "yaxis": {"title": "P(Z <= Z)"},
                },
            },
            "ltest": {
                "data": [
                    {
                        "x": list(truncateDfDict["twinLtest"]["|L|"]),
                        "y": list(truncateDfDict["twinLtest"]["N(L)"]),
                        "type": "scatter",
                        "name": "Observed",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(255, 127, 14)"},
                    },
                    {
                        "x": list(truncateDfDict["twinLtest"]["|L|"]),
                        "y": list(truncateDfDict["twinLtest"]["Untwinned"]),
                        "type": "scatter",
                        "name": "Untwinned",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dashdot"},
                    },
                    {
                        "x": list(truncateDfDict["twinLtest"]["|L|"]),
                        "y": list(truncateDfDict["twinLtest"]["Twinned"]),
                        "type": "scatter",
                        "name": "Perfect Twin",
                        "opacity": 0.8,
                        "mode": "lines",
                        "line": {"color": "rgb(31, 119, 180)", "dash": "dot"},
                    },
                ],
                "layout": {
                    "title": f"L-test (Padilla and Yeates)<br><sup>data from CTRUNCATE</sup>",
                    "xaxis": {
                        "title": "|L|",
                    },
                    "yaxis": {"title": "P(L <= L)"},
                },
                "help": """
    The Cumulative |L| plot for acentric data, where L = (I1-I2)/(I1+I2). This depends
    on the local difference in intensities. The difference operators used link to the
    neighbouring reflections taking into account possible tNCS operators.
    Note that this estimate is not as reliable as obtained via the H-test or ML Britton
    test if twin laws are available. However, it is less prone to the effects of anisotropy 
    than the H-test.

    Reference: Padilla, Yeates. A statistic for local intensity differences: robustness to 
    anisotropy and pseudo-centering and utility for detecting twinning. Acta Cryst. D59, 
    1124-30, 2003.
    """,
            },
        }
