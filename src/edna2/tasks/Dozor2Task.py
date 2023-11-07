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

__author__ = "Aaron Finke"
__contact__ = "aaron.finke@maxiv.lu.se"
__copyright__ = "MAX IV Laboratory"
__updated__ = "2023-11-03"

import os
import time
import socket
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

from edna2.tasks.AbstractTask import AbstractTask
from edna2.utils import UtilsConfig
from edna2.utils import UtilsImage


from edna2.utils import UtilsLogging
from edna2.utils import UtilsIspyb


logger = UtilsLogging.getLogger()


class Dozor2Task(AbstractTask):  # pylint: disable=too-many-instance-attributes
    """
    The Dozor2Task is responsible for executing the 'dozor2' program.
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "nx": {"type": "integer"},
                "ny": {"type": "integer"},
                "pixel": {"type": "number"},
                "exposure": {"type": "number"},
                "detector_distance": {"type": "number"},
                "X-ray_wavelength": {"type": "number"},
                "fraction_polarization": {"type": "number"},
                "orgx": {"type": "number"},
                "orgy": {"type": "number"},
                "oscillation_range": {"type": "number"},
                "starting_angle": {"type": "number"},
                "number_images": {"type": "integer"},
                "first_image_number": {"type": "integer"},
                "name_template_image": {"type": "string"},
                "library": {"type": "string"},
                "pixel_min": {"type": "number"},
                "pixel_max": {"type": "number"},
                "ix_max": {"type": "number"},
                "iy_max": {"type": "number"},
                "ix_min": {"type": "number"},
                "iy_min": {"type": "number"},
                "spot_size": {"type": "integer"},
                "bad_zona": {"type": "string"},
                "wedge_number": {"type": "integer"},
                "spot_level": {"type": "integer"},
                "dataCollectionId": {"type": "integer"},
                "masterFilePath": {"type": "string"},
            },
        }

    def getOutDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "logPath": {"type": "string"},
                "workingDirectory": {"type": "string"},
            },
        }

    def run(self, inData):
        time_start = time.perf_counter()
        logger.info("Dozor spotfinding started")
        if os.environ.get("SLURM_JOB_ID"):
            logger.info(f"SLURM job id: {os.environ.get('SLURM_JOB_ID')}")
        logger.info(f"Running on {socket.gethostname()}")
        omega = None
        outData = {}
        self.masterFilePath = inData.get("masterFilePath")
        self.dataCollectionId = inData.get("dataCollectionId", 0)

        dozorInData = self.getInDataSchema()["properties"]
        for k, v in dozorInData.items():
            dozorInData[k] = inData.get(k)
        del dozorInData["masterFilePath"]
        del dozorInData["dataCollectionId"]
        if self.dataCollectionId and not self.masterFilePath:
            self.masterFilePath = UtilsIspyb.getXDSMasterFilePath(self.dataCollectionId)
        if self.masterFilePath and Path(self.masterFilePath).exists():
            dozorMasterFileInData, omega = self.parseDozorInDataFromMasterFile(
                self.masterFilePath
            )
            for k, v in dozorMasterFileInData.items():
                dozorInData[k] = v if dozorInData[k] is None else dozorInData[k]
            dozorInData["name_template_image"] = str(self.masterFilePath).replace(
                "master", "??????"
            )
        dozorInData["library"] = UtilsConfig.get(self, "library")
        logger.debug(f"dozorInData: {dozorInData}")
        dozorInputFileData = self.generateDozorInputFileData(dozorInData)
        dozorDat = self.writeDozorInputFileData(
            dozorInputFileData, self.getWorkingDirectory()
        )

        dozorSetup = UtilsConfig.get(self, "dozorSetup")
        cmd = ""
        cmd += f"source {dozorSetup}\n" if dozorSetup is not None else ""
        cmd += "dozor dozor.dat"

        self.runCommandLine(cmd)
        time_dozor = time.perf_counter()
        logger.info(f"Dozor finished. Process Time: {time_dozor - time_start} s")
        dozorLogFile = self.getLogPath()

        dozor_df = self.parseDozorLogFile(dozorLogFile)
        if omega is None:
            starting_angle = dozorInData["starting_angle"]
            oscillation_range = dozorInData["oscillation_range"]
            number_images = dozorInData["number_images"]
            omega = [
                (x * oscillation_range) + starting_angle
                for x in range(0, number_images)
            ]
        dozor_df["angle"] = omega

        self.saveDozorDataFrame(
            dozor_df, self.getWorkingDirectory(), self.dataCollectionId
        )
        self.generateDozorPlot(
            dozor_df,
            self.masterFilePath,
            self.getWorkingDirectory(),
            self.dataCollectionId,
        )
        time_end = time.perf_counter()
        logger.info(f"Process Time: {time_end - time_start} s")
        return outData

    @staticmethod
    def parseDozorInDataFromMasterFile(masterFile):
        imageHeader = UtilsImage.getImageMetadataFromH5MasterFile(masterFile)
        if not imageHeader:
            logger.error("Could not read image header")
            return

        dozorInData = {
            "nx": imageHeader["nx"],
            "ny": imageHeader["ny"],
            "pixel": imageHeader["x_pixel_size"],
            "exposure": imageHeader["count_time"],
            "detector_distance": imageHeader["detector_distance"],
            "X-ray_wavelength": imageHeader["wavelength"],
            "fraction_polarization": UtilsConfig.get(
                "Dozor2Task", "polarization_fraction", 1.0
            ),
            "orgx": imageHeader["beam_center_x"],
            "orgy": imageHeader["beam_center_y"],
            "oscillation_range": imageHeader["omega_range_average"],
            "starting_angle": imageHeader["starting_angle"],
            "number_images": imageHeader["num_images"],
            "first_image_number": 1,
            "pixel_min": 0,
            "pixel_max": imageHeader["saturation_value"],
            "ix_max": 0,
            "iy_max": 0,
            "ix_min": 0,
            "iy_min": 0,
        }
        omega = imageHeader["omega"]
        return dozorInData, omega

    @staticmethod
    def generateDozorInputFileData(dozorInData):
        inputData = ""
        for k, v in dozorInData.items():
            inputData += f"{k} {v}\n" if v is not None else ""
        return inputData

    @staticmethod
    def writeDozorInputFileData(dozorInputFileData, dirf):
        dirf = Path(dirf)
        with open(dirf / "dozor.dat", "w") as fp:
            fp.write(dozorInputFileData)
        return dirf / "dozor.dat"

    @staticmethod
    def parseDozorLogFile(dozorLogFile):
        """
        This is neggia 1.2.0 (Copyright Dectris 2020)
         Program dozor /A.Popov & G.Bourenkov/
         Version 2.0.2 //  20.03.2020
         Copyright 2014-2020 by Alexander Popov and Gleb Bourenkov
         [generic_data_plugin] - INFO - generic_open
               + library          = </mxn/groups/sw/mxsw/xds_related/dectris-neggia.so>
               + template_name    = </gpfs/offline1/visitors/biomax/20170251/20220927/raw/Tau/Tau-BrD1/repack-Tau-BrD1-1-14200_master.h5>
               + dll_filename     = </mxn/groups/sw/mxsw/xds_related/dectris-neggia.so>
               + image_data_filename   = </gpfs/offline1/visitors/biomax/20170251/20220927/raw/Tau/Tau-BrD1/repack-Tau-BrD1-1-14200_master.h5>
         [generic_data_plugin] - INFO - generic_get_header
         N    | SPOTS     Main     Visible
        image | num.of   Score    Resolution
        ------------------------------------
            1 |   864     36.71    1.69
            2 |   824     38.11    1.69
            3 |   861     36.69    1.66
            4 |   889     37.85    1.69
        ...
        """
        skiprow = 0
        with open(dozorLogFile, "r") as fp:
            for n, line in enumerate(fp):
                if "-------------" in line:
                    skiprow = n + 1
                    break

        colspecs = [(1, 5), (6, 7), (7, 17), (17, 27), (27, -1)]
        df = pd.read_fwf(
            dozorLogFile,
            skiprows=skiprow,
            names=["NImage", "|", "NumSpots", "DozorScore", "Resolution"],
            skipfooter=1,
            colspecs=colspecs,
        )
        df = df.drop("|", axis=1)
        return df

    @staticmethod
    def saveDozorDataFrame(df, workingDirectory, dataCollectionId):
        try:
            df.to_csv(workingDirectory / f"dozor_{dataCollectionId:06d}.csv")
        except Exception as e:
            logger.error(f"Couldn't save CSV file: {e}")

    @staticmethod
    def generateDozorPlot(df, masterFile, workingDirectory, dataCollectionId=0):
        df["DozorScore"] = df["DozorScore"] * 10
        fig, ax1 = plt.subplots(figsize=(7, 5), layout="constrained")
        fig.suptitle(masterFile, fontsize=8)
        ax2 = ax1.twinx()
        ax2.invert_yaxis()
        ax3 = ax1.twiny()
        nf1 = df.plot.scatter(
            x="NImage",
            y="DozorScore",
            ax=ax1,
            s=2,
            c="orange",
            label="Dozor Score",
            legend=False,
        )
        nf2 = df.plot.scatter(
            x="NImage",
            y="NumSpots",
            ax=ax1,
            s=2,
            c="blue",
            label="Number of Spots",
            legend=False,
        )
        nf3 = df.plot.scatter(
            x="NImage",
            y="Resolution",
            ax=ax2,
            s=2,
            c="purple",
            label="Visible Resolution",
            legend=False,
        )
        ax3.scatter(df["angle"], df["DozorScore"], s=0)

        minResolution = max(df["Resolution"])
        maxResolution = min(df["Resolution"])
        maxImageNumber = max(df["NImage"])
        minImageNumber = min(df["NImage"])
        minAngle = min(df["angle"])
        maxAngle = max(df["angle"])
        maxDozorValue = max(df["DozorScore"])
        minDozorValue = min(df["DozorScore"])

        if maxDozorValue < 0.001 and minDozorValue < 0.001:
            maxDozorValue = 0.5
            minDozorValue = -0.5

        if maxResolution is None or maxResolution > 0.8:
            maxResolution = 0.8
        else:
            maxResolution = int(maxResolution * 10.0) / 10.0

        if minResolution is None or minResolution < 4.5:
            minResolution = 4.5
        else:
            minResolution = int(minResolution * 10.0) / 10.0 + 1

        ax2.set_ylim(minResolution, maxResolution)
        ax1.set_xlim(minImageNumber, maxImageNumber)
        ax1.set_ylabel("Number of Spots/Dozor Score(*10)", fontsize=12)
        ax2.set_ylabel("Resolution (Å)", fontsize=12)
        ax3.set_xlabel("Angle (degrees)", fontsize=12, labelpad=8)
        ax1.set_xlabel("Image Number", fontsize=12, labelpad=8)
        ax2.set_axisbelow(True)
        ax3.set_axisbelow(True)
        ax2.yaxis.grid(color="gray", linestyle="--")
        ax3.xaxis.grid(color="gray", linestyle="--")
        ax1.tick_params(axis="both", which="major", labelsize=12, pad=4, direction="in")
        ax3.tick_params(axis="both", which="major", labelsize=12, pad=4, direction="in")
        ax2.tick_params(axis="both", which="major", labelsize=12, pad=4, direction="in")

        lns = [nf1, nf2, nf3]
        labs = [l.get_label() for l in lns]
        blue_line = mpatches.Patch(color="blue", label="Number of Spots")
        orange_line = mpatches.Patch(color="orange", label="Dozor Score")
        purple_line = mpatches.Patch(color="purple", label="Visible Resolution")

        fig.legend(
            handles=[blue_line, orange_line, purple_line],
            ncols=3,
            loc="outside lower center",
            frameon=False,
        )
        try:
            plt.savefig(workingDirectory / f"dozor_{dataCollectionId:06d}.png")
        except Exception as e:
            logger.error(f"Couldn't save png file: {e}")
