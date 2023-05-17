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
__date__ = "06/02/2023"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx

# mxPluginExec/plugins/EDPluginGroupXDS-v1.0/plugins/EDPluginXDSv1_0.py
# mxPluginExec/plugins/EDPluginGroupXDS-v1.0/plugins/EDPluginXDSIndexingv1_0.py

import os
import math
import shutil
import numpy as np
from pathlib import Path
import sys

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsImage
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsDetector
from edna2.utils import UtilsSymmetry


logger = UtilsLogging.getLogger()


R2D = 180 / math.pi


class XSCALETask(AbstractTask):
    """
    Runs XSCALE for merging statistics
    """

    def run(self, inData):
        commandLine = "/mxn/groups/sw/mxsw/XDS/xscale_par"
        self.isAnom = inData["isAnom"]
        self.merge = inData["merge"]
        listXSCALE_INP = self.generateXSCALE_INP(inData=inData, isAnom=self.isAnom, merge=self.merge)
        self.writeXSCALE_INP(listXSCALE_INP, self.getWorkingDirectory())
        self.setLogFileName("xscale.log")
        self.runCommandLine(commandLine, listCommand=[])
        # Work in progress!
        outData = self.parseXSCALEOutput()
        return outData

    def generateXSCALE_INP(self, inData, isAnom, merge):

        friedels_law = False if isAnom else True
        xdsAsciiPath = Path(inData["xdsAsciiPath_anom"]) if isAnom else Path(inData["xdsAsciiPath_noAnom"])
        Path(self.getWorkingDirectory() / xdsAsciiPath.name).symlink_to(xdsAsciiPath)
        bins = " ".join([str(x) for x in inData["bins"]])
        sgNumber = inData["sgNumber"]
        a,b,c,alpha,beta,gamma = inData["cell"].values()
        res = inData["res"]
        output_file_name = f"{'' if merge else 'un'}merged_{'anom' if isAnom else 'noanom'}_XSCALE.hkl"
        
        listXSCALE_INP = [
            f"OUTPUT_FILE= {output_file_name}",
            f"FRIEDEL'S_LAW= {'TRUE' if friedels_law else 'FALSE'}",
            f"MERGE={'TRUE' if merge else 'FALSE'}",
            f"INPUT_FILE= {xdsAsciiPath.name}",
            f"UNIT_CELL_CONSTANTS= {a} {b} {c} {alpha} {beta} {gamma}",
            f"SPACE_GROUP_NUMBER= {sgNumber}",
            f"RESOLUTION_SHELLS= {bins}"
        ]

        return listXSCALE_INP

    def writeXSCALE_INP(self, listXSCALE_INP, workingDirectory):
        fileName = "XSCALE.INP"
        filePath = workingDirectory / fileName
        with open(str(filePath), "w") as f:
            for line in listXSCALE_INP:
                f.write(line + '\n')

        
    def parseXSCALEOutput(self):
        outData = {}
        xscaleLPPath = self.getWorkingDirectory() / "XSCALE.LP"
        xscaleInpPath = self.getWorkingDirectory() / "XSCALE.INP"
        outData["xscaleLp"] = str(xscaleLPPath)
        outData["xscaleInp"] = str(xscaleInpPath)
        try:
            with open(Path(xscaleLPPath),'r') as fp:
                lines =  [l.strip('\n') for l in fp.readlines()]
        except IOError:
            logger.error("Could not open the specified XSCALE output file for reading: {0}".format(xscaleLPPath))
            return None

        completeness_entry_begin = [i for i,s in enumerate(lines) if 'LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected' in s][-1]
        completeness_entry_end = [i for i,s in enumerate(lines[completeness_entry_begin:]) if 'total' in s][0]
        completenessEntries = self._extractCompletenessEntries(lines[completeness_entry_begin+1:completeness_entry_begin+completeness_entry_end+1])
        outData['completenessEntries'] = completenessEntries
        return outData

    def _extractCompletenessEntries(self, lines):
        """
        Since the latest XDS version there's no guarantee the fields
        will be separated by whitespace. What's fixed is the size of
        each field. So we'll now use fixed offsets to extract the
        fields.

        The Fortran code uses this format statement:
        1130  FORMAT(F9.2,I12,I8,I10,F11.1,'%',F10.1,'%',F9.1,'%',I9,F8.2,      &
                F8.1,'%',F8.1,A1,I6,A1,F8.3,I8)
        """
        logger.debug(f'extracting completeness entries...')
        outData = {"completenessEntries": []}
        offsets = {
        'res': (0,9),
        'observed': (9, 21),
        'unique': (21, 29),
        'possible': (29, 39),
        'complete': (39, 50),
        'rfactor': (51, 61),
        'isig': (81, 89),
        'rmeas': (89, 98),
        'half_dataset_correlation': (99, 107),
        }

        for line in lines:
            #deal with blank lines and total completeness
            if line.strip() == '':
                continue
            if 'total' in line:
                res_dict = {}
                for (name, (start, end)) in offsets.items():
                    value = float(line[start:end]) if not 'total' in line[start:end] else 'total'
                    res_dict[name] = value
                res_dict.pop('res', None)
                outData['total_completeness'] = res_dict
            else:
                res_dict = {}
                for (name, (start, end)) in offsets.items():
                    value = float(line[start:end])
                    res_dict[name] = value
                outData['completenessEntries'].append(res_dict)
            for res in outData['completenessEntries']:
                res["multiplicity"] = round(res["observed"] / res["possible"], 2)
        return outData
