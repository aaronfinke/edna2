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
__date__ = "02/05/2024"

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.CCP4Tasks import DimpleTask
from edna2.utils.UtilsPDB import pdbQuery, generatePdbSearchQuery, fetchPdbFileFromAccessionCode
from edna2.utils import UtilsCCTBX
from edna2.utils import UtilsLogging

from cctbx.crystal import symmetry
from cctbx import sgtbx

import json

logger = UtilsLogging.getLogger()


class CellToMRTask(AbstractTask):
    """
    given a unit cell and space group,
    this will search the PDB for best matches
    and use them as MR search models for a given
    list of reflections.
    "tolerance" here refers to the percent error
    for the PDB search parameters, default 0.01 (1%) for
    cell edges and 0.03 (3%) for cell angles.
    """

    def run(self, inData):
        self.inputMtz = inData["inputMtz"]
        self.unitCell = inData["unitCell"]
        self.spaceGroup = inData["spaceGroup"]


        self.tolerance = inData.get("tolerance",(0.01,0.03))
        self.numResults = inData.get("numResults",5)

        self.pdbQuery = generatePdbSearchQuery(unitCell=self.unitCell, spaceGroup=self.spaceGroup, dev=self.tolerance)
        if self.pdbQuery is None:
            logger.error("Error in generating search query")
            return
        
        with open(self.getWorkingDirectory() / "pdbQuery.json",'w') as fp:
            json.dump(json.loads(self.pdbQuery), fp, indent=2)
        
        self.queryResults = pdbQuery(self.pdbQuery)
        if not self.queryResults:
            logger.error("no results found. Stopping.")
            return
        
        with open(self.getWorkingDirectory() / "pdbQueryResults.json",'w') as fp:
            json.dump(self.queryResults, fp, indent=2)

        self.searchModels = [ x['identifier'] for x in self.queryResults if x['score'] == 1.0 ]
        if len(self.searchModels) > self.numResults:
            self.searchModels = self.searchModels[0:self.numResults]
        
        dimpleRuns = []
        for pdbCode in self.searchModels:
            inDataDimple = {
                "inputMtz": self.inputMtz,
                "inputPdb": pdbCode,
                "doSubmit": True
            }
            dimpleRun = DimpleTask(inData=inDataDimple, workingDirectorySuffix=pdbCode)
            dimpleRuns.append(dimpleRun)
            dimpleRun.start()
        index = 0

        #XXX: Probably shouldn't use an infinite while loop here, maybe timeout?
        while True:
            for run in dimpleRuns:
                if not run._process.is_alive():
                    run.join()
            if all(not t._process.is_alive() for t in dimpleRuns):
                break

        outData = self.queryResults
        return outData

            








        
