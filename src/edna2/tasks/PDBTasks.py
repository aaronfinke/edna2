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
from edna2.utils.UtilsPDB import pdbQuery
from cctbx.crystal import symmetry
from cctbx import sgtbx

class PDBAPIUnitCellQueryTask(AbstractTask):
    def run(self,inData):
        self.unitCell = inData["unitCell"]
        self.spaceGroup = inData['spaceGroup']
        self.deviation = inData.get('deviation',(0.01,0.03))

        pass

    def generatePdbSearchQuery(self):
        unitCell = self.unitCell
        spaceGroup = self.spaceGroup
        dev = self.deviation

        
        
