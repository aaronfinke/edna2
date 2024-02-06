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
from edna2.utils import UtilsCCTBX
from cctbx.crystal import symmetry
from cctbx import sgtbx

class PDBAPIUnitCellQueryTask(AbstractTask):
    def run(self,inData):
        self.unitCell = inData["unitCell"]
        self.spaceGroup = inData['spaceGroup']
        self.deviation = inData.get('deviation',(0.01,0.03))

        pass

    def cellAngleEquals(self,angle,deg):
        assert angle.lower in ["alpha","beta","gamma"]
        assert isinstance(deg,int)
        return {
                "type": "terminal",
                "service": "text",
                "parameters": {
                "attribute": f"cell.angle_{angle.lower}",
                "operator": "equals",
                "negation": False,
                "value": deg
                }
                }
    def cellAngleRange(self,angle,aRange):
        assert angle.lower in ["alpha","beta","gamma"]
        return {
                "type": "terminal",
                "service": "text",
                "parameters": {
                "attribute": f"cell.angle_{angle}",
                "operator": "range",
                "negation": False,
                "value": {
                    "from": aRange[0],
                    "to": aRange[1],
                    "include_lower": True,
                    "include_upper": True
                }
                }
            }


    def generatePdbSearchQuery(self):
        unitCell = self.unitCell
        spaceGroup = sgtbx.space_group_info(self.spaceGroup)
        dev = self.deviation
        edgeDev = dev[0]
        angleDev = dev[1]
        if isinstance(unitCell,str):
            unitCell = UtilsCCTBX.parseUnitCell(unitCell)

        unitCell_str = UtilsCCTBX.parseUnitCell_str(unitCell)
        
        cell_a = unitCell["cell_a"]
        cell_a_lims = (cell_a - (cell_a*edgeDev),cell_a+(cell_a+edgeDev))
        cell_b = unitCell["cell_b"]
        cell_b_lims = (cell_b - (cell_b*edgeDev),cell_b+(cell_b+edgeDev))
        cell_c = unitCell["cell_c"]
        cell_a_lims = (cell_c - (cell_c*edgeDev),cell_c+(cell_c+edgeDev))
        cell_alpha = unitCell["cell_alpha"]
        cell_beta = unitCell["cell_beta"]
        cell_gamma = unitCell["cell_gamma"]
        
        symm = symmetry(unit_cell=unitCell_str,space_group_info=spaceGroup)
        
        crystalSystem = spaceGroup.group().crystal_system()
        if crystalSystem == 'Cubic':
            cell_alpha_node = self.cellAngleEquals('alpha',90)
            cell_beta_node = self.cellAngleEquals('beta',90)
            cell_gamma_node = self.cellAngleEquals('gamma',90)
        elif crystalSystem == 'Hexagonal':
            cell_alpha_node = self.cellAngleEquals('alpha',90)
            cell_beta_node = self.cellAngleEquals('beta',90)
            cell_gamma_node = self.cellAngleEquals('gamma',120)
        elif crystalSystem == 'Tetragonal' or crystalSystem == 'Orthorhombic':
            cell_alpha_node = self.cellAngleEquals('alpha',90)
            cell_beta_node = self.cellAngleEquals('beta',90)
            cell_gamma_node = self.cellAngleEquals('gamma',90)
        elif crystalSystem == 'Monoclinic':
            cell_alpha_lims = (cell_alpha - (cell_alpha*angleDev),cell_alpha+(cell_alpha+angleDev))
            cell_alpha_node = self.cellAngleRange('alpha',cell_alpha_lims)
            cell_beta_node = self.cellAngleEquals('beta',90)
            cell_gamma_lims = (cell_gamma - (cell_gamma*angleDev),cell_gamma+(cell_gamma+angleDev))
            cell_gamma_node = self.cellAngleRange('gamma',cell_gamma_lims)
        elif crystalSystem == 'Triclinic':
            cell_alpha_lims = (cell_alpha - (cell_alpha*angleDev),cell_alpha+(cell_alpha+angleDev))
            cell_alpha_node = self.cellAngleRange('alpha',cell_alpha_lims)
            cell_beta_lims = (cell_beta - (cell_beta*angleDev),cell_beta+(cell_beta+angleDev))
            cell_beta_node = self.cellAngleRange('beta',cell_alpha_lims)
            cell_gamma_lims = (cell_gamma - (cell_gamma*angleDev),cell_gamma+(cell_gamma+angleDev))
            cell_gamma_node = self.cellAngleRange('gamma',cell_gamma_lims)




            








        
