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
__date__ = "21/04/2019"

import os
import pathlib
import configparser
from cctbx.sgtbx import space_group_info

from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()

def parseSpaceGroup(spaceGroup):
    """parses space group and returns the space
    group number and string."""
    if not spaceGroup:
        logger.info("No space group supplied")
        return 0, ""
    try:
        spaceGroupInfo = space_group_info(spaceGroup).symbol_and_number()
        spaceGroupString = spaceGroupInfo.split("No. ")[0][:-2]
        spaceGroupNumber = int(spaceGroupInfo.split("No. ")[1][:-1])
        logger.debug("Supplied space group is {}, number {}".format(spaceGroupString, spaceGroupNumber))
    except:
        logger.debug("Could not parse space group")
        spaceGroupNumber = 0
        spaceGroupString = ""
    return spaceGroupNumber, spaceGroupString

def parseUnitCell(unitCell: str):
    """parse unit cell and return as a dict
    assumes a string with constants separated by commas"""
    try:
        unitCellList = [float(x) for x in unitCell.split(',')]
        #if there are zeroes parsed in, need to deal with it
        if 0.0 in unitCellList:
            raise Exception
        unitCell = {"cell_a": unitCellList[0],
                    "cell_b": unitCellList[1],
                    "cell_c": unitCellList[2],
                    "cell_alpha": unitCellList[3],
                    "cell_beta": unitCellList[4],
                    "cell_gamma": unitCellList[5]}
        logger.info("Supplied unit cell is {cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**unitCell))
    except:
        logger.debug("could not parse unit cell")
        unitCell = None
    return unitCell

def parseUnitCell_str(unitCell: str):
    """parse unit cell and return as a string
    assumes a string with constants separated by commas"""
    try:
        unitCellList = [float(x) for x in unitCell.split(',')]
        #if there are zeroes parsed in, need to deal with it
        if 0.0 in unitCellList:
            raise Exception
        unitCell = {"cell_a": unitCellList[0],
                    "cell_b": unitCellList[1],
                    "cell_c": unitCellList[2],
                    "cell_alpha": unitCellList[3],
                    "cell_beta": unitCellList[4],
                    "cell_gamma": unitCellList[5]}
        logger.info("Supplied unit cell is {cell_a} {cell_b} {cell_c} {cell_alpha} {cell_beta} {cell_gamma}".format(**unitCell))
    except:
        logger.debug("could not parse unit cell")
        unitCell = None
        
    return "{cell_a},{cell_b},{cell_c},{cell_alpha},{cell_beta},{cell_gamma}".format(**unitCell) if unitCell else None
