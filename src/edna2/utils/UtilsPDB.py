#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the 'Software'), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ["A. Finke"]
__license__ = "MIT"
__date__ = "13/06/2023"

import os
import time
import pathlib
import tempfile
import shutil

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
from edna2.utils import UtilsCCTBX
from edna2.utils import UtilsPath
from cctbx import sgtbx
from cctbx.crystal import symmetry
import tempfile
import requests,json
from pathlib import Path

logger = UtilsLogging.getLogger()

def fetchPdbFileFromAccessionCode(pdbCode, targetDirectory=Path.cwd()):
    """
    downloads the PDB file of the code, and returns the 
    path of the downloaded file. Saves the file to a temp
    directory first, then copies it to target directory.
    Returns final path to file. Uses get_pdb from iotbx.
    """
    from iotbx.pdb.fetch import get_pdb

    cwd = Path.cwd()

    mirror = UtilsConfig.get("PDB","pdbMirror", None)
    if mirror is None:
        logger.warning("No PDB mirror indicated- using RCSB as default")
        mirror = "rcsb"

    tempdir = tempfile.TemporaryDirectory()
    os.chdir(tempdir.name)
    try:
        logger.info(f"Downloading PDB File {pdbCode}.pdb from mirror {mirror}...")
        pdbFile = get_pdb(pdbCode, data_type="pdb",mirror=mirror, log=open(os.devnull,'w'), format="pdb")
        pdbPath = Path(pdbFile)
        UtilsPath.systemCopyFile(pdbPath,targetDirectory / pdbPath.name)
    except Exception as e:
        logger.error(f"Could not download and copy pdb file: {e}")
        os.chdir(cwd)
        tempdir.cleanup()
        return None
    
    os.chdir(cwd)
    tempdir.cleanup()
    return targetDirectory / pdbPath.name


def pdbQuery(query):
    """
    send a query to pdb API
    query should be a dict or JSON string
    see: https://search.rcsb.org/#search-api
    """
    queryURL = UtilsConfig.get("PDB","pdbQuery", None)
    if isinstance(query,dict):
        query_string = json.dumps(query)
    else:
        query_string = query

    #make sure it is valid json
    try:
        test = json.loads(query_string)
    except ValueError:
        logger.error('query is not JSON formattable')
        return
    
    query = f"{queryURL}?json={query_string}"
    results = requests.get(query)

    # any reqiest with a status code > 300 indicates something went wrong
    if results.status_code > 300:
        logger.error("PDB API query failed")
        return []
    result_json = json.loads(results.text)
    if result_json["total_count"] < 1:
        return []
    return result_json["result_set"]
    
def generatePdbSearchQuery(unitCell, spaceGroup, dev=(0.01,0.03)) -> str:
    """
    generate a search query JSON string 
    for the PDB API
    requires unit cell and SG
    """
    edgeDev = dev[0]
    angleDev = dev[1]
    unitCell_str = UtilsCCTBX.parseUnitCell_str(unitCell=unitCell)
    unitCell = UtilsCCTBX.parseUnitCell(unitCell=unitCell)
    
    cell_a = unitCell["cell_a"]
    cell_a_lims = (cell_a - (cell_a*edgeDev),cell_a+(cell_a*edgeDev))
    cell_b = unitCell["cell_b"]
    cell_b_lims = (cell_b - (cell_b*edgeDev),cell_b+(cell_b*edgeDev))
    cell_c = unitCell["cell_c"]
    cell_c_lims = (cell_c - (cell_c*edgeDev),cell_c+(cell_c*edgeDev))
    cell_alpha = unitCell["cell_alpha"]
    cell_beta = unitCell["cell_beta"]
    cell_gamma = unitCell["cell_gamma"]
    logger.debug(f"cell_a = {cell_a_lims}, cell_b={cell_b_lims},cell_c={cell_c_lims}")
    spaceGroupInfo = sgtbx.space_group_info(spaceGroup)    
    #make sure unit cell/SG are sensible
    try:
        symm = symmetry(unit_cell=unitCell_str,space_group_info=spaceGroupInfo)
    except Exception as e:
        logger.error(f"unit cell and SG not accepted: {e}")
        return
    
    spaceGroupString = spaceGroupInfo.type().lookup_symbol()
    spaceGroupMatch = spaceGroupMatches(spaceGroupString)

    cell_a_node = cellEdgeRange('a',cell_a_lims)
    cell_b_node = cellEdgeRange('b',cell_b_lims)
    cell_c_node = cellEdgeRange('c',cell_c_lims)

    crystalSystem = spaceGroupInfo.group().crystal_system()
    if crystalSystem == 'Cubic':
        cell_alpha_node = cellAngleEquals('alpha',90)
        cell_beta_node = cellAngleEquals('beta',90)
        cell_gamma_node = cellAngleEquals('gamma',90)
    elif crystalSystem == 'Hexagonal':
        cell_alpha_node = cellAngleEquals('alpha',90)
        cell_beta_node = cellAngleEquals('beta',90)
        cell_gamma_node = cellAngleEquals('gamma',120)
    elif crystalSystem == 'Tetragonal' or crystalSystem == 'Orthorhombic':
        cell_alpha_node = cellAngleEquals('alpha',90)
        cell_beta_node = cellAngleEquals('beta',90)
        cell_gamma_node = cellAngleEquals('gamma',90)
    elif crystalSystem == 'Monoclinic':
        cell_alpha_node = cellAngleEquals('alpha',90)
        cell_beta_lims = (cell_beta - (cell_beta*angleDev),cell_beta+(cell_beta+angleDev))
        cell_beta_node = cellAngleRange('beta',cell_beta_lims)
        cell_gamma_node = cellAngleEquals('gamma',90)
    elif crystalSystem == 'Triclinic':
        cell_alpha_lims = (cell_alpha - (cell_alpha*angleDev),cell_alpha+(cell_alpha+angleDev))
        cell_alpha_node = cellAngleRange('alpha',cell_alpha_lims)
        cell_beta_lims = (cell_beta - (cell_beta*angleDev),cell_beta+(cell_beta+angleDev))
        cell_beta_node = cellAngleRange('beta',cell_beta_lims)
        cell_gamma_lims = (cell_gamma - (cell_gamma*angleDev),cell_gamma+(cell_gamma+angleDev))
        cell_gamma_node = cellAngleRange('gamma',cell_gamma_lims)
    
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
            cell_a_node,
            cell_b_node,
            cell_c_node,
            cell_alpha_node,
            cell_beta_node,
            cell_gamma_node,
            spaceGroupMatch
            ],
            "label": "text"
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
            "start": 0,
            "rows": 25
            },
            "results_content_type": [
            "experimental"
            ],
            "sort": [
            {
                "sort_by": "score",
                "direction": "desc"
            }
            ],
            "scoring_strategy": "combined"
        }
        }
    return json.dumps(query)

def cellAngleEquals(angle,deg):
    assert angle.lower() in ["alpha","beta","gamma"]
    assert isinstance(deg,int)
    return {
            "type": "terminal",
            "service": "text",
            "parameters": {
            "attribute": f"cell.angle_{angle}",
            "operator": "equals",
            "negation": False,
            "value": deg
            }
            }
def cellAngleRange(angle,aRange):
    assert angle.lower() in ["alpha","beta","gamma"]
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

def cellEdgeRange(edge,arange):
    assert edge.lower() in ["a","b","c"]
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": f"cell.length_{edge}",
          "operator": "range",
          "negation": False,
          "value": {
            "from": arange[0],
            "to": arange[1],
            "include_lower": True,
            "include_upper": True
          }
        }
      }

def spaceGroupMatches(spaceGroup):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "symmetry.space_group_name_H_M",
          "operator": "exact_match",
          "negation": False,
          "value": spaceGroup
        }
      }
