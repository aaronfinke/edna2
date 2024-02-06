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

import requests,json

logger = UtilsLogging.getLogger()

def fetchPdbFileFromAccessionCode(pdbCode):
    """
    downloads the PDB file of the code, and returns the 
    path of the downloaded file.
    """
    from iotbx.pdb.fetch import fetch, get_pdb
    mirror = UtilsConfig.get("PDB","pdbMirror", None)
    if mirror is None:
        logger.warning("No PDB mirror indicated- using RCSB as default")
        mirror = "rcsb"
    try:
        pdbFile = get_pdb(pdbCode, data_type="pdb",mirror=mirror, log=logger, format="pdb")
        return pdbFile
    except:
        logger.error("Could not download pdb file!")
    return None

def pdbQuery(query):
    """
    send a query to pdb API
    query should be a dict
    see: https://search.rcsb.org/#search-api
    """
    queryURL = UtilsConfig.get("PDB","pdbQuery", None)
    query_string = json.dumps(query)
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
    
