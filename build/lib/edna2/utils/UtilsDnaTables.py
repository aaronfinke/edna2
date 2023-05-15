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

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# kernel/src/EDUtilsTable.py

import json
import xmltodict

from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


def getDict(dnaTablesPath):
    with open(str(dnaTablesPath)) as f:
        dnaTables = f.read()
    if "</dna_tables>" not in dnaTables:
        # Fix for bug in MOSFLM
        dnaTables += "</dna_tables>"
    orderedDictDnaTables = xmltodict.parse(dnaTables)
    dictDnaTables = json.loads(json.dumps(orderedDictDnaTables))
    return dictDnaTables


def getTables(dictDnaTables, tableName):
    listTables = []
    listTable = dictDnaTables["dna_tables"]["table"]
    for table in listTable:
        if tableName == table["@name"]:
            listTables.append(table)
    return listTables


def getListParam(table):
    if isinstance(table["list"], list):
        listParam = table["list"]
    else:
        listParam = [table["list"]]
    return listParam


def getItemValue(dictParameter, key):
    value = None
    if isinstance(dictParameter["item"], list):
        listItem = dictParameter["item"]
    else:
        listItem = [dictParameter["item"]]
    for item in listItem:
        if item["@name"] == key:
            value = item["#text"]
    return _convertFromString(value)


def _convertFromString(value):
    if value is not None:
        try:
            if "." in value:
                value = float(value)
            else:
                value = int(value)
        except ValueError:
            # The value is returned as a string...
            pass
    return value


def getListValue(listParameter, key1, key2):
    value = None
    for dictParameter in listParameter:
        if dictParameter["@name"] == key1:
            if isinstance(dictParameter["item"], list):
                for item in dictParameter["item"]:
                    if item["@name"] == key2:
                        value = item["#text"]
            else:
                value = dictParameter["item"]["#text"]
    return _convertFromString(value)
