#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the 'Software')), to deal in
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
__date__ = "12/02/2024"


import unittest

from edna2.utils import UtilsPDB


class UtilsPDBUnitTest(unittest.TestCase):
    def test_searchQuery(self):
        unitCell = "58.3,58.3,150.65,90,90,90"
        sg = "P 41 21 2"
        query = UtilsPDB.generatePdbSearchQuery(unitCell=unitCell, spaceGroup=sg)
        result = UtilsPDB.pdbQuery(query)
        self.assertTrue(result)

    def test_fetchPDB(self):
        pdbCode = '13PK'
        pdb = UtilsPDB.fetchPdbFileFromAccessionCode(pdbCode=pdbCode)
        assert pdb.exists()


if __name__ == '__main__':
    unittest.main()
