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
__date__ = "14/04/2020"


import numpy as np

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.SubWedgeAssembly import SubWedgeAssembly
from edna2.tasks.XDSTasks import XDSIndexAndIntegration

from edna2.utils import UtilsImage

class RadiationDamageProcessing(AbstractTask):
    """
    This task receives a list of images or data collection ids and
    returns result of indexing
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {
                "dataCollectionId": {"type": "integer"},
                "imagePath": {
                    "type": "array",
                    "items": {
                        "type": "string",
                    },
                },
            },
        }

    def run(self, inData):
        outData = {}
        # First get the list of subWedges
        if "subWedge" in inData:
            list_sub_wedge = inData["subWedge"]
            directory_prefix = ""
        else:
            first_image = inData["imagePath"][0]
            last_image = inData["imagePath"][-1]
            prefix = UtilsImage.getPrefix(first_image)
            first_image_number = UtilsImage.getImageNumber((first_image))
            last_image_number = UtilsImage.getImageNumber((last_image))
            directory_prefix = "{0}_{1}_{2}".format(
                prefix, first_image_number, last_image_number
            )
            list_sub_wedge = self.getListSubWedge(inData)
        indata_xds_integration = {
            "subWedge": list_sub_wedge
        }
        xds_integration = XDSIndexAndIntegration(
            inData=indata_xds_integration,
            workingDirectorySuffix=directory_prefix,
        )
        xds_integration.execute()
        outData = xds_integration.outData
        return outData




    def getListSubWedge(self, inData):
        list_sub_wedges = None
        sub_wedge_assembly = SubWedgeAssembly(inData=inData)
        sub_wedge_assembly.execute()
        if sub_wedge_assembly.isSuccess():
            list_sub_wedges = sub_wedge_assembly.outData["subWedge"]
        return list_sub_wedges
