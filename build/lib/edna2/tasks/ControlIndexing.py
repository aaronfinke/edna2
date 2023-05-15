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

import os
import numpy as np

from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.SubWedgeAssembly import SubWedgeAssembly
from edna2.tasks.ControlDozor import ControlDozor
from edna2.tasks.XDSTasks import XDSIndexing
from edna2.utils import UtilsImage


class ControlIndexing(AbstractTask):
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

    def run(self, in_data):
        outData = {}
        # First get the list of subWedges
        if "subWedge" in in_data:
            list_sub_wedge = in_data["subWedge"]
        else:
            list_sub_wedge = self.getListSubWedge(in_data)
        # # Get list of spots from Dozor
        # listOutDataControlDozor = self.runControlDozor(listSubWedge)
        # listDozorSpotFile = []
        # for outDataControlDozor in listOutDataControlDozor:
        #     if "dozorSpotFile" in outDataControlDozor["imageQualityIndicators"][0]:
        #         dozorSpotFile = outDataControlDozor["imageQualityIndicators"][0][
        #             "dozorSpotFile"
        #         ]
        #         listDozorSpotFile.append(dozorSpotFile)
        # imageDict = listSubWedge[0]
        # Run XDS indexing
        xds_indexin_in_data = {
            "subWedge": list_sub_wedge
            # "dozorSpotFile": listDozorSpotFile,
        }
        xds_indexing_task = XDSIndexing(
            inData=xds_indexin_in_data
            # workingDirectorySuffix=UtilsImage.getPrefix(imageDict["image"][0]["path"]),
        )
        xds_indexing_task.execute()
        resultIndexing = None
        xparm_path = None
        spot_path = None
        if xds_indexing_task.isSuccess():
            xdsIndexingOutData = xds_indexing_task.outData
            if os.path.exists(xdsIndexingOutData["xparmXdsPath"]):
                xparm_path = xdsIndexingOutData["xparmXdsPath"]
            if os.path.exists(xdsIndexingOutData["spotXdsPath"]):
                spot_path = xdsIndexingOutData["spotXdsPath"]
            resultIndexing = ControlIndexing.getResultIndexingFromXds(
                xdsIndexingOutData
            )
        outData = {
            "resultIndexing": resultIndexing,
            "xparmXdsPath": xparm_path,
            "spotXdsPath": spot_path
        }
        return outData

    @staticmethod
    def getListSubWedge(in_data):
        list_sub_wedge = None
        # First check if we have data collection ids or image list
        # if "dataCollectionId" in inData:
        #     # TODO: get list of data collections from ISPyB
        #         logger.warning("Not implemented!")
        # el
        if "imagePath" in in_data or "fastCharacterisation":
            # Read the header(s)
            sub_wedge_assembly = SubWedgeAssembly(
                inData=in_data
                # workingDirectorySuffix=UtilsImage.getPrefix(listImagePath[0]),
            )
            sub_wedge_assembly.execute()
            list_sub_wedge = sub_wedge_assembly.outData["subWedge"]
        else:
            raise RuntimeError("No dataCollectionId or imagePath in inData")
        return list_sub_wedge

    @staticmethod
    def readImageHeaders(in_data):
        # Read the header(s)
        sub_wedge_assembly = SubWedgeAssembly(
            inData=in_data
            # workingDirectorySuffix=UtilsImage.getPrefix(listImagePath[0]),
        )
        sub_wedge_assembly.execute()
        list_sub_wedge = sub_wedge_assembly.outData["subWedge"]
        return list_sub_wedge

    @staticmethod
    def runControlDozor(listSubWedge):
        listControlDozor = []
        listOutDataControlDozor = []
        for subWedge in listSubWedge:
            listSubWedgeImage = subWedge["image"]
            for image in listSubWedgeImage:
                # listImage.append(image['path'])
                inDataControlDozor = {"image": [image["path"]]}
                controlDozor = ControlDozor(
                    inData=inDataControlDozor,
                    workingDirectorySuffix=UtilsImage.getPrefixNumber(image["path"]),
                )
                listControlDozor.append(controlDozor)
                controlDozor.start()
        for controlDozor in listControlDozor:
            controlDozor.join()
            if controlDozor.isSuccess():
                listOutDataControlDozor.append(controlDozor.outData)
        return listOutDataControlDozor

    @staticmethod
    def getResultIndexingFromXds(xdsIndexingOutData):
        idxref = xdsIndexingOutData["idxref"]
        xparamDict = xdsIndexingOutData["xparm"]
        if "A" in xparamDict:
            # Calculate MOSFLM UB matrix
            A = np.array(xparamDict["A"])
            B = np.array(xparamDict["B"])
            C = np.array(xparamDict["C"])

            volum = np.cross(A, B).dot(C)
            Ar = np.cross(B, C) / volum
            Br = np.cross(C, A) / volum
            Cr = np.cross(A, B) / volum
            UBxds = np.array([Ar, Br, Cr]).transpose()

            BEAM = np.array(xparamDict["beam"])
            ROT = np.array(xparamDict["rot"])
            wavelength = 1 / np.linalg.norm(BEAM)

            xparamDict["cell_volum"] = volum
            xparamDict["wavelength"] = wavelength
            xparamDict["Ar"] = Ar.tolist()
            xparamDict["Br"] = Br.tolist()
            xparamDict["Cr"] = Cr.tolist()
            xparamDict["UB"] = UBxds.tolist()

            normROT = float(np.linalg.norm(ROT))
            CAMERA_z = np.true_divide(ROT, normROT)
            CAMERA_y = np.cross(CAMERA_z, BEAM)
            normCAMERA_y = float(np.linalg.norm(CAMERA_y))
            CAMERA_y = np.true_divide(CAMERA_y, normCAMERA_y)
            CAMERA_x = np.cross(CAMERA_y, CAMERA_z)
            CAMERA = np.transpose(np.array([CAMERA_x, CAMERA_y, CAMERA_z]))

            mosflmUB = CAMERA.dot(UBxds) * xparamDict["wavelength"]
            # mosflmUB = UBxds*xparamDict["wavelength"]
            # xparamDict["mosflmUB"] = mosflmUB.tolist()

            reciprocCell = XDSIndexing.reciprocal(xparamDict["cell"])
            B = XDSIndexing.BusingLevy(reciprocCell)
            mosflmU = np.dot(mosflmUB, np.linalg.inv(B)) / xparamDict["wavelength"]
            # xparamDict[

            resultIndexing = {
                "spaceGroupNumber": idxref["spaceGroupNumber"],
                "cell": {
                    "a": idxref["a"],
                    "b": idxref["b"],
                    "c": idxref["c"],
                    "alpha": idxref["alpha"],
                    "beta": idxref["beta"],
                    "gamma": idxref["gamma"],
                },
                "xBeam": idxref["xBeam"],
                "yBeam": idxref["yBeam"],
                "distance": idxref["distance"],
                "qualityOfFit": idxref["qualityOfFit"],
                "mosaicity": idxref["mosaicity"],
                "XDS_xparm": xparamDict,
                "mosflmB": mosflmU.tolist(),
                "mosflmU": mosflmU.tolist(),
            }
        else:
            resultIndexing = {}
        return resultIndexing
