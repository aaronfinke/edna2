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
__date__ = "16/10/2020"


from edna2.tasks.AbstractTask import AbstractTask
from edna2.tasks.ReadImageHeader import ReadImageHeader
from edna2.tasks.ControlDozor import ControlDozor
from edna2.tasks.XDSTasks import XDSIndexing
from edna2.tasks.MosflmTasks import MosflmGeneratePredictionTask


class ControlPredictionTask(AbstractTask):
    """
    This task receives a list of images or data collection ids and
    returns result of indexing
    """

    def getInDataSchema(self):
        return {
            "type": "object",
            "properties": {  }
        }

    def run(self, inData):
        outData = {}
        xdsIndexingOutData = inData["indexingSolution"]
        # Run MOSFLM prediction
        mosflmUB = xdsIndexingOutData["xparm"]["mosflmUB"]
        mosflmU = xdsIndexingOutData["xparm"]["mosflmU"]
        matrix = {
            "matrixA": mosflmUB,
            "missettingsAngles": [
                0.0,
                0.0,
                0.0
            ],
            "matrixU": mosflmU,
            "cell": {
                "a": xdsIndexingOutData["idxref"]["a"],
                "b": xdsIndexingOutData["idxref"]["b"],
                "c": xdsIndexingOutData["idxref"]["c"],
                "alpha": xdsIndexingOutData["idxref"]["alpha"],
                "beta": xdsIndexingOutData["idxref"]["beta"],
                "gamma": xdsIndexingOutData["idxref"]["gamma"]
            },
        }
        mosflmInData = MosflmGeneratePredictionTask.generateMOSFLMInData(inData={"subWedge": listSubWedge})
        mosflmInData["matrix"] = matrix
        mosflmInData.update({
            "mosaicityEstimation": xdsIndexingOutData["idxref"]["mosaicity"],
            "deviationAngular": 1.0,
            "refinedDistance": xdsIndexingOutData["idxref"]["distance"],
            "symmetry": "P1"
        })
        mosflmGeneratePredictionTask =MosflmGeneratePredictionTask(inData=mosflmInData)
        mosflmGeneratePredictionTask.execute()
        pass
        return outData

