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
# kernel/src/EDUtilsImage.py

import os
import re
import fabio
import pathlib


def __compileAndMatchRegexpTemplate(pathToImage):
    listResult = []
    if not isinstance(pathToImage, pathlib.Path):
        pathToImage = pathlib.Path(str(pathToImage))
    baseImageName = pathToImage.name
    regexp = re.compile(r"(.*)([^0^1^2^3^4^5^6^7^8^9])([0-9]*)\.(.*)")
    match = regexp.match(baseImageName)
    if match is not None:
        listResult = [
            match.group(0),
            match.group(1),
            match.group(2),
            match.group(3),
            match.group(4),
        ]
    return listResult


def getImageNumber(pathToImage):
    iImageNumber = None
    listResult = __compileAndMatchRegexpTemplate(pathToImage)
    if listResult is not None:
        iImageNumber = int(listResult[3])
    return iImageNumber


def getTemplate(pathToImage, symbol="#"):
    template = None
    listResult = __compileAndMatchRegexpTemplate(pathToImage)
    if listResult is not None:
        prefix = listResult[1]
        separator = listResult[2]
        imageNumber = listResult[3]
        suffix = listResult[4]
        hashes = ""
        for i in range(len(imageNumber)):
            hashes += symbol
        template = prefix + separator + hashes + "." + suffix
    return template


def getPrefix(pathToImage):
    prefix = None
    listResult = __compileAndMatchRegexpTemplate(pathToImage)
    if listResult is not None:
        prefix = listResult[1]
    return prefix


def getSuffix(pathToImage):
    suffix = None
    listResult = __compileAndMatchRegexpTemplate(pathToImage)
    if listResult is not None:
        suffix = listResult[4]
    return suffix


def getPrefixNumber(pathToImage):
    prefix = getPrefix(pathToImage)
    number = getImageNumber(pathToImage)
    prefixNumber = "{0}_{1:04d}".format(prefix, number)
    return prefixNumber


def splitPrefixRunNumber(path_to_image):
    file_name = pathlib.Path(path_to_image).name
    list_parts = str(file_name).split("_")
    pre_prefix = "_".join(list_parts[0:-2])
    run_number = int(list_parts[-2])
    return pre_prefix, run_number


def getH5FilePath(filePath, batchSize=100, hasOverlap=False, isFastMesh=False):
    if type(filePath) == str:
        filePath = pathlib.Path(filePath)
    imageNumber = getImageNumber(filePath)
    prefix = getPrefix(filePath)
    if hasOverlap or filePath.name.startswith("ref-"):
        h5ImageNumber = 1
        h5FileNumber = imageNumber
    elif (
        isFastMesh
        or filePath.name.startswith("mesh-")
        or filePath.name.startswith("line-")
    ):
        h5ImageNumber = int((imageNumber - 1) / 100) + 1
        h5FileNumber = 1
    else:
        h5ImageNumber = 1
        h5FileNumber = int((imageNumber - 1) / batchSize) * batchSize + 1
    h5MasterFileName = "{prefix}_{h5FileNumber}_master.h5".format(
        prefix=prefix, h5FileNumber=h5FileNumber
    )
    h5MasterFilePath = filePath.parent / h5MasterFileName
    h5DataFileName = "{prefix}_{h5FileNumber}_data_{h5ImageNumber:06d}.h5".format(
        prefix=prefix, h5FileNumber=h5FileNumber, h5ImageNumber=h5ImageNumber
    )
    h5DataFilePath = filePath.parent / h5DataFileName
    return h5MasterFilePath, h5DataFilePath, h5FileNumber


def mergeCbf(listPath, outputPath):
    firstImage = None
    no_images = len(listPath)
    for imagePath in listPath:
        if firstImage is None:
            firstImage = fabio.open(imagePath)
        else:
            image = fabio.open(imagePath)
            firstImage.data += image.data
    header_contents = firstImage.header["_array_data.header_contents"]
    list_header_contents = header_contents.split("\n")
    index = 0
    done = False
    while not done:
        line = list_header_contents[index]
        if line.startswith("# Angle_increment"):
            list_angle_increment = line.split(" ")
            angle_increment = float(list_angle_increment[2])
            new_angle_increment = angle_increment * no_images
            list_angle_increment[2] = str(new_angle_increment)
            new_line = " ".join(list_angle_increment)
            list_header_contents[index] = new_line
            done = True
        else:
            index += 1
    new_header_contents = "\n".join(list_header_contents)
    firstImage.header["_array_data.header_contents"] = new_header_contents
    firstImage.write(outputPath)
    return


def mergeCbfInDirectory(cbfDirectory, prefix=None, newPrefix=None):
    path_to_dir = pathlib.Path(cbfDirectory)
    index = None
    list_dir = [str(path) for path in list(path_to_dir.glob("*.cbf"))]
    list_dir.sort()
    list_of_image_lists = []
    listImage = None
    for cbf_file in list_dir:
        if not os.path.basename(cbf_file).startswith("ref-"):
            image_no = getImageNumber(cbf_file)
            if index is None or index != image_no:
                print("Starting image: {0}".format(cbf_file))
                index = image_no
                if listImage is not None:
                    list_of_image_lists.append(listImage)
                listImage = [cbf_file]
            else:
                print("Adding image no {0}".format(image_no))
                listImage.append(cbf_file)
            index += 1
    list_of_image_lists.append(listImage)
    image_number = 1
    firstImage = list_of_image_lists[0][0]
    directory = os.path.dirname(firstImage)
    old_run_number = None
    for list_image in list_of_image_lists:
        pre_prefix, run_number = splitPrefixRunNumber(list_image[0])
        new_run_Number = run_number + 10
        if old_run_number is None or old_run_number != new_run_Number:
            old_run_number = new_run_Number
            image_number = 1
        new_cbf_path = "{0}/{1}_{2}_{3:04d}.cbf".format(
            directory, pre_prefix, new_run_Number, image_number
        )
        mergeCbf(list_image, new_cbf_path)
        image_number += 1
