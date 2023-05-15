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
__date__ = "28/03/2022"

import json

from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


def subWedgeMerge(inData):
    list_result_sub_wedge_merge = []
    # 1. Sort the incoming subwedges in groupes of similar experimental conditions
    list_sub_wedge_sorted = sortSubWedgesOnExperimentalCondition(inData)

    # 2. Merge subwedges which are adjascent with respect to the orientation matrix
    for list_sub_wedge in list_sub_wedge_sorted:
        list_sub_wedge_merged = mergeListOfSubWedgesWithAdjascentRotationAxis(
            list_sub_wedge
        )
        for index, sub_wedge in enumerate(list_sub_wedge_merged):
            sub_wedge["subWedgeNumber"] = index + 1
        list_result_sub_wedge_merge += list_sub_wedge_merged
    return list_result_sub_wedge_merge




def compareTwoValues(value1, value2, tolerance=0.001):
    """
    This method compares two values which can be of type double, integer or string.
    The types of the two values two be compared must be identical.
    If the value is of type double, an optional tolerance can be passed as an
    argument to the method. If both values are None the methods returns True.
    """
    return_value = False
    if value1 is None and value2 is None:
        return_value = True
    elif type(value1) != type(value2):
        error_message = "Types of values different : value1=%r, value2=%r" % (
            type(value1),
            type(value2),
        )
        logger.error(error_message)
        raise RuntimeError(error_message)
    elif type(value1) == float:
        if abs(value1 - value2) < tolerance:
            return_value = True
    elif type(value1) == int:
        if value1 == value2:
            return_value = True
    elif type(value1) == str:
        if value1 == value2:
            return_value = True
    else:
        error_message = "Unknown value type : %r for value %r" % (type(value1), value1)
        logger.error(error_message)
        raise RuntimeError(error_message)
    return return_value


def compareTwoParameters(exp_cond_1, exp_cond_2, type1, type2, tolerance=0.001):
    return_value = True
    if type1 in exp_cond_1 and type1 in exp_cond_2:
        if type2 in exp_cond_1[type1] and type2 in exp_cond_2[type1]:
            value1 = exp_cond_1[type1][type2]
            value2 = exp_cond_2[type1][type2]
            return_value = compareTwoValues(value1, value2, tolerance=tolerance)
    return return_value


def isSameExperimentalCondition(exp_cond_1, exp_cond_2):
    """
    This method compares two experimental condition objects in order to verify if
    they can belong to the same sub wedge. The following parameters are checked:
    beam.exposureTime [s], tolerance 0.001
    beam.wavelength [A], tolerance 0.001
    detector.beamPositionX [mm], tolerance 0.1
    detector.beamPositionY [mm], tolerance 0.1
    detector.distance [mm], tolerance 0.1
    detector.name [string]
    detector.numberPixelX [int]
    detector.numberPixelY [int]
    detector.serialNumber [string]
    detector.twoTheta [degrees], tolerance 0.1
    goniostat.oscillationWidth [degrees], tolerance 0.001
    goniostat.rotationAxis [string]
    """
    return_value = True

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "beam", "exposureTime"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "beam", "wavelength"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "detector", "distance"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "detector", "name"):
        return_value = False

    if not compareTwoParameters(
        exp_cond_1, exp_cond_2, "detector", "beamPositionX", tolerance=0.1
    ):
        return_value = False

    if not compareTwoParameters(
        exp_cond_1, exp_cond_2, "detector", "beamPositionY", tolerance=0.1
    ):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "detector", "numberPixelX"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "detector", "numberPixelY"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "detector", "serialNumber"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "detector", "twoTheta"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "goniostat", "oscillationWidth"):
        return_value = False

    if not compareTwoParameters(exp_cond_1, exp_cond_2, "goniostat", "rotationAxis"):
        return_value = False

    return return_value


def isSameExperimentalConditionInSubWedge(subWedge1, subWedge2):
    """
    This method compares two XSDataSubWedge objects in order to verify if
    they can belong to the same sub wedge. The two experimental condition
    objects are compared using the method "isSameExperimentalCondition"
    """
    exp_cond_1 = subWedge1["experimentalCondition"]
    exp_cond_2 = subWedge2["experimentalCondition"]
    return isSameExperimentalCondition(exp_cond_1, exp_cond_2)


def sortIdenticalObjects(list_objects, method_comparison):
    """
    This method takes as input a list of objects with the same type.
    It returns a list containing list of identical objects. If an object
    is not identical to another it is returned in a list alone.
    """
    list_results = []
    number_of_objects = len(list_objects)
    if number_of_objects == 0:
        pass
    elif number_of_objects == 1:
        list_results.append(list_objects)
    else:
        # More than 1 object
        list_remaining_objects = list_objects
        while len(list_remaining_objects) > 0:
            index = 0
            current_object = list_remaining_objects[0]
            list_remaining_objects = list_remaining_objects[1:]
            list_sorted_objects = []
            list_sorted_objects.append(current_object)
            # print "Before loop: ", current_object, list_remaining_objects, iNumberOfRemainingObjects
            while index < len(list_remaining_objects):
                compare_object = list_remaining_objects[index]
                are_equal = method_comparison(current_object, compare_object)
                # print "   In the loop:", index, compare_object, are_equal, list_sorted_objects, list_remaining_objects
                if are_equal:
                    list_sorted_objects.append(compare_object)
                    list_remaining_objects_tmp = list_remaining_objects[0:index]
                    if index < len(list_remaining_objects):
                        list_remaining_objects_tmp.extend(
                            list_remaining_objects[index + 1 :]
                        )
                    list_remaining_objects = list_remaining_objects_tmp
                else:
                    index += 1
            list_results.append(list_sorted_objects)
    return list_results


def sortSubWedgesOnExperimentalCondition(listSubWedge):
    """
    This method sorts a list of sub wedges into a new list containing lists of
    sub wegdes with identical experimental conditions.
    """
    # Sort it
    listSubWedgeSorted = sortIdenticalObjects(
        listSubWedge, isSameExperimentalConditionInSubWedge
    )
    return listSubWedgeSorted


def mergeTwoSubWedgesAdjascentInRotationAxis(sub_wedge_1, sub_wedge_2):
    """
    This method takes as input two sub wedges and merges them to an unique subwedge, if possible,
    and returns the resulting merged sub wedge. If the merge is not possible a None is returned.
    """
    sub_wedge_merged = None
    # First check that the two sub wedges have identical experimental conditions
    if isSameExperimentalConditionInSubWedge(sub_wedge_1, sub_wedge_2):
        # Check if sub wedges are adjascent:
        rotation_axis_end_1 = sub_wedge_1["experimentalCondition"]["goniostat"][
            "rotationAxisEnd"
        ]
        rotation_axis_start_2 = sub_wedge_2["experimentalCondition"]["goniostat"][
            "rotationAxisStart"
        ]
        if compareTwoValues(rotation_axis_end_1, rotation_axis_start_2, 0.001):
            # Same sub wedge! Let's merge them
            sub_wedge_merged = json.loads(json.dumps(sub_wedge_1))
            sub_wedge_2 = json.loads(json.dumps(sub_wedge_2))
            rotation_axis_end_2 = sub_wedge_2["experimentalCondition"]["goniostat"][
                "rotationAxisEnd"
            ]
            sub_wedge_merged["experimentalCondition"]["goniostat"][
                "rotationAxisEnd"
            ] = rotation_axis_end_2
            for image in sub_wedge_2["image"]:
                sub_wedge_merged["image"].append(image)
    return sub_wedge_merged


def mergeListOfSubWedgesWithAdjascentRotationAxis(list_of_sub_wegdes):
    """
    This method merges sub wedges in a list if they are adjascent in phi.
    """
    # Copy the incoming list to a new list
    list_of_sub_wegdes_copy = json.loads(json.dumps(list_of_sub_wegdes))
    list_of_merged_sub_wedges = []
    if len(list_of_sub_wegdes_copy) == 0:
        pass
    elif len(list_of_sub_wegdes_copy) == 1:
        list_of_merged_sub_wedges = list_of_sub_wegdes_copy
    else:
        # First sort the list as function of rotation axis start
        list_of_sub_wegdes_copy.sort(
            key=lambda x: x["experimentalCondition"]["goniostat"]["rotationAxisStart"]
        )
        # Then loop through the subwedges and merge them if possible
        list_of_remaining_sub_wedges = list_of_sub_wegdes_copy
        current_sub_wedge = list_of_remaining_sub_wedges[0]
        list_of_remaining_sub_wedges = list_of_remaining_sub_wedges[1:]
        while len(list_of_remaining_sub_wedges) > 0:
            next_sub_wedge = list_of_remaining_sub_wedges[0]
            list_of_remaining_sub_wedges = list_of_remaining_sub_wedges[1:]
            sub_wedge_merged = mergeTwoSubWedgesAdjascentInRotationAxis(
                current_sub_wedge, next_sub_wedge
            )
            if sub_wedge_merged is None:
                list_of_merged_sub_wedges.append(current_sub_wedge)
                current_sub_wedge = next_sub_wedge
            else:
                current_sub_wedge = sub_wedge_merged
        list_of_merged_sub_wedges.append(current_sub_wedge)
    return list_of_merged_sub_wedges
