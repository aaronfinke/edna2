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
# mxv1/src/EDHandlerESRFPyarchv1_0.py

import os
import subprocess
import time
import pathlib
import tempfile

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()

DEFAULT_TIMEOUT = 120  # s


def getWorkingDirectory(task, inData, workingDirectorySuffix=None):
    parentDirectory = inData.get("workingDirectory", None)
    if parentDirectory is None:
        parentDirectory = os.getcwd()
    parentDirectory = pathlib.Path(parentDirectory)
    if workingDirectorySuffix is None:
        # Create unique directory
        workingDirectory = tempfile.mkdtemp(
            prefix=task.__class__.__name__ + "_", dir=parentDirectory
        )
        os.chmod(workingDirectory, 0o755)
        workingDirectory = pathlib.Path(workingDirectory)
    else:
        # Here we assume that the user knows what he is doing and there's no
        # race condition for creating the working directory!
        workingDirectoryName = (
            task.__class__.__name__ + "_" + str(workingDirectorySuffix)
        )
        workingDirectory = parentDirectory / workingDirectoryName
        index = 1
        while workingDirectory.exists():
            workingDirectoryName = (
                task.__class__.__name__
                + "_"
                + str(workingDirectorySuffix)
                + "_{0:02d}".format(index)
            )
            workingDirectory = parentDirectory / workingDirectoryName
            index += 1
        workingDirectory.mkdir(mode=0o775, parents=True, exist_ok=False)
    workingDirectory = stripDataDirectoryPrefix(workingDirectory)
    return workingDirectory


def createPyarchFilePath(file_path):
    """
    This method translates from an ESRF "visitor" path to a "pyarch" path:
    /data/visitor/mx415/id14eh1/20100209 -> /data/pyarch/2010/id14eh1/mx415/20100209
    """
    pyarch_file_path = None
    if isinstance(file_path, str):
        file_path = pathlib.Path(file_path)
    list_of_directories = list(file_path.parts)
    if UtilsConfig.isEMBL():
        if "p13" in list_of_directories[0:3] or "P13" in list_of_directories[0:3]:
            pyarch_file_path = os.path.join("/data/ispyb/p13", *list_of_directories[4:])
        else:
            pyarch_file_path = os.path.join("/data/ispyb/p14", *list_of_directories[4:])
        return pyarch_file_path
    list_beamlines = [
        "bm07",
        "id23eh1",
        "id23eh2",
        "id29",
        "id30a1",
        "id30a2",
        "id30a3",
        "id30b",
    ]

    if (
        "data" in list_of_directories
        and len(list_of_directories) > 5
        and list_of_directories[1] != "data"
    ):
        while list_of_directories[1] != "data" and len(list_of_directories) > 5:
            del list_of_directories[1]

    # Check that we have at least four levels of directories:
    if len(list_of_directories) > 5:
        data_directory = list_of_directories[1]
        second_directory = list_of_directories[2]
        third_directory = list_of_directories[3]
        fourth_directory = list_of_directories[4]
        fifth_directory = list_of_directories[5]
        year = fifth_directory[0:4]
        proposal = None
        beamline = None
        if data_directory == "data" and second_directory == "gz":
            if third_directory == "visitor":
                proposal = fourth_directory
                beamline = fifth_directory
            elif fourth_directory == "inhouse":
                proposal = fifth_directory
                beamline = third_directory
            else:
                raise RuntimeError(
                    "Illegal path for UtilsPath.createPyarchFilePath: "
                    + "{0}".format(file_path)
                )
            list_of_remaining_directories = list_of_directories[6:]
        elif data_directory == "data" and second_directory == "visitor":
            proposal = list_of_directories[3]
            beamline = list_of_directories[4]
            list_of_remaining_directories = list_of_directories[5:]
        elif data_directory == "data" and second_directory in list_beamlines:
            beamline = second_directory
            proposal = list_of_directories[4]
            list_of_remaining_directories = list_of_directories[5:]
        if proposal is not None and beamline is not None:
            pyarch_file_path = pathlib.Path("/data/pyarch") / year / beamline
            pyarch_file_path = pyarch_file_path / proposal
            for directory in list_of_remaining_directories:
                pyarch_file_path = pyarch_file_path / directory
    if pyarch_file_path is None:
        logger.warning(
            "UtilsPath.createPyarchFilePath: path not converted for"
            + " pyarch: %s " % file_path
        )
    return pyarch_file_path


def waitForFile(file, expectedSize=None, timeOut=DEFAULT_TIMEOUT):
    file_path = pathlib.Path(file)
    final_size = None
    has_timed_out = False
    should_continue = True
    file_dir = file_path.parent
    if os.name != "nt" and file_dir.exists():
        # Patch provided by Sebastien 2018/02/09 for forcing NFS cache:
        # logger.debug("NFS cache clear, doing os.fstat on directory {0}".format(fileDir))
        fd = os.open(file_dir.as_posix(), os.O_DIRECTORY)
        stat_result = os.fstat(fd)
        os.close(fd)
        # logger.debug("Results of os.fstat: {0}".format(statResult))
    # Check if file is there
    if file_path.exists():
        file_size = file_path.stat().st_size
        if expectedSize is not None:
            # Check size
            if file_size > expectedSize:
                should_continue = False
        final_size = file_size
    if should_continue:
        logger.info("Waiting for file %s" % file_path)
        #
        time_start = time.time()
        while should_continue and not has_timed_out:
            if os.name != "nt" and file_dir.exists():
                # Patch provided by Sebastien 2018/02/09 for forcing NFS cache:
                # logger.debug("NFS cache clear, doing os.fstat on directory {0}".format(fileDir))
                fd = os.open(file_dir.as_posix(), os.O_DIRECTORY)
                stat_result = os.fstat(fd)  # noqa F841
                os.close(fd)
                # logger.debug("Results of os.fstat: {0}".format(statResult))
            time_elapsed = time.time() - time_start
            # Check if time out
            if time_elapsed > timeOut:
                has_timed_out = True
                str_warning = "Timeout while waiting for file %s" % file_path
                logger.warning(str_warning)
            else:
                # Check if file is there
                if file_path.exists():
                    file_size = file_path.stat().st_size
                    if expectedSize is not None:
                        # Check that it has right size
                        if file_size > expectedSize:
                            should_continue = False
                    else:
                        should_continue = False
                    final_size = file_size
            if should_continue:
                # Sleep 1 s
                time.sleep(1)
    return has_timed_out, final_size


def stripDataDirectoryPrefix(data_directory):
    """ Removes any paths before /data/..., e.g. /gpfs/easy/data/..."""
    list_paths = str(data_directory).split(os.sep)
    if "data" in list_paths:
        while list_paths[1] != "data":
            list_paths = [list_paths[0]] + list_paths[2:]
        new_data_directory = os.sep.join(list_paths)
    else:
        new_data_directory = data_directory
    return pathlib.Path(new_data_directory)

def systemCopyFile(from_path, to_path):
    p = subprocess.Popen(["cp", from_path, to_path])
    p.wait()

def systemRmTree(treePath, ignore_errors=False):
    try:
        if ignore_errors:
            subprocess.check_call(f"rm -rf {treePath}", shell=True)
        else:
            subprocess.check_call(f"rm -r {treePath} 2>&1 > /dev/null", shell=True)
    except subprocess.CalledProcessError:
        if not ignore_errors:
            raise

def systemCopyTree(from_path, to_path, dirs_exists_ok=False):
    if os.path.exists(to_path):
        if dirs_exists_ok:
            systemRmTree(to_path)
        else:
            raise FileExistsError(to_path)
    p = subprocess.Popen(["cp", "-r", from_path, to_path])
    p.wait()
