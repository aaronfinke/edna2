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

import os
import pathlib
import configparser


def getConfigDir():
    """
    Returns the directory used for configuration.
    If EDNA2_CONFIG is defined the value is returned,
    otherwise <project_base>/config.
    """
    if "EDNA2_CONFIG" in os.environ:
        configDir = pathlib.Path(os.environ["EDNA2_CONFIG"])
    else:
        pathFile = pathlib.Path(__file__)
        pathProjectBase = pathFile.parents[3]
        configDir = pathProjectBase / "config"
    return configDir


def getSite():
    """
    Returns the EDNA2_SITE variable from environment.
    If not defined returns 'Default'.
    """
    site = "Default"
    if "EDNA2_SITE" in os.environ:
        site = os.environ["EDNA2_SITE"]
    return site


def setSite(site):
    """
    Sets the EDNA2_SITE variable.
    """
    os.environ["EDNA2_SITE"] = site


def getConfig(site=None):
    config = configparser.ConfigParser()
    if site is None:
        site = getSite()
    configFile = site + ".ini"
    configDir = getConfigDir()
    configPath = configDir / configFile.lower()
    if configPath.exists():
        config.read(configPath.as_posix())
    return config


def getTaskConfig(taskName, site=None):
    dictConfig = {}
    config = getConfig(site)
    sections = config.sections()
    # First search in included configs
    if "Include" in sections:
        for site in config["Include"]:
            dictConfig.update(getTaskConfig(taskName, site))
    # Then update with the current config
    if taskName in sections:
        dictConfig.update(dict(config[taskName]))
    # Substitute ${} from os.environ
    for key in dictConfig:
        dictConfig[key] = os.path.expandvars(dictConfig[key])
    return dictConfig


def isEMBL():
    """
    Returns true if EMBL config
    """
    return getSite().lower().startswith("embl")


def isESRF():
    """
    Returns true if ESRF config
    """
    return getSite().lower().startswith("esrf")


def get(task, parameterName, defaultValue=None):
    if isinstance(task, str):
        taskConfig = getTaskConfig(task)
    else:
        taskConfig = getTaskConfig(task.__class__.__name__)
    value = taskConfig.get(parameterName.lower(), defaultValue)
    return value
