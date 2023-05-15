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
import time
import graypy
import logging
import logging.handlers

from edna2.utils import UtilsConfig


def addGrayLogHandler(logger):
    server = UtilsConfig.get("Logging", "graylog_server")
    port = UtilsConfig.get("Logging", "graylog_port")
    if server is not None and port is not None:
        graylogHandler = graypy.GELFUDPHandler(server, int(port))
        logger.addHandler(graylogHandler)


def addStreamHandler(logger):
    streamHandler = logging.StreamHandler()
    logFileFormat = "%(asctime)s %(levelname)-8s %(message)s"
    formatter = logging.Formatter(logFileFormat)
    streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)


def addFileHandler(logger):
    logPath = UtilsConfig.get("Logging", "log_file_path")
    if logPath is not None:
        if "DATE" in logPath:
            logPath = logPath.replace(
                "DATE", time.strftime("%Y-%m-%d", time.localtime(time.time()))
            )
        if not os.path.exists(os.path.dirname(logPath)):
            os.makedirs(os.path.dirname(logPath))
        maxBytes = int(UtilsConfig.get("Logging", "log_file_maxbytes", 1e7))
        backupCount = int(UtilsConfig.get("Logging", "log_file_backupCount", 0))
        fileHandler = logging.handlers.RotatingFileHandler(
            logPath, maxBytes=maxBytes, backupCount=backupCount
        )
        logFileFormat = UtilsConfig.get("Logging", "log_file_format")
        if logFileFormat is None:
            logFileFormat = (
                "%(asctime)s %(module)-20s %(funcName)-15s %(levelname)-8s %(message)s"
            )
        formatter = logging.Formatter(logFileFormat)
        fileHandler.setFormatter(formatter)
        logger.addHandler(fileHandler)


def setLoggingLevel(logger, level):
    if level is None:
        level = UtilsConfig.get("Logging", "level")
    if level is None:
        level = "INFO"
    if level == "DEBUG":
        loggingLevel = logging.DEBUG
    elif level == "INFO":
        loggingLevel = logging.INFO
    elif level == "WARNING":
        loggingLevel = logging.WARNING
    elif level == "ERROR":
        loggingLevel = logging.ERROR
    elif level == "CRITICAL":
        loggingLevel = logging.CRITICAL
    elif level == "FATAL":
        loggingLevel = logging.FATAL
    else:
        raise RuntimeError('Unknown logging level: "{0}"'.format(level))
    logger.setLevel(loggingLevel)


def getLogger(level=None):
    logger = logging.getLogger("edna2")
    # Check if handlers already added:
    hasGraylogHandler = False
    hasStreamHandler = False
    hasFileHandler = False
    for handler in logger.handlers:
        if isinstance(handler, graypy.GELFUDPHandler):
            hasGraylogHandler = True
        elif isinstance(handler, logging.handlers.RotatingFileHandler):
            hasFileHandler = True
        elif isinstance(handler, logging.StreamHandler):
            hasStreamHandler = True
    if not hasGraylogHandler:
        addGrayLogHandler(logger)
    if not hasStreamHandler:
        addStreamHandler(logger)
    if not hasFileHandler:
        addFileHandler(logger)
    # Set level
    setLoggingLevel(logger, level)
    return logger
