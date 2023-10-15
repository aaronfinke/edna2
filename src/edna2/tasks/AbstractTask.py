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
import signal
import psutil
import json
import pathlib
import billiard
import traceback
import jsonschema
import subprocess
import socket

from edna2.utils import UtilsPath
from edna2.utils import UtilsLogging
from edna2.utils import UtilsConfig

logger = UtilsLogging.getLogger()


class EDNA2Process(billiard.Process):
    """
    See https://stackoverflow.com/a/33599967.
    """

    def __init__(self, *args, **kwargs):
        billiard.Process.__init__(self, *args, **kwargs)
        self._pconn, self._cconn = billiard.Pipe()
        self._exception = None

    def run(self):
        try:
            billiard.Process.run(self)
            self._cconn.send(None)
        except BaseException as e:
            tb = traceback.format_exc()
            self._cconn.send((e, tb))
    
    def timeOut(self):
        raise billiard.TimeoutError

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception


class AbstractTask():  # noqa R0904
    """
    Parent task to all EDNA2 tasks.
    """

    def __init__(self, inData, workingDirectorySuffix=None):
        self._dictInOut = billiard.Manager().dict()
        self._dictInOut["inData"] = json.dumps(inData, default=str)
        self._dictInOut["outData"] = json.dumps({})
        self._dictInOut["isFailure"] = False
        self._dictInOut["timeOut"] = None
        self._dictInOut["timeoutExit"] = False
        self._process = EDNA2Process(target=self.executeRun, args=())
        self._workingDirectorySuffix = workingDirectorySuffix
        self._workingDirectory = None
        self._logFileName = None
        self._errorLogFileName = None
        self._slurmLogFileName = None
        self._slurmErrorLogFileName = None
        self._schemaPath = pathlib.Path(__file__).parents[1] / "schema"
        self._persistInOutData = True
        self._oldDir = os.getcwd()
        self._slurmId = None
        self._slurmHostname = None
        self._jobName = type(self).__name__

    def getSchemaUrl(self, schemaName):
        return "file://" + str(self._schemaPath / schemaName)

    def executeRun(self):
        inData = self.getInData()
        hasValidInDataSchema = False
        hasValidOutDataSchema = False
        if self.getInDataSchema() is not None:
            instance = inData
            schema = self.getInDataSchema()
            try:
                jsonschema.validate(instance=instance, schema=schema)
                hasValidInDataSchema = True
            except Exception as e:
                logger.exception(e)
        else:
            hasValidInDataSchema = True
        if hasValidInDataSchema:
            self._workingDirectory = UtilsPath.getWorkingDirectory(
                self, inData, workingDirectorySuffix=self._workingDirectorySuffix
            )
            self.writeInputData(inData)
            self._oldDir = os.getcwd()
            os.chdir(str(self._workingDirectory))
            outData = self.run(inData)
            os.chdir(self._oldDir)
        else:
            raise RuntimeError("Schema validation error for inData")
        if self.getOutDataSchema() is not None:
            instance = outData
            schema = self.getOutDataSchema()
            try:
                jsonschema.validate(instance=instance, schema=schema)
                hasValidOutDataSchema = True
            except Exception as e:
                logger.exception(e)
        else:
            hasValidOutDataSchema = True
        if hasValidOutDataSchema:
            self.writeOutputData(outData)
        else:
            raise RuntimeError("Schema validation error for outData")
        if not os.listdir(str(self._workingDirectory)):
            os.rmdir(str(self._workingDirectory))

    def getInData(self):
        return json.loads(self._dictInOut["inData"])

    def setInData(self, inData):
        self._dictInOut["inData"] = json.dumps(inData, default=str)

    inData = property(getInData, setInData)

    def getOutData(self):
        return json.loads(self._dictInOut["outData"])

    def setOutData(self, outData):
        self._dictInOut["outData"] = json.dumps(outData, default=str)

    outData = property(getOutData, setOutData)

    def writeInputData(self, inData):
        # Write input data
        if self._persistInOutData and self._workingDirectory is not None:
            jsonName = "inData" + self.__class__.__name__ + ".json"
            with open(str(self._workingDirectory / jsonName), "w") as f:
                f.write(json.dumps(inData, default=str, indent=4))

    def writeOutputData(self, outData):
        self.setOutData(outData)
        if self._persistInOutData and self._workingDirectory is not None:
            jsonName = "outData" + self.__class__.__name__ + ".json"
            with open(str(self._workingDirectory / jsonName), "w") as f:
                f.write(json.dumps(outData, default=str, indent=4))
    
    def setTimeout(self,timeOut):
        if timeOut is not None:
            timeOut = float(timeOut)
        self._dictInOut["timeOut"] = timeOut

    def getTimeout(self):
        return self._dictInOut.get("timeOut", None)

    def getLogPath(self):
        if self._logFileName is None:
            self._logFileName = self.__class__.__name__ + ".log.txt"
        logPath = self._workingDirectory / self._logFileName
        return logPath
    
    def getSlurmLogPath(self):
        if self._slurmLogFileName is None:
            return None
        slurmLogPath = self._workingDirectory / self._slurmLogFileName
        return slurmLogPath
    
    def getSlurmErrorLogPath(self):
        if self._slurmErrorLogFileName is None:
            return None
        slurmErrorLogPath = self._workingDirectory / self._slurmErrorLogFileName
        return slurmErrorLogPath

    def setLogFileName(self, logFileName):
        self._logFileName = logFileName

    def getLogFileName(self):
        return self._logFileName

    def getErrorLogPath(self):
        if self._errorLogFileName is None:
            self._errorLogFileName = self.__class__.__name__ + ".error.txt"
        errorLogPath = self._workingDirectory / self._errorLogFileName
        return errorLogPath
    
    def getSlurmLogFileName(self):
        return self._slurmLogFileName
    
    def getSlurmErrorLogFileName(self):
        return self._slurmErrorLogFileName
    
    def setSlurmLogFileName(self, slurmLogFileName):
        self._slurmLogFileName = slurmLogFileName

    def setSlurmErrorLogFileName(self, slurmErrorLogFileName):
        self._slurmErrorLogFileName = slurmErrorLogFileName

    def setErrorLogFileName(self, errorLogFileName):
        self._errorLogFileName = errorLogFileName

    def getErrorLogFileName(self):
        return self._errorLogFileName
        
    def getLog(self):
        with open(str(self.getLogPath())) as f:
            log = f.read()
        return log

    def getErrorLog(self):
        with open(str(self.getErrorLogPath())) as f:
            errorLog = f.read()
        return errorLog

    def getSlurmLog(self):
        if self._slurmLogFileName is None:
            return None
        with open(self.getSlurmLogPath()) as f:
            slurmLog = f.read()
        return slurmLog
    
    def getSlurmErrorLog(self):
        if self._slurmErrorLogFileName is None:
            return None
        with open(self.getSlurmErrorLogPath()) as f:
            slurmLog = f.read()
        return slurmLog
    
    def getSlurmId(self):
        return self._slurmId
    
    def getSlurmHostname(self):
        return self._slurmHostname
    
    def setSlurmHostname(self,hostname):
        self._slurmHostname = hostname


    def submitCommandLine(self, commandLine, ignoreErrors, timeout=UtilsConfig.get("Slurm","time","01:00:00"), partition=None, jobName="EDNA2"):
        jobName = "EDNA2_" + self._jobName
        exclusive = UtilsConfig.get("Slurm","is_exclusive",False)
        nodes = UtilsConfig.get("Slurm","nodes",1)
        core = UtilsConfig.get("Slurm","cores",10)
        mem = UtilsConfig.get("Slurm","mem",4000)
        workingDir = str(self._workingDirectory)
        if workingDir.startswith("/mntdirect/_users"):
            workingDir = workingDir.replace("/mntdirect/_users", "/home/esrf")
        script = "#!/bin/bash\n"
        script += '#SBATCH --job-name="{0}"\n'.format(jobName)
        script += "#SBATCH --partition={0}\n".format(partition) if partition else ""
        script += "#SBATCH --exclusive\n" if exclusive else ""
        script += "#SBATCH --mem={0}\n".format(mem) if not exclusive else ""
        script += "#SBATCH --nodes={0}\n".format(nodes)
        # script += "#SBATCH --nodes=1\n"  # Necessary for not splitting jobs! See ATF-57
        script += "#SBATCH --cpus-per-task={0}\n".format(core) if not exclusive else ""
        script += "#SBATCH --time={0}\n".format(timeout) if timeout else "0"
        script += "#SBATCH --chdir={0}\n".format(workingDir)
        script += f"#SBATCH --output={jobName}_%j.out\n"
        script += f"#SBATCH --error={jobName}_%j.err\n"
        script += commandLine + "\n"
        shellFile = self._workingDirectory / (jobName + "_slurm.sh")
        with open(str(shellFile), "w") as f:
            f.write(script)
            f.close()
        shellFile.chmod(0o755)
        slurmCommandLine = "sbatch --wait {0}".format(shellFile)
        pipes = subprocess.Popen(
            slurmCommandLine,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            start_new_session=True,
            cwd=str(self._workingDirectory),
        )
        while True:
            line = pipes.stdout.readline()
            if line:
                break
        line = line.decode('utf-8')
        if "Submitted batch job" in line:
            self._slurmJobId = int(line.split()[-1])
        if self._slurmJobId:
            self.setSlurmLogFileName(f"{jobName}_{self._slurmJobId}.out")
            self.setSlurmErrorLogFileName(f"{jobName}_{self._slurmJobId}.err")
            slurm_hostname = ""
            try:
                squeue = subprocess.run(["squeue","-j",str(self._slurmJobId)], capture_output=True)
                slurm_hostname = squeue.stdout.decode('utf-8').strip('\n').split()[-1]
                while slurm_hostname == "(None)":
                    squeue = subprocess.run(["squeue","-j",str(self._slurmJobId)], capture_output=True)
                    slurm_hostname = squeue.stdout.decode('utf-8').strip('\n').split()[-1]
            except Exception as e:
                logger.error(f"Failed to get slurm host: {e}")
            if slurm_hostname:
                self.setSlurmHostname(slurm_hostname)
            logger.debug(f"Job {jobName} submitted to slurm on host {self._slurmHostname} with slurmJobId {self._slurmJobId}")


        stdout, stderr = pipes.communicate()
        # slurmLogPath = self._workingDirectory / (jobName + "_slurm.log")
        # slurmErrorLogPath = self._workingDirectory / (jobName + "_slurm.error.log")
        slurmLogPath = self.getLogPath()
        slurmErrorLogPath = self.getErrorLogPath()
        if len(stdout) > 0:
            log = str(stdout, "utf-8")
            with open(str(slurmLogPath), "w") as f:
                f.write(log)
        if len(stderr) > 0:
            if not ignoreErrors:
                logger.warning(
                    "Error messages from command {0}".format(commandLine.split(" ")[0])
                )
            with open(str(slurmErrorLogPath), "w") as f:
                f.write(str(stderr, "utf-8"))
        if pipes.returncode != 0:
            # Error!
            warningMessage = "{0}, code {1}".format(stderr, pipes.returncode)
            logger.warning(warningMessage)
            # raise RuntimeError(errorMessage)
        return pipes.returncode

    def runCommandLine(
        self,
        commandLine,
        logPath=None,
        listCommand=None,
        ignoreErrors=False,
        doSubmit=False,
        partition=None,
    ):
        if logPath is None:
            logPath = self.getLogPath()
        jobName = type(self).__name__
        logFileName = os.path.basename(logPath)
        errorLogPath = self.getErrorLogPath()
        errorLogFileName = os.path.basename(errorLogPath)
        commandLine += " 1>{0} 2>{1}".format(logFileName, errorLogFileName)
        if listCommand:
            commandLine += " << EOF-EDNA2\n"
            for command in listCommand:
                commandLine += command + "\n"
            commandLine += "EOF-EDNA2"
        commandLogFileName = jobName + ".commandLine.txt"
        commandLinePath = self._workingDirectory / commandLogFileName
        with open(str(commandLinePath), "w") as f:
            f.write(commandLine)
        if doSubmit:
            self.submitCommandLine(commandLine, partition, ignoreErrors)
        else:
            pipes = subprocess.Popen(
                commandLine,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                close_fds=True,
                start_new_session=True,
                cwd=str(self._workingDirectory),
            )
            logger.debug(f"pid = {pipes.pid}")
            stdout, stderr = pipes.communicate()
            if len(stdout) > 0:
                log = str(stdout, "utf-8")
                with open(str(logPath), "w") as f:
                    f.write(log)
            if len(stderr) > 0:
                if not ignoreErrors:
                    logger.warning(
                        "Error messages from command {0}".format(
                            commandLine.split(" ")[0]
                        )
                    )
                errorLogPath = self._workingDirectory / errorLogFileName
                with open(str(errorLogPath), "w") as f:
                    f.write(str(stderr, "utf-8"))
            if pipes.returncode != 0:
                # Error!
                errorMessage = "{0}, code {1}".format(stderr, pipes.returncode)
                raise RuntimeError(errorMessage)

    def onError(self):
        pass

    def start(self):
        self._process.start()

    def join(self):
        timeOut = self.getTimeout()
        if timeOut is not None:
            logger.debug(f"timeout for {self.__class__.__name__}: {timeOut}")
        self._process.join(timeout=timeOut)
        # deal with timeouts
        if self._process.exitcode is None:
            # to ensure all subprocesses generated by process terminate cleanly
            # logger.debug(f"exitcode for {self} = {self._process.exitcode}")
            # logger.debug(f"pid for {self}: {self._process.pid}")
            current_process = psutil.Process(self._process.pid)
            children = current_process.children(recursive=True)
            for child in children:
                # logger.debug('Child pid for child {} is {}'.format(child,child.pid))
                try:
                    child.send_signal(signal.SIGTERM)
                except:
                    pass
            self._process.timeOut()
        if self._process.exception:
            error, trace = self._process.exception
            logger.error(error)
            logger.error(trace)
            self._dictInOut["isFailure"] = True
            if "TimeoutError" in trace.split()[-1]:
                self._dictInOut["timeoutExit"] = True
            self.onError()

    def execute(self):
        self.start()
        self.join()
    
    def setTimeoutExit(self):
        self._dictInOut["timeoutExit"] = True
    
    def isTimeoutExit(self):
        return self._dictInOut["timeoutExit"]

    def setFailure(self):
        self._dictInOut["isFailure"] = True

    def isFailure(self):
        return self._dictInOut["isFailure"]

    def isSuccess(self):
        return not self.isFailure()

    def getWorkingDirectory(self):
        return self._workingDirectory

    def setWorkingDirectory(self, inData):
        self._workingDirectory = UtilsPath.getWorkingDirectory(self, inData)

    def getInDataSchema(self):
        return None

    def getOutDataSchema(self):
        return None

    def setPersistInOutData(self, value):
        self._persistInOutData = value
