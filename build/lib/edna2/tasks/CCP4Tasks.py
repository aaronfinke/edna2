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
__date__ = "10/05/2019"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxPluginExec/plugins/EDPluginGroupCCP4-v1.0/plugins/EDPluginExecAimlessv1_0.py
# mxPluginExec/plugins/EDPluginGroupCCP4-v1.0/plugins/EDPluginExecPointlessv1_0.py

import os
import re

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()


class AimlessTask(AbstractTask):
    """
    Execution of CCP4 aimless
    """

    def run(self, inData):
        outData = {}
        input_file = inData['input_file']
        output_file = inData['output_file']
        symoplib = UtilsConfig.get('CCP4', 'symoplib')
        # ccp4setup = UtilsConfig.get('CCP4', 'ccp4setup')
        # if ccp4setup is None:
        #     commandLine = ""
        # else:
        #     commandLine = ". " + ccp4setup + '\n'
        commandLine = ""
        if symoplib is None:
            commandLine += 'aimless HKLIN {0} HKLOUT {1}'.format(
                input_file, output_file)
        else:
            commandLine += 'aimless HKLIN {0} HKLOUT {1} SYMINFO {2}'.format(
                input_file, output_file, symoplib)
        logger.info("Command line: {0}".format(commandLine))
        start_image = inData['start_image']
        end_image = inData['end_image']
        projectName = inData.get('dataCollectionID', 'EDNA_proc')
        resolution = inData.get('res', 0.0)
        anom = inData['anom']
        listCommand = [
            'bins 15',
            'run 1 batch {0} to {1}'.format(start_image, end_image),
            'name run 1 project {0} crystal DEFAULT dataset NATIVE'.format(projectName),
            'scales constant',
            'resolution 50 {0}'.format(resolution),
            'cycles 100',
            'anomalous {0}'.format('ON' if anom else 'OFF'),
            'output MERGED UNMERGED',
            'END'
            ]
        self.setLogFileName('aimless.log')
        self.runCommandLine(commandLine, listCommand=listCommand)
        outData['isSuccess'] = os.path.exists(output_file)
        return outData

class PointlessTask(AbstractTask):
    """
    Executes the CCP4 program pointless
    """

    def run(self, inData):
        input_file = inData['input_file']
        output_file = inData['output_file']
        commandLine = 'pointless'
        if UtilsConfig.isEMBL():
            commandLine += ' -c'
        commandLine += " xdsin {0} hklout {1}".format(input_file, output_file)
        listCommand = ['setting symmetry-based']
        if 'choose_spacegroup' in inData:
            listCommand += 'choose spacegroup {0}'.format(
                inData['choose_spacegroup'])
        self.setLogFileName('pointless.log')
        self.runCommandLine(commandLine, listCommand=listCommand)
        outData = self.parsePointlessOutput(self.getLogPath())
        return outData

    @classmethod
    def parsePointlessOutput(cls, logPath):
        sgre = re.compile(""" \* Space group = '(?P<sgstr>.*)' \(number\s+(?P<sgnumber>\d+)\)""")
        outData = {'isSuccess': False}
        if logPath.exists():
            with open(str(logPath)) as f:
                log = f.read()
            m = sgre.search(log)
            if m is not None:
                d = m.groupdict()
                sgnumber = d['sgnumber']
                sgstr = d['sgstr']

                outData['sgnumber'] = int(sgnumber)
                outData['sgstr'] = sgstr
                outData['isSuccess'] = True
                # Search first for unit cell after the Laue group...
                unitCellRe = re.compile("""  Laue group confidence.+\\n\\n\s+Unit cell:(.+)""")
                m2 = unitCellRe.search(log)
                if m2 is None:
                    # Then search it from the end...
                    unitCellRe = re.compile(""" \* Cell Dimensions : \(obsolete \- refer to dataset cell dimensions above\)\\n\\n(.+)""")
                    m2 = unitCellRe.search(log)
                if m2 is not None:
                    listCell = m2.groups()[0].split()
                    cell = {
                        'length_a': float(listCell[0]),
                        'length_b': float(listCell[1]),
                        'length_c': float(listCell[2]),
                        'angle_alpha': float(listCell[3]),
                        'angle_beta': float(listCell[4]),
                        'angle_gamma': float(listCell[5])
                    }
                    outData['cell'] = cell
        return outData

