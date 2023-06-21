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

__authors__ = ["D. Fastus"]
__license__ = "MIT"
__date__ = "12/06/2023"

import pathlib
import sys
import pickle
import os
import json
import shutil

from edna2.tasks.AbstractTask import AbstractTask
# from edna2.tasks.PhenixTasks import ProcPredTask
# from edna2.tasks.CCP4Tasks import DimpleTask
from edna2.utils import UtilsLogging

# import edna2.utils.UtilsPDB as UtilsPDB
from edna2.utils import UtilsPDB


logger = UtilsLogging.getLogger()

class AlphaFoldTask(AbstractTask):
    """
    Runs an AlphaFold2 prediction
    """

    def run(self, inData):
        fasta_path = inData.get("fasta_path")
        output_Dir = self._workingDirectory 
        # ALPHAFOLD_DATA_DIR = os.environ.get('ALPHAFOLD_DATA_DIR', None)
    
        try:
            with open(fasta_path, mode="r") as file:
                line = file.readline()

                if not line.startswith('>'):
                    logger.error("The input is not a fasta file!")
                    sys.exit(1)
                else:
                    line = line.strip()
                    fasta_name = line[1:5]

        except Exception as e:
            logger.error(f"{fasta_path} can not be open or does not exist!", exc_info=True)
            sys.exit(1)

        commandLine = 'module purge \n'
        commandLine += 'module add fosscuda/2020b AlphaFold \n'
        commandLine += 'export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b \n\n'

        commandLine += """alphafold \\
        --fasta_paths={0} \\
        --max_template_date=2021-11-01 \\
        --output_dir={1} \\
        --data_dir=$ALPHAFOLD_DATA_DIR""".format(fasta_path, output_Dir)

        logger.info("Command line: {0}".format(commandLine))
        self.runCommandLine(commandLine, ignoreErrors=True)
        # self.submitCommandLine(commandLine, jobName=f"{fasta_name}", ignoreErrors=True, mem=0, partition="v100", time="01-00:00")
        # self.monitorCommandLine(job=f"{fasta_name}_slurm.sh", name=f"AlphaFold prediction of {fasta_name}")

        # logPath = self.getWorkingDirectory() / 'AlphaFold.log'
        # outData = self.parseXtriageLogFile(logPath)
        outData = {}
        outData['isSuccess']  = True # self.check_out(out_dir=output_Dir)

        return outData

        # check if the output are complete and have the right format
    def check_out(self, out_dir):
        """
        Check if all output files are created (first one or ranked_0 is the most relevant in the pipeline)
        """
        files_to_check = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb', "ranking_debug.json"]

        # Check if all files are present in the directory
        for file in files_to_check:
            file_path = os.path.join(out_dir, file)

            if not os.path.isfile(file_path):
                logger.error("The files from the AlphaFold prediction are not successfully generated, the {0} is missing...".format(file), exc_info=True)
                return False
            
        logger.info("The files from the AlphaFold prediciton are successfully generated in {0}".format(out_dir))
        return True
    
    def parseAlphafoldPredictionLogFile(self, output_dir):
        """
        Extract AlphaFold prediction information
        """

            #  dict_keys (['aatype'
            # 'between_segment_residues',
            # 'domain_name',
            # 'residue index'
            # 'sed_length'
            # 'sequence',
            # 'deletion matrix int', 'msa', 'num alignments',
            # 'msa_species_identifiers',
            # 'template aatype',
            # 'template_all_atom_masks', 'template_all_atom_positions',
            # 'template_ domain_names'
            # 'template sequence',
            # 'template_sum_probs'])
        AlphaFoldResults = {
            "overall/features" : {
                "featuresTime" : None,
                "DomainName" : None,
                "Template_ domain_names": None
            },
            # dict_keys (['distogram', 'experimentally_ resolved', 'masked msa',
            # 'predicted_aligned_error', "predicted dot'. structure module"
            # 'olddt'. 'alianed confidence probs'
            # 'max_predicted_aligned_error', 'ptm',
            # 'iptm', 'ranking_confidence'])
            "ProteinModels": {
                "result_model_1(ranked_0)" : None,
                "result_model_2(ranked_1)" : None,
                "result_model_3(ranked_2)" : None,
                "result_model_4(ranked_3)" : None,
                "result_model_5(ranked_5)" : None,
            },
        }

        extract = []
        # open features.pkl
        try:
            feature_dict = pickle.load(open(f"{output_dir}/features.pkl", "rb"))
        except:
            logger.error("features.pkl file could not be parsed")
            return None
        
        AlphaFoldResults["overall/features"]["DomainName"] = feature_dict["domain_name"]
        AlphaFoldResults["overall/features"]["featuyy"] = feature_dict["domain_name"]
        AlphaFoldResults["overall/features"]["DomainName"] = feature_dict["domain_name"]
        AlphaFoldResults["overall/features"]["featuyy"] = feature_dict["domain_name"]
        AlphaFoldResults["overall/features"]["DomainName"] = feature_dict["domain_name"]
        AlphaFoldResults["overall/features"]["featuyy"] = feature_dict["domain_name"]

        # open ranking_debug and model information
        try:
            with open(f"{output_dir}/ranking_debug.json","r") as file:
                ranking_dict = json.load(file)
        except:
            logger.error("ranking_debug.json log file could not be parsed")
            return None
        
        try:
            with open(f"{output_dir}/timings.json","r") as file:
                timings_dict = json.load(file)
        except:
            logger.error("timings.json log file could not be parsed")
            return None
        
        for num in range(1,6):
            try:
                result_model = pickle.load(open(f"{output_dir}/result_model_{num}.pkl"),"rb")
            except:
                logger.error(f"result_model_{num}.pkl file could not be parsed")
                return None
        
        return AlphaFoldResults
        
class visualizeProteinPrediction(AbstractTask):
    """
    Visualize and interpret predicted proteins
    """
    pass
