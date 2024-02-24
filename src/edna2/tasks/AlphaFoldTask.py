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
from edna2.utils import UtilsConfig

# import edna2.utils.UtilsPDB as UtilsPDB
from edna2.utils import UtilsPDB


logger = UtilsLogging.getLogger()


class AlphaFoldTask(AbstractTask):
    """
    Runs an AlphaFold2 prediction
    """

    def run(self, inData):
        fasta_path = inData.get("FASTA_file")
        output_Dir = self._workingDirectory

        try:
            with open(fasta_path, mode="r") as file:
                line = file.read()

                if not line.startswith(">"):
                    logger.error("The input is not a fasta file!")
                    sys.exit(1)
                else:
                    if line.count(">") == 1:
                        monomer = True
                    else:
                        monomer = False

                    line = line.strip()
                    fasta_name = line[1:5]

        except Exception as e:
            logger.error(f"{fasta_path} can not be open or does not exist!", exc_info=True)
            sys.exit(1)

        commandLine = UtilsConfig.get("AlphaFold","alphaFold_env") if UtilsConfig.get("AlphaFold","alphaFold_env") else ""
        commandLine += " "
        commandLine += "alphafold "
        commandLine += f"--fasta_paths={fasta_path} "
        commandLine += f"--max_template_date=2021-11-01 "
        if monomer:
            commandLine += f"--model_preset=monomer "
        else:
            commandLine += f"--model_preset=multimer "
        commandLine += f"--output_dir={output_Dir} "
        commandLine += "--data_dir=$ALPHAFOLD_DATA_DIR "

        logger.info("Command line: {0}".format(commandLine))
        logPath = self.getWorkingDirectory() / "AlphaFold.log"
        self.runCommandLine(commandLine, ignoreErrors=True, logPath=logPath)
        # self.submitCommandLine(commandLine, jobName=f"{fasta_name}", ignoreErrors=True, mem=0, partition="v100", time="01-00:00")
        # self.monitorCommandLine(job=f"{fasta_name}_slurm.sh", name=f"AlphaFold prediction of {fasta_name}")

        outputDir = str(pathlib.Path(output_Dir).resolve()) + f"/{fasta_name}"
        outData = {}
        outData["outputDir"] = outputDir

        outData["AlphaFoldResults"] = self.parseAlphafoldResultFiles(outputDir)
        outData["isSuccess"] = self.check_out(outputDir)

        return outData

        # check if the output are complete and have the right format

    def check_out(self, output_dir):
        """
        Check if all output files are created (first one or ranked_0 is the most relevant in the pipeline)
        """
        files_to_check = [
            "ranked_0.pdb",
            "relaxed_model_1.pdb",
            "result_model_1.pkl",
            "unrelaxed_model_1.pdb",
            "ranking_debug.json",
        ]

        # Check if all files are present in the directory
        for file in files_to_check:
            file_path = os.path.join(output_dir, file)

            if not os.path.isfile(file_path):
                logger.error(
                    "The files from the AlphaFold prediction are not successfully generated, the {0} is missing...".format(
                        file
                    ),
                    exc_info=True,
                )
                return False

        logger.info("The files from the AlphaFold prediciton are successfully generated in {0}".format(output_dir))
        return True

    def parseAlphafoldLogFile(self, pathToLogFile: pathlib.Path):
        outData = {"logPath": str(pathToLogFile)}

        if pathToLogFile.is_file():
            with open(pathToLogFile, "r") as logfile:
                pass

    def parseAlphafoldResultFiles(self, output_dir):
        """
        Extract AlphaFold prediction information
        """
        # overall/features information
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

        # ProteinModels / model speficic information
        # dict_keys (['distogram', 'experimentally_ resolved', 'masked msa',
        # 'predicted_aligned_error', "predicted dot'. structure module"
        # 'olddt'. 'alianed confidence probs'
        # 'max_predicted_aligned_error', 'ptm',
        # 'iptm', 'ranking_confidence'])

        AlphaFoldResults = {
            "overall/features": None,
            "ProteinModels": None,
        }

        # open features.pkl
        try:
            feature_dict = pickle.load(open(f"{output_dir}/features.pkl", "rb"))
        except:
            logger.error("features.pkl file could not be parsed")
            return None

        AlphaFoldResults["overall/features"]["DomainName"] = feature_dict["domain_name"]
        AlphaFoldResults["overall/features"]["SequenceLength"] = feature_dict["seq_length"][0]
        AlphaFoldResults["overall/features"]["Sequence"] = feature_dict["sequence"]
        AlphaFoldResults["overall/features"]["NumberOfAlignments"] = feature_dict["num_alignments"]

        AlphaFoldResults["overall/features"]["TemplateDomains"] = {}

        AlphaFoldResults["overall/features"]["TemplateDomains"]["NumberOfTemplateDomains"] = feature_dict[
            "template_sum_probs"
        ]

        for name, sequence in zip(feature_dict["template_domain_names"], feature_dict["template_sequence"]):
            AlphaFoldResults["overall/features"]["TemplateDomains"][name] = sequence

        # open ranking_debug and model information
        try:
            with open(f"{output_dir}/ranking_debug.json", "r") as file:
                ranking_dict = json.load(file)

                for model in ranking_dict["order"]:
                    if model in ranking_dict["plddts"].keys():
                        key = f"result_{model}(ranked_{ranking_dict['order'].index(model)})"
                        AlphaFoldResults["ProteinModels"][key] = {"plddt": ranking_dict["plddts"][model]}
        except:
            logger.error("ranking_debug.json log file could not be parsed")
            return None

        try:
            with open(f"{output_dir}/timings.json", "r") as file:
                timings_dict = json.load(file)

                AlphaFoldResults["overall/features"]["featuresTime"] = timings_dict["features"]
        except:
            logger.error("timings.json log file could not be parsed")
            return None

        for num in range(1, 6):
            try:
                result_model = pickle.load(open(f"{output_dir}/result_model_{num}.pkl"), "rb")
            except:
                logger.error(f"result_model_{num}.pkl file could not be parsed")
                return None

        return AlphaFoldResults
