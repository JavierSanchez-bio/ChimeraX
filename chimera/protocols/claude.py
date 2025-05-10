# **************************************************************************
# *
# * Authors:   Javier Sanchez
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import re
import subprocess
import pyworkflow.utils as pwutils
from pyworkflow.protocol import MultiPointerParam
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects import SetOfAtomStructs, AtomStruct
from pyworkflow.protocol import String

from chimera.protocols.protocol_base import ChimeraProtBase


class ChimeraProtDiscrepancies(EMProtocol):
    """
    Protocol to find atom discrepancies of all atomic models versus all of the rest.
    """
    _label = 'Find discrepancies'


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('structures', MultiPointerParam, pointerClass="AtomStruct",
                      label='Atomic structure', important=True,
                      help='Select the set of atomic structures to be aligned and analyzed.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertStep)
        self._insertFunctionStep(self.create_chimerax_script)
        self._insertFunctionStep(self.run_chimerax_script_step)
        self._insertFunctionStep(self.analyze_alignment_regions)
        self._insertFunctionStep(self.convertStep)

    # --------------------------- STEPS functions -----------------------------
    def convertStep(self):
        self.extra_files = []
        for i, atomstruct in enumerate(self.structures):
            ori_file = atomstruct.get().getFileName()
            dest_file = os.path.basename(ori_file)
            dest_path = self._getExtraPath(dest_file)
            pwutils.createLink(ori_file, dest_path)
            if not os.path.exists(dest_path):
                raise Exception(f"Failed to create link for {ori_file} to {dest_path}")
            self.extra_files.append(dest_file)
        print(f"Conversion step completed: {self.extra_files}")

    def create_chimerax_script(self):
        project_path = self.getProject().getPath()
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        os.makedirs(output_path, exist_ok=True)
        save_path = os.path.join(project_path, output_path)


        # Generate the ChimeraX command script as a string
        chimerax_script = ""

        # Open all models
        for i, model_file in enumerate(self.extra_files):
            chimerax_script += f"open {model_file}\n"

        # Align every model with every other model and save RMSD values
        rmsd_counter = 1
        for i in range(len(self.extra_files)):
            for j in range(i + 1, len(self.extra_files)):
                model1 = os.path.splitext(os.path.basename(self.extra_files[i]))[0]
                model2 = os.path.splitext(os.path.basename(self.extra_files[j]))[0]
                chimerax_script += f"matchmaker #{i + 1} to #{j + 1} showAlignment true\n"
                chimerax_script += f"sequence header {rmsd_counter} rmsd save {save_path}/rmsd_{model1}_{model2}.txt\n"
                rmsd_counter += 1

        # Save all aligned models
        for i, model_file in enumerate(self.extra_files):
            original_name = os.path.splitext(os.path.basename(model_file))[0]
            chimerax_script += f"save {save_path}/align_{original_name}.cif models #{i + 1}\n"

        # Exit ChimeraX
        chimerax_script += "exit\n"

        # Save the script to a file for debugging purposes
        script_path = os.path.join(output_path, 'chimerax_script.cxc')
        with open(script_path, 'w') as script_file:
            script_file.write(chimerax_script)

        print(f"ChimeraX script created at {script_path}")
        return script_path

    def run_chimerax_script(self, script_path, output_log_path):
        if not os.path.exists(script_path):
            raise Exception(f"ChimeraX script not found at {script_path}")

        # Get the path for the Chimera executable - Version 1.6.1
        self.scipion_path = os.environ.get('SCIPION_HOME', None)
        if self.scipion_path is None:
            raise Exception(
                "SCIPION_HOME environment variable is not set. Please set it to the Scipion installation directory.")
        self.chimerax_executable = os.path.join(self.scipion_path, 'software/em/chimerax-1.6.1/bin/ChimeraX')

        # Run the ChimeraX script and capture the output
        print(f"Chimera: {self.chimerax_executable}")
        result = subprocess.run([self.chimerax_executable, '--nogui', script_path], capture_output=True, text=True)

        # Save the output log to a file
        with open(output_log_path, 'w') as log_file:
            log_file.write(result.stdout)

        if result.returncode != 0:
            with open(output_log_path, 'r') as log_file:
                log_contents = log_file.read()
            raise Exception(f"ChimeraX script failed with return code {result.returncode}. See log for details:\n{log_contents}")

        print(f"ChimeraX script executed successfully. Log saved at {output_log_path}")

    def run_chimerax_script_step(self):
        script_path = os.path.join(self.getWorkingDir(), 'extra', 'chimerax_script.cxc')
        output_log_path = os.path.join(self.getWorkingDir(), 'extra', 'chimerax_output.log')

        self.run_chimerax_script(script_path, output_log_path)
        print(f"ChimeraX script executed. Log saved at {output_log_path}")

        # New methods for alignment region detection

    def detect_alignment_regions(self, rmsd_file, threshold=None):

        """

        Detect well-aligned regions in the RMSD file

        """

        if threshold is None:
            threshold = 5.0

        regions = []

        with open(rmsd_file, 'r') as f:

            lines = f.readlines()[1:]  # Skip header

        current_region = None

        for line in lines:

            try:

                pos, rmsd = line.strip().split(':')

                pos = int(pos)

                rmsd = float(rmsd) if rmsd != 'None' else None

                if rmsd is not None and rmsd < threshold:

                    if current_region is None:

                        current_region = {'start': pos, 'end': pos}

                    else:

                        current_region['end'] = pos

                elif current_region is not None:

                    regions.append(current_region)

                    current_region = None

            except:

                continue

        if current_region is not None:
            regions.append(current_region)

        return regions

    def map_rmsd_to_structure(self, rmsd_file, structure_file):
        """
        Map RMSD indices to actual structural regions
        """
        # 1. Read RMSD file to get alignment indices
        rmsd_regions = self.detect_alignment_regions(rmsd_file)

        # 2. Parse structure file to get chain and residue information
        def parse_structure_file(structure_file):
            chains = {}
            current_chain = None
            residue_count = 0

            with open(structure_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain_id = line[21]
                        residue_num_str = line[22:26].strip()

                        try:
                            residue_num = int(re.search(r'\d+', residue_num_str).group())
                        except (ValueError, AttributeError):
                            print(f"Warning: Skipping invalid residue number '{residue_num_str}' in line: {line}")
                            continue

                        if chain_id not in chains:
                            chains[chain_id] = {
                                'start_residue': residue_num,
                                'residues': []
                            }

                        chains[chain_id]['residues'].append(residue_num)

            # Finalize chain information
            for chain_id, chain_data in chains.items():
                chain_data['end_residue'] = max(chain_data['residues'])

            return chains

        # 3. Map RMSD indices to structural regions
        structure_chains = parse_structure_file(structure_file)

        mapped_regions = []
        for region in rmsd_regions:
            # Determine which chain this region corresponds to
            for chain_id, chain_data in structure_chains.items():
                if (region['start'] >= 1 and
                        region['end'] <= len(chain_data['residues'])):
                    mapped_region = {
                        'chain': chain_id,
                        'start_residue': chain_data['start_residue'] + region['start'] - 1,
                        'end_residue': chain_data['start_residue'] + region['end'] - 1
                    }
                    mapped_regions.append(mapped_region)

        return mapped_regions

    def analyze_alignment_regions(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        rmsd_files = [f for f in os.listdir(output_path) if f.startswith('rmsd_') and f.endswith('.txt')]

        alignment_report = "Alignment Regions Report:\n"
        for rmsd_file in rmsd_files:
            full_rmsd_path = os.path.join(output_path, rmsd_file)

            # Extract model names to find corresponding structure file
            model_names = re.findall(r'rmsd_(.+)\.txt', rmsd_file)[0]
            structure_file = os.path.join(output_path, model_names.split('_')[1] + '.cif')
            print(f"Structure: {structure_file}")

            # Map RMSD regions to actual structure
            mapped_regions = self.map_rmsd_to_structure(full_rmsd_path, structure_file)

            alignment_report += f"\nAlignment for {model_names}:\n"
            for region in mapped_regions:
                alignment_report += (
                    f"  Aligned region in chain {region['chain']}: "
                    f"Residues {region['start_residue']} to {region['end_residue']}\n"
                )

        report_path = os.path.join(output_path, 'alignment_regions_report.txt')

        with open(report_path, 'w') as f:

            f.write(alignment_report)

        print("Alignment regions analysis complete.")