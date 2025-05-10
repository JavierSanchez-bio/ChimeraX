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
import subprocess
import re
import shutil
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
        self._insertFunctionStep(self.create_folders)
        self._insertFunctionStep(self.add_rmsd)
        self._insertFunctionStep(self.final_models)

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
                chimerax_script += "setattr a occupancy 1111.11\n"  # Arbitrary number for later changes
                chimerax_script += f"sequence header {rmsd_counter} rmsd save {save_path}/rmsd_{model1}_{model2}.txt\n"
                chimerax_script += f"save {save_path}/fasta_{model1}_{model2}.fasta format fasta alignment {rmsd_counter}\n"
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
            raise Exception(
                f"ChimeraX script failed with return code {result.returncode}. See log for details:\n{log_contents}")

        print(f"ChimeraX script executed successfully. Log saved at {output_log_path}")

    def run_chimerax_script_step(self):
        script_path = os.path.join(self.getWorkingDir(), 'extra', 'chimerax_script.cxc')
        output_log_path = os.path.join(self.getWorkingDir(), 'extra', 'chimerax_output.log')

        self.run_chimerax_script(script_path, output_log_path)
        print(f"ChimeraX script executed. Log saved at {output_log_path}")

    def create_folders(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        for fasta_file in os.listdir(output_path):
            if fasta_file.startswith('fasta_') and fasta_file.endswith('.fasta'):
                model1, model2 = fasta_file[6:-6].split('_')
                folder_name = f"{model1}_{model2}"
                folder_path = os.path.join(output_path, folder_name)
                os.makedirs(folder_path, exist_ok=True)

                # Link the FASTA / RMSD / CIF files
                src_fasta_file = os.path.join(output_path, fasta_file)
                dest_fasta_file = os.path.join(folder_path, fasta_file)
                pwutils.createLink(src_fasta_file, dest_fasta_file)

                rmsd_file = f"rmsd_{model1}_{model2}.txt"
                src_rmsd_file = os.path.join(output_path, rmsd_file)
                dest_rmsd_file = os.path.join(folder_path, rmsd_file)
                pwutils.createLink(src_rmsd_file, dest_rmsd_file)

                cif_file = f"align_{model1}.cif"
                src_cif_file = os.path.join(output_path, cif_file)
                dest_cif_file = os.path.join(folder_path, cif_file)
                pwutils.createLink(src_cif_file, dest_cif_file)
                out_model1_file = os.path.join(folder_path, f"out_{model1}.cif")
                pwutils.copyFile(src_cif_file, out_model1_file)
                cif_file = f"align_{model2}.cif"
                src_cif_file = os.path.join(output_path, cif_file)
                dest_cif_file = os.path.join(folder_path, cif_file)
                pwutils.createLink(src_cif_file, dest_cif_file)
                out_model2_file = os.path.join(folder_path, f"out_{model2}.cif")
                pwutils.copyFile(src_cif_file, out_model2_file)

        for folder in os.listdir(output_path):
            folder_path = os.path.join(output_path, folder)
            if os.path.isdir(folder_path) and "_" in folder:
                model1, model2 = folder.split('_')
                out_model1_file = os.path.join(folder_path, f"out_{model1}.cif")
                out_model2_file = os.path.join(folder_path, f"out_{model2}.cif")

                text_to_append = "\nloop_\n_scipion_attributes.name\n_scipion_attributes.recipient\n_scipion_attributes.specifier\n_scipion_attributes.value\n"

                with open(out_model1_file, 'a') as file1:
                    file1.write(text_to_append)

                with open(out_model2_file, 'a') as file2:
                    file2.write(text_to_append)

    def add_rmsd(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        self.aa_rmsd = []  # Top-level array that stores every aminoacid rmsd
        self.occ_position = [[], []]  # Stores the position of occupancy column in each model
        for folder in os.listdir(output_path):
            folder_path = os.path.join(output_path, folder)
            if os.path.isdir(folder_path) and "_" in folder:
                model1, model2 = folder.split('_')
                fasta_file = os.path.join(folder_path, f"fasta_{model1}_{model2}.fasta")

                mat1, mat2 = [], []
                with open(fasta_file, 'r') as f:
                    lines = f.readlines()
                    reading_model1, reading_model2 = False, False
                    for line in lines:
                        if line.startswith(f">{model1}"):
                            reading_model1 = True
                            reading_model2 = False
                            continue
                        elif line.startswith(f">{model2}"):
                            reading_model1 = False
                            reading_model2 = True
                            continue
                        elif line.startswith(">"):
                            reading_model1, reading_model2 = False, False

                        if reading_model1:
                            mat1.extend([0 if char == '.' else 1 for char in line.strip()])
                        elif reading_model2:
                            mat2.extend([0 if char == '.' else 1 for char in line.strip()])

                rmsd_file = os.path.join(folder_path, f"rmsd_{model1}_{model2}.txt")

                if not os.path.exists(rmsd_file):
                    raise Exception(f"RMSD file not found: {rmsd_file}")

                mat_rmsd = []

                with open(rmsd_file, 'r') as f:

                    lines = f.readlines()

                    for line in lines[1:]:  # Skip the first line

                        parts = line.split(":")

                        if len(parts) != 2:
                            continue

                        try:

                            Aa = int(parts[0].strip())

                            rmsd_value = float(parts[1].strip()) if parts[1].strip() != "None" else 250.00

                            while len(mat_rmsd) <= Aa:
                                mat_rmsd.append(0)

                            mat_rmsd[Aa] += rmsd_value

                        except ValueError:

                            continue

                # mat_pos_ori stores the number of the amino acids in the model and mat_pos_align the position in the alignment
                mat_pos1_ori = []
                mat_pos1_align = []

                pos = 0
                for i, value in enumerate(mat1):
                    if value != 0:
                        mat_pos1_ori.append(pos)
                        mat_pos1_align.append(i)
                        pos += 1

                mat_pos2_ori = []
                mat_pos2_align = []
                pos = 0
                for i, value in enumerate(mat2):
                    if value != 0:
                        mat_pos2_ori.append(pos)
                        mat_pos2_align.append(i)
                        pos += 1

                # Write RMSD residues at the end of each file
                out_model1_file = os.path.join(folder_path, f"out_{model1}.cif")
                out_model2_file = os.path.join(folder_path, f"out_{model2}.cif")
                occupancy1 = []
                occupancy2 = []

                aa_occupancy1 = [[], []]  # First line: amino acid number / Second line: rmsd value
                with open(out_model1_file, 'a') as file1:
                    for i in range(len(mat_pos1_ori)):
                        rmsd_value = mat_rmsd[mat_pos1_align[i] + 1]
                        file1.write(f"RMSD residues A:{i + 1}   {rmsd_value}\n")
                        occupancy1.append(rmsd_value)
                        aa_occupancy1[0].append(i + 1)
                        aa_occupancy1[1].append(rmsd_value)
                self.aa_rmsd.append({model1: aa_occupancy1})

                aa_occupancy2 = [[], []]  # First line: amino acid number / Second line: rmsd value
                with open(out_model2_file, 'a') as file2:
                    for i in range(len(mat_pos2_ori)):
                        rmsd_value = mat_rmsd[mat_pos2_align[i] + 1]
                        file2.write(f"RMSD residues A:{i + 1}   {rmsd_value}\n")
                        occupancy2.append(rmsd_value)
                        aa_occupancy2[0].append(i + 1)
                        aa_occupancy2[1].append(rmsd_value)
                self.aa_rmsd.append({model2: aa_occupancy2})

                # Add occupancy to each model
                with open(out_model1_file, 'r') as file1, open(out_model2_file, 'r') as file2:
                    lines1 = file1.readlines()
                    lines2 = file2.readlines()

                subs = []
                with open(out_model1_file, 'w') as mod_file1, open(out_model2_file, 'w') as mod_file2:
                    for line in lines1:
                        if line.startswith('ATOM'):
                            atom_line = line
                            clean_line = re.sub(r'\s+', ' ', atom_line).strip()
                            columns = clean_line.split()

                            if '1111.11' in columns:
                                strip_position = columns.index('1111.11')
                                subs.append(strip_position)

                            index = int(columns[8])
                            new_value = round(occupancy1[index - 1], 2)
                            modified_line = re.sub(r'1111\.11', str(new_value), atom_line)
                            mod_file1.write(modified_line)
                        else:
                            mod_file1.write(line)

                    for line in lines2:
                        if line.startswith('ATOM'):
                            atom_line = line
                            clean_line = re.sub(r'\s+', ' ', atom_line).strip()
                            columns = clean_line.split()

                            if '1111.11' in columns:
                                strip_position = columns.index('1111.11')
                                subs.append(strip_position)

                            index = int(columns[8])
                            new_value = round(occupancy2[index - 1], 2)
                            modified_line = re.sub(r'1111\.11', str(new_value), atom_line)
                            mod_file2.write(modified_line)
                        else:
                            mod_file2.write(line)

                    subs = list(set(subs))  # Remove duplicates
                    self.occ_position[0].append(model1)
                    self.occ_position[0].append(model2)
                    if len(subs) > 1:
                        self.occ_position[1].append(subs[0])
                        self.occ_position[1].append(subs[1])
                    else:
                        self.occ_position[1].append(subs[0])
                        self.occ_position[1].append(subs[0])

                    unique = {}
                    for k, v in zip(*self.occ_position):
                        if k not in unique:
                            unique[k] = v

                    self.occ_position = [list(unique.keys()), list(unique.values())]

    def final_models(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        final_output_path = os.path.join(output_path, 'FINAL-OUTPUTS')
        print(f"Position of occupancy:", self.occ_position)
        os.makedirs(final_output_path, exist_ok=True)

        cif_counter = {}

        for folder in os.listdir(output_path):
            folder_path = os.path.join(output_path, folder)
            if os.path.isdir(folder_path) and folder != 'FINAL-OUTPUTS':
                for file_name in os.listdir(folder_path):
                    if file_name.startswith('out_') and file_name.endswith('.cif'):
                        src_file = os.path.join(folder_path, file_name)
                        base_name = os.path.basename(src_file)
                        if base_name in cif_counter:
                            cif_counter[base_name] += 1
                            new_file_name = f"{os.path.splitext(base_name)[0]}_{cif_counter[base_name]}.cif"
                        else:
                            cif_counter[base_name] = 1
                            new_file_name = base_name
                        dest_file = os.path.join(final_output_path, new_file_name)
                        shutil.copy(src_file, dest_file)

                rmsd_dict = {}
                for i, model in enumerate(self.occ_position[0]):
                    occ = self.occ_position[1][i]
                    rmsd_dict[model] = []

                    model_files = [
                        file_name for file_name in os.listdir(final_output_path)
                        if file_name.startswith(f"out_{model}")
                    ]

                    for model_file in model_files:
                        rmsd_matrix = []
                        model_file_path = os.path.join(final_output_path, model_file)

                        # Open the file and process lines starting with "ATOM"
                        with open(model_file_path, 'r') as file:
                            for line in file:
                                if line.startswith('ATOM'):
                                    clean_line = re.sub(r'\s+', ' ', line).strip()
                                    columns = clean_line.split()
                                    if len(columns) > occ:
                                        rmsd_value = float(columns[occ])
                                        rmsd_matrix.append(rmsd_value)

                        rmsd_dict[model].append(rmsd_matrix)

            averaged_rmsd_dict = {}
            for model, matrices in rmsd_dict.items():
                if matrices:
                    averaged_matrix = [
                        sum(values) / len(values) for values in zip(*matrices)
                    ]

                    averaged_rmsd_dict[model] = averaged_matrix

            for i, (model, averaged_matrix) in enumerate(averaged_rmsd_dict.items()):
                occ = self.occ_position[1][i]  # Get the occupancy position for the current model

                model_files = [
                    file_name for file_name in os.listdir(final_output_path)
                    if file_name.startswith(f"out_{model}_")
                ]
                for model_file in model_files:
                    file_path = os.path.join(final_output_path, model_file)
                    os.remove(file_path)

                remaining_file = os.path.join(final_output_path, f"out_{model}.cif")

                updated_lines = []
                counter = 0
                with open(remaining_file, 'r') as file:
                    for line in file:
                        if line.startswith('ATOM'):
                            clean_line = re.sub(r'\s+', ' ', line).strip()
                            columns = clean_line.split()
                            if len(columns) > occ:
                                columns[occ] = f"{averaged_matrix[counter]:.2f}"
                                counter += 1
                            updated_line = ' '.join(columns) + '\n'
                            updated_lines.append(updated_line)
                        else:
                            updated_lines.append(line)

                with open(remaining_file, 'w') as file:
                    file.writelines(updated_lines)

        # Modify amino acids RMSD:
        final_aa_rmsd = []
        grouped_rmsd = {}

        for entry in self.aa_rmsd:
            for model_name, matrix in entry.items():
                if model_name not in grouped_rmsd:
                    grouped_rmsd[model_name] = []
                grouped_rmsd[model_name].append(matrix)

        # Compute the mean matrix for each model and store it in final_aa_rmsd
        for model_name, matrices in grouped_rmsd.items():
            mean_matrix = [[], []]

            num_matrices = len(matrices)
            for line_index in range(2):
                summed_values = [0] * len(matrices[0][line_index])
                for matrix in matrices:
                    for i, value in enumerate(matrix[line_index]):
                        summed_values[i] += value
                mean_matrix[line_index] = [val / num_matrices for val in summed_values]

            final_aa_rmsd.append({model_name: mean_matrix})

        # Add RMSD to the end of the final models:
        for file_name in os.listdir(final_output_path):
            if file_name.startswith("out_") and file_name.endswith(".cif"):
                model = file_name[4:-4]  # Extract the model name from the file name (e.g., "out_model.cif" -> "model")
                file_path = os.path.join(final_output_path, file_name)

                # Find the corresponding matrix in final_aa_rmsd
                matching_entry = next((entry for entry in final_aa_rmsd if model in entry), None)
                if not matching_entry:
                    raise Exception(f"No matching RMSD matrix found for model {model} in final_aa_rmsd")

                rmsd_matrix = matching_entry[model][1]  # Second line of the matrix contains the RMSD values

                # Read the file and update the RMSD values
                updated_lines = []
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                    modify_rmsd = False  # Flag to start modifying RMSD lines

                    for line in lines:
                        if line.startswith("_scipion_attributes.value"):
                            modify_rmsd = True
                            updated_lines.append(line)
                        elif modify_rmsd and line.startswith("RMSD residues A:"):
                            parts = line.split()
                            position = int(parts[2].strip("A:"))
                            rmsd_value = round(rmsd_matrix[position - 1], 6)
                            updated_line = f"RMSD residues A:{position}   {rmsd_value}\n"
                            updated_lines.append(updated_line)
                        else:
                            updated_lines.append(line)

                # Write the updated lines back to the file
                with open(file_path, 'w') as file:
                    file.writelines(updated_lines)


        # Define outputs directly within final_models
        for file_name in os.listdir(final_output_path):
            file_path = os.path.join(final_output_path, file_name)
            if os.path.isfile(file_path):
                # Create an AtomStruct object for each file
                output = AtomStruct(filename=file_path)

                # Use the file name (without extension) as the output name
                output_name = os.path.splitext(file_name)[0]

                # Define the output dynamically
                self._defineOutputs(**{output_name: output})

                # Relate the output to the input structures
                self._defineSourceRelation(self.structures, output)

