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

    # --------------------------- STEPS functions -----------------------------
    def convertStep(self):
        self.extra_files = []
        for i, atomstruct in self.structures:
            ori_file = atomstruct.get().getFileName()
            dest_file = str(i) + '.cif'  # maintain .cif extension
            pwutils.createLink(ori_file, self._getExtraPath(dest_file))
            self.extra_files.append(dest_file)
        print(f"Conversion step completed: {self.extra_files}")
