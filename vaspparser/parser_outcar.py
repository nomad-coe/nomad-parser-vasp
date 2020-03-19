# Copyright 2019-2018 Markus Scheidgen
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import sys
import numpy as np

from ase.io import read as ase_read

from nomadcore.simple_parser import SimpleMatcher
from nomadcore.baseclasses import ParserInterface, MainHierarchicalParser

from nomad.parsing import Backend


"""
A very basic VASP OUTCAR parser. It is only used to get nomad repository
metadata in absence of the VASP .xml file.
"""


class VaspOutcarParser(ParserInterface):

    def get_metainfo_filename(self):
        return 'vasp.nomadmetainfo.json'

    def get_parser_info(self):
        return {
            'name': 'vaspoutcar_parser',
            'version': '1.0.0'
        }

    def setup_version(self):
        self.setup_main_parser(None)

    def setup_main_parser(self, _):
        self.main_parser = MainParser(self.parser_context)


class MainParser(MainHierarchicalParser):
    def __init__(self, parser_context, *args, **kwargs):
        super().__init__(parser_context, *args, **kwargs)
        self.lattice_vectors = []

        self.root_matcher = SimpleMatcher(
            name='root',
            startReStr=r'vasp.',
            weak=True,
            sections=['section_run', 'section_system', 'section_method'],
            subMatchers=[
                SimpleMatcher(r'\svasp.(?P<program_version>\d+.\d+.\d+)\s.*')
            ]
        )

    def parse(self, mainfile, *args, **kwargs):
        self.ase = ase_read(mainfile, format='vasp-out')
        super().parse(mainfile, *args, **kwargs)

    def add_lattice_vector(self, backend, data):
        self.lattice_vectors.append(list(float(x.strip()) for x in data[1].split(',')))

    def onClose_section_method(self, backend, *args, **kwargs):
        backend.openNonOverlappingSection('section_XC_functionals')
        backend.addValue('XC_functional_name', 'GGA_X_PBE')
        backend.closeNonOverlappingSection('section_XC_functionals')
        backend.addValue('electronic_structure_method', 'DFT')

    def onClose_section_system(self, backend, *args, **kwargs):
        backend.addArrayValues('atom_labels', np.array(self.ase.get_chemical_symbols()))
        backend.addArrayValues('atom_positions', self.ase.get_positions() * 1e-10)
        backend.addArrayValues('lattice_vectors', self.ase.get_cell() * 1e-10)
        backend.addArrayValues('configuration_periodic_dimensions', self.ase.get_pbc())

    def onClose_section_run(self, backend, *args, **kwargs):
        backend.addValue('program_name', 'VASP')
        backend.addValue('program_basis_set_type', 'plane waves')


if __name__ == "__main__":
    parser = VaspOutcarParser(backend=Backend)
    parser.parse(sys.argv[1])
    print(parser.parser_context.super_backend)
