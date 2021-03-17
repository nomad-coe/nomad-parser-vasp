#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import pint
import numpy as np

from nomad.datamodel import EntryArchive
from vaspparser.vasp_parser import VASPParser


@pytest.fixture(scope='module')
def parser():
    return VASPParser()


def test_vasprunxml_static(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.static', archive, None)

    assert len(archive.section_run) == 1

    sec_run = archive.section_run[0]

    assert sec_run.program_compilation_datetime == pint.Quantity(1366564273.0, 's')

    sec_method = sec_run.section_method[0]
    assert len(sec_method.x_vasp_incar_in) == 27
    assert len(sec_method.x_vasp_incar_out) == 112
    assert sec_method.x_vasp_incar_in['LCHARG']
    assert len(sec_method.k_mesh_points) == 888
    assert pytest.approx(sec_method.k_mesh_points[-1][2], 0.48437500)
    assert pytest.approx(sec_method.k_mesh_weights[100], 0.00073242)
    assert pytest.approx(sec_method.x_vasp_tetrahedrons_list[0], [4, 1, 2, 2, 18])
    assert pytest.approx(sec_method.section_method_atom_kind[0].method_atom_kind_mass.magnitude, 24.305)
    assert len(sec_method.section_XC_functionals) == 2

    sec_system = sec_run.section_system[-1]
    assert len(sec_system.atom_labels) == 1
    assert pytest.approx(sec_system.lattice_vectors[0][0].magnitude, -1.78559323e-10)

    sec_scc = sec_run.section_single_configuration_calculation[-1]
    assert pytest.approx(sec_scc.energy_total.magnitude, 2.3264377e-19)
    assert np.shape(sec_scc.atom_forces) == (1, 3)
    assert pytest.approx(sec_scc.stress_tensor[2][2].magnitude, -2.78384438e+08)
    assert len(sec_scc.section_dos[0].dos_energies) == 5000
    assert len(sec_scc.section_dos[0].dos_lm) == 9
    assert pytest.approx(sec_scc.section_dos[0].dos_values_lm[0][0][0][-1], 3.40162245e+17)
    assert np.shape(sec_scc.section_eigenvalues[0].eigenvalues_values) == (1, 888, 37)
    assert pytest.approx(sec_scc.section_scf_iteration[2].energy_total_T0_scf_iteration.magnitude, -2.27580485e-19,)


def test_vasprunxml_relax(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.relax', archive, None)

    assert len(archive.section_run[0].section_method) == 1

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 3
    assert pytest.approx(sec_systems[1].atom_positions[1][0].magnitude, 3.2771907e-10)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 3
    assert [len(scc.section_scf_iteration) for scc in sec_sccs] == [12, 10, 6]
    assert pytest.approx(sec_sccs[0].energy_free.magnitude, -1.14352735e-18)
    assert np.mean(sec_sccs[1].atom_forces.magnitude) == 0.0
    assert pytest.approx(sec_sccs[2].stress_tensor[2][2].magnitude, -2.02429105e+08)
    assert pytest.approx(sec_sccs[2].energy_reference_lowest_unoccupied[0].magnitude, 7.93718304e-19)
    assert pytest.approx(sec_sccs[2].energy_reference_highest_occupied[1].magnitude, 7.93702283e-19)
    assert [len(scc.section_eigenvalues) for scc in sec_sccs] == [0, 0, 1]
    assert [len(scc.section_dos) for scc in sec_sccs] == [0, 0, 1]


def test_vasprunxml_bands(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.bands', archive, None)

    sec_k_band = archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    assert len(sec_k_band.section_k_band_segment) == 6
    assert np.shape(sec_k_band.section_k_band_segment[0].band_energies.magnitude) == (1, 128, 37)
    assert pytest.approx(sec_k_band.section_k_band_segment[1].band_energies[0][1][1].magnitude, -6.27128785e-18)
    assert sec_k_band.section_k_band_segment[5].band_occupations[0][127][5] == 0.0


def test_outcar(parser):
    archive = EntryArchive()
    parser.parse('tests/data/OUTCAR', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_version == '5.3.2 13Sep12 complex serial LinuxIFC'
    assert pytest.approx(sec_run.section_basis_set_cell_dependent[0].basis_set_planewave_cutoff.magnitude, 8.3313185e-17)

    sec_method = sec_run.section_method[0]
    assert sec_method.section_XC_functionals[0].XC_functional_name == 'GGA_X_PBE'
    assert len(sec_method.section_method_atom_kind) == 2
    assert sec_method.section_method_atom_kind[1].method_atom_kind_pseudopotential_name == 'Ag'

    sec_system = sec_run.section_system[0]
    assert pytest.approx(sec_system.lattice_vectors[1][0].magnitude, 3.141538e-10)
    assert pytest.approx(sec_system.atom_positions[1][0].magnitude, 3.14154e-10)

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert pytest.approx(sec_scc.energy_total.magnitude, -2.23390885e-18)
    assert sec_scc.atom_forces[0][0].magnitude == 0.0
    assert pytest.approx(sec_scc.stress_tensor[0][0].magnitude, 7.060258e+09)
    assert pytest.approx(sec_scc.energy_reference_lowest_unoccupied[0].magnitude, 9.40461662e-19)
    assert pytest.approx(sec_scc.energy_reference_highest_occupied[0].magnitude, 9.51212268e-19)
    sec_eigs = sec_scc.section_eigenvalues[0]
    assert np.shape(sec_eigs.eigenvalues_values) == (1, 145, 15)
    assert pytest.approx(sec_eigs.eigenvalues_values[0][9][14].magnitude, 1.41810256e-18)
    assert sec_eigs.eigenvalues_occupation[0][49][9] == 2.0
    sec_dos = sec_scc.section_dos[0]
    assert len(sec_dos.dos_energies) == 301
    assert sec_dos.dos_integrated_values[0][-1] == 30.0
    assert pytest.approx(sec_dos.dos_energies[2].magnitude, -5.05662967e-18)
    assert len(sec_dos.dos_lm) == 16
    assert pytest.approx(sec_dos.dos_values_lm[5][0][0][-15], 1.51481425e+17)
    assert pytest.approx(sec_dos.dos_values_lm[2][0][1][-16], 1.71267009e+16)
