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
import numpy as np

from nomad.units import ureg
from nomad.datamodel import EntryArchive
from vaspparser.vasp_parser import VASPParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return VASPParser()


@pytest.fixture(scope='module')
def silicon_dos(parser):
    archive = EntryArchive()
    parser.parse('tests/data/dos_si.xml', archive, None)
    return archive


@pytest.fixture(scope='module')
def silicon_band(parser):
    archive = EntryArchive()
    parser.parse('tests/data/band_si.xml', archive, None)
    return archive


def test_vasprunxml_static(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.static', archive, None)

    assert len(archive.section_run) == 1

    sec_run = archive.section_run[0]

    assert sec_run.program_compilation_datetime.magnitude == 1366564273.0

    sec_method = sec_run.section_method[0]
    assert len(sec_method.x_vasp_incar_in) == 27
    assert len(sec_method.x_vasp_incar_out) == 112
    assert sec_method.x_vasp_incar_in['LCHARG']
    assert len(sec_method.k_mesh_points) == 888
    assert sec_method.k_mesh_points[-1][2] == approx(0.48437500)
    assert sec_method.k_mesh_weights[100] == approx(0.00073242)
    assert (sec_method.x_vasp_tetrahedrons_list[0] == [4, 1, 2, 2, 18]).all()
    assert sec_method.section_method_atom_kind[0].method_atom_kind_mass.magnitude == approx(24.305)
    assert len(sec_method.section_XC_functionals) == 2

    sec_system = sec_run.section_system[-1]
    assert len(sec_system.atom_labels) == 1
    assert sec_system.lattice_vectors[0][0].magnitude == approx(-1.78559323e-10)

    sec_scc = sec_run.section_single_configuration_calculation[-1]
    assert sec_scc.energy_total.value.magnitude == approx(-2.3264377e-19)
    assert np.shape(sec_scc.forces_total.value) == (1, 3)
    assert sec_scc.stress_total.value[2][2].magnitude == approx(-2.78384438e+08)
    assert len(sec_scc.dos_electronic[0].energies) == 5000
    assert sec_scc.dos_electronic[0].total[0].value[1838].magnitude == approx(1.94581094e+19)
    assert len(sec_scc.dos_electronic[0].atom_projected) == 9
    assert sec_scc.dos_electronic[0].atom_projected[0].value[-1].magnitude == approx(3.40162245e+17)
    assert np.shape(sec_scc.eigenvalues[0].band_energies[887].value) == (37,)
    assert sec_scc.section_scf_iteration[2].energy_total_T0_scf_iteration.magnitude == approx(-2.27580485e-19,)


def test_vasprunxml_relax(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.relax', archive, None)

    assert len(archive.section_run[0].section_method) == 1

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 3
    assert sec_systems[1].atom_positions[1][0].magnitude == approx(3.2771907e-10)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 3
    assert [len(scc.section_scf_iteration) for scc in sec_sccs] == [12, 10, 6]
    assert sec_sccs[0].energy_free.value.magnitude == approx(-1.14352735e-18)
    assert np.mean(sec_sccs[1].forces_total.value.magnitude) == 0.0
    assert sec_sccs[2].stress_total.value[2][2].magnitude == approx(-2.02429105e+08)
    assert sec_sccs[2].energy_reference_lowest_unoccupied[0].magnitude == approx(7.93718304e-19)
    assert sec_sccs[2].energy_reference_highest_occupied[1].magnitude == approx(7.93702283e-19)
    assert [len(scc.eigenvalues) for scc in sec_sccs] == [0, 0, 1]
    assert [len(scc.dos_electronic) for scc in sec_sccs] == [0, 0, 1]


def test_vasprunxml_bands(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.bands', archive, None)

    sec_k_band = archive.section_run[0].section_single_configuration_calculation[0].band_structure_electronic[0]
    assert len(sec_k_band.band_structure_segment) == 6
    assert np.shape(sec_k_band.band_structure_segment[0].band_energies[127].value.magnitude) == (37,)
    assert sec_k_band.band_structure_segment[1].band_energies[1].value[1].magnitude == approx(-6.27128785e-18)
    assert sec_k_band.band_structure_segment[5].band_energies[127].occupations[5] == 0.0


def test_band_silicon(silicon_band):
    """Tests that the band structure of silicon is parsed correctly.
    """
    scc = silicon_band.section_run[-1].section_single_configuration_calculation[0]
    band = scc.band_structure_electronic[-1]
    segments = band.band_structure_segment
    energies = []
    for segment in segments:
        energies.append([b.value.to(ureg.electron_volt).magnitude for b in segment.band_energies])
    energies = np.array(energies)

    # Check that an energy reference is reported
    energy_reference = scc.energy_reference_fermi
    if energy_reference is None:
        energy_reference = scc.energy_reference_highest_occupied
    assert energy_reference is not None
    energy_reference = energy_reference.to(ureg.electron_volt).magnitude

    # Check that an approporiately sized band gap is found at the given
    # reference energy
    energies = energies.flatten()
    energies.sort()
    lowest_unoccupied_index = np.searchsorted(energies, energy_reference, "right")[0]
    highest_occupied_index = lowest_unoccupied_index - 1
    gap = energies[lowest_unoccupied_index] - energies[highest_occupied_index]
    assert gap == approx(0.5091)


def test_dos_silicon(silicon_dos):
    """Tests that the DOS of silicon is parsed correctly.
    """
    scc = silicon_dos.section_run[-1].section_single_configuration_calculation[0]
    dos = scc.dos_electronic[-1]
    energies = dos.energies.to(ureg.electron_volt).magnitude
    values = np.array([d.value for d in dos.total])

    # Check that an energy reference is reported
    energy_reference = scc.energy_reference_fermi
    if energy_reference is None:
        energy_reference = scc.energy_reference_highest_occupied
    assert energy_reference is not None
    energy_reference = energy_reference.to(ureg.electron_volt).magnitude

    # Check that an appropriately sized band gap is found at the given
    # reference energy
    nonzero = np.unique(values.nonzero())
    energies = energies[nonzero]
    energies.sort()
    lowest_unoccupied_index = np.searchsorted(energies, energy_reference, "right")[0]
    highest_occupied_index = lowest_unoccupied_index - 1
    gap = energies[lowest_unoccupied_index] - energies[highest_occupied_index]
    assert gap == approx(0.83140)


def test_outcar(parser):
    archive = EntryArchive()
    parser.parse('tests/data/OUTCAR', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_version == '5.3.2 13Sep12 complex serial LinuxIFC'
    assert sec_run.section_basis_set_cell_dependent[0].basis_set_planewave_cutoff.magnitude == approx(8.3313185e-17)

    sec_method = sec_run.section_method[0]
    assert sec_method.section_XC_functionals[0].XC_functional_name == 'GGA_X_PBE'
    assert len(sec_method.section_method_atom_kind) == 2
    assert sec_method.section_method_atom_kind[1].method_atom_kind_pseudopotential_name == 'Ag'

    sec_system = sec_run.section_system[0]
    assert sec_system.lattice_vectors[1][0].magnitude == approx(3.141538e-10)
    assert sec_system.atom_positions[1][0].magnitude == approx(3.14154e-10)

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert sec_scc.energy_total.value.magnitude == approx(-2.23390885e-18)
    assert sec_scc.forces_total.value[0][0].magnitude == 0.0
    assert sec_scc.stress_total.value[0][0].magnitude == approx(7.060258e+09)
    assert sec_scc.energy_reference_lowest_unoccupied[0].magnitude == approx(9.40461662e-19)
    assert sec_scc.energy_reference_highest_occupied[0].magnitude == approx(9.51212268e-19)
    sec_scfs = sec_scc.section_scf_iteration
    assert len(sec_scfs) == 11
    assert sec_scfs[4].energy_total_scf_iteration.magnitude == approx(-1.20437432e-18)
    assert sec_scfs[7].energy_sum_eigenvalues_scf_iteration.magnitude == approx(-1.61008452e-17)
    sec_eigs = sec_scc.eigenvalues[0]
    assert np.shape(sec_eigs.band_energies[144].value) == (15,)
    assert sec_eigs.band_energies[9].value[14].magnitude == approx(1.41810256e-18)
    assert sec_eigs.band_energies[49].occupations[9] == 2.0
    sec_dos = sec_scc.dos_electronic[0]
    assert len(sec_dos.energies) == 301
    assert sec_dos.total[0].value[282].magnitude == approx(9.84995713e+20)
    assert sec_dos.total[0].value_integrated[-1] == 30.0
    assert sec_dos.atom_projected[10].value[-15].magnitude == approx(1.51481425e+17)
    assert sec_dos.atom_projected[5].value[-16].magnitude == approx(1.71267009e+16)


def test_broken_xml(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasprun.xml.broken', archive, None)
