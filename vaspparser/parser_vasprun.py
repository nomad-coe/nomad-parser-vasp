# Copyright 2016-2018 Fawzi Mohamed, Lauri Himanen, Danio Brambila, Ankit Kariryaa, Henning Glawe
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

from __future__ import division
from builtins import range
from builtins import object
import xml.etree.ElementTree
from xml.etree.ElementTree import ParseError
import bisect
from datetime import datetime
import re
import traceback
import numpy as np
import ase.geometry
import ase.data
import xml.etree.ElementTree as ET
import logging

from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.unit_conversion.unit_conversion import convert_unit

from nomad.metainfo import Quantity
from nomad.datamodel.metainfo.public import section_run as msection_run
from nomad.datamodel.metainfo.public import section_method as msection_method
from nomad.datamodel.metainfo.public import section_sampling_method as msection_sampling_method
from nomad.datamodel.metainfo.public import section_frame_sequence as msection_frame_sequence
from nomad.datamodel.metainfo.public import section_single_configuration_calculation as msection_single_configuration_calculation
from nomad.datamodel.metainfo.public import section_basis_set as msection_basis_set
from nomad.datamodel.metainfo.public import section_XC_functionals as msection_XC_functionals
from nomad.datamodel.metainfo.public import section_system as msection_system
from nomad.datamodel.metainfo.public import section_k_band as msection_k_band
from nomad.datamodel.metainfo.public import section_k_band_segment as msection_k_band_segment
from nomad.datamodel.metainfo.public import section_k_band_normalized as msection_k_band_normalized
from nomad.datamodel.metainfo.public import section_k_band_segment_normalized as msection_k_band_segment_normalized
from nomad.datamodel.metainfo.public import section_eigenvalues as msection_eigenvalues
from nomad.datamodel.metainfo.public import section_method_atom_kind as msection_method_atom_kind
from nomad.datamodel.metainfo.public import section_basis_set_cell_dependent as msection_basis_set_cell_dependent
from nomad.datamodel.metainfo.public import section_dos as msection_dos
from nomad.datamodel.metainfo.common import section_method_basis_set as msection_method_basis_set

eV2J = convert_unit_function("eV", "J")
eV2JV = np.vectorize(eV2J)

vasp_to_metainfo_type_mapping = {
    'string': ['C'],
    'int': ['i'],
    'logical': ['b', 'C'],
    'float': ['f']}

special_paths = {
    'cubic': 'ΓXMΓRX,MR',
    'fcc': 'ΓXWKΓLUWLK,UX',
    'bcc': 'ΓHNΓPH,PN',
    'tetragonal': 'ΓXMΓZRAZXR,MA',
    'orthorhombic': 'ΓXSYΓZURTZ,YT,UX,SR',
    'hexagonal': 'ΓMKΓALHA,LM,KH',
    'monoclinic': 'ΓYHCEM1AXH1,MDZ,YD'}


def seconds_from_epoch(date):
    epoch = datetime(1970, 1, 1)
    ts = date - epoch
    return ts.seconds + ts.microseconds / 1000.0


trueRe = re.compile(
    r"\s*(?:\.?[Tt](?:[Rr][Uu][Ee])?\.?|1|[Yy](?:[Ee][Ss])?|[Jj][Aa]?)\s*$")
falseRe = re.compile(
    r"\s*(?:\.?[fF](?:[Aa][Ll][Ss][Ee])?\.?|0|[Nn](?:[Oo]|[Ee][Ii][Nn])?)\s*$")


def to_bool(value):
    if falseRe.match(value):
        return False
    elif trueRe.match(value):
        return True
    else:
        return None


meta_type_transformers = {
    'C': lambda x: x.strip(),
    'i': lambda x: int(float(x.strip())),
    'f': lambda x: float(x.strip()),
    'b': to_bool,
}


class MyXMLParser(ET.XMLParser):

    rx = re.compile("&#([0-9]+);|&#x([0-9a-fA-F]+);")

    def feed(self, data):
        m = self.rx.search(data)
        if m is not None:
            target = m.group(1)
            if target:
                num = int(target)
            else:
                num = int(m.group(2), 16)
            if not(
                num in (0x9, 0xA, 0xD) or 0x20 <= num <= 0xD7FF or 0xE000 <= num <= 0xFFFD
                or 0x10000 <= num <= 0x10FFFF):
                # is invalid xml character, cut it out of the stream
                mstart, mend = m.span()
                mydata = data[:mstart] + data[mend:]
        else:
            mydata = data
        super(MyXMLParser, self).feed(mydata)


def transform2(y):
    if '**' in y:
        return float('nan')
    else:
        return y


def get_vector(el, transform=float, field="v"):
    """ returns the vasp style vector contained in the element el (using field v).
    single elements are converted using the function convert"""
#
#    for x in el.findall(field):
#        for y in re.split(r"\s+", x.text.strip()):
    return [[transform(transform2(y)) for y in re.split(r"\s+", x.text.strip())] for x in el.findall(field)]


class VasprunContext(object):

    def __init__(self, logger=None):
        if logger is None:
            logger = logging.getLogger(__name__)
        self.logger = logger

        self._parser = None
        self.bands = None
        self.kpoints = None
        self.weights = None
        self.tetrahedrons = None
        self.tetrahedronVolume = None
        self.ispin = None
        self.ibrion = None
        self.labels = None
        self.vbTopE = None
        self.ebMinE = None
        self.e_fermi = None
        self.cell = None
        self.angstrom_cell = None
        self.unknown_incars = {}
        self._energies = []

    section_map = {
        "modeling": ["section_run", "section_method"],
        "structure": ["section_system"],
        "calculation": ["section_single_configuration_calculation"]
    }

    @property
    def parser(self):
        return self._parser

    @parser.setter
    def parser(self, name):
        self._parser = name

    def on_end_generator(self, element, path_str):
        run = self.parser.run
        program_name = g(element, "i/[@name='program']")
        if program_name.strip().upper() == 'VASP':
            run.program_name = 'VASP'
        else:
            raise Exception('unexpected program name: %s' % program_name)

        version = (
            g(element, "i/[@name='version']", "") + " " +
            g(element, "i/[@name='subversion']", "") + " " +
            g(element, "i/[@name='platform']", ""))

        if not version.isspace():
            run.program_version = version

        run.program_basis_set_type = "plane waves"
        date = g(element, "i/[@name='date']")
        pdate = None
        time = g(element, "i/[@name='time']")

        if date:
            pdate = datetime.strptime(date.strip(), "%Y %m %d")
        if pdate and time:
            pdate = datetime.combine(pdate.date(), datetime.strptime(
                time.strip(), "%H:%M:%S").timetz())
        if pdate:
            run.program_compilation_datetime = seconds_from_epoch(pdate)
        for i in element:
            if i.tag != "i" or not i.attrib.get("name") in set(
                [
                    'program', 'version', 'subversion', 'platform', 'program_version',
                    'date', 'time']):
                pass

    def on_end_incar(self, element, path_str):
        run = self.parser.run
        m_env = self.parser.metainfo_env
        section_method = run.m_create(msection_method)
        dft_plus_u = False
        ibrion = None
        nsw = 0
        for el in element:
            if el.tag == "v":
                name = el.attrib.get("name", None)
                meta = m_env.resolve_definition('x_vasp_incar_' + name, Quantity)
                if not meta:
                    self.logger.warn(
                        "Unknown INCAR parameter (not registered in the meta data): %s %s %r" % (
                            el.tag, el.attrib, el.text))
                    continue
                vector_val = np.asarray(get_vector(el))
                setattr(section_method, meta.get('name'), vector_val)
            elif el.tag == "i":
                name = el.attrib.get("name", None)
                meta = m_env.resolve_definition('x_vasp_incar_' + name, Quantity)
                val_type = el.attrib.get("type")
                if not meta:
                    self.logger.warn(
                        "Unknown INCAR parameter (not registered in the meta data): %s %s %r" % (
                            el.tag, el.attrib, el.text))
                elif val_type:
                    expected_meta_type = {
                        'string': ['C'],
                        'int': ['i'],
                        'logical': ['b', 'C']
                    }.get(val_type)
                    metainfo_type = self._metainfo_type(meta)
                    if not expected_meta_type:
                        self.logger.warn("Unknown value type %s encountered in INCAR: %s %s %r" % (
                            val_type, el.tag, el.attrib, el.text))
                    elif not metainfo_type in expected_meta_type:
                        self.logger.warn("type mismatch between meta data %s and INCAR type %s for %s %s %r" % (
                            metainfo_type, val_type, el.tag, el.attrib, el.text))
                    else:
                        shape = meta.get("shape")
                        converter = meta_type_transformers.get(metainfo_type)
                        if not converter:
                            self.logger.warn(
                                "could not find converter for %s when handling meta info %s" % (
                                    metainfo_type, meta))
                        elif shape:
                            vals = re.split(r"\s+", el.text.strip())
                            setattr(section_method, meta['name'], [converter(x) for x in vals])
                        else:
                            setattr(section_method, meta["name"], converter(el.text))
                    if name == 'GGA':
                        # FIXME tmk: many options are not coded yet. See
                        # https://www.vasp.at/wiki/index.php/GGA
                        f_map = {
                            '91': ['GGA_X_PW91', 'GGA_C_PW91'],
                            'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
                            'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
                            'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
                            'MK': ['GGA_X_OPTB86_VDW']
                        }
                        functs = f_map.get(el.text.strip(), None)
                        if not functs:
                            self.logger.warn("Unknown XC functional %s" % el.text.strip())
                        else:
                            for f in functs:
                                section_XC_functionals = section_method.m_create(msection_XC_functionals)
                                section_XC_functionals.XC_functional_name = f
                    elif name == "ISPIN":
                        self.ispin = int(el.text.strip())
                    elif name == "LDAU":
                        if re.match(".?[Tt](?:[rR][uU][eE])?.?|[yY](?:[eE][sS])?|1", el.text.strip()):
                            dft_plus_u = True
                    elif name == "IBRION":
                        ibrion = int(el.text.strip())
                    elif name == "NSW":
                        nsw = int(el.text.strip())
            else:
                self.logger.warn(
                    "unexpected tag %s %s %r in incar" % (el.tag, el.attrib, el.text))

        if ibrion is None:
            ibrion = -1 if nsw == 0 or nsw == 1 else 0
        if nsw == 0:
            ibrion = -1
        self.ibrion = ibrion
        section_method.electronic_structure_method = 'DFT+U' if dft_plus_u else 'DFT'

    def on_end_kpoints(self, element, path_str):
        run = self.parser.run
        self.bands = None
        self.kpoints = None
        self.weights = None
        for el in element:
            if el.tag == "generation":
                param = el.attrib.get("param", None)  # eg. listgenerated, Monkhorst-Pack, Gamma
                if param:
                    run.section_method[-1].x_vasp_k_points_generation_method = param
                if param == "listgenerated":
                    # This implies a path on k-space, potentially a bandstructure calculation
                    # Save k-path info into a dictionary
                    self.bands = {
                        "divisions": g(el, "i/[@name='divisions']", None),
                        "points": get_vector(el)
                    }

                elif param in ["Monkhorst-Pack", "Gamma"]:
                    # This implies a (2D|3D) mesh on k-space, i.e., not a badstructure calculation
                    # Hence, do nothing: k-points will be stored in the `varray` if-block
                    pass
                else:
                    self.logger.warn("Unknown k point generation method '%s'" % param)
            elif el.tag == "varray":
                name = el.attrib.get("name", None)
                if name == "kpointlist":
                    self.kpoints = np.asarray(get_vector(el))
                    run.section_method[-1].k_mesh_points = self.kpoints
                elif name == "weights":
                    self.weights = np.asarray(get_vector(el))
                    run.section_method[-1].k_mesh_weights = self.weights.flatten()
                elif name == "tetrahedronlist":
                    self.tetrahedrons = np.asarray(get_vector(el), dtype=np.int)
                    run.section_method[-1].x_vasp_tetrahedrons_list = self.tetrahedrons
                else:
                    self.logger.warn("Unknown array %s in kpoints" % name)
            elif el.tag == "i":
                name = el.attrib.get("name", None)
                if name == "volumeweight":
                    ang2m = convert_unit_function("angstrom", "m")

                    # get volume and transform to meters^3
                    vol_cubic_angs = float(el.text.strip())
                    vol_cubic_meters = ang2m(ang2m(ang2m(vol_cubic_angs)))

                    run.section_method[-1].x_vasp_tetrahedron_volume = vol_cubic_meters
            else:
                self.logger.pwarn("Unknown tag %s in kpoints" % el.tag)

    def on_end_structure(self, element, path_str):
        run = self.parser.run
        section_system = run.m_create(msection_system)
        self.cell = None
        for el in element:
            if (el.tag == "crystal"):
                for cell_el in el:
                    if cell_el.tag == "varray":
                        name = cell_el.attrib.get("name", None)
                        if name == "basis":
                            conv = convert_unit_function("angstrom", "m")
                            self.cell = get_vector(
                                cell_el, lambda x: conv(float(x)))
                            self.angstrom_cell = np.array(get_vector(cell_el))
                            section_system.simulation_cell = np.asarray(self.cell)
                            section_system.configuration_periodic_dimensions = np.ones(3, dtype=bool)
                        elif name == "rec_basis":
                            pass
                        else:
                            self.logger.warn(
                                "Unexpected varray %s in crystal" % name)
                    elif cell_el.tag == "i":
                        if cell_el.attrib.get("name") != "volume":
                            self.logger.warn(
                                "Unexpected i value %s in crystal" % cell_el.attrib)
                    else:
                        self.logger.warn("Unexpected tag %s %s %r in crystal" % (
                            cell_el.tag, cell_el.attrib, cell_el.text))
            elif el.tag == "varray":
                name = el.attrib.get("name", None)
                if name == "positions":
                    pos = get_vector(el)
                    section_system.atom_positions = np.dot(np.asarray(pos), self.cell)
                elif name == "selective":
                    atom_sel = get_vector(el, transform=lambda item: item == 'T')
                    section_system.x_vasp_selective_dynamics = np.asarray(atom_sel, dtype=np.bool)
                else:
                    self.logger.warn(
                        "Unexpected varray in structure %s" % el.attrib)
            elif el.tag == "nose":
                nose = get_vector(el)
                section_system.x_vasp_nose_thermostat = nose
            else:
                self.logger.warn(
                    "Unexpected tag in structure %s %s %r" % (el.tag, el.attrib, el.text))
        if self.labels is not None:
            section_system.atom_labels = self.labels

    def on_end_eigenvalues(self, element, path_str):
        if path_str != "modeling/calculation/eigenvalues":
            return True
        run = self.parser.run
        eigenvalues = None
        occupation = None
        for el in element:
            if el.tag == "array":
                for arrray_el in el:
                    if arrray_el.tag == "dimension":
                        pass
                    elif arrray_el.tag == "field":
                        pass
                    elif arrray_el.tag == "set":
                        isp = -1
                        for spin_el in arrray_el:
                            if spin_el.tag == "set":
                                ik = -1
                                isp += 1
                                for kEl in spin_el:
                                    if kEl.tag == "set":
                                        ik += 1
                                        bands = np.asarray(
                                            get_vector(kEl, field="r"))
                                        if eigenvalues is None:
                                            eigenvalues = np.zeros(
                                                (self.ispin, self.kpoints.shape[0], bands.shape[0]), dtype=float)
                                            occupation = np.zeros(
                                                (self.ispin, self.kpoints.shape[0], bands.shape[0]), dtype=float)
                                        eigenvalues[isp, ik] = bands[:, 0]
                                        occupation[isp, ik] = bands[:, 1]
                                    else:
                                        self.logger.warn(
                                            "unexpected tag %s in k array of the eigenvalues" % kEl.tag)
                            else:
                                self.logger.warn(
                                    "unexpected tag %s in spin array of the eigenvalues" % spin_el.tag)
                    else:
                        self.logger.warn(
                            "unexpected tag %s in array of the eigenvalues" % arrray_el.tag)
                if eigenvalues is not None:

                    ev = eV2JV(eigenvalues)
                    vbTopE = []
                    ebMinE = []
                    for ispin in range(occupation.shape[0]):
                        vbTopE.append(float('-inf'))
                        ebMinE.append(float('inf'))
                        for ik in range(occupation.shape[1]):
                            ebIndex = bisect.bisect_right(
                                -occupation[ispin, ik, :], -0.5) - 1
                            vbTopIndex = ebIndex - 1
                            if vbTopIndex >= 0:
                                vbTopK = ev[ispin, ik, vbTopIndex]
                                if vbTopK > vbTopE[ispin]:
                                    vbTopE[ispin] = vbTopK
                            if ebIndex < ev.shape[2]:
                                ebMinK = ev[ispin, ik, ebIndex]
                                if ebMinK < ebMinE[ispin]:
                                    ebMinE[ispin] = ebMinK
                    self.vbTopE = vbTopE
                    self.ebMinE = ebMinE
                    run.section_single_configuration_calculation[-1].energy_reference_highest_occupied = np.array(vbTopE)
                    run.section_single_configuration_calculation[-1].energy_reference_lowest_unoccupied = np.array(ebMinE)
                    if self.bands:
                        divisions = int(self.bands['divisions'])
                        section_k_band = run.section_single_configuration_calculation[-1].m_create(msection_k_band)
                        nsegments = self.kpoints.shape[0] // divisions
                        kpt = np.reshape(
                            self.kpoints, (nsegments, divisions, 3))
                        energies = np.reshape(
                            ev, (self.ispin, nsegments, divisions, bands.shape[0]))
                        occ = np.reshape(
                            occupation, (self.ispin, nsegments, divisions, bands.shape[0]))
                        for isegment in range(nsegments):
                            section_k_band_segment = section_k_band.m_create(msection_k_band_segment)
                            section_k_band_segment.band_energies = energies[:, isegment, :, :]
                            section_k_band_segment.band_occupations = occ[:, isegment, :, :]
                            section_k_band_segment.band_k_points = kpt[isegment]
                            # "band_segm_labels"
                            section_k_band_segment.band_segm_start_end = np.asarray(
                                [kpt[isegment, 0], kpt[isegment, divisions - 1]])

                        section_k_band_normalized = run.section_single_configuration_calculation[-1].m_create(msection_k_band_normalized)
                        for isegment in range(nsegments):
                            section_k_band_segment_normalized = section_k_band_normalized.m_create(msection_k_band_segment_normalized)
                            section_k_band_segment_normalized.band_energies_normalized = energies[:, isegment, :, :] - max(self.vbTopE)
                            section_k_band_segment_normalized.band_occupations_normalized = occ[:, isegment, :, :]
                            section_k_band_segment_normalized.band_k_points_normalized = kpt[isegment]
                            section_k_band_segment_normalized.band_segm_start_end_normalized = np.asarray(
                                [kpt[isegment, 0], kpt[isegment, divisions - 1]])

                    else:
                        section_eigenvalues = run.section_single_configuration_calculation[-1].m_create(msection_eigenvalues)
                        section_eigenvalues.eigenvalues_values = ev
                        section_eigenvalues.eigenvalues_occupation = occupation
            else:
                self.logger.warn("unexpected tag %s in the eigenvalues" % el.tag)

    def on_start_calculation(self, element, path_str):
        run = self.parser.run
        sscc = run.m_create(msection_single_configuration_calculation)
        if self.waveCut:
            section_basis_set = sscc.m_create(msection_basis_set)
            section_basis_set.mapping_section_basis_set_cell_dependent = self.waveCut

    def on_end_modeling(self, element, path_str):
        run = self.parser.run
        run.section_method[-1].x_vasp_unknown_incars = self.unknown_incars
        if self.ibrion is None or self.ibrion == -1:
            return
        section_sampling_method = run.m_create(msection_sampling_method)
        if self.ibrion == 0:
            sampling_method = "molecular_dynamics"
        else:
            sampling_method = "geometry_optimization"
        section_sampling_method.sampling_method = sampling_method
        section_frame_sequence = run.m_create(msection_frame_sequence)
        section_frame_sequence.frame_sequence_to_sampling_ref = section_sampling_method
        section_frame_sequence.frame_sequence_local_frames_ref = run.section_single_configuration_calculation[-1]

    def on_end_calculation(self, element, path_str):
        e_conv = eV2J
        f_conv = convert_unit_function("eV/angstrom", "N")
        p_conv = convert_unit_function("eV/angstrom^3", "Pa")
        run = self.parser.run
        sscc = run.section_single_configuration_calculation[-1]
        sscc.single_configuration_calculation_to_system_ref = run.section_system[-1]
        sscc.single_configuration_to_calculation_method_ref = run.section_method[-1]
        for el in element:
            if el.tag == "energy":
                for en_el in el:
                    if en_el.tag == "i":
                        name = en_el.attrib.get("name", None)
                        if name == "e_fr_energy":
                            value = e_conv(float(en_el.text.strip()))
                            sscc.energy_free = value
                        elif name == "e_wo_entrp":
                            value = e_conv(float(en_el.text.strip()))
                            sscc.energy_total = value
                        elif name == "e_0_energy":
                            value = e_conv(float(en_el.text.strip()))
                            sscc.energy_total_T0 = value
                            self._energies.append(value)
                    elif en_el.tag == "varray":
                        pass

            elif el.tag == 'varray':
                name = el.attrib.get("name", None)
                if name == "forces":
                    f = get_vector(el, lambda x: f_conv(float(x)))
                    sscc.atom_forces = f
                elif name == 'stress':
                    f = get_vector(el, lambda x: p_conv(float(x)))
                    sscc.stress_tensor = f

    def on_end_atominfo(self, element, path_str):
        run = self.parser.run
        atom_types = []
        labels = []
        labels2 = None
        atom_types_desc = []
        for el in element:
            if el.tag in ['atoms', 'types']:
                pass

            elif el.tag == "array":
                name = el.attrib.get("name", None)
                if name == "atoms":
                    for atoms_el in el:
                        if atoms_el.tag == "dimension":
                            pass
                        elif atoms_el.tag == "field":
                            pass
                        elif atoms_el.tag == "set":
                            for atoms_line in atoms_el:
                                if atoms_line.tag != "rc":
                                    self.logger.warn(
                                        "unexpected tag %s in atoms array in atominfo" % atoms_line.tag)
                                else:
                                    line = atoms_line.findall("c")
                                    labels.append(line[0].text.strip())
                                    atom_types.append(int(line[1].text.strip()))
                        else:
                            self.logger.warn(
                                "unexpected tag %s in atoms array in atominfo" % atoms_el.tag)
                elif name == "atomtypes":
                    keys = []
                    field_types = []
                    for atoms_el in el:
                        if atoms_el.tag == "dimension":
                            pass
                        elif atoms_el.tag == "field":
                            keys.append(atoms_el.text.strip())
                            field_types.append(
                                atoms_el.attrib.get("type", "float"))
                        elif atoms_el.tag == "set":
                            expected_keys = [
                                "atomspertype", "element", "mass", "valence", "pseudopotential"]
                            if keys != expected_keys:
                                self.logger.warn(
                                    "unexpected fields in atomtype: %s vs %s" % (keys, expected_keys))
                            for atoms_line in atoms_el:
                                if atoms_line.tag != "rc":
                                    self.logger.warn(
                                        "unexpected tag %s in atoms array in atominfo" % atoms_line.tag)
                                else:
                                    line = atoms_line.findall("c")
                                    type_desc = {}
                                    for i, k in enumerate(keys):
                                        field_type = field_types[i]
                                        value = line[i].text
                                        if field_type == "float":
                                            value = float(value)
                                        elif field_type == "int":
                                            value = int(value)
                                        else:
                                            pass
                                        type_desc[k] = value
                                    atom_types_desc.append(type_desc)
                        else:
                            self.logger.warn(
                                "unexpected tag %s in atomtypes array in atominfo" % atoms_el.tag)
                    n_el = {}
                    kind_labels = []
                    for atom_desc in atom_types_desc:
                        section_method_atom_kind = run.section_method[-1].m_create(msection_method_atom_kind)
                        if 'element' in atom_desc:
                            elName = atom_desc['element'].strip()
                            try:
                                elNr = ase.data.chemical_symbols.index(elName)
                                section_method_atom_kind.method_atom_kind_atom_number = elNr
                            except Exception as e:
                                self.logger.error(
                                    "error finding element number for %r" % atom_desc['element'].strip(),
                                    exc_info=e)
                            n_el_now = 1 + n_el.get(elName, 0)
                            n_el[elName] = n_el_now
                            el_label = elName + \
                                (str(n_el_now) if n_el_now > 1 else "")
                            kind_labels.append(el_label)
                            section_method_atom_kind.method_atom_kind_label = el_label
                            if "mass" in atom_desc:
                                section_method_atom_kind.method_atom_kind_mass = atom_desc["mass"]
                            if "valence" in atom_desc:
                                section_method_atom_kind.method_atom_kind_explicit_electrons = atom_desc["valence"]
                            if "pseudopotential" in atom_desc:
                                section_method_atom_kind.method_atom_kind_pseudopotential_name = atom_desc["pseudopotential"]
                    run.section_method[-1].x_vasp_atom_kind_refs = run.section_method[-1].section_method_atom_kind
                    labels2 = [kind_labels[i - 1] for i in atom_types]
                else:
                    self.logger.warn(
                        "unexpected array named %s in atominfo" % name)
            else:
                self.logger.warn("unexpected tag %s in atominfo" % el.tag)
        self.labels = np.asarray(labels2) if labels2 else np.asarray(labels)

    def _metainfo_type(self, meta):
        dtype = meta.type
        if dtype == str:
            return 'C'
        elif dtype == bool:
            return 'b'
        elif dtype == int or type(dtype) == np.dtype:
            return 'i'
        elif dtype == float:
            return 'f'
        else:
            return

    def _incar_out_tag(self, el):
        if (el.tag != "i"):
            self.logger.warn("unexpected tag %s %s %r in incar" % (el.tag, el.attrib, el.text))

        else:
            run = self.parser.run
            name = el.attrib.get("name", None)
            val_type = el.attrib.get("type")
            try:
                meta = self.parser.metainfo_env.resolve_definition('x_vasp_incarOut_' + name, Quantity)
            except Exception:
                return

            if not meta:
                # Unknown_Incars_Begin: storage into a dictionary
                if not val_type:
                    # On vasp's xml files, val_type *could* be absent if incar value is float
                    val_type = 'float'

                # map vasp's datatype to nomad's datatype [b, f, i, C, D, R]
                nomad_dtype_str = vasp_to_metainfo_type_mapping[val_type][0]

                converter = meta_type_transformers.get(nomad_dtype_str)
                text_value = el.text.strip()  # text representation of incar value
                try:
                    pyvalue = converter(text_value)  # python data type
                except Exception:
                    pyvalue = text_value

                # save (name, pyvalue) into a dict
                self.unknown_incars[name] = pyvalue
                # Unknown_Incars_end
            else:
                if not val_type:
                    val_type = 'float'

                vasp_metainfo_type = vasp_to_metainfo_type_mapping.get(val_type)[0]
                metainfo_type = self._metainfo_type(meta)
                if not vasp_metainfo_type:
                    self.logger.warn("Unknown value type %s encountered in INCAR out: %s %s %r" % (
                        val_type, el.tag, el.attrib, el.text))

                elif metainfo_type != vasp_metainfo_type:
                    if (metainfo_type == 'C' and vasp_metainfo_type == 'b'):
                        pass
                    elif (metainfo_type == 'i' and vasp_metainfo_type == 'f'):
                        pass
                    else:
                        self.logger.warn(
                            "Data type mismatch: %s. Vasp_type: %s, metainfo_type: %s " % (
                                name, vasp_metainfo_type, metainfo_type))
                try:
                    shape = meta.get("shape")
                    converter = meta_type_transformers.get(metainfo_type)
                    if not converter:
                        self.logger.warn(
                            "could not find converter for %s when handling meta info %s" %
                            (metainfo_type, meta))
                    elif shape:
                        vals = re.split(r"\s+", el.text.strip())
                        setattr(run.section_method[-1], meta["name"], [converter(x) for x in vals])
                    else:
                        # If-block to handle incars without value
                        if not el.text:
                            el.text = ''
                        setattr(run.section_method[-1], meta["name"], converter(el.text))

                except Exception:
                    self.logger.warn("Exception trying to handle incarOut %s: %s" % (
                        name, traceback.format_exc()))

                if name == 'ENMAX' or name == 'PREC':
                    if name == 'ENMAX':
                        self.enmax = converter(el.text)
                    if name == 'PREC':
                        if 'acc' in converter(el.text):
                            self.prec = 1.3
                        else:
                            self.prec = 1.0
                if name == 'GGA':
                    f_map = {
                        '91': ['GGA_X_PW91', 'GGA_C_PW91'],
                        'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
                        'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
                        'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
                        'MK': ['GGA_X_OPTB86_VDW'],
                        '--': ['GGA_X_PBE', 'GGA_C_PBE']  # should check potcar
                    }
                    functs = f_map.get(el.text.strip(), None)
                    if not functs:
                        self.logger.warn("Unknown XC functional %s" % el.text.strip())
                    else:
                        for f in functs:
                            section_XC_functionals = run.section_method[-1].m_create(msection_XC_functionals)
                            section_XC_functionals.XC_functional_name = f
                elif name == "ISPIN":
                    self.ispin = int(el.text.strip())

    def _separator_scan(self, element, depth=0):
        for separators in element:
            if separators.tag == "separator":
                for el in separators:
                    if el.tag == "i":
                        self._incar_out_tag(el)
                    elif el.tag == "separator":
                        self._separator_scan(el, depth + 1)
                    else:
                        pass
            elif separators.tag == "i":
                self._incar_out_tag(separators)
            else:
                pass

    def on_end_parameters(self, element, path_str):
        self._separator_scan(element)
        run = self.parser.run
        try:
            self.prec
            try:
                self.enmax
                self.waveCut = run.m_create(msection_basis_set_cell_dependent)
                self.waveCut.basis_set_planewave_cutoff = eV2J(self.enmax * self.prec)

                section_method_basis_set = run.section_method[-1].m_create(msection_method_basis_set)
                section_method_basis_set.mapping_section_method_basis_set_cell_associated = self.waveCut
            except AttributeError:
                import traceback
                traceback.print_exc()
                self.logger.warn(
                    "Missing ENMAX for calculating plane wave basis cut off ")
        except AttributeError:
            self.logger.pwarn(
                "Missing PREC for calculating plane wave basis cut off ")

    def on_end_dos(self, element, path_str):
        "density of states"
        run = self.parser.run
        section_dos = run.section_single_configuration_calculation[-1].m_create(msection_dos)
        for el in element:
            if el.tag == "i":
                if el.attrib.get("name") == "efermi":
                    self.e_fermi = eV2J(float(el.text.strip()))
                    section_dos.energy_reference_fermi = np.array([self.e_fermi] * self.ispin)
                else:
                    self.logger.warn("unexpected tag %s %s in dos" % (el.tag, el.attrib))
            elif el.tag == "total":
                for el1 in el:
                    if el1.tag == "array":
                        for el2 in el1:
                            if el2.tag == "dimension" or el2.tag == "field":
                                pass
                            elif el2.tag == "set":
                                dosL = []
                                for spin_component in el2:
                                    if spin_component.tag == "set":
                                        dosL.append(
                                            get_vector(spin_component, field="r"))
                                    else:
                                        self.logger.warn("unexpected tag %s %s in dos total array set" % (
                                            spin_component.tag, spin_component.attrib))
                                dos_a = np.asarray(dosL)
                                if len(dos_a.shape) != 3:
                                    raise Exception("unexpected shape %s (%s) for total dos (ragged arrays?)" % (
                                        dos_a.shape), dos_a.dtype)
                                dos_e = eV2JV(dos_a[0, :, 0])
                                dos_i = dos_a[:, :, 2]
                                dos_v = dos_a[:, :, 1]

                                # Convert the DOS values to SI. VASP uses the
                                # following units in the output:
                                # states/eV/cell. This means that the volume
                                # dependence has been introduced by multiplying
                                # by the cell volume
                                # the integrated dos value is the number of electrons until that energy level
                                # and thus not directly energy dependent anymore
                                joule_in_ev = convert_unit(1, "eV", "J")
                                dos_v = dos_v / joule_in_ev

                                section_dos.dos_energies = dos_e
                                cell_volume = np.abs(np.linalg.det(self.cell))
                                section_dos.dos_values = dos_v * cell_volume
                                section_dos.dos_integrated_values = dos_i
                            else:
                                self.logger.warn("unexpected tag %s %s in dos total array" % (
                                    el2.tag, el2.attrib))
                    else:
                        self.logger.warn("unexpected tag %s %s in dos total" % (
                            el2.tag, el2.attrib))
            elif el.tag == "partial":
                for el1 in el:
                    if el1.tag == "array":
                        lm = []
                        for el2 in el1:
                            if el2.tag == "dimension":
                                pass
                            elif el2.tag == "field":
                                if el2.text.strip() == "energy":
                                    pass
                                else:
                                    strLm = {
                                        "s": [0, 0],
                                        "p": [1, -1],
                                        "px": [1, 0],
                                        "py": [1, 1],
                                        "pz": [1, 2],
                                        "d": [2, -1],
                                        "dx2": [2, 0],
                                        "dxy": [2, 1],
                                        "dxz": [2, 2],
                                        "dy2": [2, 3],
                                        "dyz": [2, 4],
                                        "dz2": [2, 5]
                                    }
                                    lm.append(
                                        strLm.get(el2.text.strip(), [-1, -1]))
                            elif el2.tag == "set":
                                dosL = []
                                for atom in el2:
                                    if atom.tag == "set":
                                        atomL = []
                                        dosL.append(atomL)
                                        for spin_component in atom:
                                            if spin_component.tag == "set":
                                                atomL.append(
                                                    get_vector(spin_component, field="r"))
                                            else:
                                                self.logger.warn("unexpected tag %s %s in dos partial array set set" % (
                                                    spin_component.tag, spin_component.attrib))
                                    else:
                                        self.logger.warn("unexpected tag %s %s in dos partial array set" % (
                                            spin_component.tag, spin_component.attrib))
                                dosLM = np.asarray(dosL)
                                assert len(
                                    dosLM.shape) == 4, "invalid shape dimension in projected dos (ragged arrays?)"
                                section_dos.dos_values_lm = dosLM[:, :, :, 1:]
                            else:
                                self.logger.warn("unexpected tag %s %s in dos total array" % (
                                    el2.tag, el2.attrib))
                        section_dos.dos_lm = np.asarray(lm)
                        section_dos.dos_m_kind = 'polynomial'
                    else:
                        self.logger.warn("unexpected tag %s %s in dos total" % (
                            el2.tag, el2.attrib))
            else:
                self.logger.warn("unexpected tag %s %s in dos" % (el2.tag, el2.attrib))

    def on_end_projected(self, element, path_str):
        return None


class XmlParser(object):
    @staticmethod
    def extract_callbacks(obj):
        """extracts all callbacks from the object obj

        triggers should start with onStart_ or onEnd__ and then have a valid section name.
        They will be called with this object, the event and current element
        """
        triggers = {}
        for attr in dir(obj):
            if attr.startswith("on_start_"):
                triggers[attr] = getattr(obj, attr)
            elif attr.startswith("on_end_"):
                triggers[attr] = getattr(obj, attr)
        return triggers

    @staticmethod
    def maybe_get(el, path, default=None):
        i = el.findall(path)
        if i:
            return i.pop().text
        else:
            return default

    def __init__(self, parser_info, super_context, metainfo_env=None, callbacks=None, section_map=None):
        self.f_in = None
        self.parser_info = parser_info
        self.super_context = super_context
        self.callbacks = callbacks if callbacks is not None else XmlParser.extract_callbacks(
            super_context)
        self.section_map = section_map if section_map is not None else super_context.section_map
        self.path = []
        self.tagSections = {}
        self.metainfo_env = metainfo_env
        self.run = None
        self.archive = None

    def parse(self, main_file_uri, f_in, archive):
        if self.path:
            raise Exception(
                "Parse of %s called with non empty path, parse already in progress?" % main_file_uri)
        self.main_file_uri = main_file_uri
        self.f_in = f_in
        self.archive = archive
        self.run = archive.m_create(msection_run)
        self.super_context.parser = self
        # there are invalid characters like esc in the files, we do not want to crash on them
        xml_parser = MyXMLParser()
        try:
            for event, el in xml.etree.ElementTree.iterparse(self.f_in, events=["start", "end"], parser=xml_parser):
                if event == 'start':
                    path_str = "/".join([x.tag for x in self.path]) + "/" + el.tag
                    callback = self.callbacks.get("on_start_" + el.tag, None)

                    if callback:
                        callback(el, path_str)
                    self.path.append(el)
                elif event == 'end':
                    last_el = self.path.pop()
                    if last_el != el:
                        raise Exception(
                            "mismatched path at end, got %s expected %s" % (last_el, el))
                    tag = el.tag
                    path_str = "/".join([x.tag for x in self.path]) + "/" + tag
                    callback = self.callbacks.get("on_end_" + tag, None)
                    if callback:
                        if not callback(el, path_str):
                            # if callback does not return True then assume that the current element has been processed
                            # and can be removed
                            el.clear()
                            if self.path:
                                self.path[-1].remove(el)
                    elif len(self.path) == 1:
                        el.clear()
                        self.path[-1].remove(el)
                else:
                    raise Exception("Unexpected event %s" % event)
        except ParseError as e:
            self.super_context.logger.warn("Could not complete parsing: %s" % e, exc_info=e)
        except Exception as e:
            import traceback
            traceback.print_exc()
        else:
            pass


g = XmlParser.maybe_get

parser_info = {
    "name": "parser_vasprun",
    "version": "1.0"
}
