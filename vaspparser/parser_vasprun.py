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
import sys
import bisect
from datetime import datetime
import os
import re
import traceback
import numpy as np
import ase.geometry
import ase.data
from math import pi
import xml.etree.ElementTree as ET
import logging

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.unit_conversion.unit_conversion import convert_unit

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


def secondsFromEpoch(date):
    epoch = datetime(1970, 1, 1)
    ts = date-epoch
    return ts.seconds + ts.microseconds/1000.0


trueRe = re.compile(
    r"\s*(?:\.?[Tt](?:[Rr][Uu][Ee])?\.?|1|[Yy](?:[Ee][Ss])?|[Jj][Aa]?)\s*$")
falseRe = re.compile(
    r"\s*(?:\.?[fF](?:[Aa][Ll][Ss][Ee])?\.?|0|[Nn](?:[Oo]|[Ee][Ii][Nn])?)\s*$")


def toBool(value):
    if falseRe.match(value):
        return False
    elif trueRe.match(value):
        return True
    else:
        backend.pwarn("Unexpected value for boolean field: %s" % (value))
        return None

metaTypeTransformers = {
    'C': lambda x: x.strip(),
    'i': lambda x: int(float(x.strip())),
    'f': lambda x: float(x.strip()),
    'b': toBool,
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
            if not(num in (0x9, 0xA, 0xD) or 0x20 <= num <= 0xD7FF
                   or 0xE000 <= num <= 0xFFFD or 0x10000 <= num <= 0x10FFFF):
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


def getVector(el, transform=float, field="v"):
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

        self.parser = None
        self.bands = None
        self.kpoints = None
        self.weights = None
        self.tetrahedrons = None
        self.tetrahedronVolume = None
        self.ispin = None
        self.ibrion = None
        self.lastSystemDescription = None
        self.labels = None
        self.singleConfCalcs = []
        self.vbTopE = None
        self.ebMinE = None
        self.eFermi = None
        self.cell = None
        self.angstrom_cell = None
        self.unknown_incars = {}

    sectionMap = {
        "modeling": ["section_run", "section_method"],
        "structure": ["section_system"],
        "calculation": ["section_single_configuration_calculation"]
    }

    def startedParsing(self, parser):
        self.parser = parser

    def onEnd_generator(self, parser, event, element, pathStr):
        backend = parser.backend
        program_name = g(element, "i/[@name='program']")
        if program_name.strip().upper() == "VASP":
            backend.addValue("program_name", "VASP")
        else:
            raise Exception("unexpected program name: %s" % program_name)
        version = (g(element, "i/[@name='version']", "") + " " +
                   g(element, "i/[@name='subversion']", "") + " " +
                   g(element, "i/[@name='platform']", ""))
        if not version.isspace():
            backend.addValue("program_version", version)
        backend.addValue("program_basis_set_type", "plane waves")
        date = g(element, "i/[@name='date']")
        pdate = None
        time = g(element, "i/[@name='time']")
        if date:
            pdate = datetime.strptime(date.strip(), "%Y %m %d")
        if pdate and time:
            pdate = datetime.combine(pdate.date(), datetime.strptime(
                time.strip(), "%H:%M:%S").timetz())
        if pdate:
            backend.addValue("program_compilation_datetime",
                             secondsFromEpoch(pdate))
        for i in element:
            if i.tag != "i" or not i.attrib.get("name") in set(["program", "version", "subversion", "platform", "program_version", "date", "time"]):
                backend.pwarn("unexpected tag %s %s %r in generator" %
                              (i.tag, i.attrib, i.text))

    def onEnd_incar(self, parser, event, element, pathStr):
        backend = parser.backend
        metaEnv = parser.backend.metaInfoEnv()
        dft_plus_u = False
        ibrion = None
        nsw = 0
        for el in element:
            if el.tag == "v":
                name = el.attrib.get("name", None)
                meta = metaEnv['x_vasp_incar_' + name]
                if not meta:
                    backend.pwarn("Unknown INCAR parameter (not registered in the meta data): %s %s %r" % (
                        el.tag, el.attrib, el.text))
                    continue
                #- -
                vector_val = np.asarray(getVector(el))
                backend.addArrayValues(meta.get('name'), vector_val)
            elif el.tag == "i":
                name = el.attrib.get("name", None)
                meta = metaEnv['x_vasp_incar_' + name]
                valType = el.attrib.get("type")
                if not meta:
                    backend.pwarn("Unknown INCAR parameter (not registered in the meta data): %s %s %r" % (
                        el.tag, el.attrib, el.text))
                elif valType:
                    expectedMetaType = {
                        'string': ['C'],
                        'int': ['i'],
                        'logical': ['b', 'C']
                    }.get(valType)
                    if not expectedMetaType:
                        backend.pwarn("Unknown value type %s encountered in INCAR: %s %s %r" % (
                            valType, el.tag, el.attrib, el.text))
                    elif not meta.get('dtypeStr') in expectedMetaType:
                        backend.pwarn("type mismatch between meta data %s and INCAR type %s for %s %s %r" % (
                            meta.get('dtypeStr'), valType, el.tag, el.attrib, el.text))
                    else:
                        shape = meta.get("shape", None)
                        dtypeStr = meta.get("dtypeStr", None)
                        converter = metaTypeTransformers.get(dtypeStr)
                        if not converter:
                            backend.pwarn(
                                "could not find converter for dtypeStr %s when handling meta info %s" % (dtypeStr, ))
                        elif shape:
                            vals = re.split("\s+", el.text.strip())
                            backend.addValue(
                                meta["name"], [converter(x) for x in vals])
                        else:
                            backend.addValue(meta["name"], converter(el.text))
                    if name == 'GGA':
                        # FIXME tmk: many options are not coded yet. See
                        # https://www.vasp.at/wiki/index.php/GGA
                        fMap = {
                            '91': ['GGA_X_PW91', 'GGA_C_PW91'],
                            'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
                            'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
                            'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
                            'MK': ['GGA_X_OPTB86_VDW']
                        }
                        functs = fMap.get(el.text.strip(), None)
                        if not functs:
                            backend.pwarn("Unknown XC functional %s" %
                                          el.text.strip())
                        else:
                            for f in functs:
                                backend.openNonOverlappingSection(
                                    "section_XC_functionals")
                                backend.addValue("XC_functional_name", f)
                                backend.closeNonOverlappingSection(
                                    "section_XC_functionals")
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
                backend.pwarn("unexpected tag %s %s %r in incar" %
                              (el.tag, el.attrib, el.text))
        if ibrion is None:
            ibrion = -1 if nsw == 0 or nsw == 1 else 0
        if nsw == 0:
            ibrion = -1
        self.ibrion = ibrion
        if dft_plus_u:
            backend.addValue("electronic_structure_method", "DFT+U")
        else:
            backend.addValue("electronic_structure_method", "DFT")

    def onEnd_kpoints(self, parser, event, element, pathStr):
        backend = parser.backend
        self.bands = None
        self.kpoints = None
        self.weights = None
        for el in element:
            if el.tag == "generation":
                param = el.attrib.get("param", None) # eg. listgenerated, Monkhorst-Pack, Gamma
                if param:
                    backend.addValue(
                        "x_vasp_k_points_generation_method", param)
                if param == "listgenerated":
                    # This implies a path on k-space, potentially a bandstructure calculation
                    # Save k-path info into a dictionary
                    self.bands = {
                        "divisions": g(el, "i/[@name='divisions']", None),
                        "points": getVector(el)
                    }

                elif param in ["Monkhorst-Pack", "Gamma"]:
                    # This implies a (2D|3D) mesh on k-space, i.e., not a badstructure calculation
                    # Hence, do nothing: k-points will be stored in the `varray` if-block
                    pass
                else:
                    backend.pwarn("Unknown k point generation method '%s'" %(param))
            elif el.tag == "varray":
                name = el.attrib.get("name", None)
                if name == "kpointlist":
                    self.kpoints = np.asarray(getVector(el))
                    backend.addArrayValues("k_mesh_points", self.kpoints)
                elif name == "weights":
                    self.weights = np.asarray(getVector(el))
                    backend.addArrayValues(
                        "k_mesh_weights", self.weights.flatten())
                elif name == "tetrahedronlist":
                    self.tetrahedrons = np.asarray(getVector(el), dtype=np.int)
                    backend.addArrayValues(
                        "x_vasp_tetrahedrons_list", self.tetrahedrons)
                else:
                    backend.pwarn("Unknown array %s in kpoints" % name)
            elif el.tag == "i":
                name = el.attrib.get("name", None)
                if name == "volumeweight":
                    ang2m = convert_unit_function("angstrom", "m")

                    # get volume and transform to meters^3
                    vol_cubic_angs = float(el.text.strip())
                    vol_cubic_meters = ang2m(ang2m(ang2m(vol_cubic_angs)))

                    backend.addArrayValues("x_vasp_tetrahedron_volume",
                        vol_cubic_meters)
            else:
                backend.pwarn("Unknown tag %s in kpoints" % el.tag)


    def onEnd_structure(self, parser, event, element, pathStr):
        backend = parser.backend
        gIndexes = parser.tagSections[pathStr]
        self.lastSystemDescription = gIndexes["section_system"]
        self.cell = None
        for el in element:
            if (el.tag == "crystal"):
                for cellEl in el:
                    if cellEl.tag == "varray":
                        name = cellEl.attrib.get("name", None)
                        if name == "basis":
                            conv = convert_unit_function("angstrom", "m")
                            self.cell = getVector(
                                cellEl, lambda x: conv(float(x)))
                            self.angstrom_cell = np.array(getVector(cellEl))
                            backend.addArrayValues(
                                "simulation_cell", np.asarray(self.cell))
                            backend.addArrayValues(
                                "configuration_periodic_dimensions", np.ones(3, dtype=bool))
                        elif name == "rec_basis":
                            pass
                        else:
                            backend.pwarn(
                                "Unexpected varray %s in crystal" % name)
                    elif cellEl.tag == "i":
                        if cellEl.attrib.get("name") != "volume":
                            backend.pwarn(
                                "Unexpected i value %s in crystal" % cellEl.attrib)
                    else:
                        backend.pwarn("Unexpected tag %s %s %r in crystal" % (
                            cellEl.tag, cellEl.attrib, cellEl.text))
            elif el.tag == "varray":
                name = el.attrib.get("name", None)
                if name == "positions":
                    pos = getVector(el)
                    backend.addArrayValues(
                        "atom_positions", np.dot(np.asarray(pos), self.cell))
                elif name == "selective":
                    atom_sel = getVector(el, transform=lambda item: item == 'T')
                    backend.addArrayValues(
                        "x_vasp_selective_dynamics", np.asarray(atom_sel, dtype=np.bool))
                else:
                    backend.pwarn(
                        "Unexpected varray in structure %s" % el.attrib)
            elif el.tag == "nose":
                nose = getVector(el)
                backend.addArrayValues("x_vasp_nose_thermostat", nose)
            else:
                backend.pwarn("Unexpected tag in structure %s %s %r" %
                              (el.tag, el.attrib, el.text))
        if self.labels is not None:
            backend.addArrayValues("atom_labels", self.labels)

    def onEnd_eigenvalues(self, parser, event, element, pathStr):
        if pathStr != "modeling/calculation/eigenvalues":
            return True
        backend = parser.backend
        eigenvalues = None
        occupation = None
        for el in element:
            if el.tag == "array":
                for arrEl in el:
                    if arrEl.tag == "dimension":
                        pass
                    elif arrEl.tag == "field":
                        pass
                    elif arrEl.tag == "set":
                        isp = -1
                        for spinEl in arrEl:
                            if spinEl.tag == "set":
                                ik = -1
                                isp += 1
                                for kEl in spinEl:
                                    if kEl.tag == "set":
                                        ik += 1
                                        bands = np.asarray(
                                            getVector(kEl, field="r"))
                                        if eigenvalues is None:
                                            eigenvalues = np.zeros(
                                                (self.ispin, self.kpoints.shape[0],  bands.shape[0]), dtype=float)
                                            occupation = np.zeros(
                                                (self.ispin, self.kpoints.shape[0],  bands.shape[0]), dtype=float)
                                        eigenvalues[isp, ik] = bands[:, 0]
                                        occupation[isp, ik] = bands[:, 1]
                                    else:
                                        backend.pwarn(
                                            "unexpected tag %s in k array of the eigenvalues" % kEl.tag)
                            else:
                                backend.pwarn(
                                    "unexpected tag %s in spin array of the eigenvalues" % spinEl.tag)
                    else:
                        backend.pwarn(
                            "unexpected tag %s in array of the eigenvalues" % arrEl.tag)
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
                    backend.addArrayValues(
                        "energy_reference_highest_occupied", np.array(vbTopE))
                    backend.addArrayValues(
                        "energy_reference_lowest_unoccupied", np.array(ebMinE))
                    if self.bands:
                        divisions = int(self.bands['divisions'])
                        backend.openNonOverlappingSection("section_k_band")
                        nsegments = self.kpoints.shape[0] // divisions
                        kpt = np.reshape(
                            self.kpoints, (nsegments, divisions, 3))
                        energies = np.reshape(
                            ev, (self.ispin, nsegments, divisions, bands.shape[0]))
                        occ = np.reshape(
                            occupation, (self.ispin, nsegments, divisions, bands.shape[0]))
                        for isegment in range(nsegments):
                            backend.openNonOverlappingSection(
                                "section_k_band_segment")
                            backend.addArrayValues(
                                "band_energies", energies[:, isegment, :, :])
                            backend.addArrayValues(
                                "band_occupations", occ[:, isegment, :, :])
                            backend.addArrayValues(
                                "band_k_points", kpt[isegment])
                            # "band_segm_labels"
                            backend.addArrayValues("band_segm_start_end", np.asarray(
                                [kpt[isegment, 0], kpt[isegment, divisions - 1]]))
                            backend.closeNonOverlappingSection(
                                "section_k_band_segment")
                        backend.closeNonOverlappingSection("section_k_band")
                        backend.openNonOverlappingSection(
                            "section_k_band_normalized")
                        for isegment in range(nsegments):
                            backend.openNonOverlappingSection(
                                "section_k_band_segment_normalized")
                            backend.addArrayValues(
                                "band_energies_normalized", energies[:, isegment, :, :] - max(self.vbTopE))
                            backend.addArrayValues(
                                "band_occupations_normalized", occ[:, isegment, :, :])
                            backend.addArrayValues(
                                "band_k_points_normalized", kpt[isegment])
                            backend.addArrayValues("band_segm_start_end_normalized", np.asarray(
                                [kpt[isegment, 0], kpt[isegment, divisions - 1]]))
                            backend.closeNonOverlappingSection(
                                "section_k_band_segment_normalized")

                        backend.closeNonOverlappingSection(
                            "section_k_band_normalized")
                    else:
                        backend.openNonOverlappingSection(
                            "section_eigenvalues")
                        backend.addArrayValues("eigenvalues_values", ev)
                        backend.addArrayValues(
                            "eigenvalues_occupation", occupation)
                        backend.closeNonOverlappingSection(
                            "section_eigenvalues")
            else:
                backend.pwarn("unexpected tag %s in the eigenvalues" % el.tag)

    def onEnd_scstep(self, parser, event, element, pathStr):
        pass

    def onStart_calculation(self, parser, event, element, pathStr):
        backend = parser.backend
        gIndexes = parser.tagSections[pathStr]
        self.singleConfCalcs.append(
            gIndexes["section_single_configuration_calculation"])
        if self.waveCut:
            backend.openNonOverlappingSection("section_basis_set")
            backend.addValue(
                "mapping_section_basis_set_cell_dependent", self.waveCut)
            backend.closeNonOverlappingSection("section_basis_set")

    def onEnd_modeling(self, parser, event, element, pathStr):
        backend = parser.backend
        backend.addValue("x_vasp_unknown_incars", self.unknown_incars)
        if self.ibrion is None or self.ibrion == -1:
            return
        samplingGIndex = backend.openSection("section_sampling_method")
        if self.ibrion == 0:
            sampling_method = "molecular_dynamics"
        else:
            sampling_method = "geometry_optimization"
        backend.addValue("sampling_method", sampling_method)
        backend.closeSection("section_sampling_method", samplingGIndex)
        frameSequenceGIndex = backend.openSection("section_frame_sequence")
        backend.addValue("frame_sequence_to_sampling_ref", samplingGIndex)
        backend.addArrayValues(
            "frame_sequence_local_frames_ref", np.asarray(self.singleConfCalcs))
        backend.closeSection("section_frame_sequence", frameSequenceGIndex)


    def onEnd_calculation(self, parser, event, element, pathStr):
        eConv = eV2J
        fConv = convert_unit_function("eV/angstrom", "N")
        pConv = convert_unit_function("eV/angstrom^3", "Pa")
        backend = parser.backend
        backend.addValue(
            "single_configuration_calculation_to_system_ref", self.lastSystemDescription)
        gIndexes = parser.tagSections["/modeling"]
        backend.addValue(
            "single_configuration_to_calculation_method_ref", gIndexes["section_method"])
        for el in element:
            if el.tag == "energy":
                for enEl in el:
                    if enEl.tag == "i":
                        name = enEl.attrib.get("name", None)
                        if name == "e_fr_energy":
                            value = eConv(float(enEl.text.strip()))
                            backend.addValue("energy_free", value)
                        elif name == "e_wo_entrp":
                            value = eConv(float(enEl.text.strip()))
                            backend.addValue("energy_total", value)
                        elif name == "e_0_energy":
                            value = eConv(float(enEl.text.strip()))
                            backend.addValue("energy_total_T0", value)
                        else:
                            backend.pwarn(
                                "Unexpected i tag with name %s in energy section" % name)
                    elif enEl.tag == "varray":
                        name = enEl.attrib.get("name", None)
                        if name == "forces":
                            f = getVector(enEl, lambda x: fConv(float(x)))
                            backend.addValue("atom_forces", f)
                        elif name == 'stress':
                            f = getVector(enEl, lambda x: pConv(float(x)))
                            backend.addValue("stress_tensor", f)

    def onEnd_atominfo(self, parser, event, element, pathStr):
        nAtoms = None
        nAtomTypes = None
        atomTypes = []
        labels = []
        labels2 = None
        atomTypesDesc = []
        backend = parser.backend
        for el in element:
            if el.tag == "atoms":
                nAtoms = int(el.text.strip())
            elif el.tag == "types":
                nAtomTypes = int(el.text.strip())
            elif el.tag == "array":
                name = el.attrib.get("name", None)
                if name == "atoms":
                    for atomsEl in el:
                        if atomsEl.tag == "dimension":
                            pass
                        elif atomsEl.tag == "field":
                            pass
                        elif atomsEl.tag == "set":
                            for atomsLine in atomsEl:
                                if atomsLine.tag != "rc":
                                    backend.pwarn(
                                        "unexpected tag %s in atoms array in atominfo" % atomsLine.tag)
                                else:
                                    line = atomsLine.findall("c")
                                    labels.append(line[0].text.strip())
                                    atomTypes.append(int(line[1].text.strip()))
                        else:
                            backend.pwarn(
                                "unexpected tag %s in atoms array in atominfo" % atomsEl.tag)
                elif name == "atomtypes":
                    keys = []
                    fieldTypes = []
                    for atomsEl in el:
                        if atomsEl.tag == "dimension":
                            pass
                        elif atomsEl.tag == "field":
                            keys.append(atomsEl.text.strip())
                            fieldTypes.append(
                                atomsEl.attrib.get("type", "float"))
                        elif atomsEl.tag == "set":
                            expectedKeys = ["atomspertype", "element",
                                            "mass", "valence", "pseudopotential"]
                            if keys != expectedKeys:
                                backend.pwarn(
                                    "unexpected fields in atomtype: %s vs %s" % (keys, expectedKeys))
                            for atomsLine in atomsEl:
                                if atomsLine.tag != "rc":
                                    backend.pwarn(
                                        "unexpected tag %s in atoms array in atominfo" % atomsLine.tag)
                                else:
                                    line = atomsLine.findall("c")
                                    typeDesc = {}
                                    for i, k in enumerate(keys):
                                        fieldType = fieldTypes[i]
                                        value = line[i].text
                                        if fieldType == "float":
                                            value = float(value)
                                        elif fieldType == "int":
                                            value = int(value)
                                        else:
                                            pass
                                        typeDesc[k] = value
                                    atomTypesDesc.append(typeDesc)
                        else:
                            backend.pwarn(
                                "unexpected tag %s in atomtypes array in atominfo" % atomsEl.tag)
                    kindIds = []
                    nEl = {}
                    kindLabels = []
                    for atomDesc in atomTypesDesc:
                        kindId = backend.openSection(
                            "section_method_atom_kind")
                        if 'element' in atomDesc:
                            elName = atomDesc['element'].strip()
                            try:
                                elNr = ase.data.chemical_symbols.index(elName)
                                backend.addValue(
                                    "method_atom_kind_atom_number", elNr)
                            except Exception as e:
                                self.logger.error(
                                    "error finding element number for %r" % atomDesc['element'].strip(),
                                    exc_info=e)
                            nElNow = 1 + nEl.get(elName, 0)
                            nEl[elName] = nElNow
                            elLabel = elName + \
                                (str(nElNow) if nElNow > 1 else "")
                            kindLabels.append(elLabel)
                            backend.addValue("method_atom_kind_label", elLabel)
                            if "mass" in atomDesc:
                                backend.addValue(
                                    "method_atom_kind_mass", atomDesc["mass"])
                            if "valence" in atomDesc:
                                backend.addValue(
                                    "method_atom_kind_explicit_electrons", atomDesc["valence"])
                            if "pseudopotential" in atomDesc:
                                backend.addValue(
                                    "method_atom_kind_pseudopotential_name", atomDesc["pseudopotential"])
                        kindIds.append(kindId)
                        backend.closeSection(
                            "section_method_atom_kind", kindId)
                    backend.addArrayValues("x_vasp_atom_kind_refs", np.asarray(
                        [kindIds[i-1] for i in atomTypes]))
                    labels2 = [kindLabels[i-1] for i in atomTypes]
                else:
                    backend.pwarn(
                        "unexpected array named %s in atominfo" % name)
            else:
                backend.pwarn("unexpected tag %s in atominfo" % el.tag)
        self.labels = np.asarray(labels2) if labels2 else np.asarray(labels)

    def incarOutTag(self, el):
        backend = self.parser.backend
        metaEnv = self.parser.backend.metaInfoEnv()
        if (el.tag != "i"):
            backend.pwarn("unexpected tag %s %s %r in incar" %
                          (el.tag, el.attrib, el.text))
        else:
            name    = el.attrib.get("name", None)
            valType = el.attrib.get("type")
            meta = metaEnv['x_vasp_incarOut_' + name]


            if not meta:
                # Unknown_Incars_Begin: storage into a dictionary
                if not valType:
                    # On vasp's xml files, valType *could* be absent if incar value is float
                    valType = 'float'

                # map vasp's datatype to nomad's datatype [b, f, i, C, D, R]
                nomad_dtypeStr = vasp_to_metainfo_type_mapping[valType][0]

                converter = metaTypeTransformers.get(nomad_dtypeStr)
                text_value = el.text.strip() # text representation of incar value
                try:
                    pyvalue = converter(text_value) # python data type
                except Exception:
                    pyvalue = text_value

                # save (name, pyvalue) into a dict
                self.unknown_incars[name] = pyvalue
                # Unknown_Incars_end
            else:
                if not valType:
                    valType = 'float'

                vasp_metainfo_type = vasp_to_metainfo_type_mapping.get(valType)[0]
                metainfo_type = meta.get('dtypeStr')
                if not vasp_metainfo_type:
                    backend.pwarn("Unknown value type %s encountered in INCAR out: %s %s %r" % (
                        valType, el.tag, el.attrib, el.text))

                elif metainfo_type != vasp_metainfo_type:
                    if  (metainfo_type == 'C' and vasp_metainfo_type == 'b'):
                        pass
                    elif  (metainfo_type == 'i' and vasp_metainfo_type == 'f'):
                        pass
                    else:
                        backend.pwarn("Data type mismatch: %s. Vasp_type: %s, metainfo_type: %s " %
                        (name, vasp_metainfo_type, metainfo_type))
                try:
                    shape = meta.get("shape", None)
                    converter = metaTypeTransformers.get(metainfo_type)
                    if not converter:
                        backend.pwarn(
                            "could not find converter for dtypeStr %s when handling meta info %s" %
                            (metainfo_type, meta ))
                    elif shape:
                        vals = re.split("\s+", el.text.strip())
                        backend.addValue(
                            meta["name"], [converter(x) for x in vals])
                    else:
                        # If-block to handle incars without value
                        if el.text == None:
                            el.text = ''
                        backend.addValue(meta["name"], converter(el.text))

                except:
                    backend.pwarn("Exception trying to handle incarOut %s: %s" % (
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
                    fMap = {
                        '91': ['GGA_X_PW91', 'GGA_C_PW91'],
                        'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
                        'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
                        'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
                        'MK': ['GGA_X_OPTB86_VDW'],
                        '--': ['GGA_X_PBE', 'GGA_C_PBE']  # should check potcar
                    }
                    functs = fMap.get(el.text.strip(), None)
                    if not functs:
                        backend.pwarn("Unknown XC functional %s" %
                                      el.text.strip())
                    else:
                        for f in functs:
                            backend.openNonOverlappingSection(
                                "section_XC_functionals")
                            backend.addValue("XC_functional_name", f)
                            backend.closeNonOverlappingSection(
                                "section_XC_functionals")
                elif name == "ISPIN":
                    self.ispin = int(el.text.strip())


    def separatorScan(self, element, backend, depth=0):
        for separators in element:
            if separators.tag == "separator":
                separatorName = separators.attrib.get("name")
                for el in separators:
                    if el.tag == "i":
                        self.incarOutTag(el)
                    elif el.tag == "separator":
                        self.separatorScan(el, backend, depth + 1)
                    else:
                        # backend.pwarn("unexpected tag %s %s in parameters separator %s at depth %d" % (
                        #     el.tag, el.attrib, separatorName, depth))
                        pass
            elif separators.tag == "i":
                self.incarOutTag(separators)
            else:
                # backend.pwarn("unexpected tag %s %s in parameters at depth %d" % (
                #     separators.tag, separators.attrib, depth))
                pass

    def onEnd_parameters(self, parser, event, element, pathStr):
        self.separatorScan(element, parser.backend)
        backend = parser.backend
        try:
            self.prec
            try:
                self.enmax
                self.waveCut = backend.openNonOverlappingSection(
                    "section_basis_set_cell_dependent")
                backend.addValue("basis_set_planewave_cutoff",
                                 eV2J(self.enmax*self.prec))
                backend.closeNonOverlappingSection(
                    "section_basis_set_cell_dependent")
                backend.openNonOverlappingSection("section_method_basis_set")
                backend.addValue(
                    "mapping_section_method_basis_set_cell_associated", self.waveCut)
                backend.closeNonOverlappingSection("section_method_basis_set")
            except AttributeError:
                import traceback
                traceback.print_exc()
                backend.pwarn(
                    "Missing ENMAX for calculating plane wave basis cut off ")
        except AttributeError:
            backend.pwarn(
                "Missing PREC for calculating plane wave basis cut off ")

    def onEnd_dos(self, parser, event, element, pathStr):
        "density of states"
        backend = parser.backend
        backend.openNonOverlappingSection("section_dos")
        for el in element:
            if el.tag == "i":
                if el.attrib.get("name") == "efermi":
                    self.eFermi = eV2J(float(el.text.strip()))
                    backend.addValue("dos_fermi_energy", self.eFermi)
                    backend.addArrayValues(
                        "energy_reference_fermi", np.array([self.eFermi]*self.ispin))
                else:
                    backend.pwarn("unexpected tag %s %s in dos" %
                                  (el.tag, el.attrib))
            elif el.tag == "total":
                for el1 in el:
                    if el1.tag == "array":
                        for el2 in el1:
                            if el2.tag == "dimension" or el2.tag == "field":
                                pass
                            elif el2.tag == "set":
                                dosL = []
                                for spinComponent in el2:
                                    if spinComponent.tag == "set":
                                        dosL.append(
                                            getVector(spinComponent, field="r"))
                                    else:
                                        backend.pwarn("unexpected tag %s %s in dos total array set" % (
                                            spinComponent.tag, spinComponent.attrib))
                                dosA = np.asarray(dosL)
                                if len(dosA.shape) != 3:
                                    raise Exception("unexpected shape %s (%s) for total dos (ragged arrays?)" % (
                                        dosA.shape), dosA.dtype)
                                dosE = eV2JV(dosA[0, :, 0])
                                dosI = dosA[:, :, 2]
                                dosV = dosA[:, :, 1]

                                # Convert the DOS values to SI. VASP uses the
                                # following units in the output:
                                # states/eV/cell. This means that the volume
                                # dependence has been introduced by multiplying
                                # by the cell volume
                                # the integrated dos value is the number of electrons until that energy level
                                # and thus not directly energy dependent anymore
                                joule_in_ev = convert_unit(1, "eV", "J")
                                dosV = dosV / joule_in_ev

                                if self.vbTopE:
                                    eRef = max(self.vbTopE)
                                else:
                                    eRef = self.eFermi
                                backend.addArrayValues("dos_energies", dosE)
                                backend.addArrayValues(
                                    "dos_energies_normalized", dosE - eRef)
                                backend.addArrayValues("dos_values", dosV)
                                backend.addArrayValues(
                                    "dos_integrated_values", dosI)
                            else:
                                backend.pwarn("unexpected tag %s %s in dos total array" % (
                                    el2.tag, el2.attrib))
                    else:
                        backend.pwarn("unexpected tag %s %s in dos total" % (
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
                                        for spinComponent in atom:
                                            if spinComponent.tag == "set":
                                                atomL.append(
                                                    getVector(spinComponent, field="r"))
                                            else:
                                                backend.pwarn("unexpected tag %s %s in dos partial array set set" % (
                                                    spinComponent.tag, spinComponent.attrib))
                                    else:
                                        backend.pwarn("unexpected tag %s %s in dos partial array set" % (
                                            spinComponent.tag, spinComponent.attrib))
                                dosLM = np.asarray(dosL)
                                assert len(
                                    dosLM.shape) == 4, "invalid shape dimension in projected dos (ragged arrays?)"
                                backend.addArrayValues(
                                    "dos_values_lm", dosLM[:, :, :, 1:])
                            else:
                                backend.pwarn("unexpected tag %s %s in dos total array" % (
                                    el2.tag, el2.attrib))
                        backend.addArrayValues("dos_lm", np.asarray(lm))
                        backend.addValue("dos_m_kind", "polynomial")
                    else:
                        backend.pwarn("unexpected tag %s %s in dos total" % (
                            el2.tag, el2.attrib))
            else:
                backend.pwarn("unexpected tag %s %s in dos" %
                              (el2.tag, el2.attrib))
        backend.closeNonOverlappingSection("section_dos")

    def onEnd_projected(self, parser, event, element, pathStr):
        "projected eigenvalues"
        return None


class XmlParser(object):
    @staticmethod
    def extractCallbacks(obj):
        """extracts all callbacks from the object obj

        triggers should start with onStart_ or onEnd__ and then have a valid section name.
        They will be called with this object, the event and current element
        """
        triggers = {}
        for attr in dir(obj):
            if attr.startswith("onStart_"):
                triggers[attr] = getattr(obj, attr)
            elif attr.startswith("onEnd_"):
                triggers[attr] = getattr(obj, attr)
        return triggers

    @staticmethod
    def maybeGet(el, path, default=None):
        i = el.findall(path)
        if i:
            return i.pop().text
        else:
            return default

    def __init__(self, parserInfo, superContext, callbacks=None, sectionMap=None):
        self.fIn = None
        self.parserInfo = parserInfo
        self.superContext = superContext
        self.callbacks = callbacks if callbacks is not None else XmlParser.extractCallbacks(
            superContext)
        self.sectionMap = sectionMap if sectionMap is not None else superContext.sectionMap
        self.path = []
        self.tagSections = {}

    def parse(self, mainFileUri, fIn, backend):
        if self.path:
            raise Exception(
                "Parse of %s called with non empty path, parse already in progress?" % mainFileUri)
        self.mainFileUri = mainFileUri
        self.fIn = fIn
        self.backend = backend
        backend.startedParsingSession(
            mainFileUri=mainFileUri,
            parserInfo=self.parserInfo)
        self.superContext.startedParsing(self)
        # there are invalid characters like esc in the files, we do not want to crash on them
        xmlParser = MyXMLParser()
        try:
            for event, el in xml.etree.ElementTree.iterparse(self.fIn, events=["start", "end"], parser=xmlParser):
                if event == 'start':
                    pathStr = "/".join([x.tag for x in self.path]
                                       ) + "/" + el.tag
                    sectionsToOpen = self.sectionMap.get(el.tag, None)
                    if sectionsToOpen:
                        gIndexes = {}
                        for sect in sectionsToOpen:
                            gIndexes[sect] = backend.openSection(sect)
                        self.tagSections[pathStr] = gIndexes
                    callback = self.callbacks.get("onStart_" + el.tag, None)
                    if callback:
                        callback(self, event, el, pathStr)
                    self.path.append(el)
                elif event == 'end':
                    lastEl = self.path.pop()
                    if lastEl != el:
                        raise Exception(
                            "mismatched path at end, got %s expected %s" % (lastEl, el))
                    tag = el.tag
                    pathStr = "/".join([x.tag for x in self.path]) + "/" + tag
                    callback = self.callbacks.get("onEnd_" + tag, None)
                    if callback:
                        if not callback(self, event, el, pathStr):
                            # if callback does not return True then assume that the current element has been processed
                            # and can be removed
                            el.clear()
                            if self.path:
                                self.path[-1].remove(el)
                    elif len(self.path) == 1:
                        self.backend.pwarn("Skipping level 1 tag %s" % tag)
                        el.clear()
                        self.path[-1].remove(el)
                    sectionsToClose = self.sectionMap.get(tag, None)
                    if sectionsToClose:
                        gIndexes = self.tagSections[pathStr]
                        del self.tagSections[pathStr]
                        for sect in reversed(sectionsToClose):
                            self.backend.closeSection(sect, gIndexes[sect])
                        self.tagSections[pathStr] = gIndexes
                else:
                    raise Exception("Unexpected event %s" % event)
        except ParseError as e:
            self.superContext.logger.warn("Could not complete parsing: %s" % e, exc_info=e)
            backend.finishedParsingSession(
                parserStatus="ParseSuccess",
                parserErrors=["exception: %s" % e]
            )
        except Exception as e:
            import traceback
            traceback.print_exc()
            backend.finishedParsingSession(
                parserStatus="ParseFailure",
                parserErrors=["exception: %s" % e]
            )
        else:
            backend.finishedParsingSession(
                parserStatus="ParseSuccess",
                parserErrors=None
            )


g = XmlParser.maybeGet

parserInfo = {
    "name": "parser_vasprun",
    "version": "1.0"
}
