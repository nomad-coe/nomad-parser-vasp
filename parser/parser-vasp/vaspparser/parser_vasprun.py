# coding=utf-8
from __future__ import division
from builtins import map
from builtins import range
from builtins import object
import xml.etree.ElementTree
import logging, sys, bisect
import setup_paths
from datetime import datetime
import os, logging, re, traceback
from nomadcore.parser_backend import JsonParseEventsWriterBackend
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import numpy as np
import setup_paths
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
from nomadcore.unit_conversion.unit_conversion import convert_unit
import ase.geometry
import ase.data
from math import pi

eV2J = convert_unit_function("eV", "J")
eV2JV = np.vectorize(eV2J)

def crystal_structure_from_cell(cell, eps=1e-4):
    """Return the crystal structure as a string calculated from the cell.
    """
    cellpar = ase.geometry.cell_to_cellpar(cell=cell)
    abc = cellpar[:3]
    angles = cellpar[3:] / 180 * pi
    a, b, c = abc
    alpha, beta, gamma = angles
    if abc.ptp() < eps and abs(angles - pi / 2).max() < eps:
        return 'cubic'
    elif abc.ptp() < eps and abs(angles - pi / 3).max() < eps:
        return 'fcc'
    elif abc.ptp() < eps and abs(angles - np.arccos(-1 / 3)).max() < eps:
        return 'bcc'
    elif abs(a - b) < eps and abs(angles - pi / 2).max() < eps:
        return 'tetragonal'
    elif abs(angles - pi / 2).max() < eps:
        return 'orthorhombic'
    elif (abs(a - b) < eps and
          abs(gamma - pi / 3 * 2) < eps and
          abs(angles[:2] - pi / 2).max() < eps):
        return 'hexagonal'
    elif (c >= a and c >= b and alpha < pi / 2 and
          abs(angles[1:] - pi / 2).max() < eps):
        return 'monoclinic'
    else:
       raise ValueError('Cannot find crystal structure')


special_points = {
    'cubic': {'Γ': [0, 0, 0],
              'M': [1 / 2, 1 / 2, 0],
              'R': [1 / 2, 1 / 2, 1 / 2],
              'X': [0, 1 / 2, 0]},
    'fcc': {'Γ': [0, 0, 0],
            'K': [3 / 8, 3 / 8, 3 / 4],
            'L': [1 / 2, 1 / 2, 1 / 2],
            'U': [5 / 8, 1 / 4, 5 / 8],
            'W': [1 / 2, 1 / 4, 3 / 4],
            'X': [1 / 2, 0, 1 / 2]},
    'bcc': {'Γ': [0, 0, 0],
            'H': [1 / 2, -1 / 2, 1 / 2],
            'P': [1 / 4, 1 / 4, 1 / 4],
            'N': [0, 0, 1 / 2]},
    'tetragonal': {'Γ': [0, 0, 0],
                   'A': [1 / 2, 1 / 2, 1 / 2],
                   'M': [1 / 2, 1 / 2, 0],
                   'R': [0, 1 / 2, 1 / 2],
                   'X': [0, 1 / 2, 0],
                   'Z': [0, 0, 1 / 2]},
    'orthorhombic': {'Γ': [0, 0, 0],
                     'R': [1 / 2, 1 / 2, 1 / 2],
                     'S': [1 / 2, 1 / 2, 0],
                     'T': [0, 1 / 2, 1 / 2],
                     'U': [1 / 2, 0, 1 / 2],
                     'X': [1 / 2, 0, 0],
                     'Y': [0, 1 / 2, 0],
                     'Z': [0, 0, 1 / 2]},
    'hexagonal': {'Γ': [0, 0, 0],
                  'A': [0, 0, 1 / 2],
                  'H': [1 / 3, 1 / 3, 1 / 2],
                  'K': [1 / 3, 1 / 3, 0],
                  'L': [1 / 2, 0, 1 / 2],
                  'M': [1 / 2, 0, 0]}}


special_paths = {
    'cubic': 'ΓXMΓRX,MR',
    'fcc': 'ΓXWKΓLUWLK,UX',
    'bcc': 'ΓHNΓPH,PN',
    'tetragonal': 'ΓXMΓZRAZXR,MA',
    'orthorhombic': 'ΓXSYΓZURTZ,YT,UX,SR',
    'hexagonal': 'ΓMKΓALHA,LM,KH',
    'monoclinic': 'ΓYHCEM1AXH1,MDZ,YD'}


def get_special_points(cell, eps=1e-4):
    """Return dict of special points.

    The definitions are from a paper by Wahyu Setyawana and Stefano
    Curtarolo::

        http://dx.doi.org/10.1016/j.commatsci.2010.05.010

    lattice: str
        One of the following: cubic, fcc, bcc, orthorhombic, tetragonal,
        hexagonal or monoclinic.
    cell: 3x3 ndarray
        Unit cell.
    eps: float
        Tolerance for cell-check.
    """

    lattice = crystal_structure_from_cell(cell)

    cellpar = ase.geometry.cell_to_cellpar(cell=cell)
    abc = cellpar[:3]
    angles = cellpar[3:] / 180 * pi
    a, b, c = abc
    alpha, beta, gamma = angles

    # Check that the unit-cells are as in the Setyawana-Curtarolo paper:
    if lattice == 'cubic':
        assert abc.ptp() < eps and abs(angles - pi / 2).max() < eps
    elif lattice == 'fcc':
        assert abc.ptp() < eps and abs(angles - pi / 3).max() < eps
    elif lattice == 'bcc':
        angle = np.arccos(-1 / 3)
        assert abc.ptp() < eps and abs(angles - angle).max() < eps
    elif lattice == 'tetragonal':
        assert abs(a - b) < eps and abs(angles - pi / 2).max() < eps
    elif lattice == 'orthorhombic':
        assert abs(angles - pi / 2).max() < eps
    elif lattice == 'hexagonal':
        assert abs(a - b) < eps
        assert abs(gamma - pi / 3 * 2) < eps
        assert abs(angles[:2] - pi / 2).max() < eps
    elif lattice == 'monoclinic':
        assert c >= a and c >= b
        assert alpha < pi / 2 and abs(angles[1:] - pi / 2).max() < eps
    if lattice != 'monoclinic':
        return special_points[lattice]

    # Here, we need the cell:
    eta = (1 - b * cos(alpha) / c) / (2 * sin(alpha)**2)
    nu = 1 / 2 - eta * c * cos(alpha) / b
    return {'Γ': [0, 0, 0],
            'A': [1 / 2, 1 / 2, 0],
            'C': [0, 1 / 2, 1 / 2],
            'D': [1 / 2, 0, 1 / 2],
            'D1': [1 / 2, 0, -1 / 2],
            'E': [1 / 2, 1 / 2, 1 / 2],
            'H': [0, eta, 1 - nu],
            'H1': [0, 1 - eta, nu],
            'H2': [0, eta, -nu],
            'M': [1 / 2, eta, 1 - nu],
            'M1': [1 / 2, 1 - eta, nu],
            'M2': [1 / 2, eta, -nu],
            'X': [0, 1 / 2, 0],
            'Y': [0, 0, 1 / 2],
            'Y1': [0, 0, -1 / 2],
            'Z': [1 / 2, 0, 0]}

def findLabel(labels, value):
    for k, v in labels.items():
        if np.all(np.abs(v-value) < 1.e-5):
            return k
    return "?"

def secondsFromEpoch(date):
    epoch = datetime(1970,1,1)
    ts=date-epoch
    return ts.seconds + ts.microseconds/1000.0

trueRe = re.compile(r"\s*(?:\.?[Tt](?:[Rr][Uu][Ee])?\.?|1|[Yy](?:[Ee][Ss])?|[Jj][Aa]?)\s*$")
falseRe = re.compile(r"\s*(?:\.?[fF](?:[Aa][Ll][Ss][Ee])?\.?|0|[Nn](?:[Oo]|[Ee][Ii][Nn])?)\s*$")
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
    'i': lambda x: int(x.strip()),
    'f': lambda x: float(x.strip()),
    'b': toBool,
}


import xml.etree.ElementTree as ET

class MyXMLParser(ET.XMLParser):

    rx = re.compile("&#([0-9]+);|&#x([0-9a-fA-F]+);")

    def feed(self,data):
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
        super(MyXMLParser,self).feed(mydata)
def transform2(y):
  if '**' in y: return float('nan')
  else: return y

def getVector(el, transform = float, field = "v"):
    """ returns the vasp style vector contained in the element el (using field v).
    single elements are converted using the function convert"""
#
#    for x in el.findall(field):
#        for y in re.split(r"\s+", x.text.strip()):
    return [[transform(transform2(y)) for y in re.split(r"\s+", x.text.strip())] for x in el.findall(field)]

class VasprunContext(object):
    def __init__(self):
        self.parser = None
        self.bands = None
        self.kpoints = None
        self.weights = None
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
            pdate = datetime.combine(pdate.date(), datetime.strptime(time.strip(), "%H:%M:%S").timetz())
        if pdate:
            backend.addValue("program_compilation_datetime", secondsFromEpoch(pdate))
        for i in element:
            if i.tag != "i" or not i.attrib.get("name") in set(["program", "version", "subversion", "platform", "program_version", "date", "time"]):
                backend.pwarn("unexpected tag %s %s %r in generator" % (i.tag, i.attrib, i.text))

    def onEnd_incar(self, parser, event, element, pathStr):
        backend = parser.backend
        metaEnv = parser.backend.metaInfoEnv()
        dft_plus_u = False
        ibrion = None
        nsw = 0
        for el in element:
            if (el.tag != "i"):
                backend.pwarn("unexpected tag %s %s %r in incar" % (el.tag, el.attrib, el.text))
            else:
                name = el.attrib.get("name", None)
                meta = metaEnv['x_vasp_incar_' + name]
                valType = el.attrib.get("type")
                if not meta:
                    backend.pwarn("Unknown INCAR parameter (not registered in the meta data): %s %s %r" % (el.tag, el.attrib, el.text))
                elif valType:
                    expectedMetaType = {
                        'string': ['C'],
                        'int': ['i'],
                        'logical': ['b','C']
                    }.get(valType)
                    if not expectedMetaType:
                        backend.pwarn("Unknown value type %s encountered in INCAR: %s %s %r" % (valType, el.tag, el.attrib, el.text))
                    elif not meta.get('dtypeStr') in expectedMetaType:
                        backend.pwarn("type mismatch between meta data %s and INCAR type %s for %s %s %r" % ( meta.get('dtypeStr'), valType, el.tag, el.attrib, el.text))
                    else:
                        shape = meta.get("shape", None)
                        dtypeStr = meta.get("dtypeStr", None)
                        converter = metaTypeTransformers.get(dtypeStr)
                        if not converter:
                            backend.pwarn("could not find converter for dtypeStr %s when handling meta info %s" % (dtypeStr, ))
                        elif shape:
                            vals = re.split("\s+", el.text.strip())
                            backend.addValue(meta["name"], [converter(x) for x in vals])
                        else:
                            backend.addValue(meta["name"], converter(el.text))
                    if name == 'GGA':
                        fMap = {
                            '91': ['GGA_X_PW91', 'GGA_C_PW91'],
                            'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
                            'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
                            'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL']
                        }
                        functs = fMap.get(el.text.strip(), None)
                        if not functs:
                            backend.pwarn("Unknown XC functional %s" % el.text.strip())
                        else:
                            for f in functs:
                                backend.openNonOverlappingSection("section_XC_functionals")
                                backend.addValue("XC_functional_name", f)
                                backend.closeNonOverlappingSection("section_XC_functionals")
                    elif name == "ISPIN":
                        self.ispin = int(el.text.strip())
                    elif name == "LDAU":
                        if re.match(".?[Tt](?:[rR][uU][eE])?.?|[yY](?:[eE][sS])?|1", el.text.strip()):
                            dft_plus_u = True
                    elif name == "IBRION":
                        ibrion = int(el.text.strip())
                    elif name == "NSW":
                        nsw = int(el.text.strip())
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
                param = el.attrib.get("param", None)
                if param:
                    backend.addValue("x_vasp_k_points_generation_method", param)
                if param == "listgenerated":
                    self.bands = {
                        "divisions": g(el,"i/[@name='divisions']", None),
                        "points": getVector(el)
                    }
                elif param == "Monkhorst-Pack":
                    pass
                else:
                    backend.pwarn("unknown k point generation method")
            elif el.tag == "varray":
                name = el.attrib.get("name", None)
                if name == "kpointlist":
                    self.kpoints = np.asarray(getVector(el))
                    backend.addArrayValues("k_mesh_points", self.kpoints)
                elif name == "weights":
                    self.weights = np.asarray(getVector(el))
                    backend.addArrayValues("k_mesh_weights", self.weights.flatten())
                else:
                    backend.pwarn("Unknown array %s in kpoints" % name)
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
                            self.cell = getVector(cellEl, lambda x: conv(float(x)))
                            self.angstrom_cell = np.array(getVector(cellEl))
                            backend.addArrayValues("simulation_cell", np.asarray(self.cell))
                            backend.addArrayValues("configuration_periodic_dimensions", np.ones(3, dtype=bool))
                        elif name == "rec_basis":
                            pass
                        else:
                            backend.pwarn("Unexpected varray %s in crystal" % name)
                    elif cellEl.tag == "i":
                        if cellEl.attrib.get("name") != "volume":
                            backend.pwarn("Unexpected i value %s in crystal" % cellEl.attrib)
                    else:
                        backend.pwarn("Unexpected tag %s %s %r in crystal" % (cellEl.tag, cellEl.attrib, cellEl.text))
            elif el.tag == "varray":
                name = el.attrib.get("name", None)
                if name == "positions":
                    pos = getVector(el)
                    backend.addArrayValues("atom_positions", np.dot(np.asarray(pos), self.cell))
                else:
                    backend.pwarn("Unexpected varray in structure %s" % el.attrib)
            else:
                backend.pwarn("Unexpected tag in structure %s %s %r" % el.tag, el.attrib, el.text)
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
                                        bands = np.asarray(getVector(kEl, field = "r"))
                                        if eigenvalues is None:
                                            eigenvalues = np.zeros((self.ispin, self.kpoints.shape[0],  bands.shape[0]), dtype = float)
                                            occupation = np.zeros((self.ispin, self.kpoints.shape[0],  bands.shape[0]), dtype = float)
                                        eigenvalues[isp, ik] = bands[:,0]
                                        occupation[isp, ik] = bands[:,1]
                                    else:
                                        backend.pwarn("unexpected tag %s in k array of the eigenvalues" % kEl.tag)
                            else:
                                backend.pwarn("unexpected tag %s in spin array of the eigenvalues" % spinEl.tag)
                    else:
                        backend.pwarn("unexpected tag %s in array of the eigenvalues" % arrEl.tag)
                if eigenvalues is not None:

                    ev = eV2JV(eigenvalues)
                    vbTopE = []
                    ebMinE = []
                    for ispin in range(occupation.shape[0]):
                        vbTopE.append(float('-inf'))
                        ebMinE.append(float('inf'))
                        for ik in range(occupation.shape[1]):
                            ebIndex = bisect.bisect_right(-occupation[ispin, ik, :], -0.5) - 1
                            vbTopIndex = ebIndex -1
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
                    backend.addArrayValues("energy_reference_highest_occupied", np.array(vbTopE))
                    backend.addArrayValues("energy_reference_lowest_unoccupied", np.array(ebMinE))
                    if self.bands:
                        divisions = int(self.bands['divisions'])
                        backend.openNonOverlappingSection("section_k_band")
                        nsegments = self.kpoints.shape[0] // divisions
                        kpt = np.reshape(self.kpoints, (nsegments, divisions, 3))
                        energies = np.reshape(ev, (self.ispin, nsegments, divisions ,  bands.shape[0]))
                        occ = np.reshape(occupation, (self.ispin, nsegments, divisions, bands.shape[0]))
                        for isegment in range(nsegments):
                            backend.openNonOverlappingSection("section_k_band_segment")
                            backend.addArrayValues("band_energies", energies[:, isegment, :, :])
                            backend.addArrayValues("band_occupations", occ[:, isegment, :, :])
                            backend.addArrayValues("band_k_points", kpt[isegment])
                            # "band_segm_labels"
                            backend.addArrayValues("band_segm_start_end", np.asarray([kpt[isegment, 0], kpt[isegment, divisions - 1]]))
                            backend.closeNonOverlappingSection("section_k_band_segment")
                        backend.closeNonOverlappingSection("section_k_band")
                        backend.openNonOverlappingSection("section_k_band_normalized")
                        specialPoints = {}
                        try:
                            specialPoints = get_special_points(convert_unit_function("m","angstrom")(self.cell))
                        except:
                            logging.exception("failed to get special points")
                        for isegment in range(nsegments):
                            backend.openNonOverlappingSection("section_k_band_segment_normalized")
                            backend.addArrayValues("band_energies_normalized", energies[:, isegment, :, :] - max(self.vbTopE))
                            backend.addArrayValues("band_occupations_normalized", occ[:, isegment, :, :])
                            backend.addArrayValues("band_k_points_normalized", kpt[isegment])
                            backend.addArrayValues("band_segm_start_end_normalized", np.asarray([kpt[isegment, 0], kpt[isegment, divisions - 1]]))
                            backend.addValue("band_segm_labels_normalized",
                                             [findLabel(specialPoints, kpt[isegment, 0]),
                                              findLabel(specialPoints, kpt[isegment, divisions - 1])])
                            backend.closeNonOverlappingSection("section_k_band_segment_normalized")

                        backend.closeNonOverlappingSection("section_k_band_normalized")
                    else:
                        backend.openNonOverlappingSection("section_eigenvalues")
                        backend.addArrayValues("eigenvalues_values", ev)
                        backend.addArrayValues("eigenvalues_occupation", occupation)
                        backend.closeNonOverlappingSection("section_eigenvalues")
            else:
                backend.pwarn("unexpected tag %s in the eigenvalues" % el.tag)

    def onEnd_scstep(self, parser, event, element, pathStr):
        pass

    def onStart_calculation(self, parser, event, element, pathStr):
        gIndexes = parser.tagSections[pathStr]
        self.singleConfCalcs.append(gIndexes["section_single_configuration_calculation"])
        if self.waveCut:
            backend.openNonOverlappingSection("section_basis_set")
            backend.addValue("mapping_section_basis_set_cell_dependent", self.waveCut)
            backend.closeNonOverlappingSection("section_basis_set")

    def onEnd_modeling(self, parser, event, element, pathStr):
        backend = parser.backend
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
        backend.addArrayValues("frame_sequence_local_frames_ref", np.asarray(self.singleConfCalcs))
        backend.closeSection("section_frame_sequence", frameSequenceGIndex)

    def onEnd_calculation(self, parser, event, element, pathStr):
        eConv = eV2J
        fConv = convert_unit_function("eV/angstrom", "N")
        pConv = convert_unit_function("eV/angstrom^3", "Pa")
        backend = parser.backend
        backend.addValue("single_configuration_calculation_to_system_ref", self.lastSystemDescription)
        gIndexes = parser.tagSections["/modeling"]
        backend.addValue("single_configuration_to_calculation_method_ref", gIndexes["section_method"])
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
                            backend.pwarn("Unexpected i tag with name %s in energy section" % name)
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
                                    backend.pwarn("unexpected tag %s in atoms array in atominfo" % atomsLine.tag)
                                else:
                                    line = atomsLine.findall("c")
                                    labels.append(line[0].text.strip())
                                    atomTypes.append(int(line[1].text.strip()))
                        else:
                            backend.pwarn("unexpected tag %s in atoms array in atominfo" % atomsEl.tag)
                elif name == "atomtypes":
                    keys = []
                    fieldTypes = []
                    for atomsEl in el:
                        if atomsEl.tag == "dimension":
                            pass
                        elif atomsEl.tag == "field":
                            keys.append(atomsEl.text.strip())
                            fieldTypes.append(atomsEl.attrib.get("type", "float"))
                        elif atomsEl.tag == "set":
                            expectedKeys = ["atomspertype", "element", "mass", "valence", "pseudopotential"]
                            if keys != expectedKeys:
                                backend.pwarn("unexpected fields in atomtype: %s vs %s" % (keys, expectedKeys))
                            for atomsLine in atomsEl:
                                if atomsLine.tag != "rc":
                                    backend.pwarn("unexpected tag %s in atoms array in atominfo" % atomsLine.tag)
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
                            backend.pwarn("unexpected tag %s in atomtypes array in atominfo" % atomsEl.tag)
                    kindIds = []
                    nEl={}
                    kindLabels = []
                    for atomDesc in atomTypesDesc:
                        kindId = backend.openSection("section_method_atom_kind")
                        if 'element' in atomDesc:
                            elName = atomDesc['element'].strip()
                            try:
                                elNr = ase.data.chemical_symbols.index(elName)
                                backend.addValue("method_atom_kind_atom_number", elNr)
                            except:
                                logging.exception("error finding element number for %r",atomDesc['element'].strip())
                            nElNow = 1 + nEl.get(elName,0)
                            nEl[elName] = nElNow
                            elLabel = elName + (str(nElNow) if nElNow > 1 else  "")
                            kindLabels.append(elLabel)
                            backend.addValue("method_atom_kind_label", elLabel)
                            if "mass" in atomDesc:
                                backend.addValue("method_atom_kind_mass", atomDesc["mass"])
                            if "valence" in atomDesc:
                                backend.addValue("method_atom_kind_explicit_electrons", atomDesc["valence"])
                            if "pseudopotential" in atomDesc:
                                backend.addValue("method_atom_kind_pseudopotential_name", atomDesc["pseudopotential"])
                        kindIds.append(kindId)
                        backend.closeSection("section_method_atom_kind", kindId)
                    backend.addArrayValues("x_vasp_atom_kind_refs", np.asarray([kindIds[i-1] for i in atomTypes]))
                    labels2 = [kindLabels[i-1] for i in atomTypes]
                else:
                    backend.pwarn("unexpected array named %s in atominfo" % name)
            else:
                backend.pwarn("unexpected tag %s in atominfo" % el.tag)
        self.labels = np.asarray(labels2) if labels2 else np.asarray(labels)

    def incarOutTag(self, el):
        backend = self.parser.backend
        metaEnv = self.parser.backend.metaInfoEnv()
        if (el.tag != "i"):
            backend.pwarn("unexpected tag %s %s %r in incar" % (el.tag, el.attrib, el.text))
        else:
            name = el.attrib.get("name", None)
            meta = metaEnv['x_vasp_incarOut_' + name]
            valType = el.attrib.get("type")
            if not meta:
                backend.pwarn("Unknown INCAR out parameter (not registered in the meta data): %s %s %r" % (el.tag, el.attrib, el.text))
            else:
                if valType:
                    expectedMetaType = {
                        'string': ['C'],
                        'int': ['i'],
                        'logical': ['b','C']
                    }.get(valType)
                    if not expectedMetaType:
                        backend.pwarn("Unknown value type %s encountered in INCAR out: %s %s %r" % (valType, el.tag, el.attrib, el.text))
                    elif not meta.get('dtypeStr') in expectedMetaType:
                        backend.pwarn("type mismatch between meta data %s and INCAR type %s for %s %s %r" % ( meta.get('dtypeStr'), valType, el.tag, el.attrib, el.text))
                try:
                    shape = meta.get("shape", None)
                    dtypeStr = meta.get("dtypeStr", None)
                    converter = metaTypeTransformers.get(dtypeStr)
                    if not converter:
                        backend.pwarn("could not find converter for dtypeStr %s when handling meta info %s" % (dtypeStr, ))
                    elif shape:
                        vals = re.split("\s+", el.text.strip())
                        backend.addValue(meta["name"], [converter(x) for x in vals])
                    else:
                        backend.addValue(meta["name"], converter(el.text))
                except:
                    backend.pwarn("Exception trying to handle incarOut %s: %s" % (name, traceback.format_exc()))
                if name == 'ENMAX' or name == 'PREC':
                    if name =='ENMAX': self.enmax=converter(el.text)
                    if name =='PREC' :
                      if 'acc' in converter(el.text):
                        self.prec=1.3
                      else:
                        self.prec=1.0
                if name == 'GGA':
                    fMap = {
                        '91': ['GGA_X_PW91', 'GGA_C_PW91'],
                        'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
                        'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
                        'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
                        '--': ['GGA_X_PBE', 'GGA_C_PBE'] # should check potcar
                    }
                    functs = fMap.get(el.text.strip(), None)
                    if not functs:
                        backend.pwarn("Unknown XC functional %s" % el.text.strip())
                    else:
                        for f in functs:
                            backend.openNonOverlappingSection("section_XC_functionals")
                            backend.addValue("XC_functional_name", f)
                            backend.closeNonOverlappingSection("section_XC_functionals")
                elif name == "ISPIN":
                    self.ispin = int(el.text.strip())

    def separatorScan(self, element, backend, depth = 0):
        for separators in element:
            if separators.tag == "separator":
                separatorName = separators.attrib.get("name")
                for el in separators:
                    if el.tag == "i":
                        self.incarOutTag(el)
                    elif el.tag == "separator":
                        self.separatorScan(el, backend, depth + 1)
                    else:
                        backend.pwarn("unexpected tag %s %s in parameters separator %s at depth %d" % (el.tag, el.attrib, separatorName, depth))
            elif separators.tag == "i":
                self.incarOutTag(separators)
            else:
                backend.pwarn("unexpected tag %s %s in parameters at depth %d" % (separators.tag, separators.attrib, depth))

    def onEnd_parameters(self, parser, event, element, pathStr):
        self.separatorScan(element, parser.backend)
        backend = parser.backend
        try:
           self.prec
           try:
              self.enmax
              self.waveCut = backend.openNonOverlappingSection("section_basis_set_cell_dependent")
              backend.addValue("basis_set_planewave_cutoff", eV2J(self.enmax*self.prec))
              backend.closeNonOverlappingSection("section_basis_set_cell_dependent")
              backend.openNonOverlappingSection("section_method_basis_set")
              backend.addValue("mapping_section_method_basis_set_cell_associated", self.waveCut)
              backend.closeNonOverlappingSection("section_method_basis_set")
           except AttributeError:
              backend.pwarn("Missing ENMAX for calculating plane wave basis cut off ")
        except AttributeError:
           backend.pwarn("Missing PREC for calculating plane wave basis cut off ")

    def onEnd_dos(self, parser, event, element, pathStr):
        "density of states"
        backend = parser.backend
        backend.openNonOverlappingSection("section_dos")
        for el in element:
            if el.tag == "i":
                if el.attrib.get("name") == "efermi":
                    self.eFermi = eV2J(float(el.text.strip()))
                    backend.addValue("dos_fermi_energy", self.eFermi)
                    backend.addArrayValues("energy_reference_fermi", np.array([self.eFermi]*self.ispin))
                else:
                    backend.pwarn("unexpected tag %s %s in dos" % (el.tag, el.attrib))
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
                                        dosL.append(getVector(spinComponent, field = "r"))
                                    else:
                                        backend.pwarn("unexpected tag %s %s in dos total array set" % (spinComponent.tag, spinComponent.attrib))
                                dosA = np.asarray(dosL)
                                if len(dosA.shape) != 3:
                                    raise Exception("unexpected shape %s (%s) for total dos (ragged arrays?)" % (dosA.shape), dosA.dtype)
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
                                backend.addArrayValues("dos_energies_normalized", dosE - eRef)
                                backend.addArrayValues("dos_values", dosV)
                                backend.addArrayValues("dos_integrated_values", dosI)
                            else:
                                backend.pwarn("unexpected tag %s %s in dos total array" % (el2.tag, el2.attrib))
                    else:
                        backend.pwarn("unexpected tag %s %s in dos total" % (el2.tag, el2.attrib))
            elif el.tag == "partial":
                for el1 in el:
                    if el1.tag == "array":
                        lm=[]
                        for el2 in el1:
                            if el2.tag == "dimension":
                                pass
                            elif el2.tag == "field":
                                if el2.text.strip() == "energy":
                                    pass
                                else:
                                    strLm = {
                                        "s": [0,0],
                                        "p": [1,-1],
                                        "px":[1,0],
                                        "py":[1,1],
                                        "pz":[1,2],
                                        "d": [2,-1],
                                        "dx2":[2,0],
                                        "dxy":[2,1],
                                        "dxz":[2,2],
                                        "dy2":[2,3],
                                        "dyz":[2,4],
                                        "dz2":[2,5]
                                    }
                                    lm.append(strLm.get(el2.text.strip(), [-1,-1]))
                            elif el2.tag == "set":
                                dosL = []
                                for atom in el2:
                                    if atom.tag == "set":
                                        atomL = []
                                        dosL.append(atomL)
                                        for spinComponent in atom:
                                            if spinComponent.tag == "set":
                                                atomL.append(getVector(spinComponent, field = "r"))
                                            else:
                                                backend.pwarn("unexpected tag %s %s in dos partial array set set" % (spinComponent.tag, spinComponent.attrib))
                                    else:
                                        backend.pwarn("unexpected tag %s %s in dos partial array set" % (spinComponent.tag, spinComponent.attrib))
                                dosLM = np.asarray(dosL)
                                assert len(dosLM.shape) == 4, "invalid shape dimension in projected dos (ragged arrays?)"
                                backend.addArrayValues("dos_values_lm", dosLM[:,:,:,1:])
                            else:
                                backend.pwarn("unexpected tag %s %s in dos total array" % (el2.tag, el2.attrib))
                        backend.addArrayValues("dos_lm", np.asarray(lm))
                        backend.addValue("dos_m_kind", "polynomial")
                    else:
                        backend.pwarn("unexpected tag %s %s in dos total" % (el2.tag, el2.attrib))
            else:
                backend.pwarn("unexpected tag %s %s in dos" % (el2.tag, el2.attrib))
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
    def maybeGet(el, path, default = None):
        i = el.findall(path)
        if i:
            return i.pop().text
        else:
            return default

    def __init__(self, parserInfo, superContext, callbacks = None, sectionMap = None):
        self.fIn = None
        self.parserInfo = parserInfo
        self.superContext = superContext
        self.callbacks = callbacks if callbacks is not None else XmlParser.extractCallbacks(superContext)
        self.sectionMap = sectionMap if sectionMap is not None else superContext.sectionMap
        self.path = []
        self.tagSections = {}

    def parse(self, mainFileUri, fIn, backend):
        if self.path:
            raise Exception("Parse of %s called with non empty path, parse already in progress?" % mainFileUri)
        self.mainFileUri = mainFileUri
        self.fIn = fIn
        self.backend = backend
        backend.startedParsingSession(
            mainFileUri = mainFileUri,
            parserInfo = self.parserInfo)
        self.superContext.startedParsing(self)
        # there are invalid characters like esc in the files, we do not want to crash on them
        xmlParser = MyXMLParser()
        try:
            for event, el in xml.etree.ElementTree.iterparse(self.fIn, events=["start","end"], parser = xmlParser):
                if event == 'start':
                    pathStr = "/".join([x.tag for x in self.path]) + "/" + el.tag
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
                        raise Exception("mismatched path at end, got %s expected %s" % (lastEl, el))
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
        except:
            logging.exception("failure when parsing %s", self.mainFileUri)
            backend.finishedParsingSession(
                parserStatus = "ParseFailure",
                parserErrors = ["exception: %s" % sys.exc_info()[1]]
            )
        else:
            backend.finishedParsingSession(
                parserStatus = "ParseSuccess",
                parserErrors = None
            )

g = XmlParser.maybeGet

parserInfo = {
  "name": "parser_vasprun",
  "version": "1.0"
}

if __name__ == "__main__":
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../../nomad-meta-info/meta_info/nomad_meta_info/vasp.nomadmetainfo.json"))
    metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)
    superContext =  VasprunContext()
    parser = XmlParser(parserInfo, superContext)
    backend = JsonParseEventsWriterBackend(metaInfoEnv, sys.stdout)
    parser.parse(sys.argv[1], open(sys.argv[2]), backend)
