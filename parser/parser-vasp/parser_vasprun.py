import xml.etree.ElementTree
import logging, sys
import setup_paths
from datetime import datetime
import os, logging, re
from nomadcore.utils import goInteractive
from nomadcore.parser_backend import JsonParseEventsWriterBackend
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import numpy as np
import setup_paths
from nomadcore.unit_conversion.unit_conversion import convert_unit_function

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

def getVector(el, transform = float, field = "v"):
    """ returns the vasp style vector contained in the element el (using field v).
    single elements are converted using the function convert"""
    return map(lambda x: map(transform, re.split(r"\s+", x.text.strip())), el.findall(field))

class VasprunContext(object):
    def __init__(self):
        self.parser = None
        self.bands = None
        self.kpoints = None
        self.weights = None
        self.ispin = None
        self.lastSystemDescription = None

    sectionMap = {
        "modeling": ["section_run", "section_method"],
        "structure": ["section_system_description"],
        "calculation": ["single_configuration_calculation"]
    }

    def startedParsing(self, parser):
        self.parser = parser

    def onEnd_generator(self, parser, event, element):
        backend = parser.backend
        program_name = g(element, "i/[@name='name']")
        if program_name:
            backend.addValue("program_name", program_name)
        version = (g(element, "i/[@name='version']", "") + " " +
                   g(element, "i/[@name='subversion']", "") + " " +
                   g(element, "i/[@name='platform']", ""))
        if not version.isspace():
            backend.addValue("program_version", version)
        date = g(element, "i/[@name='date']")
        pdate = None
        time = g(element, "i/[@name='time']")
        if date:
            pdate = datetime.strptime(date, "%Y %m %d")
        if pdate and time:
            pdate = datetime.combine(pdate.date(), datetime.strptime(time, "%H:%M:%S").timetz())
        if pdate:
            backend.addValue("program_compilation_datetime", secondsFromEpoch(pdate))
        for i in element:
            if i.tag != "i" or not i.attrib.get("name") in set(["name", "version", "subversion", "platform", "program_version", "date", "time"]):
                backend.pwarn("unexpected tag %s %s %r in generator" % (i.tag, i.attrib, i.text))

    def onEnd_incar(self, parser, event, element):
        backend = parser.backend
        metaEnv = parser.backend.metaInfoEnv()
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
                            backend.addValue(meta["name"], map(converter, vals))
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


    def onEnd_kpoints(self, parser, event, element):
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
                elif name == "weights":
                    self.weights = np.asarray(getVector(el))
                else:
                    backend.pwarn("Unknown array %s in kpoints" % name)
            else:
                backend.pwarn("Unknown tag %s in kpoints" % el.tag)


    def onEnd_structure(self, parser, event, element):
        backend = parser.backend
        gIndexes = parser.tagSections["/".join(map(lambda x: x.tag, parser.path[:])) + "/" + element.tag]
        self.lastSystemDescription = gIndexes["section_system_description"]
        for el in element:
            if (el.tag == "crystal"):
                cell = None
                for cellEl in el:
                    if cellEl.tag == "varray":
                        name = cellEl.attrib.get("name", None)
                        if name == "basis":
                            conv = convert_unit_function("angstrom","m")
                            cell = getVector(cellEl, lambda x: conv(float(x)))
                        elif name =="rec_basis":
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
                    backend.addArrayValues("atom_position", np.asarray(pos))
                else:
                    backend.pwarn("Unexpected varray in structure %s" % el.attrib)
            else:
                backend.pwarn("Unexpected tag in structure %s %s %r" % el.tag, el.attrib, el.text)

    def onEnd_eigenvalues(self, parser, event, element):
        backend = parser.backend
        eConv = convert_unit_function("eV", "J")
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
                        ik = -1
                        isp = -1
                        for spinEl in arrEl:
                            if spinEl.tag == "set":
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
                    ev = map(eConv, eigenvalues)
                    if self.bands:
                        divisions = int(self.bands['divisions'])
                        backend.openNonOverlappingSection("section_k_band")
                        backend.addArrayValues("band_energies", np.reshape(ev, (self.kpoints.shape[0]/divisions, self.ispin, divisions ,  bands.shape[0])))
                        backend.closeNonOverlappingSection("section_k_band")
                    else:
                        backend.openNonOverlappingSection("section_eigenvalues_group")
                        for isp in range(self.ispin):
                            backend.openNonOverlappingSection("section_eigenvalues")
                            backend.addArrayValues("eigenvalues_eigenvalues", ev[isp])
                            backend.closeNonOverlappingSection("section_eigenvalues")
                        backend.closeNonOverlappingSection("section_eigenvalues_group")

            else:
                backend.pwarn("unexpected tag %s in the eigenvalues" % el.tag)

    def onEnd_scstep(self, parser, event, element):
        pass

    def onEnd_calculation(self, parser, event, element):
        eConv = convert_unit_function("eV", "J")
        fConv = convert_unit_function("eV/angstrom", "N")
        pConv = convert_unit_function("eV/angstrom^3", "Pa")
        backend = parser.backend
        backend.addValue("single_configuration_calculation_to_system_description_ref", self.lastSystemDescription)
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
                        elif name == "e_wo_entr":
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
        try:
            for event, el in xml.etree.ElementTree.iterparse(self.fIn, events=["start","end"]):
                if event == 'start':
                    sectionsToOpen = self.sectionMap.get(el.tag, None)
                    if sectionsToOpen:
                        pathStr = "/".join(map(lambda x: x.tag, self.path)) + "/" + el.tag
                        gIndexes = {}
                        for sect in sectionsToOpen:
                            gIndexes[sect] = backend.openSection(sect)
                        self.tagSections[pathStr] = gIndexes
                    callback = self.callbacks.get("onStart_" + el.tag, None)
                    if callback:
                        callback(self, event, el)
                    self.path.append(el)
                elif event == 'end':
                    lastEl = self.path.pop()
                    if lastEl != el:
                        raise Exception("mismatched path at end, got %s expected %s" % (lastEl, el))
                    tag = el.tag
                    callback = self.callbacks.get("onEnd_" + tag, None)
                    if callback:
                        if not callback(self, event, el):
                            el.clear()
                            if self.path:
                                self.path[-1].remove(el)
                    elif len(self.path) == 1:
                        self.backend.pwarn("Skipping level 1 tag %s" % tag)
                        el.clear()
                        self.path[-1].remove(el)
                    sectionsToClose = self.sectionMap.get(tag, None)
                    if sectionsToClose:
                        pathStr = "/".join(map(lambda x: x.tag, self.path)) + "/" + tag
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
                parserErrors = ["exception: %s" % sys.exc_value]
            )
        else:
            backend.finishedParsingSession(
                parserStatus = "ParseSuccess",
                parserErrors = []
            )

g = XmlParser.maybeGet


metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/vasp.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)

parserInfo = {
  "name": "parser_vasprun",
  "version": "1.0"
}

if __name__ == "__main__":
    superContext =  VasprunContext()
    parser = XmlParser(parserInfo, superContext)
    backend = JsonParseEventsWriterBackend(metaInfoEnv, sys.stdout)
    parser.parse(sys.argv[1], open(sys.argv[2]), backend)
