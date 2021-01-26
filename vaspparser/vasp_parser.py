import os
import numpy as np
import logging
import pint
from datetime import datetime
import ase
import re

from .metainfo import m_env
from nomad.parsing import FairdiParser
from nomad.parsing.file_parser.xml_parser import XMLParser
from nomad.parsing.file_parser.text_parser import TextParser, Quantity
from nomad.datamodel.metainfo.common_dft import Run, Method, XCFunctionals,\
    SingleConfigurationCalculation, ScfIteration, MethodAtomKind, System, Eigenvalues,\
    KBand, KBandSegment, Dos, BasisSetCellDependent, MethodBasisSet


def get_key_values(val_in):
    val = [v for v in val_in.split('\n') if '=' in v]
    data = {}
    pattern = re.compile(r'([A-Z]+)\s*=\s*([\w\-]+\s{0,3}[\d\. ]*)')

    def convert(v):
        if isinstance(v, list):
            v = [convert(vi) for vi in v]
        elif isinstance(v, str):
            try:
                v = float(v) if '.' in v else int(v)
            except Exception:
                pass
        else:
            pass
        return v

    for v in val:
        res = pattern.findall(v)
        for resi in res:
            vi = resi[1].split()
            vi = vi[0] if len(vi) == 1 else vi
            vi = vi.strip() if isinstance(vi, str) else vi
            vi = vi == 'T' if vi in ['T', 'F'] else vi
            data[resi[0]] = convert(vi)
    return data


class OutcarParser(TextParser):
    def __init__(self):
        self._chemical_symbols = None

        super().__init__(None)

    def init_quantities(self):
        self._quantities = []

        def str_to_array(val_in):
            val = [v.split() for v in val_in.strip().split('\n') if '--' not in v]
            return np.array([v[0:3] for v in val], float), np.array([v[3:6] for v in val], float)

        def str_to_stress(val_in):
            val = [float(v) for v in val_in.strip().split()]
            stress = np.zeros((3, 3))
            stress[0][0] = val[0]
            stress[1][1] = val[1]
            stress[2][2] = val[2]
            stress[0][1] = stress[1][0] = val[3]
            stress[1][2] = stress[2][1] = val[4]
            stress[0][2] = stress[2][0] = val[5]
            return stress

        def str_to_header(val_in):
            version, build_date, build_type, platform, date, time, parallel = val_in.split()
            parallel = 'parallel' if parallel == 'running' else parallel
            sub_version = '%s %s %s' % (build_date, build_type, parallel)
            date = date.replace('.', ' ')
            return dict(version=version, sub_version=sub_version, platform=platform, date=date, time=time)

        scf_iteration = [
            Quantity(
                'total_energy', r'free energy\s*TOTEN\s*=\s*([\d\.]+)\s*eV',
                repeats=False, dtype=float),
            Quantity(
                'energy_entropy0', r'energy without entropy\s*=\s*([\d\.]+)',
                repeats=False, dtype=float),
            Quantity(
                'energy_T0', r'energy\(sigma\->0\)\s*=\s*([\d\.]+)',
                repeats=False, dtype=float),
            Quantity(
                'energy_components',
                r'Free energy of the ion-electron system \(eV\)\s*\-+([\s\S]+?)\-{100}',
                str_operation=get_key_values, convert=False)
        ]

        self._quantities.append(Quantity(
            'calculation',
            r'(\-+\s*Iteration\s*\d+\(\s*1\)\s*\-+[\s\S]+?FREE ENERGIE OF THE ION-ELECTRON SYSTEM \(eV\))'
            r'([\s\S]+?\-{100})',
            repeats=True, sub_parser=TextParser(quantities=[
                Quantity(
                    'scf_iteration',
                    r'Iteration\s*\d+\(\s*\d+\)([\s\S]+?energy\(sigma\->0\)\s*=\s*[\d\.]+)',
                    repeats=True, sub_parser=TextParser(quantities=scf_iteration)),
                Quantity(
                    'total_energy',
                    r'FREE ENERGIE[\s\S]+?free\s*energy\s*TOTEN\s*=\s*([\d\.]+)',
                    repeats=False, dtype=float),
                Quantity(
                    'energy_entropy0',
                    r'FREE ENERGIE[\s\S]+?energy\s*without\s*entropy\s*=\s*([\d\.]+)',
                    repeats=False, dtype=float),
                Quantity(
                    'energy_T0',
                    r'FREE ENERGIE[\s\S]+?energy\(sigma\->0\)\s*=\s*([\d\.]+)',
                    repeats=False, dtype=float),
                Quantity(
                    'stress',
                    r'in kB\s*([\-\d\. E]+)', str_operation=str_to_stress, convert=False),
                Quantity(
                    'positions_forces',
                    r'POSITION\s*TOTAL\-FORCE \(eV/Angst\)\s*\-+\s*([\d\.\s\-E]+)',
                    str_operation=str_to_array, convert=False),
                Quantity(
                    'lattice_vectors',
                    r'direct lattice vectors\s*reciprocal lattice vectors([\d\.\s\-E]+)',
                    str_operation=str_to_array, convert=False),
                Quantity(
                    'converged',
                    r'aborting loop because (EDIFF is reached)', repeats=False,
                    dtype=str, convert=False),
                Quantity(
                    'fermi_energy', r'E\-fermi :\s*([\d\.]+)', dtype=str, repeats=False),
                Quantity(
                    'eigenvalues',
                    r'band No\.\s*band energies\s*occupation\s*([\d\.\s\-]+?)(?:k\-point|spin|\-{10})',
                    repeats=True, dtype=float)
            ])))

        self._quantities.append(Quantity(
            'header',
            r'vasp\.([\d\.]+)\s*(\w+)\s*[\s\S]+?\)\s*(\w+)\s*'
            r'executed on\s*(\w+)\s*date\s*([\d\.]+)\s*([\d\:]+)\s*(\w+)',
            repeats=False, str_operation=str_to_header, convert=False))

        self._quantities.append(Quantity(
            'parameters', r'Startparameter for this run:([\s\S]+?)\-{100}',
            str_operation=get_key_values, repeats=False, convert=False))

        self._quantities.append(Quantity(
            'ions_per_type', r'ions per type =([ \d]+)', dtype=int, repeats=False))

        self._quantities.append(Quantity(
            'species', r'TITEL\s*=\s*(\w+) ([A-Z][a-z]*)', dtype=str, repeats=True))

        self._quantities.append(Quantity(
            'mass_valence', r'POMASS\s*=\s*([\d\.]+);\s*ZVAL\s*=\s*([\d\.]+)\s*mass and valenz',
            dtype=float, repeats=True
        ))

        self._quantities.append(Quantity(
            'kpoints',
            r'k-points in reciprocal lattice and weights:[\s\S]+?\n([\d\.\s\-]+)',
            repeats=False, dtype=float))

        self._quantities.append(Quantity(
            'nbands', r'NBANDS\s*=\s*(\d+)', dtype=int, repeats=False))


class VASPParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/vasp', code_name='VASP', code_homepage='https://www.vasp.at/',
            mainfile_mime_re=r'(application/.*)|(text/.*)',
            mainfile_name_re=r'xml(\.gz|\.bz|\.xz)?',
            mainfile_contents_re=(
                r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
                r'?\s*<modeling>'
                r'?\s*<generator>'
                r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
                r'?|^\svasp[\.\d]+\s*\w+\s*\(build'),
            supported_compressions=['gz', 'bz2', 'xz'], mainfile_alternative=True)

        self._metainfo_env = m_env
        self.vasprun_parser = XMLParser()
        self.outcar_parser = OutcarParser()
        self.parser = None
        self._incar = None
        self._kpoints_info = None
        self._atom_info = None
        self._header = None

        self.metainfo_mapping = {
            'e_fr_energy': 'energy_free', 'e_wo_entrp': 'energy_total',
            'e_0_energy': 'energy_total_T0', 'hartreedc': 'energy_hartree_error',
            'XCdc': 'energy_XC', 'forces': 'atom_forces', 'stress': 'stress_tensor',
            'energy_total': 'energy_free', 'energy_T0': 'energy_total_T0',
            'energy_entropy0': 'energy_total', 'DENC': 'energy_hartree_error',
            'EXHF': 'energy_X', 'EBANDS': 'energy_sum_eigenvalues'}

        self.xc_functional_mapping = {
            '91': ['GGA_X_PW91', 'GGA_C_PW91'], 'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
            'RP': ['GGA_X_RPBE', 'GGA_C_PBE'], 'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
            'MK': ['GGA_X_OPTB86_VDW'], '--': ['GGA_X_PBE', 'GGA_C_PBE']}

    def init_parser(self, mainfile, logger):
        if 'OUTCAR' in mainfile:
            self.parser = 'outcar'
            self.outcar_parser.mainfile = mainfile
            self.outcar_parser.logger = logger
        else:
            self.parser = 'vasprun'
            self.vasprun_parser.mainfile = mainfile
            self.vasprun_parser.logger = logger
        self._incar = None
        self._kpoints_info = None
        self._atom_info = None
        self._header = None

    def reuse_parser(self, parser):
        self.outcar_parser.quantities = parser.outcar_parser.quantities

    def _get_key_values(self, path):
        if self.parser == 'outcar':

            if not os.path.isfile(path):
                return dict()
            with self.outcar_parser.open(path) as f:
                text = f.read()
                text = text.decode() if isinstance(text, bytes) else text
                data = get_key_values(text)
            return data

        elif self.parser == 'vasprun':
            if self.vasprun_parser.root is None:
                return dict()

            elements = self.vasprun_parser.root.findall(os.path.join('./', path))

            dtypes = {'string': str, 'int': int, 'logical': bool, '': float}
            names = [e.attrib.get('name', None) for e in elements]
            types = [e.attrib.get('type', '') for e in elements]
            values = [e.text if e.text else '' for e in elements]

            def convert(val, dtype):
                if isinstance(val, list):
                    return [convert(v, dtype) for v in val]
                else:
                    return dtype(val)

            for i in range(len(values)):
                if names[i] is None:
                    continue
                dtype = dtypes.get(types[i], str)
                values[i] = values[i].split() if dtype != str else values[i].strip()
                values[i] = values[i][0] if len(values[i]) == 1 else values[i]
                array = isinstance(values[i], list) or path.endswith('/v/')

                if dtype == bool:
                    values[i] = [v == 'T' for v in values[i]] if array else values[i] == 'T'
                if dtype == float:
                    # prevent nan printouts
                    if array:
                        values[i] = ['nan' if '*' in v else v for v in values[i]]
                    else:
                        values[i] = 'nan' if '*' in values[i] else values[i]
                # using numpy array does not work
                values[i] = convert(values[i], dtype)

            return dict(zip(names, values))

        return dict()

    @property
    def header(self):
        if self._header is None:
            self._header = dict()
            if self.parser == 'outcar':
                self._header = self.outcar_parser.get('header', {})
                self._header['program'] = 'vasp'
            elif self.parser == 'vasprun':
                self._header = self._get_key_values('generator/i')
            for key, val in self._header.items():
                if not isinstance(val, str):
                    self._header[key] = ' '.join(val)
        return self._header

    # why make a distinction between incar_in and incar_out?
    @property
    def incar(self):
        if self._incar is None:
            self._incar = dict(incar=None, incar_out=None)
        if self._incar['incar'] is None:
            self.get_incar()
        if self._incar['incar_out'] is None:
            self.get_incar_out()

        incar = self._incar['incar']
        incar.update(self._incar['incar_out'])

        return incar

    def _fix_incar(self, incar):
        # fix for LORBIT, list is read
        lorbit = incar.get('LORBIT', None)
        if isinstance(lorbit, list):
            incar['LORBIT'] = lorbit[0]

    def get_incar(self):
        if self._incar is not None and self._incar.get('incar', None) is not None:
            return self._incar.get('incar')

        incar = dict()
        if self._incar is None:
            self._incar = dict(incar=None, incar_out=None)
        if self.parser == 'outcar':
            path = os.path.join(self.outcar_parser.maindir, 'INCAR%s' % os.path.basename(
                self.outcar_parser.mainfile).strip('OUTCAR'))
            path = path if os.path.isfile(path) else os.path.join(
                self.outcar_parser.maindir, 'INCAR')
            incar = self._get_key_values(path)
        elif self.parser == 'vasprun':
            incar = self._get_key_values('incar/i')
        self._fix_incar(incar)
        self._incar['incar'] = incar
        return incar

    def get_incar_out(self):
        if self._incar is not None and self._incar.get('incar_out', None) is not None:
            return self._incar.get('incar_out')

        incar = dict()
        if self._incar is None:
            self._incar = dict(incar=None, incar_out=None)
        if self.parser == 'outcar':
            incar = self.outcar_parser.get('parameters', {})
        elif self.parser == 'vasprun':
            for tag in ['i', 'v']:
                incar.update(self._get_key_values('parameters//%s' % tag))

        self._incar['incar_out'] = incar
        self._fix_incar(incar)
        return incar

    @property
    def ispin(self):
        return self.incar.get('ISPIN', 1)

    @property
    def ibrion(self):
        val = self.incar.get('IBRION,', None)
        if val is None:
            val = -1 if self.incar.get('NSW', 0) in [0, 1] else 0
        return val

    @property
    def n_calculations(self):
        calculations = []
        if self.parser == 'outcar':
            calculations = self.outcar_parser.get('calculation')
        elif self.parser == 'vasprun':
            calculations = self.vasprun_parser.get('calculation/')

        if isinstance(calculations, dict):
            return 1
        elif isinstance(calculations, list):
            return len(calculations)
        else:
            return 0

    @property
    def kpoints_info(self):
        if self._kpoints_info is None:
            self._kpoints_info = dict()
            if self.parser == 'outcar':
                kpts_occs = self.outcar_parser.get('kpoints')
                if kpts_occs is not None:
                    kpts_occs = np.reshape(kpts_occs, (len(kpts_occs) // 4, 4)).T
                    self._kpoints_info['k_mesh_points'] = kpts_occs[0:3].T
                    self._kpoints_info['k_mesh_weights'] = kpts_occs[3].T

            elif self.parser == 'vasprun':
                self._kpoints_info['x_vasp_k_points_generation_method'] = self.vasprun_parser.get(
                    'kpoints/generation/param')
                self._kpoints_info['k_mesh_points'] = self.vasprun_parser.get(
                    'kpoints/varray/[@name="kpointlist"]/v')
                self._kpoints_info['k_mesh_weights'] = self.vasprun_parser.get(
                    'kpoints/varray/[@name="weights"]/v')
                self._kpoints_info['x_vasp_tetrahedrons_list'] = self.vasprun_parser.get(
                    'kpoints/varray/[@name="tetrahedronlist"]/v')
                self._kpoints_info['divisions'] = self.vasprun_parser.get(
                    'kpoints/generation/i/[@name="divisions"]')
                self._kpoints_info['points'] = self.vasprun_parser.get('kpoints/generation/v')
                volumeweight = self.vasprun_parser.get('kpoints/i/[@name="volumeweight"]')
                if volumeweight is not None:
                    volumeweight = pint.Quantity(volumeweight, 'angstrom ** 3').to('m**3')
                    # TODO set propert unit in metainfo
                    self._kpoints_info['x_vasp_tetrahedron_volume'] = volumeweight.magnitude

        return self._kpoints_info

    @property
    def n_bands(self):
        for n in range(self.n_calculations):
            if self.parser == 'outcar':
                val = self.outcar_parser.get('calculation', [{}] * (n + 1))[n].get(
                    'eigenvalues', [None])[0]
                if val is not None:
                    return len(val) // 3
            elif self.parser == 'vasprun':
                val = self.vasprun_parser.get(
                    'calculation[%d]/eigenvalues/array/set[1]/set[1]/set[1]/r' % (n + 1))
                if val is not None:
                    return len(val)
        # check consistency with eigenvalues
        return self.incar.get('NBANDS', 0)

    @property
    def n_dos(self):
        if self.parser == 'outcar':
            path = os.path.join(self.outcar_parser.maindir, 'DOSCAR%s' % os.path.basename(
                self.outcar_parser.mainfile).strip('OUTCAR'))
            path = path if os.path.isfile(path) else os.path.join(
                self.outcar_parser.maindir, 'DOSCAR')
            with open(path) as f:
                for n in range(6):
                    line = f.readline()
                return int(line.split()[2])

        elif self.parser == 'vasprun':
            for n in range(self.n_calculations):
                val = self.vasprun_parser.get(
                    'calculation[%d]/dos/total/array/set[1]/set[1]/r' % (n + 1))
                if val is not None:
                    return len(val)
        return self._incar.get('NEDOS', 0)

    @property
    def atom_info(self):
        if self._atom_info is None:
            self._atom_info = {}

            if self.parser == 'outcar':
                ions = self.outcar_parser.get('ions_per_type', [])
                species = self.outcar_parser.get('species', [])
                ions = [ions] if isinstance(ions, int) else ions
                mass_valence = self.outcar_parser.get('mass_valence', [])
                assert len(ions) == len(species)
                self.atom_info['n_atoms'] = sum(ions)
                self.atom_info['n_types'] = len(species)

                element = []
                atomtype = []
                for n in range(len(ions)):
                    element.extend([str(species[n][1])] * ions[n])
                    atomtype.extend([(n + 1)] * ions[n])
                self.atom_info['atoms'] = dict(element=element, atomtype=atomtype)

                self.atom_info['atomtypes'] = dict(
                    atomspertype=ions, element=[s[0] for s in species],
                    mass=[m[0] for m in mass_valence], valence=[m[1] for m in mass_valence],
                    pseudopotential=[s[1] for s in species])

            elif self.parser == 'vasprun':
                self._atom_info['n_atoms'] = self.vasprun_parser.get('atominfo/atoms')
                self._atom_info['n_types'] = self.vasprun_parser.get('atominfo/types')

                array_keys = [
                    e.get('name', None) for e in self.vasprun_parser.get('atominfo/array/')]
                for key in array_keys:
                    array_info = {}
                    fields = self.vasprun_parser.get('atominfo/array/[@name="%s"]/field' % key)
                    rcs = self.vasprun_parser.get('atominfo/array/[@name="%s"]/set/rc/' % key, 1)
                    for n in range(max(len(rcs), 1)):
                        val = self.vasprun_parser.get(
                            'atominfo/array/[@name="%s"]/set/rc[%d]/c' % (key, n + 1))
                        for i in range(len(fields)):
                            array_info.setdefault(fields[i], [])
                            array_info[fields[i]].append(val[i])
                    self._atom_info[key] = array_info

        return self._atom_info

    def get_n_scf(self, n_calc):
        if self.parser == 'outcar':
            return len(self.outcar_parser.get(
                'calculation', [{}] * (n_calc + 1))[n_calc].get('scf_iteration', []))
        elif self.parser == 'vasprun':
            return len(self.vasprun_parser.get('calculation[%d]/scstep/' % (n_calc + 1)))
        return 0

    def get_structure(self, n_calc):
        if self.parser == 'outcar':
            cell = self.outcar_parser.get(
                'calculation', [{}] * (n_calc + 1))[n_calc].get('lattice_vectors', [None])[0]
            positions = self.outcar_parser.get(
                'calculation', [{}] * (n_calc + 1))[n_calc].get('positions_forces', [None])[0]
            selective = None
            nose = None

        elif self.parser == 'vasprun':
            cell = self.vasprun_parser.get(
                'calculation[%d]/structure/crystal/varray/[@name="basis"]/v' % (n_calc + 1))
            positions = self.vasprun_parser.get(
                'calculation[%d]/structure/varray/[@name="positions"]/v' % (n_calc + 1))
            selective = self.vasprun_parser.get(
                'calculation[%d]/structure/varray/[@name="selective"]/v' % (n_calc + 1))
            nose = self.vasprun_parser.get('calculation[%d]/structure/nose/v' % (n_calc + 1))

        else:
            cell = None
            positions = None
            selective = None
            nose = None

        if positions is not None:
            positions = pint.Quantity(np.dot(positions, cell), 'angstrom')
        if cell is not None:
            cell = pint.Quantity(cell, 'angstrom')

        return dict(cell=cell, positions=positions, selective=selective, nose=nose)

    def get_energies(self, n_calc, n_scf):
        if self.parser == 'outcar':
            energies = dict()
            if n_scf is None:
                section = self.outcar_parser.get(
                    'calculation', [{}] * (n_calc + 1))[n_calc]
            else:
                section = self.outcar_parser.get('calculation', [{}] * (
                    n_calc + 1))[n_calc].get('scf_iteration', [{}] * (n_scf + 1))[n_scf]
            for key in ['total_energy', 'energy_T0', 'energy_entropy0']:
                energies[key] = section.get(key)

            energies.update(section.get('energy_components', {}))
            return energies

        elif self.parser == 'vasprun':
            if n_scf is None:
                return self._get_key_values(
                    'calculation[%d]/energy/i' % (n_calc + 1))
            else:
                return self._get_key_values(
                    'calculation[%d]/scstep[%d]/energy/i' % (n_calc + 1, n_scf + 1))
        return dict()

    def get_forces_stress(self, n_calc):
        forces = stress = None
        if self.parser == 'outcar':
            forces = self.outcar_parser.get('calculation', [{}] * (n_calc + 1))[n_calc].get(
                'positions_forces', [None, None])[1]
            stress = self.outcar_parser.get('calculation', [{}] * (n_calc + 1))[n_calc].get(
                'stress', None)
        elif self.parser == 'vasprun':
            forces = self.vasprun_parser.get(
                'calculation[%d]/varray/[@name="forces"]/v' % (n_calc + 1))
            stress = self.vasprun_parser.get(
                'calculation[%d]/varray/[@name="stress"]/v' % (n_calc + 1))

        return forces, stress

    def get_eigenvalues(self, n_calc):
        n_kpts = len(self.kpoints_info.get('k_mesh_points', []))
        eigenvalues = None

        if self.parser == 'outcar':
            eigenvalues = self.outcar_parser.get(
                'calculation', [{}] * (n_calc + 1))[n_calc].get('eigenvalues')

            if eigenvalues is None:
                return
            n_eigs = len(eigenvalues) // (self.ispin * n_kpts)
            eigenvalues = np.reshape(eigenvalues, (n_eigs, self.ispin, n_kpts, self.n_bands, 3))
            # eigenvalues can also be printed every scf iteration but we only save the
            # last one, which corresponds to the calculation
            eigenvalues = np.transpose(eigenvalues)[1:].T[-1]

        elif self.parser == 'vasprun':
            root = 'calculation[%d]/eigenvalues/array' % (n_calc + 1)
            if self.vasprun_parser.get(root) is None:
                return
            eigenvalues = self.vasprun_parser.root.findall('./%s/set//r' % root)
            eigenvalues = np.array([e.text.split() for e in eigenvalues], dtype=float)
            eigenvalues = np.reshape(eigenvalues, (self.ispin, n_kpts, self.n_bands, 2))

        return eigenvalues

    def get_total_dos(self, n_calc):
        dos_energies = dos_values = dos_integrated = e_fermi = None

        if self.parser == 'outcar':
            if n_calc != (self.n_calculations - 1):
                return dos_energies, dos_values, dos_integrated, e_fermi
            path = os.path.join(self.outcar_parser.maindir, 'DOSCAR%s' % os.path.basename(
                self.outcar_parser.mainfile).strip('OUTCAR'))
            path = path if os.path.isfile(path) else os.path.join(
                self.outcar_parser.maindir, 'DOSCAR')
            if not os.path.isfile(path):
                return dos_energies, dos_values, dos_integrated, e_fermi

            dos = []
            n_dos = 0
            with self.outcar_parser.open(path) as f:
                for i, line in enumerate(f):
                    if i < 5:
                        continue
                    if i == 5:
                        e_fermi = float(line.split()[1])
                        n_dos = int(line.split()[2])
                    if i > 5:
                        dos.append([float(v) for v in line.split()])
                    if i >= n_dos + 5:
                        break
            if not dos:
                return dos_energies, dos_values, dos_integrated, e_fermi

            # DOSCAR fomat (spin) energy dos_up dos_down integrated_up integrated_down
            dos = np.transpose(dos)
            dos_energies = dos[0]
            dos_values = dos[1: 1 + self.ispin]
            dos_integrated = dos[1 + self.ispin: 2 * self.ispin + 1]

        elif self.parser == 'vasprun':
            root = 'calculation[%d]/dos' % (n_calc + 1)
            if self.vasprun_parser.get(root) is None:
                return dos_energies, dos_values, dos_integrated, e_fermi

            dos = self.vasprun_parser.root.findall('./%s/total/array/set//r' % root)
            dos = np.array([e.text.split() for e in dos], dtype=float)
            dos = np.reshape(dos, (self.ispin, self.n_dos, 3))

            dos = np.transpose(dos)
            dos_energies = dos[0].T[0]
            dos_values = dos[1].T
            dos_integrated = dos[2].T

            # unit of dos in vasprun is states/eV/cell
            cell = self.get_structure(n_calc)['cell']
            volume = np.abs(np.linalg.det(cell.to('m').magnitude))
            dos_values *= volume

            e_fermi = self.vasprun_parser.get('%s/i/[@name="efermi"]' % root, 0.0)

        return dos_energies, dos_values, dos_integrated, e_fermi

    def get_partial_dos(self, n_calc):
        n_atoms = self.atom_info['n_atoms']
        dos = fields = None
        if self.parser == 'outcar':
            if n_calc != (self.n_calculations - 1):
                return dos, fields
            path = os.path.join(self.outcar_parser.maindir, 'DOSCAR%s' % os.path.basename(
                self.outcar_parser.mainfile).strip('OUTCAR'))
            path = path if os.path.isfile(path) else os.path.join(
                self.outcar_parser.maindir, 'DOSCAR')
            if not os.path.isfile(path):
                return dos, fields

            dos = []
            n_dos = 0
            atom = 1
            with self.outcar_parser.open(path) as f:
                for i, line in enumerate(f):
                    if i == 0:
                        if int(line.split()[2]) != 1:
                            return None, None
                    elif i < 5:
                        continue
                    elif i == 5:
                        n_dos = int(line.split()[2])
                    if i == ((n_dos + 1) * atom) + 5:
                        atom += 1
                        continue
                    if i > (n_dos + 6):
                        dos.append([float(v) for v in line.split()])
                    if i >= ((n_dos + 1) * (n_atoms + 1) + 5):
                        break
            if len(dos) == 0:
                return None, None

            dos = np.transpose(dos)[1:]
            n_lm = len(dos) // self.ispin
            dos = np.reshape(dos, (n_lm, self.ispin, n_atoms, n_dos))

            if n_lm == 3:
                fields = ['s', 'p', 'd']
            elif n_lm == 9:
                fields = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
            elif n_lm == 16:
                fields = [
                    's', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'f-3',
                    'f-2', 'f-1', 'f0', 'f1', 'f2', 'f3']
            else:
                fields = [None] * n_lm
                self.outcar_parser.logger.warn(
                    'Cannot determine lm fields for n_lm', data=dict(n_lm=n_lm))

        elif self.parser == 'vasprun':
            root = 'calculation[%d]/dos' % (n_calc + 1)
            if self.vasprun_parser.get('%s/partial' % root) is None:
                return dos, fields

            # TODO use atomprojecteddos section
            fields = self.vasprun_parser.get('%s/partial/array/field' % root)
            dos = self.vasprun_parser.root.findall('./%s/partial/array/set//r' % root)
            dos = np.array([e.text.split() for e in dos], dtype=float)
            dos = np.reshape(dos, (n_atoms, self.ispin, self.n_dos, len(fields)))

            fields = [field for field in fields if field != 'energy']
            dos = np.transpose(dos)[1:]
            dos = np.transpose(dos, axes=(0, 2, 3, 1))

        return dos, fields

    def parse_incarsout(self):
        sec_method = self.archive.section_run[-1].section_method[-1]

        unknown_incars = dict()
        for key, val in self.get_incar_out().items():
            key = 'x_vasp_incarOut_%s' % key
            if hasattr(Method, key):
                try:
                    setattr(sec_method, key, val)
                except Exception as e:
                    self.logger.warn(e)
            else:
                self.logger.warn('Unknown incar parameter', data=dict(key=key))
                # not sure why it gives out an error when np.array
                if isinstance(val, np.ndarray):
                    val = list(val)
                unknown_incars[key] = val

        sec_method.x_vasp_unknown_incars = unknown_incars

        prec = 1.3 if 'acc' in self.incar.get('PREC', '') else 1.0
        sec_basis_set_cell_dependent = self.archive.section_run[-1].m_create(
            BasisSetCellDependent)
        sec_basis_set_cell_dependent.basis_set_planewave_cutoff = pint.Quantity(
            self.incar.get('ENMAX', 0.0) * prec, 'eV')

        sec_method_basis_set = sec_method.m_create(MethodBasisSet)
        sec_method_basis_set.mapping_section_method_basis_set_cell_associated = sec_basis_set_cell_dependent

    def parse_method(self):
        sec_method = self.archive.section_run[-1].m_create(Method)

        # input incar
        for key, val in self.get_incar().items():
            key = 'x_vasp_incar_%s' % key
            try:
                setattr(sec_method, key, val)
            except Exception:
                self.logger.warn('Error setting metainfo', data=dict(key=key))

        sec_method.electronic_structure_method = 'DFT+U' if self.incar.get(
            'LDAU', False) else 'DFT'

        # kpoints
        for key, val in self.kpoints_info.items():
            if val is not None:
                try:
                    setattr(sec_method, key, val)
                except Exception:
                    self.logger.warn('Error setting metainfo', data=dict(key=key))

        # atom properties
        atomtypes = self.atom_info.get('atomtypes', {})
        element = atomtypes.get('element', [])
        atom_counts = {e: 0 for e in element}
        for i in range(len(element)):
            sec_method_atom_kind = sec_method.m_create(MethodAtomKind)
            atom_number = ase.data.atomic_numbers.get(element[i], 0)
            sec_method_atom_kind.method_atom_kind_atom_number = atom_number
            atom_label = '%s%d' % (
                element[i], atom_counts[element[i]]) if atom_counts[element[i]] > 0 else element[i]
            sec_method_atom_kind.method_atom_kind_label = str(atom_label)
            sec_method_atom_kind.method_atom_kind_mass = pint.Quantity(
                atomtypes.get('mass', [1] * (i + 1))[i], 'amu')
            sec_method_atom_kind.method_atom_kind_explicit_electrons = atomtypes.get(
                'valence', [0] * (i + 1))[i]
            pseudopotential = atomtypes.get('pseudopotential')[i]
            pseudopotential = ' '.join(pseudopotential) if isinstance(
                pseudopotential, list) else pseudopotential
            sec_method_atom_kind.method_atom_kind_pseudopotential_name = str(pseudopotential)
            atom_counts[element[i]] += 1
        sec_method.x_vasp_atom_kind_refs = sec_method.section_method_atom_kind

        self.parse_incarsout()

        gga = self.incar.get('GGA', None)
        if gga is not None:
            xc_functionals = self.xc_functional_mapping.get(gga, [])
            for xc_functional in xc_functionals:
                sec_xc_functional = sec_method.m_create(XCFunctionals)
                sec_xc_functional.XC_functional_name = xc_functional

    def parse_configurations(self):
        sec_run = self.archive.section_run[-1]

        def parse_system(n_calc):
            sec_system = sec_run.m_create(System)

            structure = self.get_structure(n_calc)
            cell = structure.get('cell', None)
            if cell is not None:
                sec_system.simulation_cell = cell
                sec_system.lattice_vectors = cell

            sec_system.configuration_periodic_dimensions = [True] * 3
            sec_system.atom_labels = self.atom_info.get('atoms', {}).get('element', [])

            positions = structure.get('positions', None)
            if positions is not None:
                sec_system.atom_positions = positions

            selective = structure.get('selective', None)
            if selective is not None:
                selective = [[s[i] == 'T' for i in range(len(s))] for s in selective]
                sec_system.x_vasp_selective_dynamics = selective

            nose = structure.get('nose')
            if nose is not None:
                sec_system.x_vasp_nose_thermostat = nose

            return sec_system

        def parse_energy(n_calc, n_scf=None):
            if n_scf is None:
                section = sec_run.m_create(SingleConfigurationCalculation)
                ext = ''
            else:
                section = sec_run.section_single_configuration_calculation[-1].m_create(ScfIteration)
                ext = '_scf_iteration'

            energies = self.get_energies(n_calc, n_scf)
            for key, val in energies.items():
                if val is None:
                    continue
                val = pint.Quantity(val, 'eV')
                metainfo_key = self.metainfo_mapping.get(key, None)
                if metainfo_key is None:
                    continue
                try:
                    setattr(section, '%s%s' % (metainfo_key, ext), val)
                except Exception:
                    self.logger.warn('Error setting metainfo', data=dict(key=key))

            return section

        def parse_eigenvalues(n_calc):
            eigenvalues = self.get_eigenvalues(n_calc)
            if eigenvalues is None:
                return

            sec_scc = sec_run.section_single_configuration_calculation[-1]
            sec_eigenvalues = sec_scc.m_create(Eigenvalues)
            eigenvalues = np.transpose(eigenvalues)
            eigs = eigenvalues[0].T
            occs = eigenvalues[1].T

            # get valence(conduction) and maximum(minimum)
            valence_max = [max([eigs[i, o[0], o[1]] for o in np.argwhere(
                occs[i] >= 0.5)]) for i in range(len(eigs))]
            conduction_min = [min([eigs[i, o[0], o[1]] for o in np.argwhere(
                occs[i] < 0.5)]) for i in range(len(eigs))]
            sec_scc.energy_reference_highest_occupied = pint.Quantity(valence_max, 'eV')
            sec_scc.energy_reference_lowest_unoccupied = pint.Quantity(conduction_min, 'eV')
            if self.kpoints_info.get('x_vasp_k_points_generation_method', None) == 'listgenerated':
                # I removed normalization since it imho it should be done by normalizer
                sec_k_band = sec_scc.m_create(KBand)
                divisions = int(self.kpoints_info.get('divisions', None))
                kpoints = self.kpoints_info.get('k_mesh_points', [])
                n_segments = len(kpoints) // divisions
                kpoints = np.reshape(kpoints, (n_segments, divisions, 3))
                eigs = np.reshape(eigs, (
                    self.ispin, n_segments, divisions, self.n_bands))
                occs = np.reshape(occs, (
                    self.ispin, n_segments, divisions, self.n_bands))
                eigs = np.transpose(eigs, axes=(1, 0, 2, 3))
                occs = np.transpose(occs, axes=(1, 0, 2, 3))
                for n in range(n_segments):
                    sec_k_band_segment = sec_k_band.m_create(KBandSegment)
                    sec_k_band_segment.band_energies = pint.Quantity(eigs[n], 'eV')
                    sec_k_band_segment.band_occupations = occs[n]
                    sec_k_band_segment.band_k_points = kpoints[n]
                    sec_k_band_segment.band_segm_start_end = np.asarray(
                        [kpoints[n, 0], kpoints[n, divisions - 1]])
            else:
                sec_eigenvalues.eigenvalues_values = pint.Quantity(eigs, 'eV')
                sec_eigenvalues.eigenvalues_occupation = occs

        def parse_dos(n_calc):
            energies, values, integrated, e_fermi = self.get_total_dos(n_calc)

            # total dos
            if values is not None:
                sec_scc = sec_run.section_single_configuration_calculation[-1]
                sec_dos = sec_scc.m_create(Dos)
                sec_dos.dos_energies = pint.Quantity(energies, 'eV')

                sec_dos.dos_values = pint.Quantity(values, '1/eV').to('1/joule').magnitude
                sec_dos.dos_integrated_values = integrated

                sec_dos.energy_reference_fermi = pint.Quantity([e_fermi] * self.ispin, 'eV')

            # partial dos
            # TODO: I do not know how the f-orbitals are arranged
            lm_converter = {
                's': [0, 0], 'p': [1, -1], 'px': [1, 0], 'py': [1, 1], 'pz': [1, 2],
                'd': [2, -1], 'dx2': [2, 0], 'dxy': [2, 1], 'dxz': [2, 2], 'dy2': [2, 3],
                'dyz': [2, 4], 'dz2': [2, 5], 'f': [3, -1], 'f-3': [3, 0], 'f-2': [3, 1],
                'f-1': [3, 2], 'f0': [3, 3], 'f1': [3, 4], 'f2': [3, 5], 'f3': [3, 6]}

            dos, fields = self.get_partial_dos(n_calc)
            if dos is not None:
                sec_dos.dos_values_lm = pint.Quantity(dos, '1/eV').to('1/joule').magnitude
                sec_dos.dos_lm = [lm_converter.get(field, [-1, -1]) for field in fields]
                sec_dos.dos_m_kind = 'polynomial'

        for n in range(self.n_calculations):
            # energies
            sec_scc = parse_energy(n, None)
            for n_scf in range(self.get_n_scf(n)):
                parse_energy(n, n_scf)

            # forces and stress
            forces, stress = self.get_forces_stress(n)
            if forces is not None:
                sec_scc.atom_forces = pint.Quantity(forces, 'eV/angstrom')
            if stress is not None:
                # TODO verify if stress unit in xml is also kbar
                sec_scc.stress_tensor = pint.Quantity(stress, 'kbar')

            # structure
            sec_system = parse_system(n)
            sec_scc.single_configuration_calculation_to_system_ref = sec_system
            sec_scc.single_configuration_to_calculation_method_ref = sec_run.section_method[-1]

            # eigenvalues
            parse_eigenvalues(n)

            # dos
            parse_dos(n)

    def parse(self, filepath, archive, logger):
        self.filepath = filepath
        self.archive = archive
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.init_parser(filepath, logger)

        sec_run = self.archive.m_create(Run)
        program_name = self.header.get('program', '')
        if program_name.strip().upper() != 'VASP':
            self.logger.error('invalid program name', data=dict(program_name=program_name))
            return
        sec_run.program_name = 'VASP'

        version = ' '.join([self.header.get(key, '') for key in [
            'version', 'subversion', 'platform']]).strip()
        if version:
            sec_run.program_version = version

        sec_run.program_basis_set_type = 'plane waves'

        date = self.header.get('date')
        if date is not None:
            date = datetime.strptime(date.strip(), '%Y %m %d').date()
            time = self.header.get('time', '0:0:0')
            time = datetime.strptime(time.strip(), '%H:%M:%S').timetz()
            dtime = datetime.combine(date, time) - datetime.utcfromtimestamp(0)
            sec_run.program_compilation_datetime = dtime.total_seconds()

        self.parse_method()

        self.parse_configurations()