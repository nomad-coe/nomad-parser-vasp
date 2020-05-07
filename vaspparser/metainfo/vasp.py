import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from vaspparser.metainfo import vasp_incars
from nomad.datamodel.metainfo import public

m_package = Package(
    name='vasp_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='vasp.nomadmetainfo.json'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    vasp_build_date = Quantity(
        type=str,
        shape=[],
        description='''
        build date as string
        ''',
        categories=[public.program_info, public.accessory_info],
        a_legacy=LegacyDefinition(name='vasp_build_date'))

    vasp_src_date = Quantity(
        type=str,
        shape=[],
        description='''
        date of last modification of the source as string
        ''',
        categories=[public.program_info, public.accessory_info],
        a_legacy=LegacyDefinition(name='vasp_src_date'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_vasp_unknown_incars = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        INCAR variables uknown wrt to Vasp Wiki
        ''',
        a_legacy=LegacyDefinition(name='x_vasp_unknown_incars'))

    x_vasp_atom_kind_refs = Quantity(
        type=public.section_method_atom_kind,
        shape=['number_of_atoms'],
        description='''
        reference to the atom kinds of each atom
        ''',
        a_legacy=LegacyDefinition(name='x_vasp_atom_kind_refs'))

    x_vasp_k_points_generation_method = Quantity(
        type=str,
        shape=[],
        description='''
        k points generation  method
        ''',
        categories=[public.settings_potential_energy_surface, vasp_incars.x_vasp_incar_param, public.settings_k_points],
        a_legacy=LegacyDefinition(name='x_vasp_k_points_generation_method'))

    x_vasp_k_points = Quantity(
        type=int,
        shape=['x_vasp_number_of_k_points'],
        description='''
        k points
        ''',
        categories=[vasp_incars.x_vasp_incar_param],
        a_legacy=LegacyDefinition(name='x_vasp_k_points'))

    x_vasp_numer_of_magmom = Quantity(
        type=int,
        shape=[],
        description='''
        number of magnetic moments, number_of_atoms for ISPIN = 2, 3*number of atoms for
        non-collinear magnetic systems
        ''',
        categories=[vasp_incars.x_vasp_incar_param],
        a_legacy=LegacyDefinition(name='x_vasp_numer_of_magmom'))

    x_vasp_tetrahedrons_list = Quantity(
        type=np.dtype(np.int32),
        shape=['N', 5],
        description='''
        Rows of 5 elements. First the weight (symmetry degeneration), then the four corner
        points of each tetrahedron.
        ''',
        categories=[vasp_incars.x_vasp_incar_param],
        a_legacy=LegacyDefinition(name='x_vasp_tetrahedrons_list'))

    x_vasp_tetrahedron_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[1],
        description='''
        Volume weight of a single tetrahedron (all tetra's must have the same volume)
        ''',
        categories=[vasp_incars.x_vasp_incar_param],
        a_legacy=LegacyDefinition(name='x_vasp_tetrahedron_volume'))

    x_vasp_nose_thermostat = Quantity(
        type=np.dtype(np.float64),
        shape=[4],
        description='''
        Nose thermostat output
        ''',
        categories=[vasp_incars.x_vasp_incar_param],
        a_legacy=LegacyDefinition(name='x_vasp_nose_thermostat'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_vasp_selective_dynamics = Quantity(
        type=np.dtype(np.bool),
        shape=['number_of_atoms', 3],
        description='''
        Boolean array to eiter allow or forbid coordinate modifications during relaxation
        ''',
        a_legacy=LegacyDefinition(name='x_vasp_selective_dynamics'))


m_package.__init_metainfo__()
