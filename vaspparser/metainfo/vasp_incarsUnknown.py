import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from vaspparser.metainfo import vasp_incars
from nomad.datamodel.metainfo import public
from vaspparser.metainfo import vasp_incarsOut

m_package = Package(
    name='vasp_incarsUnknown_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='vasp.incarsUnknown.nomadmetainfo.json'))


class x_vasp_incarUnknown_param(MCategory):
    '''
    Unknown incar parameters.  By 'unknown' we refer to incar parameters in the OLD
    (predated Sept.2019) 'vasp.nomadmetainfo.json' that are **not** listed in the current
    (Oct.2019) VaspWiki
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_vasp_incarUnknown_param'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_vasp_incarOut_ENMAX = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximum cutoff (normally specified only in POTCAR). Value prinded out after
        evaluating the input.
        ''',
        categories=[vasp_incarsOut.x_vasp_incarOut_param],
        a_legacy=LegacyDefinition(name='x_vasp_incarOut_ENMAX'))

    x_vasp_incar_ENMAX = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximum cutoff (normally specified only in POTCAR). Value prinded out after
        evaluating the input.
        ''',
        categories=[vasp_incars.x_vasp_incar_param],
        a_legacy=LegacyDefinition(name='x_vasp_incar_ENMAX'))

    x_vasp_incarOut_LCOMPAT = Quantity(
        type=bool,
        shape=[],
        description='''
        In vasp.4.2 the augmentation charges are forced to be zero at the boundary of the
        augmentation sphere, therefore results are slightly different from vasp.3.2
        (usually differences are smaller than 0.01 meV). The old behavior can be restored
        by setting LCOMPAT = .TRUE. in the INCAR file.
        ''',
        categories=[vasp_incarsOut.x_vasp_incarOut_param],
        a_legacy=LegacyDefinition(name='x_vasp_incarOut_LCOMPAT'))


m_package.__init_metainfo__()
