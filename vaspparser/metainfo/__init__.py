import sys
from nomad.metainfo import Environment
from nomad.metainfo.legacy import LegacyMetainfoEnvironment
import vaspparser.metainfo.vasp
import nomad.datamodel.metainfo.common
import nomad.datamodel.metainfo.public
import nomad.datamodel.metainfo.general
import vaspparser.metainfo.vasp_incars
import vaspparser.metainfo.vasp_incarsOut
import vaspparser.metainfo.vasp_incarsUnknown

m_env = LegacyMetainfoEnvironment()
m_env.m_add_sub_section(Environment.packages, sys.modules['vaspparser.metainfo.vasp'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.metainfo.common'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.metainfo.public'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.metainfo.general'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['vaspparser.metainfo.vasp_incars'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['vaspparser.metainfo.vasp_incarsOut'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['vaspparser.metainfo.vasp_incarsUnknown'].m_package)  # type: ignore
