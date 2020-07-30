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

import os
import logging
import gzip
import bz2
import lzma

from nomad.datamodel import EntryArchive

from .metainfo import m_env
from vaspparser.parser_vasprun import VasprunContext, XmlParser, parser_info
from nomad.parsing.parser import FairdiParser
from vaspparser.parser_outcar import VaspOutcarParser


class VASPParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/vasp', code_name='VASP', code_homepage='https://www.vasp.at/',
            mainfile_mime_re=r'(application/.*)|(text/.*)',
            mainfile_contents_re=(
                r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
                r'?\s*<modeling>'
                r'?\s*<generator>'
                r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
                r'?'),
            supported_compressions=['gz', 'bz2', 'xz']
        )

    def parse(self, filepath, archive, logger=None):
        self._metainfo_env = m_env

        super_context = VasprunContext(logger=logger)

        parser = XmlParser(parser_info, super_context, metainfo_env=m_env)

        open_file = open
        if filepath.endswith('.gz'):
            open_file = gzip.open
        elif filepath.endswith('.bz2'):
            open_file = bz2.open
        elif filepath.endswith('.xz'):
            open_file = lzma.open

        parser.parse(os.path.abspath(filepath), open_file(filepath, 'rt'), archive)
