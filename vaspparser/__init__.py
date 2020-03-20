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

from nomadcore.baseclasses import ParserInterface
import nomadcore.baseclasses

from vaspparser.parser_vasprun import parserInfo
from vaspparser.parser_vasprun import VasprunContext, XmlParser, parserInfo
from vaspparser.parser_outcar import VaspOutcarParser

from nomad.parsing.legacy import CoEInterfaceParser

class VASPRunMainParser:
    """
    The main parser class that is called for all run types. Parses the VASP
    .xml output files.
    """
    def __init__(self, parser_context):
        self.parser_context = parser_context

    def parse(self, filepath):
        # the nomadcore.baseclasses.logger is set for each parsing run
        superContext = VasprunContext(logger=nomadcore.baseclasses.logger)
        parser = XmlParser(parserInfo, superContext)
        backend = self.parser_context.super_backend

        open_file = open
        if filepath.endswith('.gz'):
            open_file = gzip.open
        elif filepath.endswith('.bz2'):
            open_file = bz2.open
        elif filepath.endswith('.xz'):
            open_file = lzma.open

        parser.parse(os.path.abspath(filepath), open_file(filepath, 'rt'), backend)


class VASPRunParserInterface(ParserInterface):
    """
    This class handles the initial setup before any parsing can happen. It
    determines which version of BigDFT was used to generate the output and then
    sets up a correct main parser.

    After the implementation has been setup, you can parse the files with
    parse().
    """
    def __init__(
            self,
            metainfo_to_keep=None, backend=None, default_units=None,
            metainfo_units=None, debug=True, log_level=logging.ERROR, store=True):

        super(VASPRunParserInterface, self).__init__(
            metainfo_to_keep, backend, default_units, metainfo_units, debug, log_level, store)

    def setup_version(self):
        """
        Setups the version by looking at the output file and the version
        specified in it.
        """
        # Setup the root folder to the fileservice that is used to access files
        dirpath, filename = os.path.split(self.parser_context.main_file)
        dirpath = os.path.abspath(dirpath)
        self.parser_context.file_service.setup_root_folder(dirpath)
        self.parser_context.file_service.set_file_id(filename, "output")
        self.main_parser = VASPRunMainParser(self.parser_context)

    def get_metainfo_filename(self):
        return "vasp.nomadmetainfo.json"

    def get_parser_info(self):
        return parserInfo


class VASPRunParser(CoEInterfaceParser):

    def __init__(self):
        super().__init__(VASPRunParserInterface)
