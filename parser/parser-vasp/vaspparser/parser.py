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
from nomadcore.baseclasses import ParserInterface
from vaspparser.parser_vasprun import parserInfo
from vaspparser.vaspmainparser import VASPMainParser
logger = logging.getLogger("nomad")


class VASPParser(ParserInterface):
    """This class handles the initial setup before any parsing can happen. It
    determines which version of BigDFT was used to generate the output and then
    sets up a correct main parser.

    After the implementation has been setup, you can parse the files with
    parse().
    """
    def __init__(self, metainfo_to_keep=None, backend=None, default_units=None, metainfo_units=None, debug=True, log_level=logging.ERROR, store=True):
        super(VASPParser, self).__init__(metainfo_to_keep, backend, default_units, metainfo_units, debug, log_level, store)

    def setup_version(self):
        """Setups the version by looking at the output file and the version
        specified in it.
        """
        # Setup the root folder to the fileservice that is used to access files
        dirpath, filename = os.path.split(self.parser_context.main_file)
        dirpath = os.path.abspath(dirpath)
        self.parser_context.file_service.setup_root_folder(dirpath)
        self.parser_context.file_service.set_file_id(filename, "output")
        self.main_parser = VASPMainParser(self.parser_context)

    def get_metainfo_filename(self):
        return "vasp.nomadmetainfo.json"

    def get_parser_info(self):
        return parserInfo
