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

from __future__ import absolute_import
import logging
import os
from vaspparser.parser_vasprun import VasprunContext, XmlParser, parserInfo
LOGGER = logging.getLogger("nomad")


class VASPMainParser(object):
    """The main parser class that is called for all run types. Parses the VASP
    XML output files.
    """
    def __init__(self, parser_context):
        """
        """
        self.parser_context = parser_context

    def parse(self, filepath):
        superContext = VasprunContext()
        parser = XmlParser(parserInfo, superContext)
        backend = self.parser_context.super_backend
        parser.parse(os.path.abspath(filepath), open(filepath), backend)
