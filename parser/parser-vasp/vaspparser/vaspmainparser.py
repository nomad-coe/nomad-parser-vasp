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
