from builtins import object
import setup_paths
from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import os, sys, json, logging
import numpy as np

# description of the input
mainFileDescription = SM(
    name = 'root',
    weak = True,
    startReStr = "",
    subMatchers = [
        SM(name = 'newRun',
           startReStr = r"^\s*vasp.(?:[0-9.]+)\s+(?:[0-9]+[A-Za-z]+[0-9]+)\s+\(build (:[^)]+)\)\s+complex\s*",
           repeats = True,
           required = True,
           forwardMatch = True,
           sections   = ['section_run','section_method'],
           subMatchers = [
               SM(name = 'header',
                  startReStr = r"^\s*vasp.(?P<program_version>[0-9.]+)\s+(?P<vasp_src_date>[0-9]+[A-Za-z]+[0-9]+)\s+\(build (?P<vasp_build_date>[^)]+)\)\s+complex\s*",
                  subMatchers = [
                      ]
              ),
           ])
])

# loading metadata from nomad-meta-info/meta_info/nomad_meta_info/gaussian.nomadmetainfo.json
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/vasp.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)

parserInfo = {
  "name": "parser_vasp",
  "version": "1.0"
}

class VaspParserContext(object):
    """main place to keep the parser status, open ancillary files,..."""
    def __init__(self):
        pass

    def startedParsing(self, path, parser):
        self.parser = parser

# which values to cache or forward (mapping meta name -> CachingLevel)
cachingLevelForMetaName = {
}

if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = VaspParserContext())
