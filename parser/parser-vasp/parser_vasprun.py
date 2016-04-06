import xml.etree.ElementTree
import logging, sys

class VasprunContext(object):
    def __init__(self,backend):
        self.backend = backend

    sectionMap = {
        "model": ["section_run", "section_method"],
        "structure": ["section_system_description"],
        "calculation": ["single_configuration_calculation"]
    }

    def onEnd_scstep(self, parser, event, element):
        print(path, element)

class XmlParser(object):
    @staticmethod
    def extractCallbacks(obj):
        """extracts all callbacks from the object obj

        triggers should start with onStart_ or onEnd__ and then have a valid section name.
        They will be called with this object, the event and current element
        """
        triggers = {}
        for attr in dir(obj):
            if attr.startswith("onStart_"):
                triggers[attr] = getattr(obj, attr)
            elif attr.startswith("onEnd_"):
                triggers[attr] = getattr(obj, attr)
        return triggers


    def __init__(self, parserInfo, superContext, callbacks, sectionMap):
        self.fIn = fIn
        self.superContext = superContext
        self.callbacks = callbacks
        self.sectionMap = sectionMap
        self.path = []
        self.tagSections = {}

    def parse(mainFileUri, fIn, backend):
        if self.path:
            raise Exception("Parse of %s called with non empty path, parse already in progress?" % mainFileUri)
        self.mainFileUri = mainFileUri
        self.fIn = fIn
        self.backend = backend
        backend.startedParsingSession(
            mainFileUri = mainFileUri,
            parserInfo = self.parserInfo)
        try:
            for event, el in xml.etree.ElementTree.iterparse(self.fIn, events=["start","end"]):
                if event == 'start':
                    sectionsToOpen = self.sectionMap.get(el.tag, None)
                    if sectionsToOpen:
                        pathStr = "/".join(self.path) + "/" + el.tag
                        gIndexes = {}
                        for sect in sectionsToOpen:
                            gIndexs[sect] = backend.openSection(sect)
                        self.tagSections[pathStr] = gIndexes
                    callback = self.callbacks.get("onStart_" + el.tag, None)
                    if callback:
                        callback(self, event, el)
                    self.path.append(el)
                elif event == 'end':
                    lastEl = path.pop()
                    if lastEl != el:
                        raise Exception("mismatched path at end, got %s expected %s" % (lastEl, el))
                    tag = el.tag
                    callback = self.callbacks.get("onEnd_" + tag, None)
                    if callback:
                        if not callback(self, event, el):
                            el.clear()
                            if path:
                                path[-1].remove(el)
                    elif len(self.path) == 1:
                        backend.addValue("parsing_message_warning_info_run","Skipping level 1 tag %s" % tag)
                        el.clear()
                        path[-1].remove(el)
                    sectionsToClose = self.sectionMap.get(tag, None)
                    if sectionsToClose:
                        pathStr = "/".join(self.path) + "/" + tag
                        gIndexes = self.tagSections[pathStr]
                        del self.tagSections[pathStr]
                        for sect in reversed(sectionsToOpen):
                            backend.closeSection(sect, gIndexes[sect])
                        self.tagSections[pathStr] = gIndexes
                else:
                    raise Exception("Unexpected event %s" % event)
        except:
            logger.exception("failure when parsing %s", self.mainFileUri)
            backend.finishedParsingSession(
                parserStatus = "ParseFailure",
                parserErrors = ["exception: %s" % sys.exc_value]
            )
        else:
            backend.finishedParsingSession(
                parserStatus = "ParseSuccess",
                parserErrors = ["exception: %s" % sys.exc_value]
            )
