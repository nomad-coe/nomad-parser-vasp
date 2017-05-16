import sys
import os
baseDir = os.path.dirname(os.path.abspath(__file__))
commonDir = os.path.normpath(os.path.join(baseDir, "../parser/parser-vasp"))

if commonDir not in sys.path:
    sys.path.insert(0, commonDir)
