package eu.nomad_lab.parsers
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import eu.{ nomad_lab => lab }
import scala.collection.breakOut

object VaspParser extends SimpleExternalParserGenerator(
  name = "VaspParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("VaspParser")) ::
      ("parserId" -> jn.JString("VaspParser" + lab.VaspVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.VaspVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """^\s*vasp.(?<version>[0-9.]+)\s+(?<srcDate>[0-9]+[A-Za-z]+[0-9]+)\s+\(build (?<buildDate>[^)]+)\)\s+complex\s*
""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/vasp/parser/parser-vasp/vaspparser/parser_vasp.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-vasp/vasprun/parser_vasp.py",
    "parser-vasp/vasprun/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/vasp.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-vasp" -> "parsers/vasp/parser/parser-vasp",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info",
    "python" -> "python-common/common/python/nomadcore"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)

object VaspRunParser extends SimpleExternalParserGenerator(
  name = "VaspRunParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("VaspRunParser")) ::
      ("parserId" -> jn.JString("VaspRunParser" + lab.VaspVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.VaspVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("application/xml"),
  mainFileRe = """\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*
?\s*<modeling>
?\s*<generator>
?\s*<i name="program" type="string">\s*vasp\s*</i>
?""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/vasp/parser/parser-vasp/vaspparser/scalainterface.py",
    "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-vasp/vasprun/__init__.py",
    "parser-vasp/vasprun/parser.py",
    "parser-vasp/vasprun/vaspmainparser.py",
    "parser-vasp/vasprun/parser_vasprun.py",
    "parser-vasp/vasprun/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/vasp.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-vasp" -> "parsers/vasp/parser/parser-vasp",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info",
    "python" -> "python-common/common/python/nomadcore"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
