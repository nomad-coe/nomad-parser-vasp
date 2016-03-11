package eu.nomad_lab.parsers
import eu.nomad_lab.DefaultPythonInterpreter
import org.{json4s => jn}
import eu.{nomad_lab => lab}

object VaspParser extends SimpleExternalParserGenerator(
      name = "VaspParser",
      parserInfo = jn.JObject(
        ("name" -> jn.JString("VaspParser")) ::
          ("parserId" -> jn.JString("VaspParser" + lab.VaspVersionInfo.version)) ::
          ("versionInfo" -> jn.JObject(
            ("nomadCoreVersion" -> jn.JString(lab.NomadCoreVersionInfo.version)) ::
              (lab.VaspVersionInfo.toMap.map{ case (key, value) =>
                (key -> jn.JString(value.toString))
              }(breakOut): List[(String, jn.JString)])
          )) :: Nil
      ),
      mainFileTypes = Seq("text/.*"),
      mainFileRe = """^\s*vasp.(?<version>[0-9.]+)\s+(?<srcDate>[0-9]+[A-Za-z]+[0-9]+)\s+\(build (?<buildDate>[^)]+)\)\s+complex\s*
""".r,
      cmd = Seq(DefaultPythonInterpreter.python2Exe(), "${envDir}/parsers/vasp/parser/parser-vasp/parser_vasp.py",
        "--uri", "${mainFileUri}", "${mainFilePath}"),
      resList = Seq(
        "parser-vasp/parser_vasp.py",
        "parser-vasp/setup_paths.py",
        "nomad_meta_info/common.nomadmetainfo.json",
        "nomad_meta_info/meta_types.nomadmetainfo.json",
        "nomad_meta_info/vasp.nomadmetainfo.json"
      ) ++ DefaultPythonInterpreter.commonFiles(),
      dirMap = Map(
        "parser-vasp" -> "parsers/vasp/parser/parser-vasp",
        "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info",
        "python" -> "python-common/common/python/nomadcore") ++ DefaultPythonInterpreter.commonDirMapping()
)
