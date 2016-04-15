package eu.nomad_lab.parsers

import eu.nomad_lab.{ parsers, DefaultPythonInterpreter }
import org.scalacheck.Properties
import org.specs2.mutable.Specification
import org.{ json4s => jn }

object VaspParserSpec extends Specification {
  "VaspParserTest" >> {
    "test with Al.out" >> {
      "test with json-events" >> {
        ParserRun.parse(VaspParser, "test/examples/oqmd/relaxation/OUTCAR", "json-events") must_== ParseResult.ParseSuccess
      }
      "test with json" >> {
        ParserRun.parse(VaspParser, "test/examples/oqmd/relaxation/OUTCAR", "json") must_== ParseResult.ParseSuccess
      }
    }
  }
}

object VaspRunParserSpec extends Specification {
  "VaspRunParserTest" >> {
    "test with Al.out" >> {
      "test with json-events" >> {
        ParserRun.parse(VaspParser, "test/examples/oqmd/relaxation/vasprun.xml", "json-events") must_== ParseResult.ParseSuccess
      }
      "test with json" >> {
        ParserRun.parse(VaspParser, "test/examples/oqmd/relaxation/vasprun.xml", "json") must_== ParseResult.ParseSuccess
      }
    }
  }
}
