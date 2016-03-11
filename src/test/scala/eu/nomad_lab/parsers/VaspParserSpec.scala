package eu.nomad_lab.parsers

import eu.nomad_lab.{parsers, DefaultPythonInterpreter}
import org.scalacheck.Properties
import org.specs2.mutable.Specification
import org.{json4s => jn}


object VaspParserSpec extends Specification {
  "VaspParserTest" >> {
    "test with Al.out">> {
      "test with json-events" >> {
        ParserRun.parse(VaspParser,"test/examples/oqmd/relaxation/0_convergence/OUTCAR","json-events") must_== ParseResult.ParseSuccess
      }
      "test with json" >> {
        ParserRun.parse(VaspParser,"test/examples/oqmd/relaxation/0_convergence/OUTCAR","json") must_== ParseResult.ParseSuccess
      }
    }
  }
}
