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
import unittest
import logging

from vaspparser import VASPRunParserInterface


def get_result(folder, metaname=None):
    """Get the results from the calculation in the given folder. By default goes through different

    Args:
        folder: The folder relative to the directory of this script where the
            parsed calculation resides.
        metaname(str): Optional quantity to return. If not specified, returns
            the full dictionary of results.
    """
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, folder, "vasprun.xml.static")
    parser = VASPRunParserInterface(None, debug=True, log_level=logging.CRITICAL)
    results = parser.parse(filename)

    if metaname is None:
        return results
    else:
        return results[metaname]


class TestDOS(unittest.TestCase):
    """Tests that the parser can handle single point calculations.
    """
    @classmethod
    def setUpClass(cls):
        cls.results = get_result(os.path.join("examples", "dos"))

    def test_program_name(self):
        result = self.results["program_name"]
        self.assertEqual(result, "VASP")

    def test_program_basis_set_type(self):
        result = self.results["program_basis_set_type"]
        self.assertEqual(result, "plane waves")

    def test_electronic_structure_method(self):
        result = self.results["electronic_structure_method"]
        self.assertEqual(result, "DFT")

    def test_dos_energies_and_values(self):
        """Test that all the DOS related values are found and that the
        dimensions match.
        """
        dos_values = self.results["dos_values"]
        dos_energies = self.results["dos_energies"]
        dos_energies_normalized = self.results["dos_energies_normalized"]
        dos_integrated_values = self.results["dos_integrated_values"]

        self.assertEqual(len(dos_energies.shape), 1)
        self.assertEqual(len(dos_values.shape), 2)
        self.assertEqual(len(dos_integrated_values.shape), 2)
        self.assertEqual(dos_energies.shape[0], dos_values.shape[1])
        self.assertEqual(dos_energies_normalized.shape[0], dos_values.shape[1])


if __name__ == '__main__':
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDOS))

    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=0).run(alltests)
