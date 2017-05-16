"""
This is a module for unit testing the BigDFT parser. The unit tests are run with
a custom backend that outputs the results directly into native python object for
easier and faster analysis.

Each property that has an enumerable list of different possible options is
assigned a new test class, that should ideally test through all the options.

The properties that can have non-enumerable values will be tested only for one
specific case inside a test class that is designed for a certain type of run
(MD, optimization, QM/MM, etc.)
"""
import os
import unittest
import logging
import setup_paths
from parser import VASPParser


#===============================================================================
def get_results(rel_path):
    """Get the given result from the calculation in the given folder by using
    the Analyzer in the nomadtoolkit package. Tries to optimize the parsing by
    giving the metainfo_to_keep argument.

    Args:
        rel_path: The path to the output file. Relative to the directory of
            this script.
        metaname: The quantity to extract.
    """
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, rel_path)
    parser = VASPParser(filename, None, debug=True, log_level=logging.WARNING)
    results = parser.parse()
    return results


#===============================================================================
def get_result(folder, metaname, optimize=True):
    if optimize:
        results = get_results(folder, None)
    else:
        results = get_results(folder)
    result = results[metaname]
    return result


#===============================================================================
class TestDOS(unittest.TestCase):
    """Tests that the parser can handle single point calculations.
    """
    @classmethod
    def setUpClass(cls):
        cls.results = get_results("dos/vasprun.xml.static")

    def test_program_name(self):
        result = self.results["program_name"]
        self.assertEqual(result, "VASP")

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

#===============================================================================
if __name__ == '__main__':
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDOS))

    alltests = unittest.TestSuite(suites)
    unittest.TextTestRunner(verbosity=0).run(alltests)
