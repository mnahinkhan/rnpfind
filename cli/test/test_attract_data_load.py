"""
Tests the Attract data loading functions for correctness

"""
import unittest

from src.rnpfind.analysis_functions import (
    analysis_method_functions,
    analysis_methods_supported_long,
    analysis_methods_supported_short,
)

# from .attract_data_load import attract_data_load
from src.rnpfind.gene_coordinates import Chromosome


class TestCorrectBindingSites(unittest.TestCase):
    """
    Tests to ensure that the ATTRACT data load function loads binding sites
    that are consistent with that of the web interface.

    """

    def test_correct_sites(self):
        """
        Checks that the binding sites match with web interface.

        """

        test_rna_info = {}
        test_rna_info["official_name"] = "Malat1"
        test_rna_info["chr_n"] = Chromosome(11)
        test_rna_info["start_coord"] = 65497738
        test_rna_info["end_coord"] = 65506516

        # attract_binding_sites = attract_data_load(test_rna_info)
        self.assertEqual(
            len(analysis_methods_supported_short),
            len(analysis_methods_supported_long),
        )
        self.assertEqual(
            len(analysis_method_functions),
            len(analysis_methods_supported_short),
        )


if __name__ == "__main__":
    unittest.main()
