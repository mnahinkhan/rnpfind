"""
Tests the RNPFind analysis functions module for internal consistency.

"""
import inspect
import unittest

from src.rnpfind.analysis_functions import (
    analysis_method_functions,
    analysis_methods_supported_long,
    analysis_methods_supported_short,
)


class TestAnalysisMethods(unittest.TestCase):
    """
    Tests to ensure that the analysis method related variables are properly set
    up.

    """

    def test_equal_number_of_methods(self):
        """
        Checks that the number of analysis functions are equal across variables
        describing them.

        """
        self.assertEqual(
            len(analysis_methods_supported_short),
            len(analysis_methods_supported_long),
        )
        self.assertEqual(
            len(analysis_method_functions),
            len(analysis_methods_supported_short),
        )

    def test_appropriate_name(self):
        """
        Checks that the analysis functions are given appropriate 'short' names.
        """
        for short_name, long_name in zip(
            analysis_methods_supported_short, analysis_methods_supported_long
        ):
            for token in short_name.split("_"):
                self.assertIn(token.lower(), long_name.lower())

    def test_method_dict_keys(self):
        """
        Checks whether analysis_method_functions describes the functions listed
        in analysis_methods_supported_short.

        """

        # equality of size from a previous test along with injection here
        # assures bijection.
        for method in analysis_methods_supported_short:
            self.assertIn(method, analysis_method_functions)

    def test_appropriate_method_type(self):
        """
        Checks that each of the analysis methods have a consistent type
        signature (only enforced for input; return type is different depending
        on analysis type)

        """

        for method in analysis_method_functions.values():
            self.assertEqual(
                list(inspect.signature(method).parameters)[0], "big_storage"
            )

            self.assertEqual(
                list(inspect.signature(method).parameters)[1], "rna_info"
            )


if __name__ == "__main__":
    unittest.main()
