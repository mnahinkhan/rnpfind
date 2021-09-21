"""
Tests the bind analysis (BindingSites) module functions for correctness.

"""
import copy
import unittest

from src.rnpfind.binding_analysis_binding_sites import BindingSites


def list_equal(al1, al2):
    """
    Check if two lists are equal
    """
    if len(al1) != len(al2):
        return False
    cal1 = copy.deepcopy(al1)
    cal2 = copy.deepcopy(al2)
    for eal1 in cal1:
        if eal1 not in cal2:
            return False
        cal2.remove(eal1)
    return True


class TestBindingSites(unittest.TestCase):
    """
    Check if BindingSites class functions correctly

    BindingSites is implemented in two modes.
    These tests will test both overlap_mode ON and OFF.
    """

    def test_list_equal(self):
        """
        Check that list_equal helper functions correctly
        """
        self.assertTrue(list_equal([1, 2, 3, 4], [1, 2, 3, 4]))
        self.assertTrue(list_equal([4, 2, 3, 1], [1, 2, 3, 4]))
        self.assertFalse(list_equal([4, 2, 3, 4], [1, 2, 3, 4]))

    def test_is_overlap_ranges(self):
        """Check that is_overlap_ranges() functions correctly

        The ranges are assumed to be half open intervals as follows:
        [start, end)
        """
        func = BindingSites.is_overlap_ranges
        self.assertTrue(func((50, 70), (60, 90)))
        self.assertTrue(func((50, 70), (50, 90)))

        self.assertFalse(func((50, 70), (70, 90)))
        self.assertFalse(func((20, 100), (400, 600)))

        self.assertTrue(func((60, 90), (50, 70)))
        self.assertTrue(func((50, 90), (50, 70)))

        self.assertFalse(func((70, 90), (50, 70)))
        self.assertFalse(func((400, 600), (20, 100)))

        self.assertTrue(func((50, 100), (50, 100)))
        self.assertTrue(func((50, 100), (99, 99)))


if __name__ == "__main__":
    unittest.main()
