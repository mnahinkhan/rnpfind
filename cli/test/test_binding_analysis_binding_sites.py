"""
Tests the bind analysis (BindingSites) module functions for correctness.

Remember that binding sites are made of (start coord, end coord, annotation)

"""
import copy
import random
import string
import unittest

from src.rnpfind.binding_analysis_binding_sites import BindingSites


def sta(tup):
    """Strip annotation from tuple"""
    return tup[0], tup[1]


def sta_list(list_of_tuples):
    """Strip annotation from every tuple in an iterable"""
    return [sta(x) for x in list_of_tuples]


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

    def test_merge_meta(self):
        """Check that _merge_meta() works correctly"""
        annotation_list = ["orange", "green", "yellow"]
        merged_list = BindingSites._merge_meta(annotation_list)
        self.assertEqual(set(merged_list), set(annotation_list))

        custom_func = "".join
        merged_list_custom_func = BindingSites._merge_meta(
            annotation_list, user_merge_func=custom_func
        )
        self.assertEqual(merged_list_custom_func, "".join(annotation_list))

    def test_collapse(self):
        """Check that _collapse() works correctly"""
        to_collapse = BindingSites([(280, 320), (300, 400), (350, 500)])
        collapsed_list = BindingSites._collapse(to_collapse)
        self.assertEqual(sta(collapsed_list), (280, 500))

    def test_add(self):
        """Check that the add() function works correctly"""

        sites = BindingSites()
        sites.add((80, 90))
        sites.add((90, 100))
        self.assertEqual(sta_list(sites), [(80, 90), (90, 100)])

        sites.add((95, 105))
        self.assertEqual(sta_list(sites), [(80, 90), (90, 105)])

        more_sites = BindingSites()
        more_sites.add((85, 95, "orange"))
        more_sites.add((88, 99, "blue"), user_merge_func="".join)
        self.assertIn(
            list(more_sites),
            [[(85, 99, "orangeblue")], [(85, 99, "blueorange")]],
        )

    def test_remove(self):
        """Check that site removal works correctly"""
        sites = BindingSites()
        sites.add((80, 90, "orange"))
        sites.add((90, 100, "blue"))

        sites.remove((80, 90, "orange"))
        self.assertEqual(list(sites), [(90, 100, "blue")])

    def test_dist(self):
        """Check if the correlation metric works correctly"""
        first_set = BindingSites()
        second_set = BindingSites()

        first_set.add((70, 90))
        second_set.add((70, 90))

        first_set.add((100, 120))
        second_set.add((100, 120))

        self.assertEqual(first_set.dist(second_set), 1)
        self.assertEqual(first_set.dist(first_set), 1)
        self.assertEqual(second_set.dist(first_set), 1)
        self.assertEqual(second_set.dist(second_set), 1)

        first_set = BindingSites()
        second_set = BindingSites()

        first_set.add((70, 90))
        first_set.add((120, 240))
        first_set.add((360, 480))

        second_set.add((95, 110))
        second_set.add((250, 310))
        second_set.add((500, 600))

        self.assertEqual(first_set.dist(second_set, bp_threshold=0), 0)

        first_set = BindingSites()
        second_set = BindingSites()

        first_set.add((70, 90))
        second_set.add((100, 200))
        nonzero_dist = first_set.dist(second_set, bp_threshold=30)

        self.assertTrue(1 > nonzero_dist > 0)

    def test_is_overlap(self):
        """Test is_overlap() functionality"""

        sites = BindingSites()
        sites.add((70, 90))
        sites.add((100, 110))
        sites.add((140, 200))
        sites.add((230, 260))

        self.assertTrue(sites.is_overlap((80, 83)))
        self.assertTrue(sites.is_overlap((90, 101)))
        self.assertFalse(sites.is_overlap((90, 100)))
        self.assertFalse(sites.is_overlap((120, 130)))

    def test_nearest_site(self):
        """Check that nearest_site() is working correctly"""

        sites = BindingSites()
        sites.add((70, 90))
        sites.add((100, 110))
        sites.add((140, 200))
        sites.add((230, 260))

        def saffe(tup):
            """Strip annotation from first element"""
            return (tup[0][0], tup[0][1]), tup[1]

        # Within
        self.assertEqual(
            saffe(sites.nearest_site((100, 105))), ((100, 110), 0)
        )
        # Closest from left
        self.assertEqual(
            saffe(sites.nearest_site((120, 123))), ((100, 110), 10)
        )
        # Closest from right
        self.assertEqual(
            saffe(sites.nearest_site((120, 133))), ((140, 200), 7)
        )
        # Equally close
        test_site = sites.nearest_site((210, 220))
        self.assertTrue(
            saffe(test_site) in (((140, 200), 10), ((230, 260), 10))
        )

    def test_distance(self):
        """Check that distance() works correctly"""

        dist_func = BindingSites.distance
        # Overlapping ranges have dist = 0
        self.assertEqual(dist_func((70, 90), (80, 100)), 0)
        # Sites right next to each other also have a dist = 0
        self.assertEqual(dist_func((70, 90), (90, 100)), 0)
        # Separate by 10
        self.assertEqual(dist_func((70, 90), (100, 100)), 10)
        self.assertEqual(dist_func((70, 90), (100, 120)), 10)
        self.assertEqual(dist_func((100, 100), (70, 90)), 10)
        self.assertEqual(dist_func((100, 120), (70, 90)), 10)

    def test_filter_overlap(self):
        """Check that filter_overlap() works correctly"""

        sites = BindingSites()
        sites.add((70, 90))
        sites.add((100, 110))
        sites.add((140, 200))
        sites.add((230, 260))

        self.assertEqual(sta_list(sites.filter_overlap((10, 40))), [])
        self.assertEqual(
            sta_list(sites.filter_overlap((10, 140))), [(70, 90), (100, 110)]
        )
        self.assertEqual(
            sta_list(sites.filter_overlap((10, 141))),
            [(70, 90), (100, 110), (140, 200)],
        )
        self.assertEqual(
            sta_list(sites.filter_overlap((200, 241))), [(230, 260)]
        )

    def test_print_bed(self):
        """Check that the BindingSites print_bed() function works correctly"""

        sites = BindingSites()
        sites.add((70, 90))
        sites.add((100, 110))
        sites.add((140, 200))
        sites.add((230, 260))
        sites.add((250, 280))

        site_name = "".join(
            random.choices(string.ascii_uppercase + string.digits, k=6)
        )

        bed_elements = sites.print_bed(name=site_name, chr_n="X").split("\n")
        bed_elements = [el for el in bed_elements if el]

        starts = [70, 100, 140, 230]
        ends = [90, 110, 200, 280]

        for i, elements in enumerate(bed_elements):
            chr_no, start, end, name = elements.split("\t")
            self.assertEqual(chr_no, "chrX")
            self.assertEqual(int(start), starts[i])
            self.assertEqual(int(end), ends[i])
            self.assertEqual(name, site_name)

    def test_return_depth(self):
        """Check that return_depth() works correctly"""
        sites = BindingSites(overlap_mode=True)
        sites.add((7, 9))
        sites.add((10, 11))
        sites.add((14, 20))
        sites.add((23, 26))
        sites.add((25, 28))
        sites.add((15, 18))
        self.assertEqual(
            sites.return_depth(),
            [0] * 7
            + [1] * 2
            + [0] * 1
            + [1] * 1
            + [0] * 3
            + [1] * 1
            + [2] * 3
            + [1] * 2
            + [0] * 3
            + [1] * 2
            + [2] * 1
            + [1] * 2,
        )

    def test_overlap_collapse(self):
        """Check that overlap_collapse() works correctly"""
        sites = BindingSites(overlap_mode=True)
        sites.add((7, 9))
        sites.add((10, 11))
        sites.add((14, 20))
        sites.add((23, 26))
        sites.add((25, 28))
        sites.add((15, 18))
        sites.add((16, 17))

        ONE_PLUS_DEPTH_SITES = [(7, 9), (10, 11), (14, 20), (23, 28)]
        TWO_PLUS_DEPTH_SITES = [(15, 18), (25, 26)]
        THREE_PLUS_DEPTH_SITES = [(16, 17)]
        FOUR_PLUS_DEPTH_SITES = []

        # BaseCoverNumber tests
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 100)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 14)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 13)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 4)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 3)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 1)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("baseCoverNumber", 0)
        self.assertEqual(sta_list(filtered_sites), FOUR_PLUS_DEPTH_SITES)

        # TopDepthRatio tests
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 1.0)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.9)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.8)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.7)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.6)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.5)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.4)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.3)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.2)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.1)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthRatio", 0.0)
        self.assertEqual(sta_list(filtered_sites), FOUR_PLUS_DEPTH_SITES)

        # TopDepthNumber tests
        filtered_sites = sites.overlap_collapse("TopDepthNumber", 3)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthNumber", 2)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthNumber", 1)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("TopDepthNumber", 0)
        self.assertEqual(sta_list(filtered_sites), FOUR_PLUS_DEPTH_SITES)

        # MinimumDepthNumber tests
        filtered_sites = sites.overlap_collapse("MinimumDepthNumber", 3)
        self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("MinimumDepthNumber", 2)
        self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        filtered_sites = sites.overlap_collapse("MinimumDepthNumber", 1)
        self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)

        # TopSitesNumber tests
        # These modes are still unimplemented apparently
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 7)
        # self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 6)
        # self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 5)
        # self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 4)
        # self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 3)
        # self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 2)
        # self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 1)
        # self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesNumber", 0)
        # self.assertEqual(sta_list(filtered_sites), FOUR_PLUS_DEPTH_SITES)

        # TopSitesRatio tests
        # These modes are still unimplemented apparently
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 7 / 7)
        # self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 6 / 7)
        # self.assertEqual(sta_list(filtered_sites), ONE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 5 / 7)
        # self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 4 / 7)
        # self.assertEqual(sta_list(filtered_sites), TWO_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 3 / 7)
        # self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 2 / 7)
        # self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 1 / 7)
        # self.assertEqual(sta_list(filtered_sites), THREE_PLUS_DEPTH_SITES)
        # filtered_sites = sites.overlap_collapse("TopSitesRatio", 0 / 7)
        # self.assertEqual(sta_list(filtered_sites), FOUR_PLUS_DEPTH_SITES)

        # Test the in_place functionality
        filtered_sites = sites.overlap_collapse(
            "MinimumDepthNumber", 2, in_place=True
        )
        self.assertIsNone(filtered_sites)
        self.assertEqual(sta_list(sites), TWO_PLUS_DEPTH_SITES)

    def test_base_cover(self):
        """Check that base_cover() functions correctly"""
        for mode in [True, False]:
            sites = BindingSites(overlap_mode=mode)
            sites.add((7, 9))
            sites.add((10, 11))
            sites.add((14, 20))
            sites.add((23, 26))
            sites.add((25, 28))
            sites.add((15, 18))
            sites.add((16, 17))

            self.assertEqual(sites.base_cover(), 14)

    def test_print_wig(self):
        """Test that print_wig() functions correctly"""
        sites = BindingSites(overlap_mode=True)
        sites.add((7, 9))
        sites.add((10, 11))
        sites.add((14, 20))
        sites.add((23, 26))
        sites.add((25, 28))
        sites.add((15, 18))
        sites.add((16, 17))

        wig = sites.print_wig()
        wig = wig.split("\n")
        # The two lines above, we won't test (although important!)
        wig = wig[2:]
        self.assertEqual(
            [int(x) for x in wig],
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                0,
                1,
                0,
                0,
                0,
                1,
                2,
                3,
                2,
                1,
                1,
                0,
                0,
                0,
                1,
                1,
                2,
                1,
                1,
            ],
        )
