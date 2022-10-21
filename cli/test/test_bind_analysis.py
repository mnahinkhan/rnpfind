"""
Tests the bind analysis module functions for correctness.

"""
import copy
import unittest

from src.rnpfind.bind_analysis import Storage
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


class TestStorage(unittest.TestCase):
    """
    Check if Storage class functions correctly
    """

    def test_list_equal(self):
        """
        Check that list_equal helper functions correctly
        """
        self.assertTrue(list_equal([1, 2, 3, 4], [1, 2, 3, 4]))
        self.assertTrue(list_equal([4, 2, 3, 1], [1, 2, 3, 4]))
        self.assertFalse(list_equal([4, 2, 3, 4], [1, 2, 3, 4]))

    def test_get_rbps(self):
        """
        Check that RBP retrieval from storage instance works correctly

        """
        storage = Storage()
        self.assertEqual(list(storage.get_rbps()), [])
        storage["rbp1"] = BindingSites([(50, 100)])
        self.assertEqual(list(s.lower() for s in storage.get_rbps()), ["rbp1"])
        storage["rbp2"] = BindingSites([(100, 200)])
        self.assertTrue(
            list_equal(
                list(s.lower() for s in storage.get_rbps()), ["rbp1", "rbp2"]
            )
        )

    def test_items(self):
        """Check that items() is correct"""
        storage = Storage()
        self.assertEqual(list(storage.items()), [])
        storage["rbp1"] = BindingSites([(50, 100)])
        items = list(storage.items())
        rbps = [k for k, v in items]
        list_of_sites = [v for k, v in items]
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(rbps, ["RBP1"])
        self.assertEqual(list_of_sites, [[(50, 100)]])
        storage["rbp2"] = BindingSites([(120, 180)])
        items = list(storage.items())
        rbps = [k for k, v in items]
        list_of_sites = [v for k, v in items]
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertTrue(list_equal(rbps, ["RBP1", "RBP2"]))
        self.assertTrue(list_equal(list_of_sites, [[(50, 100)], [(120, 180)]]))

    def test_values(self):
        """Check that values() is correct"""
        storage = Storage()
        self.assertEqual(list(storage.values()), [])
        storage["rbp1"] = BindingSites([(50, 100)])
        list_of_sites = storage.values()
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(list_of_sites, [[(50, 100)]])
        storage["rbp2"] = BindingSites([(120, 180)])
        list_of_sites = storage.values()
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertTrue(list_equal(list_of_sites, [[(50, 100)], [(120, 180)]]))

    def test_get_item(self):
        """Check that __getitem__() works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200)])
        sites = storage["rbp1"]
        sites = [(s, e) for s, e, a in sites]
        self.assertEqual(sites, [(100, 200)])
        storage["rbp2"] = BindingSites([(300, 400)])
        storage["rbp3"] = BindingSites([(500, 600)])

        filtered_storage = storage[["rbp2", "rbp3"]]
        items = list(filtered_storage.items())
        rbps = [k for k, v in items]
        list_of_sites = [v for k, v in items]
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertTrue(list_equal(rbps, ["RBP2", "RBP3"]))
        self.assertTrue(
            list_equal(list_of_sites, [[(300, 400)], [(500, 600)]])
        )

    def test_corr_reset(self):
        """Check if corr_reset functions correctly"""
        self.assertTrue(True)
        # TODO: Do proper testing

    def test_summary(self):
        """
        Check that summary() returns correct info.
        It is supposed to return the number of RBPs bound to the RNA Storage,
        and the total number of binding sites
        """
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp3"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        summary = storage.summary()
        num_rbps, num_sites = summary
        self.assertEqual(num_rbps, 3)
        self.assertEqual(num_sites, 9)

    def test_self_analysis(self):
        """
        Check if the self correlation analysis is done correctly
        """
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(300, 400), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        sym_cor_table = storage.self_analysis(progress_feedback=False)

        for rbp_i in ["RBP1", "RBP2", "RBP3"]:
            for rbp_j in ["RBP1", "RBP2", "RBP3"]:
                cor = sym_cor_table[rbp_i][rbp_j]
                op_cor = sym_cor_table[rbp_j][rbp_i]
                self.assertEqual(cor, op_cor)
                if rbp_i == rbp_j:
                    self.assertEqual(cor, 1)

        cor12 = sym_cor_table["RBP1"]["RBP2"]
        cor23 = sym_cor_table["RBP2"]["RBP3"]
        cor13 = sym_cor_table["RBP1"]["RBP3"]

        def f_score(val_1: int, val_2: int):
            return 2 * val_1 * val_2 / (val_1 + val_2)

        self.assertEqual(cor12, f_score(1, 2 / 3))
        self.assertEqual(cor23, f_score(1, 1 / 2))
        self.assertEqual(cor13, f_score(1, 1 / 3))

    def test_bind_near(self):
        """Check that the binds_near() function works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        rbps = storage.binds_near((400, 500))
        self.assertTrue(list_equal(rbps, ["RBP1", "RBP2"]))

    def test_filter(self):
        """Check that the filter function works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        filtered_storage = storage.filter(lambda rbp: rbp == "RBP1")

        rbps = filtered_storage.get_rbps()
        list_of_sites = filtered_storage.values()
        rbps = list(rbps)
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(rbps, ["RBP1"])
        self.assertEqual(list_of_sites, [[(100, 200), (300, 400), (600, 700)]])

    def test_sites_analysis(self):
        """Check that sites_analysis() function works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        sites2storage = storage.sites_analysis("rbp1")
        sites2storage = {k[0:2]: v for k, v in sites2storage.items()}

        sites_100_200 = sites2storage[(100, 200)]
        sites_300_400 = sites2storage[(300, 400)]
        sites_600_700 = sites2storage[(600, 700)]

        # (100, 200)
        rbps = list(sites_100_200.get_rbps())
        list_of_sites = sites_100_200.values()
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(rbps, ["RBP1"])
        self.assertEqual(list_of_sites, [[(100, 200)]])

        # (300, 400)
        rbps = list(sites_300_400.get_rbps())
        list_of_sites = sites_300_400.values()
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(rbps, ["RBP1", "RBP2"])
        self.assertEqual(list_of_sites, [[(300, 400)], [(200, 450)]])

        # (600, 700)
        rbps = list(sites_600_700.get_rbps())
        list_of_sites = sites_600_700.values()
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(rbps, ["RBP1", "RBP2", "RBP3"])
        self.assertEqual(
            list_of_sites, [[(600, 700)], [(600, 700)], [(600, 700)]]
        )

    def test_all_sites_in(self):
        """Check that all_sites_in() works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        filtered_storage = storage.all_sites_in((300, 500))
        rbps = list(filtered_storage.get_rbps())
        list_of_sites = filtered_storage.values()
        list_of_sites = [
            [(s, e) for s, e, a in site] for site in list_of_sites
        ]
        self.assertEqual(rbps, ["RBP1", "RBP2"])
        self.assertEqual(list_of_sites, [[(300, 400)], [(200, 450)]])

    def test_print_bed(self):
        """Check that print_bed() works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        bed_string = storage.print_bed(chr_n="X", displacement=600)
        bed_lines = bed_string.split("\n")
        bed_lines = [line for line in bed_lines if line]

        self.assertEqual(len(bed_lines), 6)
        self.assertTrue("chrX\t700\t800\tRBP1" in bed_lines)
        self.assertTrue("chrX\t900\t1000\tRBP1" in bed_lines)
        self.assertTrue("chrX\t1200\t1300\tRBP1" in bed_lines)
        self.assertTrue("chrX\t800\t1050\tRBP2" in bed_lines)
        self.assertTrue("chrX\t1200\t1300\tRBP2" in bed_lines)
        self.assertTrue("chrX\t1200\t1300\tRBP3" in bed_lines)

    def test_sum_over_all(self):
        """Check that sum_over_all() works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])

        sites = [(s, e) for s, e, a in storage.sum_over_all()]
        self.assertTrue(
            list_equal(
                sites,
                [
                    (100, 200),
                    (300, 400),
                    (600, 700),
                    (200, 450),
                    (600, 700),
                    (600, 700),
                ],
            )
        )

    def test_print_wig(self):
        """Check that print_wig() works correctly"""
        storage = Storage()
        storage["rbp1"] = BindingSites([(100, 200), (300, 400), (600, 700)])
        storage["rbp2"] = BindingSites([(200, 450), (600, 700)])
        storage["rbp3"] = BindingSites([(600, 700)])
        depth_array = storage.print_wig(chr_no=13, include_header=False).split(
            "\n"
        )
        self.assertTrue("chr13" in depth_array[0])
        depth_array = depth_array[1:]
        ideal_depth_array = (
            [0] * 100  # 0 - 99
            + [1] * 100  # 100 - 199
            + [1] * 100  # 200 - 299
            + [2] * 100  # 300 - 399
            + [1] * 50  # 400 - 449
            + [0] * 50  # 450 - 499
            + [0] * 100  # 500 - 599
            + [3] * 100  # 600 - 699
        )
        depth_array = [int(c) for c in depth_array]
        self.assertEqual(len(depth_array), 700)
        self.assertEqual(
            ideal_depth_array, depth_array[0 : len(ideal_depth_array)]
        )

        # Now test antisense
        depth_array = storage.print_wig(
            chr_no=13, include_header=False, antisense=True
        ).split("\n")
        self.assertTrue("chr13" in depth_array[0])
        depth_array = depth_array[1:]
        ideal_depth_array = (
            [3] * 100  # 600 - 699
            + [0] * 100  # 500 - 599
            + [0] * 50  # 450 - 499
            + [1] * 50  # 400 - 449
            + [2] * 100  # 300 - 399
            + [1] * 100  # 200 - 299
            + [1] * 100  # 100 - 199
            + [0] * 100  # 0 - 99
        )
        depth_array = [int(c) for c in depth_array]
        self.assertEqual(len(depth_array), 700)
        self.assertEqual(
            ideal_depth_array, depth_array[0 : len(ideal_depth_array)]
        )


if __name__ == "__main__":
    unittest.main()
