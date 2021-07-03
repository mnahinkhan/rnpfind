"""
Tests the user input functions module for correctness.

"""
import unittest

from src.rnpfind.gene_coordinates import Chromosome
from src.rnpfind.user_input import parse_genome_coord


class TestParseGenomeCoord(unittest.TestCase):
    """
    Tests to ensure that parsing of genome coordinates is correct

    Note that this does not test the ability of the parser to detect incorrect
    chromsome number based on size of chromosome

    """

    def test_normal_success(self):
        """
        Checks whether normal coordinates are validated correctly

        """

        case_1 = "3:400-14000"
        case_2 = "1:360-5000"
        case_3 = "22:3456-23456"
        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)

        self.assertTrue(result_1["success"])
        self.assertTrue(result_2["success"])
        self.assertTrue(result_3["success"])

        self.assertEqual(result_1["chr_n"], Chromosome(3))
        self.assertEqual(result_2["chr_n"], Chromosome(1))
        self.assertEqual(result_3["chr_n"], Chromosome(22))

        self.assertEqual(result_1["start_coord"], 400)
        self.assertEqual(result_2["start_coord"], 360)
        self.assertEqual(result_3["start_coord"], 3456)

        self.assertEqual(result_1["end_coord"], 14000)
        self.assertEqual(result_2["end_coord"], 5000)
        self.assertEqual(result_3["end_coord"], 23456)

    def test_invalid_chr_no(self):
        """
        Check whether genomic coordinate outside valid range fails

        """
        case_1 = "24:100-1000"
        case_2 = "0:130-2345"
        case_3 = "-1:2345-7653"
        case_4 = "23:2345-234546"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)
        result_4 = parse_genome_coord(case_4)

        self.assertFalse(result_1["success"])
        self.assertFalse(result_2["success"])
        self.assertFalse(result_3["success"])
        self.assertFalse(result_4["success"])

    def test_invalid_range(self):
        """
        Check that wrong coordinate range fails

        """

        case_1 = "2:5000-500"
        case_2 = "4:6000-5800"
        case_3 = "X:14000-4000"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)

        self.assertFalse(result_1["success"])
        self.assertFalse(result_2["success"])
        self.assertFalse(result_3["success"])

    def test_bad_start_coord(self):
        """
        Check that bad start coord results in fail-outut

        """
        case_1 = "2:abd-50000"
        case_2 = "4:cyz-580000"
        case_3 = "X:hgff-4000000"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)

        self.assertFalse(result_1["success"])
        self.assertFalse(result_2["success"])
        self.assertFalse(result_3["success"])

    def test_bad_end_coord(self):
        """
        Check that bad end coord results in fail-outut

        """
        case_1 = "2:500-xyz"
        case_2 = "4:4567-tac"
        case_3 = "X:617-rbp"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)

        self.assertFalse(result_1["success"])
        self.assertFalse(result_2["success"])
        self.assertFalse(result_3["success"])

    def test_negative_coord(self):
        """
        Check whether parsing negative coordinates fails

        """
        case_1 = "12:-234:-54"
        case_2 = "3:-234:54"
        case_3 = "X:450--600"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)

        self.assertFalse(result_1["success"])
        self.assertFalse(result_2["success"])
        self.assertFalse(result_3["success"])

    def test_letter_success(self):
        """
        Checks whether genomic coordinates with letters parse correctly

        """
        case_1 = "M:4000-15000"
        case_2 = "M:400000-15000000"
        case_3 = "MT:4000-15000"
        case_4 = "MT:400000-15000000"
        case_5 = "X:100-300"
        case_6 = "Y:300-4000"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)
        result_3 = parse_genome_coord(case_3)
        result_4 = parse_genome_coord(case_4)
        result_5 = parse_genome_coord(case_5)
        result_6 = parse_genome_coord(case_6)

        self.assertTrue(result_1["success"])
        self.assertTrue(result_2["success"])
        self.assertTrue(result_3["success"])
        self.assertTrue(result_4["success"])
        self.assertTrue(result_5["success"])
        self.assertTrue(result_6["success"])

        self.assertEqual(result_1["chr_n"], Chromosome("M"))
        self.assertEqual(result_2["chr_n"], Chromosome("M"))
        self.assertEqual(result_3["chr_n"], Chromosome("M"))
        self.assertEqual(result_4["chr_n"], Chromosome("M"))
        self.assertEqual(result_5["chr_n"], Chromosome("X"))
        self.assertEqual(result_6["chr_n"], Chromosome("Y"))

        self.assertEqual(result_1["start_coord"], 4000)
        self.assertEqual(result_2["start_coord"], 400000)
        self.assertEqual(result_3["start_coord"], 4000)
        self.assertEqual(result_4["start_coord"], 400000)
        self.assertEqual(result_5["start_coord"], 100)
        self.assertEqual(result_6["start_coord"], 300)

        self.assertEqual(result_1["end_coord"], 15000)
        self.assertEqual(result_2["end_coord"], 15000000)
        self.assertEqual(result_3["end_coord"], 15000)
        self.assertEqual(result_4["end_coord"], 15000000)
        self.assertEqual(result_5["end_coord"], 300)
        self.assertEqual(result_6["end_coord"], 4000)

    def test_gibberish(self):
        """
        Check that random gibberish fails

        """
        case_1 = "gwfefrgr"
        case_2 = "23:234:23:123"

        result_1 = parse_genome_coord(case_1)
        result_2 = parse_genome_coord(case_2)

        self.assertFalse(result_1["success"])
        self.assertFalse(result_2["success"])


if __name__ == "__main__":
    unittest.main()
