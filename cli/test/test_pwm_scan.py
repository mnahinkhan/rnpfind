"""
Tests the pwm scan module functions for correctness.

"""
import unittest
from test.read_fasta import read_fasta

from hgfind import hgfind
from pylcs import lcs2 as lcs
from src.rnpfind.pwm_scan import get_human_seq


class TestHumanSeq(unittest.TestCase):
    """
    Check if RNA sequences are extracted correctly
    """

    def test_forward_strand(self):
        """
        Checks whether sequence can be retrieved accurately on the forward
        strand

        """
        rna_info = hgfind("MALAT1")
        assert rna_info["strand"] == "+"
        seq = get_human_seq(rna_info)
        correct_seq = read_fasta("test/test_data/malat1_transcript.fasta")
        common = lcs(seq, correct_seq)
        ratio = common / min(len(seq), len(correct_seq))
        self.assertTrue(ratio > 0.95)

    def test_reverse_strand(self):
        """
        Checks whether sequence can be retrieved accurately on the reverse
        strand

        """
        rna_info = hgfind("HNRNPD")
        assert rna_info["strand"] == "-"
        seq = get_human_seq(rna_info)
        correct_seq = read_fasta("test/test_data/hnrnpd_genomic.fasta")
        common = lcs(seq, correct_seq)
        ratio = common / min(len(seq), len(correct_seq))
        self.assertTrue(ratio > 0.95)


if __name__ == "__main__":
    unittest.main()
