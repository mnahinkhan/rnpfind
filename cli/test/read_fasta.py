#!/usr/bin/env python3
"""
Read FASTA file
"""


def read_fasta(filepath: str):
    seq = ""
    with open(filepath) as file_handle:
        for line in file_handle:
            if not line.strip() or line.strip()[0] == ">":
                continue
            seq += line.strip()
    return seq


if __name__ == "__main__":
    print(read_fasta("test_data/malat1_transcript.fasta"))
