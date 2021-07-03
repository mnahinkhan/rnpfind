"""
This module is dedicated to functions that help with converting a string
representing a gene name to its coordinates on the hg38 chromosome.

"""

from typing import Union


def is_int(in_obj):
    """
    Checks if the input represents an integer and returns true iff so

    """
    try:
        int(in_obj)
        return True
    except ValueError:
        return False


class Chromosome:
    """
    This class allows one to specify chromsomes from the human genome.
    Benefits include comparison of chromosomes that are sex chromosomes,
    autosomes, or genes included in the mitochondrial DNA.

    """

    def __init__(self, n: Union[int, str]):

        if isinstance(n, int) and 1 <= n <= 22:
            self.chr_n = n
        elif isinstance(n, str) and n.upper() in ("X", "Y", "MT", "M"):
            if n.upper() == "X":
                self.chr_n = 23
            elif n.upper() == "Y":
                self.chr_n = 24
            elif n.upper() == "MT":
                self.chr_n = 25
            elif n.upper() == "M":
                self.chr_n = 25
        elif isinstance(n, str) and is_int(n) and 1 <= int(n) <= 22:
            self.chr_n = int(n)
        else:
            raise ValueError("expected X, Y, M(T), or [1-22]")

    def __str__(self):
        return (
            str(self.chr_n)
            if self.chr_n <= 22
            else "X"
            if self.chr_n == 23
            else "Y"
            if self.chr_n == 24
            else "M"
        )

    def __lt__(self, other):
        return self.chr_n < other.chr_n

    def __eq__(self, other):
        return self.chr_n == other.chr_n

    def __hash__(self):
        return hash(self.chr_n)

    def __gt__(self, other):
        return self.chr_n > other.chr_n

    def __repr__(self):
        return self.__str__()

    def __int__(self):
        return self.chr_n


if __name__ == "__main__":
    pass
