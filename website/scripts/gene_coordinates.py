# File dedicated to functions that help with converting a string representing a gene name to its coordinates on the
# hg38 chromosome.

import os
import pickle
from typing import Union, Dict

# import statistics

path_to_bio_mart = "./website/data/biomart-gene-coordinates.txt"
path_to_pickle = "./website/output-data/pickles/biomart-gene-coordinates.pickle"


class Chromosome:
    """ """
    def __init__(self, n: Union[int, str]):

        if type(n) == int and n <= 22:
            self.n = n
        elif n.upper() == "X":
            self.n = 23
        elif n.upper() == "Y":
            self.n = 24
        elif n.upper() == "MT":
            self.n = 25
        elif n.upper() == "M":
            self.n = 25
        else:
            raise ValueError("expected X, Y, M(T), or [1-22]")

    def __str__(self):
        return str(self.n) if self.n <= 22 else "X" if self.n == 23 else "Y" if self.n == 24 else "M"

    def __lt__(self, other):
        return self.n < other.n

    def __eq__(self, other):
        return self.n == other.n

    def __hash__(self):
        return hash(self.n)

    def __gt__(self, other):
        return self.n > other.n

    def __repr__(self):
        return self.__str__()

    def __int__(self):
        return self.n


def file_to_dicts(path):
    """

    :param path: 

    """
    # This function converts the BioMart file given by the path into two dictionaries: the first one takes a gene name
    # as a key and gives as value the official symbol for the gene. The second one then takes the official symbol as
    # key and gives back the chromosome number, the start coordinate, and the end coordinate as its value.

    nameToOfficial: Dict[str, str] = {}
    officialToCoord = {}

    f = open(path)
    f.readline()  # This first line just says "Gene start (bp), Gene end (bp), Gene name, Gene Synonym,
    # WikiGene name, UniProtKB Gene Name symbol, Chromosome/Scaffold Name"
    s = f.readline()  # To get to an actual line
    while s:
        start_coord, end_coord, gene_name, gene_syn_name, wiki_name, uniprot_name, str_chr_no = \
            [wss.strip() for wss in s.split('\t')]
        names = [name for name in [gene_name, gene_syn_name, wiki_name, uniprot_name] if name != ""]

        if str_chr_no == "X" or str_chr_no == "Y" or str_chr_no == "MT":
            chr_no = Chromosome(str_chr_no)
        elif len(str_chr_no) <= 2:
            chr_no = Chromosome(int(str_chr_no))
        elif "HSCHR" in str_chr_no.upper():
            str_chr_no = str_chr_no[str_chr_no.upper().find("HSCHR"):str_chr_no.upper().find("HSCHR") + 7]
            str_chr_no = str_chr_no[5:]
            if str_chr_no[1].isdigit():
                chr_no = Chromosome(int(str_chr_no))
            elif str_chr_no[0] == "X" or str_chr_no[0] == "Y":
                chr_no = Chromosome(str_chr_no[0])
            elif str_chr_no == "MT":
                chr_no = Chromosome(str_chr_no)
            else:
                chr_no = Chromosome(int(str_chr_no[0]))

        else:
            s = f.readline()
            continue

        start_coord = int(start_coord)
        end_coord = int(end_coord)

        # Maybe one of the names is already in the dictionary. Then, link all of the names to the official title
        # for that name. Otherwise, let "gene name" be official.
        if gene_name in nameToOfficial:
            official_name = nameToOfficial[gene_name]
        else:
            official_name = gene_name

        for name in names:
            if name in nameToOfficial:
                official_name = nameToOfficial[name]
                if official_name not in officialToCoord:
                    continue
                prev_chr_no, prev_start_coord, prev_end_coord = officialToCoord[official_name]
                if prev_chr_no != chr_no:
                    nameToOfficial[name] = name
                    official_name = name
                continue

            nameToOfficial[name] = official_name

        # Now get the start and end coord of this row:
        prev_chr_no, prev_start_coord, prev_end_coord = officialToCoord.get(official_name,
                                                                            (-1, [], []))
        if int(prev_chr_no) != -1 and prev_chr_no != chr_no:
            s = f.readline()
            continue
        new_chr_no = chr_no
        new_start_coord = prev_start_coord + [start_coord]
        new_end_coord = prev_end_coord + [end_coord]

        officialToCoord[official_name] = (new_chr_no, new_start_coord, new_end_coord)

        s = f.readline()
    print("Done with going through the file...")
    for official_name, coords in officialToCoord.items():
        # print(official_name)
        # print(coords)
        chr_no, start_coord_array, end_coord_array = coords
        filtered_start_coord_array = remove_outliers(start_coord_array)
        filtered_end_coord_array = remove_outliers(end_coord_array)
        start_coord = min(filtered_start_coord_array)
        end_coord = max(filtered_end_coord_array)
        officialToCoord[official_name] = (chr_no, start_coord, end_coord)
        if len(filtered_start_coord_array) != len(start_coord_array):
            print("filtered!")
            print(filtered_start_coord_array)
            print(start_coord_array)
        if len(filtered_end_coord_array) != len(end_coord_array):
            print("filtered!")
            print(filtered_end_coord_array)
            print(end_coord_array)

    return nameToOfficial, officialToCoord


def remove_outliers(array):
    """

    :param array: 

    """
    return array
    # m = statistics.median(array)
    # d = [x - m for x in array]
    # m_dev = statistics.median(d)
    # s = [x / m_dev if m_dev else 0. for x in d]
    # return [x for i, x in enumerate(array) if s[i] < cutoff]


# TODO: Consider the potential role of picklify in this file
def gene_to_coord(gene):
    """

    :param gene: 

    """

    b = True
    if os.path.isfile(path_to_pickle):
        b = False
        with open(path_to_pickle, 'rb') as handle:
            try:
                nameToOfficial, officialToCoord = pickle.load(handle)
            except ModuleNotFoundError:
                b = True

    if b:
        nameToOfficial, officialToCoord = file_to_dicts(path_to_bio_mart)
        with open(path_to_pickle, 'wb') as handle:
            pickle.dump([nameToOfficial, officialToCoord], handle, protocol=pickle.HIGHEST_PROTOCOL)

    gene = gene.upper()
    rna_info = {}
    if gene in nameToOfficial and nameToOfficial[gene] in officialToCoord:
        rna_info['success'] = True
        rna_info['chr_n'], rna_info['start_coord'], rna_info['end_coord'] = (
            officialToCoord[nameToOfficial[gene]]
        )
        rna_info['official_name'] = nameToOfficial[gene]
    else:
        rna_info['success'] = False
    return rna_info


if __name__ == "__main__":
    pass
    # import itertools
    #
    # with open(path_to_pickle, 'rb') as handle:
    #     nameToOfficial, officialToCoord = pickle.load(handle)
    #
    #     length_wanted = 1000
    #     for g in map(''.join, itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890', repeat=5)):
    #         #print(g)
    #         if g in nameToOfficial and nameToOfficial[g] in officialToCoord:
    #             c,s,e = list(officialToCoord[nameToOfficial[g]])
    #             if 4000 > e-s > 1000:
    #                 print(e-s)
    #                 print(g)
