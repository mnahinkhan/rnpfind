"""
This module contains functions for loading data from the postar database.

"""

import bisect
import os

from .config import ANNOTATION_COLUMN_DELIMITER, POSTAR_PATH

postar_all_column_names = [
    "chrom",
    "chromStart",
    "chromEnd",
    "postarID",
    "nil",
    "strand",
    "rbpName",
    "dataSource",
    "cellType",
    "expSource",
    "postarScore",
]

postar_all_column_descriptions = [
    "chromosome number",
    "start coordinate",
    "end coordinate",
    "POSTAR database ID",
    "not sure",
    "strand",
    "RBP Name",
    "Data Source",
    "Cell type",
    "experimental source",
    "score",
]

postar_columns_of_interest = [3, 7, 8, 9, 10]
postar_default_label_index = [8]
POSTAR_DEFAULT_MOUSE_OVER_INDEX = 9

postar_column_names = [
    postar_all_column_names[i] for i in postar_columns_of_interest
]

postar_column_descriptions = [
    postar_all_column_descriptions[i] for i in postar_columns_of_interest
]


class Query:
    """
    Represents a region of interest within a genome.
    Supports comparison with a line from a file containing details on binding
    sites, in a format much like POSTAR's public database files. Since the file
    is sorted, this allows quick binary search over the file using this class's
    comparison function as a comparator against the lines in the file (to find
    the region of interst, stored in an instance of this class at __init__).

    """

    def __init__(self, query):
        self.query = query

    def __lt__(self, line):
        rna_chr_no, rna_start_chr_coord, rna_end_chr_coord = self.query
        line_parts = line.split()
        # print(line_parts)

        return (line_parts[0], int(line_parts[1]), int(line_parts[2])) > (
            "chr" + str(rna_chr_no),
            rna_start_chr_coord,
            rna_end_chr_coord,
        )

    def __str__(self):
        return str(self.query)


class FileSearcher:
    """
    Class for representing a file as an indexable element. This allows binary
    search over the file using in-built libraries (bisect.bisect)
    """

    def __init__(self, file_pointer):
        self.file_pointer = file_pointer
        self.file_pointer.seek(0, os.SEEK_END)
        self.num_bytes = self.file_pointer.tell() - 1

    def __len__(self):
        return self.num_bytes

    def __getitem__(self, i):
        # TODO: Fix the three corner case bugs inherent in a binary search
        # algorithm like this
        # See: (http://pts.github.io/pts-line-bisect/line_bisect_evolution.html)
        self.file_pointer.seek(i)
        self.file_pointer.readline()

        return self.file_pointer.readline()


def binary_search_populate(file_path, rna_info, debug=False):
    """
    Searches a file containing sorted binding sites for a region of interest.
    Returns the subset of binding sites required as a generator / iterator
    object.

    :param file_path: a file path containing sorted binding sites on the genome.
    :param rna_info: a dictionary containing chromosome number, start, and end
        coordinate to slice out of the genome-wide binding sites file.
    :param debug: prints useful information if set to True, for debugging.
        (Default value = False)
    """
    # TODO: Fix a bug here that causes genes without any data to start
    # collecting the whole genome!!

    rna_chr_no = rna_info["chr_n"]
    rna_start_chr_coord = rna_info["start_coord"]
    rna_end_chr_coord = rna_info["end_coord"]
    query = Query((rna_chr_no, rna_start_chr_coord, rna_end_chr_coord))

    postar_data_file = open(file_path)

    search_file = FileSearcher(postar_data_file)
    to_seek = bisect.bisect(search_file, query)
    postar_data_file.seek(to_seek)
    postar_line_parts = postar_data_file.readline().split()
    # print(postar_line_parts)
    is_found = False
    seen = []
    not_found_counter = 0
    while postar_line_parts:
        if (
            postar_line_parts[0] == "chr" + str(rna_chr_no)
            and int(postar_line_parts[1]) > rna_start_chr_coord
            and int(postar_line_parts[2]) < rna_end_chr_coord
        ):
            is_found = True

            if debug:
                if (postar_line_parts[7]) not in seen:
                    print(";".join(postar_line_parts))
                    seen += [postar_line_parts[7]]

            rbp = postar_line_parts[6]
            start, end = postar_line_parts[1], postar_line_parts[2]
            start = int(start) - rna_start_chr_coord
            end = int(end) - rna_start_chr_coord

            # TODO: Consider reformatting the annotation for visual appeal
            annotation = ANNOTATION_COLUMN_DELIMITER.join(
                [postar_line_parts[i] for i in postar_columns_of_interest]
            )
            yield rbp, start, end, annotation

        elif is_found:
            break
        if not is_found:
            not_found_counter += 1
            if not_found_counter >= 4:
                break
        postar_line_parts = postar_data_file.readline().split()


def postar_data_load(rna_info):
    """
    Returns a generator containing binding sites on input RNA molecule, as
    found on the POSTAR database.
    :param rna_info: dictionary containing input RNA information, such as
        chromosome number, start coordinate, and end coordinate.

    """
    file_path = f"{POSTAR_PATH}/postar-human-RBP-binding-sites-sorted.txt"
    return binary_search_populate(file_path, rna_info)


if __name__ == "__main__":
    test_rna_info = ["MALAT1", 11, 65497688, 65506516]
    postar_data_load(test_rna_info)
