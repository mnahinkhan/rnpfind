"""
This module defines functions relevant to loading RBP binding sites from the
ATTRACT database.

"""

from .config import ANNOTATION_COLUMN_DELIMITER, ATTRACT_PATH
from .picklify import picklify
from .pwm_scan import get_human_seq, pwm_scan, str_to_pwm

attract_all_column_names = [
    "Gene_name",
    "Gene_id",
    "Mutated",
    "Organism",
    "Motif",
    "Len",
    "Experiment_description",
    "Database",
    "Pubmed",
    "Experiment",
    "Family",
    "Matrix_id",
    "attractScore",
]
attract_all_column_descriptions = attract_all_column_names
attract_columns_of_interest = [2, 4, 6, 7, 8, 9, 10, 11, 12]
attract_default_label_index = [8]
ATTRACT_DEFAULT_MOUSE_OVER_INDEX = 9

attract_column_names = [
    attract_all_column_names[i] for i in attract_columns_of_interest
]
attract_column_descriptions = [
    attract_all_column_descriptions[i] for i in attract_columns_of_interest
]


def generate_matrix_to_pwm_dict():
    """
    Preprocesses ATTRACT database files. In particular, this function generates
    and returns a dictionary that maps position weight matrix IDs to the PWM
    data structures used in this program (RNPFind).
    """
    attract_pwm_file_path = f"{ATTRACT_PATH}/attract-pwm.txt"
    matrix_to_pwm_dict = {}
    with open(attract_pwm_file_path) as handle:
        line = handle.readline()
        while line:
            matrix_id = line.split()[0][1:]
            raw_pwm_str = ""
            line = handle.readline()
            while line and line[0] != ">":
                raw_pwm_str += line
                line = handle.readline()
            # Now we have the raw text, we convert it to pwm and add to
            # dictionary
            matrix_to_pwm_dict[matrix_id] = str_to_pwm(raw_pwm_str)
    return matrix_to_pwm_dict


def attract_data_load(rna_info):
    """
    Loads RBP binding sites from the ATTRACT database for a given RNA molecule,
    and returns the sites as a Generator / Iterator object.
    :param rna_info: a dictionary containing the location of the RNA on the
        hg38 genome by specifying its chromosome location, start coordinates,
        and end coordinates.

    """
    attract_protein_file_path = f"{ATTRACT_PATH}/ATtRACT_db.txt"
    rna_seq = get_human_seq(rna_info)
    matrix_to_pwm_dict = picklify(generate_matrix_to_pwm_dict)
    with open(attract_protein_file_path) as handle:
        columns = handle.readline().strip().split("\t")
        assert columns == [
            "Gene_name",
            "Gene_id",
            "Mutated",
            "Organism",
            "Motif",
            "Len",
            "Experiment_description",
            "Database",
            "Pubmed",
            "Experiment_description",
            "Family",
            "Matrix_id",
            "Score",
        ]
        protein_columns = handle.readline().replace("\n", "").split("\t")
        while protein_columns != [""]:
            # Warning: Score ends with \n here, maybe remove using strip or
            # indexing. For now, we don't care about score as it seems to be
            # about literature

            # We only care about human RBPs for now.
            if protein_columns[3] != "Homo_sapiens":
                protein_columns = (
                    handle.readline().replace("\n", "").split("\t")
                )
                continue

            annotation = ANNOTATION_COLUMN_DELIMITER.join(
                [protein_columns[i] for i in attract_columns_of_interest]
            )

            rbp = protein_columns[0]

            matrix_id = protein_columns[11]

            pwm = matrix_to_pwm_dict[matrix_id]
            sites = pwm_scan(rna_seq, pwm)
            if not sites:
                protein_columns = (
                    handle.readline().replace("\n", "").split("\t")
                )
                continue

            for start, end in sites:
                yield rbp, start, end, annotation

            protein_columns = handle.readline().replace("\n", "").split("\t")
