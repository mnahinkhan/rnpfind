"""
This file contains functions for loading data from the rbpdb database.

"""

import math

from .config import (
    ANNOTATION_COLUMN_DELIMITER,
    RBPDB_MOTIF_N_REPEAT_REQ,
    RBPDB_MOTIF_PWM_LETTER_STRENGTH,
    RBPDB_PATH,
)
from .picklify import picklify
from .pwm_scan import (
    get_human_seq,
    motif_to_pwm,
    pwm_degree_of_freedom,
    pwm_scan,
    pwm_scan_naive_brute_force,
    str_to_pwm,
)

rbpdb_all_column_names = [
    "protein_id",
    "annotation_id",
    "creation_date",
    "update_date",
    "gene_name",
    "gene_description",
    "species",
    "taxID",
    "domains",
    "aliases",
    "flag_protein",
    "flag_notes",
    "some_other_id",
    "experimental_id",
    "PUBMED_ID",
    "exp_type",
    "notes",
    "seq_motif",
    "selex_file",
    "aligned_selex_file",
    "aligned_motif_file",
    "PWM_file",
    "PFM_file",
    "logo_file",
    "secondary_structure",
    "in_vivo_notes",
    "in_vivo_file",
    "flag_experimental",
]

rbpdb_all_column_descriptions = rbpdb_all_column_names
rbpdb_columns_of_interest = [
    0,
    1,
    3,
    4,
    5,
    8,
    9,
    13,
    14,
    15,
    16,
    17,
    22,
    24,
    25,
]
rbpdb_default_label_index = [8]

RBPDB_DEFAULT_MOUSE_OVER_INDEX = 9

rbpdb_column_names = [
    rbpdb_all_column_names[i] for i in rbpdb_columns_of_interest
]
rbpdb_column_descriptions = [
    rbpdb_all_column_descriptions[i] for i in rbpdb_columns_of_interest
]


def generate_rbpdb_experiment_to_columns():
    """
    Converts RBPDB file format to a useful one.
    I don't remember what columns means, but this is a one-time preprocessing
    step needed to extract data out of the RBPDB data.

    """
    rbpdb_experiment_file_path = (
        f"{RBPDB_PATH}/RBPDB_v1.3.1_experiments_human_2012-11-21.tdt"
    )
    experiment_id_to_columns_dict = {}
    with open(rbpdb_experiment_file_path) as handle:
        line = handle.readline()
        while line:
            columns = line.split("\t")
            # Here we expect the columns to be:
            # experimental_id, PUBMED_ID, exp_type, notes, seq_motif,
            # selex_file, aligned_selex_file,
            # aligned_motif_file, PWM_file, PFM_file, logo_file,
            # secondary_structure, in_vivo_notes, in_vivo_file, flag
            if columns[14] == "1":
                # The flag means this data is unreliable, according to the RBPDB
                # Readme files
                line = handle.readline()
                continue
            experimental_id = columns[0]
            experiment_id_to_columns_dict[experimental_id] = columns
            line = handle.readline()
    return experiment_id_to_columns_dict


def generate_rbpdb_protein_to_experiment_id():
    """
    Preprocesses RBPDB files to generate a dictionary mapping protein IDs to
    the experiment ID that investigates their binding data.
    Run once in initialization of RNPFind
    """
    rbpdb_protein_experiment_file_path = (
        f"{RBPDB_PATH}/RBPDB_v1.3.1_protExp_human_2012-11-21.tdt"
    )
    protein_id_to_experimental_ids_dict = {}
    with open(rbpdb_protein_experiment_file_path) as handle:
        line = handle.readline()
        while line:
            columns = line.split("\t")
            # Here we expect the columns to be:
            # protein_id, experiment_id, homolog, unique_id
            protein_id = columns[0]
            experimental_id = columns[1]
            protein_id_to_experimental_ids_dict[
                protein_id
            ] = protein_id_to_experimental_ids_dict.get(protein_id, []) + [
                experimental_id
            ]
            line = handle.readline()
    return protein_id_to_experimental_ids_dict


def generate_rbpdb_experimental_to_pwm(letter_strength, n_repeat_req):
    """
    Processes RBPDB files to generate and return a dictionary mapping RBPDB
    experiment IDs to PWM data structures representing the binding motif of the
    RBP under investigation in the experiment of interest.
    :param letter_strength: the number of times more likely a base is considered
        to be when represented by the motif string at a particular position.
    :param n_repeat_req: specifies the assumption to use for the value of n when
        RBPDB files specify a motif such as "(ATG)n", which to us humans means
        a number of repeats of "ATG". This parameter either makes n_repeat_req
        the minimum number of motif repeats for the binding motif, or makes it
        exactly n_repeat_req. I'm not too sure.

    """
    rbpdb_experiment_file_path = (
        f"{RBPDB_PATH}/RBPDB_v1.3.1_experiments_human_2012-11-21.tdt"
    )
    rbpdb_pfm_file_directory = f"{RBPDB_PATH}/rbpdb-human-pfm-matrices/"
    experimental_to_pwm_dict = {}
    with open(rbpdb_experiment_file_path) as handle:
        line = handle.readline()
        while line:
            columns = line.split("\t")
            # Here we expect the columns to be:
            # experimental_id, PUBMED_ID, exp_type, notes, seq_motif,
            # selex_file, aligned_selex_file,
            # aligned_motif_file, PWM_file, PFM_file, logo_file,
            # secondary_structure, in_vivo_notes, in_vivo_file, flag
            if columns[14] == "1":
                # The flag means this data is unreliable, according to the RBPDB
                # Readme files
                line = handle.readline()
                continue

            experimental_id = columns[0]

            assert len(experimental_id) > 0
            pfm_file = columns[9]
            seq_motifs = columns[4]
            pwms = []
            if pfm_file != "\\N":
                pfm_file_path = rbpdb_pfm_file_directory + pfm_file
                with open(pfm_file_path) as pfm_file_handle:
                    raw_pwm_str = pfm_file_handle.read()
                pwm = str_to_pwm(raw_pwm_str, is_transpose=True)
                pwms += [pwm]
            elif seq_motifs not in ("\\N", ""):
                # This experiment still generated some useful data
                seq_motifs = seq_motifs.split(";")
                i = 0
                while i != len(seq_motifs):
                    seq_motif = seq_motifs[i]
                    while ")(" in seq_motif:
                        repeat_end = seq_motif.find(")(")
                        assert seq_motif[repeat_end] == ")"
                        repeat_start = repeat_end
                        while seq_motif[repeat_start] != "(":
                            repeat_start -= 1

                        number_start = repeat_end + 2
                        assert (
                            seq_motif[number_start].isdigit()
                            or seq_motif[number_start] == "n"
                        )
                        number_end = number_start
                        while seq_motif[number_end] != ")":
                            number_end += 1

                        # deal with cases where the number of repeats is
                        # "15-30". Take minimum to be conservative.
                        # Note that most cases would be a single number like
                        # "15".
                        num_of_repeats = min(
                            [
                                int(s)
                                if s != "n"
                                else math.ceil(
                                    n_repeat_req
                                    / (repeat_end - repeat_start - 1)
                                )
                                for s in seq_motif[
                                    number_start:number_end
                                ].split("-")
                            ]
                        )

                        seq_motif = seq_motif.replace(
                            seq_motif[repeat_start : number_end + 1],
                            seq_motif[repeat_start + 1 : repeat_end]
                            * num_of_repeats,
                        )

                    maketrans = str.maketrans
                    all_letters = "wruysn"
                    upper_map = maketrans(all_letters, all_letters.upper())
                    seq_motif = seq_motif.translate(upper_map)
                    if "/" in seq_motif:
                        bracket_start = bracket_end = middle = seq_motif.find(
                            "/"
                        )
                        while seq_motif[bracket_start] != "(":
                            bracket_start -= 1
                        while seq_motif[bracket_end] != ")":
                            bracket_end += 1
                        seq_motif_1 = seq_motif.replace(
                            seq_motif[bracket_start : bracket_end + 1],
                            seq_motif[bracket_start + 1 : middle],
                        )
                        seq_motif_2 = seq_motif.replace(
                            seq_motif[bracket_start : bracket_end + 1],
                            seq_motif[middle + 1 : bracket_end],
                        )
                        seq_motifs += [seq_motif_1, seq_motif_2]
                    else:
                        pwm = motif_to_pwm(
                            seq_motif, letter_strength=letter_strength
                        )
                        pwms += [pwm]
                    i += 1

            # Now we have the raw text, we convert it to pwm and add to
            # dictionary
            experimental_to_pwm_dict[experimental_id] = pwms
            line = handle.readline()

    return experimental_to_pwm_dict


def rbpdb_data_load(rna_info):
    """
    Returns a Generator(/Iterator?) that represent binding sites loaded from the
    RBPDB database on an RNA molecule of interest.
    :param rna_info: a dictionaru containing the chromosome number, start, and
        end coordinate of the RNA molecule of interest (in hg38).
    :param out: if specified, redirects progress status to the specified
        function instead of stdout. (Default value = None)

    """
    rbpdb_protein_file_path = (
        f"{RBPDB_PATH}/RBPDB_v1.3.1_proteins_human_2012-11-21.tdt"
    )
    letter_strength = RBPDB_MOTIF_PWM_LETTER_STRENGTH
    n_repeat_req = RBPDB_MOTIF_N_REPEAT_REQ
    rna_seq = get_human_seq(rna_info)

    experiment_id_to_pwm_dict = picklify(
        generate_rbpdb_experimental_to_pwm, letter_strength, n_repeat_req
    )
    protein_id_to_experimental_ids_dict = picklify(
        generate_rbpdb_protein_to_experiment_id
    )
    experiment_id_to_columns_dict = picklify(
        generate_rbpdb_experiment_to_columns
    )
    with open(rbpdb_protein_file_path) as handle:
        _ = handle.readline().strip().split("\t")
        # columns here is expected to have the following information in the
        # following order:
        # protein_id, annotation_id, creation_date, update_date, gene_name,
        # gene_description, species, taxID, domains, aliases, flag, flag_notes,
        # some_other_id
        protein_columns = handle.readline().replace("\n", "").split("\t")
        while protein_columns != [""]:
            assert len(protein_columns) == 13
            # We only care about human RBPs for now.
            if protein_columns[10] == "0":
                protein_columns = (
                    handle.readline().replace("\n", "").split("\t")
                )
                continue
            rbp = protein_columns[4]
            protein_id = protein_columns[0]

            if protein_id not in protein_id_to_experimental_ids_dict:
                # No experiments associated. So no data to be had
                protein_columns = (
                    handle.readline().replace("\n", "").split("\t")
                )
                continue

            for experiment_id in protein_id_to_experimental_ids_dict[
                protein_id
            ]:
                assert (
                    experiment_id in experiment_id_to_pwm_dict
                    or experiment_id == "410"
                )
                if experiment_id == "410":
                    continue
                pwms = experiment_id_to_pwm_dict[experiment_id]
                for pwm in pwms:
                    assert len(pwm["A"]) > 0
                    experimental_columns = experiment_id_to_columns_dict[
                        experiment_id
                    ]
                    assert len(experimental_columns) == 15
                    total_columns = protein_columns + experimental_columns
                    annotation = ANNOTATION_COLUMN_DELIMITER.join(
                        [total_columns[i] for i in rbpdb_columns_of_interest]
                    )

                    if pwm_degree_of_freedom(pwm) >= 2048:
                        # experimentally shown that by this point naive brute
                        # force is faster. Bound could be
                        # reduced.
                        sites = pwm_scan_naive_brute_force(rna_seq, pwm)
                    else:
                        sites = pwm_scan(rna_seq, pwm)

                    if not sites:
                        continue

                    for start, end in sites:
                        yield rbp, start, end, annotation

            protein_columns = handle.readline().replace("\n", "").split("\t")
