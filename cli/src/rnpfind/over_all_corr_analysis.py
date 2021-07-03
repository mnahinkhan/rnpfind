"""
This module is dedicated to the correlation matrix generation analysis function,
a data analysis function

"""

import sys
from pathlib import Path


def generate_csv(
    symmetric_corr_table, rna_info, data_load_sources, stringency, out_dir
):
    """
    Given a correlation table of RBPs binding to an RNA molecule of interest,
    generates and saves a CSV file containing the correlation values.
    :param symmetric_corr_table: A nested dictionary such that
        symmetric_corr_table[rbp1][rbp2] returns the binding correlation value
        between rbp1 and rbp2.
    :param rna_info: A dictionary containing 'official_name' as key, with a
        corresponding value that indicates the name of the RNA under
        investigation.
    :param data_load_sources: A list of data sources that were used, such as
        'postar', 'rbpdb', etc.
    :param stringency: numerical value representing the base stingency that was
        used in analyzing the correlation coefficient.
    :returns: a filepath to the saved CSV file.

    """
    rna = rna_info["official_name"]

    path_to_save = (
        Path(out_dir)
        / "csv"
        / f"{rna.lower()}-{'-'.join(data_load_sources)}-{stringency}.csv"
    )

    path_to_save.parent.mkdir(parents=True, exist_ok=True)
    path_to_save = str(path_to_save)

    with open(path_to_save, "w+") as csv_file:
        rbps = list(symmetric_corr_table.keys())
        csv_file.write("," + ",".join(rbps))
        csv_file.write("\n")
        for rbp1 in rbps:
            csv_file.write(rbp1 + ",")
            for rbp2 in rbps:
                csv_file.write(str(symmetric_corr_table[rbp1][rbp2]))
                csv_file.write(",")
            csv_file.write("\n")

    return path_to_save


def generate_heat_map(_, rna_info, __, ___):
    """
    Function meant to generate a heatmap using correlation values.
    Unimplemented.
    :param _:
    :param rna_info:
    :param __:
    :param ___:

    """
    del rna_info
    print("HEAT MAP (supposed to be) GENERATED")


def overall_correlation_analysis(big_storage, rna_info, configs=None):
    """
    A data analysis function (designed for the CLI variant of RNPFind) that
    allows one to calculate and generate CSV files and/or heatmaps for binding
    correlation between RBPs binding on an RNA molecule of interest.

    :param big_storage: A dictionary keyed-in by data load source name. Each
        such data load source key should correspond to a Storage instance
        containing binding sites extracted from said data source (on the RNA of
        interest).
    :param rna_info: Dictionary containing key "official_name", corresponding
        to the name of the RNA molecule of interest.
    :param config: A dictionary of additional paramters. The following keys are
        useful:
         - 'base_stringency': specifies number of bases to consider to be
           competing.
         - 'out_dir': specifies the directory to write csv files to.

    """

    threshold = configs["base_stringency"]

    symmetric_corr_tables = {}
    data_load_sources = big_storage.keys()
    for data_load_source in data_load_sources:
        # print("")
        # print("Working on ", data_load_source, "data...")

        storage = big_storage[data_load_source]

        _, symmetric_corr_table, _, _ = storage.self_analysis(
            bp_threshold=threshold, progress_feedback=False
        )
        symmetric_corr_tables[data_load_source] = symmetric_corr_table
        # print("Done!")

    # print("Creating CSV files...")
    for data_load_source in data_load_sources:
        symmetric_corr_table = symmetric_corr_tables[data_load_source]
        path_saved = generate_csv(
            symmetric_corr_table,
            rna_info,
            [data_load_source],
            threshold,
            configs["out_dir"],
        )
        print("A file was saved at", path_saved + ".", file=sys.stderr)
    # print("Done!")
