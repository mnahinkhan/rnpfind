"""
Defines functions that are involved in getting some data from the user.
Used for the command line interface variant of RNPFind.
"""
from .config import data_load_sources_supported, data_load_sources_supported_short_form
from .analysis_functions import (
    analysis_methods_supported_long,
    analysis_methods_supported_short,
)
from .gene_coordinates import gene_to_coord


def get_rna_coord(rna_gene_name):
    """
    Given a gene name, informs the user of the predicted RNA coordinates.
    Allows the user to potentially revise the coordinates, if required.

    :param rna_gene_name: input RNA gene name (string)

    """
    user_input_required = True
    deduction_attempt_made = False

    rna_info = gene_to_coord(rna_gene_name)

    is_success = rna_info["success"]
    if is_success:
        rna_chr_no = rna_info["chr_n"]
        rna_start_chr_coord = rna_info["start_coord"]
        rna_end_chr_coord = rna_info["end_coord"]

        deduction_attempt_made = True
        print(
            "We have automatically deduced that this gene lies on chromosome "
            + str(rna_chr_no)
            + " from "
            + str(rna_start_chr_coord)
            + " to "
            + str(rna_end_chr_coord)
            + " (with length "
            + str(rna_end_chr_coord - rna_start_chr_coord)
            + " bases)"
        )
        print("Are you okay with the above coordinates? [y/n]: ")
        user_input = input()
        while len(user_input) != 1 or user_input.lower() not in "yn":
            print("")
            print("Please type 'y' or 'n'")
            print("Are you okay with the above coordinates? [y/n]: ")
            user_input = input()
        print("")
        if user_input.lower() == "y":
            print("Thank you!")
            user_input_required = False
        else:
            print("Alright, please give us the coordinates: ")
            print("")
            user_input_required = True

    if user_input_required:
        if not deduction_attempt_made:
            print("")
            print(
                "Sorry, we are having trouble figuring out the location of this"
                " gene on the genome. Could you tell us?"
            )
            print("")

        print("Chromosome number (1-23): > ")
        rna_chr_no = input()
        print("")
        print(
            "Thanks! What about the start coordinate of this gene on chromosome"
            + " "
            + rna_chr_no
            + "?:"
        )
        print("")
        print("Start coordinate: > ")
        rna_start_chr_coord = input()
        print("")
        print(
            "Thanks! What about the end coordinate of this gene on chromosome"
            + " "
            + rna_chr_no
            + "?:"
        )
        print("")
        print("End coordinate: > ")
        rna_end_chr_coord = input()
    return rna_chr_no, rna_start_chr_coord, rna_end_chr_coord


def get_user_rna_preference() -> dict:
    """
    Asks the user for the RNA molecule they wish to analyze (the name of the
    gene is expected). Converts the gene name into coordinates.

    :returns: a dictionary containing the official name of the RNA, along with
              its chrosome number, start coordinate, and end coordinate on the
              human genome.
    """
    print("")
    print("")
    print("Welcome to RNPFind!")
    print("")
    print("")
    print("Which RNA transcript would you like to analyse today?")
    rna = input()  # "Neat1"

    rna_chr_no, rna_start_chr_coord, rna_end_chr_coord = get_rna_coord(rna)
    rna_info = {}
    rna_info["official_name"] = rna.upper()
    rna_info["chr_n"] = rna_chr_no
    rna_info["start_coord"] = int(rna_start_chr_coord)
    rna_info["end_coord"] = int(rna_end_chr_coord)
    return rna_info


def get_user_data_source_preference():
    """
    Asks the user for their preferred method of data source retrieval.
    The user is asked to input numbers to represent their combination of
    requested data sources
    """
    print("")
    print("")
    print(
        "Which sources of data would you like to collect RBP binding data"
        " from today?"
    )
    print("")
    for i, source in enumerate(data_load_sources_supported):
        print("[" + str(i) + "]: " + source)
    print("")

    input_digits = ""
    while not (
        input_digits.isdigit()
        and all(
            [int(c) < len(data_load_sources_supported_short_form) for c in input_digits]
        )
    ):
        print("Please choose any combination from above as you like (e.g. 124)")
        print(">")
        input_digits = input()

    print("Thank you")
    return [data_load_sources_supported_short_form[int(i)] for i in input_digits]


def get_user_analysis_preference():
    """
    Asks the user for their preferred method of data analysis on retrieved
    binding sites. The user is asked to input numbers to represent their
    combination of requested data sources.
    """
    print("")
    print("")
    print("Which method of analysis on the data would you like to employ today?")
    print("")
    for i, source in enumerate(analysis_methods_supported_long):
        print("[" + str(i) + "]: " + source)
    print("")
    print(
        "Please pick just one analysis method and write the number associated"
        " with it (e.g. 3)"
    )
    print(">")
    input_digits = input()
    while not (
        len(input_digits) == 1
        and input_digits.isdigit()
        and int(input_digits) < len(analysis_methods_supported_short)
    ):
        print(
            "Please pick just one analysis method and write the number"
            " associated with it (e.g. 3)"
        )
        print(">")
        input_digits = input()

    print("Thank you")
    return analysis_methods_supported_short[int(input_digits)]
