#!/usr/bin/env python3

"""   iY Lab
  Project: Developing a tool for exploring  RNA-Protein interactions
  Name 	: Muhammad Nahin Khan
  AndrewID  : mnk1
  File Created: 01/26/2020
  Script written by M. Nahin Khan
  mnk1@andrew.cmu.edu

RNPFind is a program that allows the user to input one RNA template molecule
and to analyse the RBPs that bind on to it.

Introduction:
RNPFind takes as input the name of one RNA-encoding gene, or the coordinates of
a region of interest on the hg38 human genome. In the future we plan on
supporting an arbitrary RNA sequence as input to the program.

The program then takes as input a preferred set of data sources from which it
should extract the RNA-RBP binding information from (e.g. experimental-based
ones like POSTAR or computational-prediction-based databases like RBPDB).
Currently, RNPFind supports automated data extraction from RBPDB, ATTRACT, and
POSTAR. Furthermore, custom-data (for example data obtained from one's own
experiments) are also supported for inputting to the program. If you wish to add
support for any other data source / algorithm of your liking, I recommend
looking at "data_load_functions.py", which outlines ways for you to contribute
your own source and integrating it into RNPFind.

After loading RBP-Protein binding data, RNPFind then does internal processing to
efficiently store the binding sites of RBPs on the template RNA. For details on
how this might work, take a look at "bind_analysis.py" and
"binding_analysis_binding_sites.py".

Finally, after internal processing, RNPFind has multiple analysis features for
getting useful information out of the RBP-binding data. Some of these analysis
methods may require additional input from the user, such as certain computation
parameters. Currently, we support two main analysis features: visualization on
UCSC genome browser (within which we have integrated viewing density plots) and
correlation matrix generation. For details on supporting other analysi
functions that you may be interested in contributing, check out
"analysis_functions.py"!"""


# Responsible for managing the analysis functions that manipulate the RNA-RBP
# interaction data to get a useful output:
from scripts.analysis_functions import analysis_method_functions

# Responsible for managing the loading of RNA-RBP interaction data:
from scripts.load_data import load_data

# Functions that help with interacting with the user to get their preference:
from scripts.user_input import (
    get_user_analysis_preference,
    get_user_data_source_preference,
    get_user_rna_preference,
)


def analysis_script(transcript):
    """
    analysis_script: runs command line version of RNPFind
    """
    # Start by getting the transcript of interest to analyze
    rna_info: dict = get_user_rna_preference(transcript)

    # what data sources does the user want to collect data from today?
    # (e.g. attract, postar, etc.)
    data_load_sources = get_user_data_source_preference()
    print(data_load_sources)
    # load RNA-RBP interaction data using the selected data sources on the RNA
    # molecule of interest big_storage stores data on binding sites of RBPs on
    # the RNA molecule from each data source. For more details on how
    # big_storage is structured, consult load_data.py!
    print("Collecting data now...")
    big_storage = load_data(data_load_sources, rna_info)
    print("complete!")

    # BIOGRID is a database that stores information on protein-protein
    # interaction evidence in the literature from experiment. In future versions
    # of RNPFind, BIOGRID data should be helpful in a variety of inquiries when
    # investigating relationships between RBPs that bind on RNA molecules.

    # Todo: Consider collecting and using BIOGRID data in a meaningful way

    # Give a quick summary on the data that has been loaded, in terms of total
    # number of RBPs and binding sites collected for the RNA molecule of
    # interest
    no_rbps = 0
    no_sites = 0
    for no_rbp, no_site in [big_storage[k].summary() for k in big_storage]:
        no_rbps += no_rbp
        no_sites += no_site

    print(
        "We have populated "
        + str(no_rbps)
        + " different RBPs with "
        + str(no_sites)
        + " different binding sites on the "
        + rna_info["official_name"]
        + " sequence across the "
        + str(rna_info["end_coord"] - rna_info["start_coord"])
        + " bases specified!"
    )

    # We now proceed to perform any number of analysis methods that the user
    # may wish to apply to the data obtained
    analysis_method = get_user_analysis_preference()
    print(analysis_method)
    analysis_method_function = analysis_method_functions[analysis_method]
    analysis_method_function(big_storage, rna_info)


if __name__ == "__main__":
    import argparse

    from scripts.data_load_functions import data_load_sources_supported_short

    out_formats = analysis_method_functions.keys()
    DEFAULT_BASE_STRINGENCY = 30  # check if available elsewhere

    parser = argparse.ArgumentParser(
        add_help=False,
        description="Get binding sites of RBPs on a given transcript",
    )
    parser.add_argument(
        "transcript",
        help="Specify with the name of a gene (e.g. 'Malat1')"
        " or as hg38 chromosome coordinates given as"
        " <chr_no>:<start_coord>-<end_coord> (e.g. 5:4000-14000)"
        ". Note that chromosome number is X, Y, M(T), or a number between"
        " 1 and 22",
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "-s",
        "--sources",
        choices=data_load_sources_supported_short,
        nargs="+",
        help="Pick source(s) for RBP binding data. Pick any non-empty subset"
        f" of values from {{{', '.join(data_load_sources_supported_short)}}}."
        " If unspecified, all sources are selected.",
        metavar=("<source 1>", "<source 2>"),
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        metavar="<dir>",
        help="Directory to write output files in. If unspecified, a folder is"
        " created and written to (if possible)",
    )
    parser.add_argument(
        "-f",
        "--out-format",
        choices=out_formats,
        nargs="+",
        metavar=("<format 1>", "<format 2>"),
        help="Choose format(s) to store binding site data in. Pick any"
        f" non-empty subset from {{{', '.join(out_formats)}}}."
        " If unspecified, all formats are created in the output directory.",
    )
    parser.add_argument(
        "-b",
        "--base-stringency",
        type=int,
        help="The number of bases beween two RBP-binding-sites before they"
        " are considered to be competing (used only for csv output format)."
        f" The default value is {DEFAULT_BASE_STRINGENCY}.",
        metavar="<N>",
        default=DEFAULT_BASE_STRINGENCY,
    )
    args = parser.parse_args()

    analysis_script(args.transcript)
