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


import argparse
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

# Responsible for managing the analysis functions that manipulate the RNA-RBP
# interaction data to get a useful output:
from .analysis_functions import (
    analysis_method_functions,
    analysis_methods_supported_short,
)
from .config import (
    DEFAULT_BASE_STRINGENCY,
    RO_DATA_PATH,
    RO_DATA_TAR_NAME,
    RO_DATA_URL,
)
from .data_load_functions import data_load_sources_supported_short

# Responsible for managing the loading of RNA-RBP interaction data:
from .load_data import load_data

# Functions that help with interacting with the user to get their preference:
from .user_input import get_user_rna_preference


def rm_folder_contents(folder):
    """
    Remove all the contents in a directory.
    From: https://stackoverflow.com/a/185941/8551394

    """
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except OSError:
            print("Failed to delete %s.")


import urllib.request

from tqdm import tqdm


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(
        unit="B", unit_scale=True, miniters=1, desc=url.split("/")[-1]
    ) as t:
        urllib.request.urlretrieve(
            url, filename=output_path, reporthook=t.update_to
        )


def download_ro_data():
    """
    Deletes the ro-data directory and downloads the contents from
    rnpfind.com

    """
    # Delete folder
    try:
        shutil.rmtree(RO_DATA_PATH)
    except FileNotFoundError:
        pass

    # Create folder
    Path(RO_DATA_PATH).mkdir(parents=True, exist_ok=True)

    # Download tar file
    print("Downloading tar file...", file=sys.stderr)
    download_url(RO_DATA_URL, RO_DATA_TAR_NAME)

    # Extract tar file
    print("Extracting tar file (this will take some time)...", file=sys.stderr)
    shutil.unpack_archive(
        RO_DATA_TAR_NAME, extract_dir=Path(RO_DATA_PATH).parent
    )

    # Delete the tar file
    Path(RO_DATA_TAR_NAME).unlink()

    # Display confirmation of completion
    print("Done!", file=sys.stderr)


def rnpfind(
    transcript,
    sources=None,
    methods=None,
    base_stringency=None,
    out_dir=None,
    is_trackhub=False,
    is_trackhub_only=False,
):
    """
    Collect binding data of RBPs on RNA.

      :param transcript: gene name or genomic location to specify transcript
      :param sources: list of data sources to limit to for data collection
      :param methods: list of output formats to restrict to
      :param base_stringency: config option for csv output method
      :param out_dir: directory to write output files in
      :param is_trackhub: whether to generate trakchub structure
      :param is_trackhub_only: wheter to delete BED files in the end
    """

    # First, check if readonly data directory exists
    if not Path(RO_DATA_PATH).is_dir():
        # We assume that if the dir exists the data is fine; otherwise
        # the data needs to be downloaded

        # In case of corrupt data one would have to call download_ro_data()
        # manually (Or if one were to wish for just the data without analysis)
        print("Downloading data necessary for rnpfind...", file=sys.stderr)
        download_ro_data()

    # Start by getting the transcript of interest to analyze
    rna_info: dict = get_user_rna_preference(transcript)

    # what data sources does the user want to collect data from today?
    # (e.g. attract, postar, etc.)
    data_load_sources = (
        sources if sources else data_load_sources_supported_short
    )
    # print(
    #     f"Collecting data from: {', '.join(data_load_sources)}",
    #     file=sys.stderr,
    # )
    # load RNA-RBP interaction data using the selected data sources on the RNA
    # molecule of interest big_storage stores data on binding sites of RBPs on
    # the RNA molecule from each data source. For more details on how
    # big_storage is structured, consult load_data.py!
    big_storage = load_data(data_load_sources, rna_info)

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
        f"Collected data for {no_rbps} RBPs with {no_sites} binding sites",
        file=sys.stderr,
    )

    # We now proceed to perform any number of analysis methods that the user
    # may wish to apply to the data obtained
    analysis_methods = methods if methods else analysis_methods_supported_short
    # print(
    #     "Producing the following output formats:"
    #     f" {', '.join(analysis_methods)}",
    #     file=sys.stderr,
    # )

    if out_dir:
        # The user specified an out directory
        # Create it if it does not exist
        # TODO: consider error handling
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    else:
        # The user did not specify an out directory
        default_path = Path.cwd() / rna_info["official_name"].lower()

        while default_path.is_dir():
            # The default path already exists, append the folder name with
            # current date and time
            time_list = datetime.now().timetuple()
            time_list = [str(x) for x in time_list]
            time_date = "-".join(time_list[0:6])  # year to seconds
            default_path = default_path.parent / (
                f"{default_path.name}-{time_date}"
            )

        default_path.mkdir(parents=True)
        out_dir = str(default_path)

    rm_folder_contents(out_dir)

    for analysis_method in analysis_methods:
        print(f"Generating {analysis_method}", file=sys.stderr)
        analysis_method_function = analysis_method_functions[analysis_method]

        configs = {"out_dir": out_dir}
        if analysis_method == "csv":
            configs["base_stringency"] = (
                base_stringency if base_stringency else DEFAULT_BASE_STRINGENCY
            )
        if analysis_method == "bed":
            configs["trackhub"] = is_trackhub
            configs["trackhub-only"] = is_trackhub_only

        analysis_method_function(big_storage, rna_info, configs=configs)

    print("Done!", file=sys.stderr)


def main():
    """
    main function responsible for parsing commandline args

    """
    out_formats = analysis_method_functions.keys()

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
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--trackhub",
        action="store_true",
        help="If specified, generates a trackhub structure using bed files"
        " (only used for 'bed' output format)",
        default=False,
    )
    group.add_argument(
        "--trackhub-only",
        action="store_true",
        help="If specified, generates a trackhub structure using bed files and"
        " deletes all original bed files (only used for 'bed' output format)",
        default=False,
    )

    args = parser.parse_args()
    rnpfind(
        args.transcript,
        args.sources,
        args.out_format,
        args.base_stringency,
        args.out_dir,
        args.trackhub,
        args.trackhub_only,
    )


if __name__ == "__main__":
    main()
