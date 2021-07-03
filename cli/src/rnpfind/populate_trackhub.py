"""
This module is dedicated to functions that help generate UCSC Genome Broswer
trackhubs and populate them with binding sites and their annotations, displayed
as simple intervals of sites or as density (dense barchart/histogram) plots.

"""

import glob
import os
from pathlib import Path
from sys import platform

import trackhub

from .binding_analysis_binding_sites import OVERLAP_CONFLICT
from .config import (
    AUTOSQL_PATH,
    GENOME_VERSION,
    UCSC_TRACK_VISIBILITY,
    UCSCTOOL_PATH,
)
from .data_load_functions import (
    column_data,
    data_load_sources_supported_long,
    data_load_sources_supported_short,
)


def populate_local_track_hub(
    overarching_path, rna_info, local_stage, rbp_no_dict, rbp_peaks, rbp=""
):
    """
    Populates a local directory with files conforming to structure required by
    the UCSC specification for displaying trackhubs on their servers.

    It is assumed that all the .bb files have already been generated and are
    stored under a directory in the following structure:
        /<source 1>/<rbp1>.bb
        /<source 1>/<rbp2>.bb
        ...
        /<source 1>/<rbpn1>.bb
        /<source 2>/<rbp1>.bb
        ...
        ...
        /<source m>/<rbpnm>.bb


    This function simply copies the files into a trackhub structure and saves it
    in a specified directory.

    :param overarching_path: Directory under which the .bb files are stored in
        the structure specified above.
    :param rna_info: a dictionary used to represent RNA of interest. It should
        include the 'official_name' and 'chr_n' (chromosome number).
    :param local_stage: the directory to which the local trackhub structured
        files should be saved.
    :param rbp_no_dict: A dictionary containing data sources as keys and the
        number of RBPs discovered (to bind to RNA of interest) by each source as
        value.
    :param rbp_peaks: A dictionary containing data sources as keys and the
        highest number of RBPs discovered to bind to the RNA at one (nucleotide)
        point (just from that data source) as value.
    :param rbp: Name of one RBP of interest, among those binding to the RNA of
        interest (useful only for taking its perspective in competitive /
        cooperative relationship investigation.)

    """

    rna = rna_info["official_name"]
    rna_chr_no = rna_info["chr_n"]

    # # threshold_config_file = open(overarching_path + "threshold_config.txt")
    # with open(
    #     Path(overarching_path) / "threshold_config.txt"
    # ) as threshold_config_file:
    #     _str = threshold_config_file.read()

    assert OVERLAP_CONFLICT == "union"
    hub_name = f"rbps-on-{rna}-{OVERLAP_CONFLICT.lower()}"
    hub, _, _, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=f"RBPs on {rna}",
        long_label=(
            "RNA binding proteins on the transcript at {rna}. Overlapping"
            " sites were merged into one."
        ),
        genome=GENOME_VERSION,
        email="mnk1@alumni.cmu.edu",
    )

    for filename in glob.iglob(f"{overarching_path}/**/*.bb", recursive=True):
        category = Path(filename).parent.name
        name = Path(filename).name
        # _, _, _, _, _, _, category, name = filename.split("/")
        rbp = name.split("-")[0]
        rbp = rbp.replace(",", "_")
        rbp = rbp.replace("*", "_mut_")
        rbp = rbp.replace("(", "")
        rbp = rbp.replace(")", "")

        visibility = UCSC_TRACK_VISIBILITY
        # TODO: consider options for this for the user (add to Config at least)

        track = trackhub.Track(
            name=rbp + "-" + category + "-binding-sites",
            short_label=rbp + "-" + category,
            long_label=(
                "Binding sites of "
                + rbp
                + " derived from "
                + data_load_sources_supported_long[
                    data_load_sources_supported_short.index(category)
                ]
            ),
            source=filename,
            tracktype="bigBed 9 +",
            itemRgb="on",
            spectrum="on",
            visibility=visibility,
            chromosomes="chr" + str(rna_chr_no),
            # labelFields="", l
            # defaultLabelFields="", l
            # mouseOverField="" l
            # labelFields=",".join(
            #     [column_data[category]["names"][i]
            #     for i in column_data[category]["interest"]]
            # ),
            # defaultLabelFields=",".join(
            #     [column_data[category]["names"][i]
            #     for i in column_data[category]["default_label"]]
            # ),
            # mouseOverField=(
            #     column_data[category]["names"][column_data[category]
            #     ["default_mouse_over"]]
            # ),
            # maxItems=25
        )

        # TODO: Add options for mouse hover views of information

        trackdb.add_tracks(track)

    for filename in glob.iglob(f"{overarching_path}/**/*.bw", recursive=True):
        data_load_source = Path(filename).parent.name
        name = Path(filename).name
        # _, _, _, _, _, _, data_load_source, name = filename.split("/")
        rbp_no = rbp_no_dict[data_load_source]
        rbp_peak = max(rbp_peaks.values())
        rbp_peak = (rbp_peak // 10 + 1) * 10
        visibility = "full"

        # TODO: consider options for this for the user (add to Config at least)

        track = trackhub.Track(
            name=rna + "-" + data_load_source + "-density-plot",
            short_label="00 " + rna + "-density",
            long_label=(
                "Density plot of "
                + str(rbp_no)
                + " RBPs on "
                + rna
                + " using data from "
                + data_load_source.upper()
            ),
            source=filename,
            tracktype="bigWig",
            color="128,0,0",  # TODO: what color ought bigWig density plots be?
            visibility=visibility,
            chromosomes="chr" + str(rna_chr_no),
            viewLimits="0:" + str(rbp_peak),
            maxHeightPixels="128:50:8",
            autoScale="on",
        )

        trackdb.add_tracks(track)

    trackhub.upload.stage_hub(hub, staging=local_stage)

    return hub_name


def prepare_auto_sql(data_load_source):
    """
    Helper function for generating .as files for UCSC visualization. These files
    specify (for each data source) the column title for data that annotates the
    binding sites generated on UCSC genome browser.

    :param data_load_source: The data source for which the .as file should be
        created.

    """
    source_columns_of_interest = range(
        len(column_data[data_load_source]["names"])
    )

    no_of_extra_fields = len(source_columns_of_interest)
    name_of_file = (
        data_load_source
        + "".join([str(c) for c in source_columns_of_interest])
        + ".as"
    )

    file_path = f"{AUTOSQL_PATH}/{name_of_file}"
    template_file_path = f"{AUTOSQL_PATH}/general_template.as"
    try:
        open(file_path, "r").close()
    except FileNotFoundError:
        with open(file_path, "w") as handle:
            # TODO: make an auto generator of auto_sql template files here
            with open(template_file_path, "r") as template_handle:
                template_string = template_handle.read()
                template_string = template_string.replace(
                    "insert_source_name_here", data_load_source
                )

            handle.write(template_string)
            column_names = [
                column_data[data_load_source]["names"][i]
                for i in source_columns_of_interest
            ]
            descriptions = [
                column_data[data_load_source]["descriptions"][i]
                for i in source_columns_of_interest
            ]
            additional_str = ""
            for column, description in zip(column_names, descriptions):
                additional_str += (
                    "\t".join(
                        ["lstring", column + ";", '"' + description + '"']
                    )
                    + "\n"
                )

            additional_str += ")"
            handle.write(additional_str)

    return no_of_extra_fields, name_of_file


def convert_bed_to_bb(overarching_path, data_load_sources):
    """
    Given a directory containing .bed files, converts them all to .bb files and
    saves them. Uses the UCSC bedToBigBed tool as a subroutine.

    :param overarching_path: Directory containing the .bed files in a layout
        similar to as described in the docstring of populate_local_track_hub.

    :param data_load_sources: List of data sources for which the convertion is
        desired.

    """
    debug = False
    if GENOME_VERSION != "hg38":
        raise ValueError("Update this function for this genome version!")

    os_folder = ""
    if platform in ("linux", "linux2"):
        # linux
        os_folder = "linux"
    elif platform == "darwin":
        # OS X
        os_folder = "mac-os"
    assert os_folder

    # starting_working_directory = os.getcwd()
    for data_load_source in data_load_sources:
        no_of_extra_fields, as_file_name = prepare_auto_sql(data_load_source)
        # os.chdir(overarching_path + data_load_source + "/")

        os.system(
            f"for file in {Path(overarching_path) / data_load_source}/*.bed;"
            f" do echo $file; {UCSCTOOL_PATH}/{os_folder}/bedToBigBed"
            f" -as={AUTOSQL_PATH}/{as_file_name}"
            + " type=bed9+"
            + str(no_of_extra_fields)
            + ' "$file" '
            + f'{UCSCTOOL_PATH}/hg38.chrom.sizes "${{file%.*}}".bb; done'
            + (" >/dev/null 2>&1" if not debug else "")
        )

        # os.chdir(starting_working_directory)

    # end of function


def upload_online(local_dir, github_dir):
    """
    Given a local directory containing trackhub-structured files, uploads them
    to the internet so that UCSC Genome Browser can access it.

    This function uploads the local directory to an s3 bucket and assumes that
    s3 credentials have been configured.

    :param local_dir: Local directory with trackhub-structured files.
    :param github_dir: Directory on s3 to save the files on.

    """
    debug = False
    terminator = " >/dev/null 2>&1" if not debug else ""

    os.system("aws s3 sync " + local_dir + " " + github_dir + terminator)


def density_plot(big_storage, rna_info, data_load_sources, overarching_path):

    """
    Generates .wig files containing density information on RBPs binding to RNA
    of interest

    :param big_storage: a dictionary containing data sources as keys and
        Storage instances with the binding sites extracted from
        the corresponding data sources as values.
    :param rna_info: a dictionary specifying the RNA of interest, by specifying
        the name of the RNA and its genomic location on hg38. Keys needed are:
            'official_name' - name of RNA
            'chr_n' - chromosome in hg38 on which it lies
            'start_coord' - the coordinate on which it starts on the chromosome.
    :param data_load_sources: a list of sources for which tthe .wig files
        should be generated, such as 'rbpdb', 'postar', etc.
    :param overarching_path: the directory in which the .wig files should be
        saved.
    :returns: a dictionary containing the number of RBPs that were discovered
        by each data source on the RNA of interest.

    """

    rna = rna_info["official_name"]
    rna_chr_no = rna_info["chr_n"]
    rna_start_chr_coord = rna_info["start_coord"]
    rbp_no_dict = {}
    for data_load_source in data_load_sources:
        storage = big_storage[data_load_source]

        # TODO: check if displacement needs to be shifted by one for all data
        # sources or just RBPDB
        wig_string = storage.print_wig(
            chr_no=rna_chr_no,
            displacement=rna_start_chr_coord - 1,
            include_name=True,
            include_description=True,
            name=rna,
            description=(
                "Density plot of " + str(len(storage)) + " RBPs on " + rna
            ),
            include_header=True,
        )

        rbp_no_dict[data_load_source] = len(storage)

        folder_path = Path(overarching_path) / data_load_source
        filepath = (
            f"{rna}-{data_load_source}-{GENOME_VERSION}-density-plot.wig"
        )

        filepath = folder_path / filepath

        folder_path.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w") as density_plot_wig_file:
            density_plot_wig_file.write(wig_string)

    return rbp_no_dict


def convert_wig_to_bw(overarching_path, data_load_sources):
    """
    Given a directory containing .wig files, converts them all to .bw files and
    saves them. Uses the UCSC wigToBigWig tool as a subroutine.

    :param overarching_path: Directory containing the .wig files in a layout
        similar to as described in the docstring of populate_local_track_hub.

    :param data_load_sources: List of data sources for which the convertion is
        desired.

    """
    if GENOME_VERSION != "hg38":
        raise ValueError("Update this function for this genome version!")

    os_folder = ""
    if platform in ("linux", "linux2"):
        # linux
        os_folder = "linux"
    elif platform == "darwin":
        # OS X
        os_folder = "mac-os"
    assert os_folder

    # starting_working_directory = os.getcwd()
    for data_load_source in data_load_sources:
        # os.chdir(overarching_path + data_load_source + "/")

        os.system(
            f"for file in {overarching_path+data_load_source}/*.wig;"
            f" do echo $file; {UCSCTOOL_PATH}/{os_folder}/wigToBigWig"
            f' "$file" {UCSCTOOL_PATH}/hg38.chrom.sizes '
            '"$file.bw"; done >/dev/null 2>&1'
        )

        # os.chdir(starting_working_directory)
