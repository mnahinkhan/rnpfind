"""
Defines the ucsc_visualize function.
This function takes care of generating binding sites, converting them to
appropriate formats (BED files, then bigBED files, wig, bigwig, etc.),
generating density plots, setting them up as a track hub structure, and finally
uploading them to a server online, so that the binding sites can be visualized.
"""

from .populate_rbp_binding_sites_script import populate_binding_sites
from .populate_trackhub import (
    convert_bed_to_bb,
    convert_wig_to_bw,
    density_plot,
    populate_local_track_hub,
    upload_online,
)
from .config import GENOME_VERSION

UCSC_START_STEP = 5

def ucsc_visualize(big_storage, rna_info, out=None, total_steps=7):
    """
    This function takes care of generating binding sites, converting them to
    appropriate formats (BED files, then bigBED files, wig, bigwig, etc.),
    generating density plots, setting them up as a track hub structure, and
    finally uploading them to a server online, so that the binding sites can be
    visualized.

    The function returns a URL where the binding sites may be viewed.

    :param big_storage: the Storage instance containing the binding sites
    :param rna_info: a dictionary specifying relevant RNA information, such as
                     the name of the RNA molecule, its location in the genome,
                     etc.
    :param out: if specified, allows for outputs to be redirected.
                All status updates are directed to the "out" function, if it is
                specified.
                (Default value = None)
    :param total_steps: don't use, will be deprecated soon.
                        (Default value = 7)

    """

    if not out:
        out = print

    # Vintage code allowing a specific RBP to be colored, will be removed soon.
    # out[0] = "For this analysis method,"
    #       + " may pick one RBP of interest to compare cooperation and
    # competition against, if you like.")
    # out[0] = "Please type the name of one RBP you are interested in (if not,
    # just press enter): ")
    # main_rbp: str = input()
    # main_rbp: str = ""


    data_load_sources = big_storage.keys()

    out(f"{UCSC_START_STEP}/{total_steps}. Populating the server with bed files...")
    overarching_path = (
        populate_binding_sites(
            big_storage, rna_info, data_load_sources
        )
    )

    out(
        f"{UCSC_START_STEP + 1}/{total_steps}."
        " Converting all the bed files to bb files..."
    )
    convert_bed_to_bb(overarching_path, data_load_sources)


    out(
        f"{UCSC_START_STEP + 2}/{total_steps}."
        " Generating density plot for RBP binding..."
    )

    rbp_no_dict = (
        density_plot(
            big_storage, rna_info, data_load_sources, overarching_path
        )
    )

    out(
        f"{UCSC_START_STEP + 3}/{total_steps}."
        " Converting all the wig files to bw files..."
    )

    convert_wig_to_bw(overarching_path, data_load_sources)

    out(
        f"{UCSC_START_STEP + 4}/{total_steps}."
        " Uploading the bigBed and bigWig files to a local track hub server..."
    )

    local_dir = "./website/output-data/ucsc/ucsc-genome-track-fake/"

    date_time_folder_name = overarching_path.split("/")[-2]
    local_stage = local_dir + date_time_folder_name + "/"
    rbp_peaks = {
        k: max(big_storage[k].sum_over_all().return_depth())
        for k in big_storage
    }

    hub_name = (
        populate_local_track_hub(
            overarching_path, rna_info, local_stage, rbp_no_dict, rbp_peaks
        )
    )

    out(
        f"{UCSC_START_STEP + 5}/{total_steps}."
        " Copying the local track files to a global server..."
    )

    online_dir = "s3://rnpfind-data/ucsc-trackhub/"
    upload_online(local_dir, online_dir)

    out("The link to your hub has been generated:")

    hub_url = (
        "https://rnpfind-data.s3-us-west-1.amazonaws.com/ucsc-trackhub/"
        + date_time_folder_name + "/" + hub_name.replace(" ", "+") + ".hub.txt"
    )

    ucsc_url = (
        "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" + GENOME_VERSION
        + "&hubUrl=" + hub_url + "&position=chr" + str(rna_info['chr_n']) + ":"
        + str(rna_info['start_coord']) + "-" + str(rna_info['end_coord'])
    )

    out(ucsc_url)
    return ucsc_url
