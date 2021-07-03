"""
Defines the ucsc_visualize function.
This function takes care of generating binding sites, converting them to
appropriate formats (BED files, then bigBED files, wig, bigwig, etc.),
generating density plots, setting them up as a track hub structure, and finally
uploading them to a server online, so that the binding sites can be visualized.
"""

import os
import shutil
import sys
from pathlib import Path

from .populate_rbp_binding_sites_script import populate_binding_sites
from .populate_trackhub import (
    convert_bed_to_bb,
    convert_wig_to_bw,
    density_plot,
    populate_local_track_hub,
)


def ucsc_visualize(big_storage, rna_info, configs=None):
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
    :param print: if specified, allows for printputs to be redirected.
                All status updates are directed to the "print" function, if it is
                specified.
                (Default value = None)
    :param total_steps: don't use, will be deprecated soon.
                        (Default value = 7)
    :param configs: gives additional configurations. Supported keys include:
        - out_dir: the directory to write files to

    """

    data_load_sources = big_storage.keys()

    overarching_path = str(Path(configs["out_dir"]) / "bed")
    populate_binding_sites(
        big_storage, rna_info, data_load_sources, overarching_path
    )

    if not configs["trackhub"] and not configs["trackhub-only"]:
        # Trackhub is not needed at all...
        return

    convert_bed_to_bb(overarching_path, data_load_sources)

    rbp_no_dict = density_plot(
        big_storage, rna_info, data_load_sources, overarching_path
    )

    convert_wig_to_bw(overarching_path, data_load_sources)

    print("Generating trackhub structure...", file=sys.stderr)

    local_dir = Path(overarching_path).parent / "trackhub"

    rbp_peaks = {
        k: max(big_storage[k].sum_over_all().return_depth())
        for k in big_storage
    }

    populate_local_track_hub(
        overarching_path, rna_info, local_dir, rbp_no_dict, rbp_peaks
    )

    # Resolving links
    shutil.copytree(local_dir, f"{local_dir}-copy")
    # os.system(f"cp -rL {local_dir} {local_dir}-copy")
    os.system(f"rm -rf {local_dir}")
    os.system(f"mv {local_dir}-copy {local_dir}")

    os.system("find . -name '*.bb' | xargs rm")

    print(f"Created directory {local_dir}...", file=sys.stderr)

    if configs["trackhub-only"]:
        print("Removing original bed files...", file=sys.stderr)
        os.system(f"rm -rf {overarching_path}")

    # print("The link to your hub has been generated:")

    # hub_url = (
    #     "https://rnpfind.com/static/ucsc-tracks/"
    #     + date_time_folder_name
    #     + "/"
    #     + hub_name.replace(" ", "+")
    #     + ".hub.txt"
    # )

    # ucsc_url = (
    #     "https://genome.ucsc.edu/cgi-bin/hgTracks?db="
    #     + GENOME_VERSION
    #     + "&hubUrl="
    #     + hub_url
    #     + "&position=chr"
    #     + str(rna_info["chr_n"])
    #     + ":"
    #     + str(rna_info["start_coord"])
    #     + "-"
    #     + str(rna_info["end_coord"])
    # )

    # print(ucsc_url)
    # return ucsc_url
