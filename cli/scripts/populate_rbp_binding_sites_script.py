"""
Defines functions useful for creating and saving RBP binding sites on a given
RNA transcript, as BED files.

This module is useful for ucsc_visualize data analysis method.
"""
import os
import sys
from pathlib import Path

from .colors import green, red
from .config import GENOME_VERSION
from .data_load_functions import data_load_source_colors
from .load_data import data_source_annotation_to_columns


def populate_binding_sites(
    big_storage, rna_info, data_load_sources, overarching_path, main_rbp="",
):
    """
    Saves binding sites as BED files
    :param big_storage: a Storage instance containing binding sites data for
        RBPs.
    :param rna_info: a dictionary containing the official_name, chr_n, and
        start_coord for the RNA molecule of interest (on the genome)
    :param data_load_sources: An iterable of data load sources from which the
        binding sites are populated in the big_storage
    :param main_rbp: one RBP of interest, to allow for deducing competitive
        and cooperative relationships against
    """

    rna_chr_no = rna_info["chr_n"]
    rna_start_chr_coord = rna_info["start_coord"]

    # TODO: check if this -1 patch is necessary for all data sources or just
    # RBPDB
    displacement = rna_start_chr_coord - 1

    comp_color = red
    coop_color = green
    competitive_threshold_bp = 15
    cooperative_threshold_bp = 56

    for data_load_source in data_load_sources:
        storage = big_storage[data_load_source]

        folder_path = str(Path(overarching_path) / data_load_source) + "/"

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            os.chmod(folder_path, 0o777)
            print(f"Directory {folder_path} created...", file=sys.stderr)
        else:
            print(
                f"Directory {folder_path} already exists...", file=sys.stderr
            )

        # with open(
        #     Path(overarching_path) / "threshold_config.txt", "w"
        # ) as threshold_config_file:
        #     threshold_config_file.write(
        #         "competitive threshold used: "
        #         + str(competitive_threshold_bp)
        #         + "\n"
        #     )
        #     threshold_config_file.write(
        #         "cooperative threshold used: "
        #         + str(cooperative_threshold_bp)
        #         + "\n"
        #     )

        default_color = data_load_source_colors[data_load_source]

        def coloring_func(
            binding_site, storage=storage, default_color=default_color
        ):
            """
            Defines a coloring function for the print_bed function of the
            Storage instance. This one in particular colors based on
            competitive and cooperative RBPs in comparison with main_rbp. Sites
            that are very close to main_rbp binding sites are colored with
            comp_color (red), and those that are far but not too far are
            colored with coop_color (green). Those that are too far are neutral
            / independent / default.

            :param binding_site: binding site under consideration

            """
            competitive = main_rbp in storage.binds_near(
                binding_site, bp_threshold=competitive_threshold_bp
            )
            cooperative = main_rbp in storage.binds_near(
                binding_site, bp_threshold=cooperative_threshold_bp
            )
            return (
                comp_color
                if competitive
                else coop_color
                if cooperative
                else default_color
            )

        for rbp in storage:
            total_sites = storage[[rbp]].print_bed(
                chr_n=rna_chr_no,
                displacement=displacement,
                end_inclusion=True,
                add_annotation=True,
                include_color=True,
                include_header=False,
                conditional_color_func=coloring_func,
                is_additional_columns=True,
                annotation_to_additional_columns=data_source_annotation_to_columns[
                    data_load_source
                ],
            )

            # filepath = (
            #     rbp
            #     + "_"
            #     + data_load_source
            #     + "_"
            #     + GENOME_VERSION
            #     + "_sites.bed"
            # )
            filepath = f"{rbp}-{data_load_source}-{GENOME_VERSION}-sites.bed"

            filepath = folder_path + filepath

            Path(folder_path).mkdir(parents=True, exist_ok=True)
            with open(filepath, "w") as bed_file:
                bed_file.write(total_sites)

    return overarching_path
