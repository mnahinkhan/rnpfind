from datetime import datetime
import os

from .colors import red, green
from .config import genome_version, dedicated_analysis
from .load_data import data_source_annotation_to_columns
from .data_load_functions import data_load_source_colors


def get_overarching_path(rna):
    if not dedicated_analysis:
        year, month, day, hour, minute, sec, x, y, z = datetime.now().timetuple()
        year, month, day, hour, minute, sec = [str(x) for x in [year, month, day, hour, minute, sec]]
        time_date = "_".join([year, month, day, hour, minute, sec])
    else:
        time_date = rna
    return "./website/output-data/ucsc/rbp_binding_sites_bed_files/" + time_date + "/"


def populate_binding_sites(big_storage, rna_info, data_load_sources, main_rbp, out=None):
    if not out:
        out = lambda s: print(s)
    RNA = rna_info['official_name']
    RNA_chr_no = rna_info['chr_n']
    RNA_start_chr_coord = rna_info['start_coord']

    # TODO: check if this -1 patch is necessary for all data sources or just RBPDB
    displacement = RNA_start_chr_coord - 1

    overarching_path = get_overarching_path(RNA)

    for data_load_source in data_load_sources:
        out(f"populating {data_load_source} binding sites!")

        storage = big_storage[data_load_source]

        folder_path = overarching_path + data_load_source + "/"

        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            os.chmod(folder_path, 0o777)
            out(f"Directory {folder_path} created!")
        else:
            out(f"Directory {folder_path} already exists...")

        competitive_threshold_bp = 15
        cooperative_threshold_bp = 56

        f = open(overarching_path + "threshold_config.txt", "w")
        f.write("competitive threshold used: " + str(competitive_threshold_bp) + "\n")
        f.write("cooperative threshold used: " + str(cooperative_threshold_bp) + "\n")
        f.close()

        default_color = data_load_source_colors[data_load_source]
        comp_color = red
        coop_color = green

        def coloring_func(binding_site):
            competitive = main_rbp in storage.binds_near(binding_site, bp_threshold=competitive_threshold_bp)
            cooperative = main_rbp in storage.binds_near(binding_site, bp_threshold=cooperative_threshold_bp)
            return comp_color if competitive else coop_color if cooperative else default_color

        for rbp in storage:
            total_sites = storage[[rbp]].printBED(chrN=RNA_chr_no, displacement=displacement, endInclusion=True,
                                                  addAnnotation=True, includeColor=True, includeHeader=False,
                                                  conditionalColor_func=coloring_func, is_additional_columns=True,
                                                  annotation_to_additional_columns=data_source_annotation_to_columns[
                                                      data_load_source])

            filepath = rbp + "_" + data_load_source + "_" + genome_version + "_sites.bed"

            filepath = folder_path + filepath

            try:
                f = open(filepath, "w")
            except FileNotFoundError:
                os.makedirs(folder_path)
                f = open(filepath, "w")

            f.write(total_sites)
            f.close()
    return overarching_path
