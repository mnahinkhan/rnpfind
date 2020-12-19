from .populate_rbp_binding_sites_script import populate_binding_sites
from .populate_trackhub import populate_local_track_hub, convert_bed_to_bb, upload_online, density_plot, \
    convert_wig_to_bw
from .config import genome_version


def ucsc_visualize(big_storage, rna_info, out=None, total_steps=7):
    """

    :param big_storage: 
    :param rna_info: 
    :param out:  (Default value = None)
    :param total_steps:  (Default value = 7)

    """
    if not out:
        out = lambda s:print(s)

    # out[0] = "For this analysis method,"
    #       + " may pick one RBP of interest to compare cooperation and competition against, if you like.")
    # out[0] = "Please type the name of one RBP you are interested in (if not, just press enter): ")
    # main_rbp: str = input()
    main_rbp: str = ""


    data_load_sources = big_storage.keys()
    # out[0] = "Thank you!"
    start_step  = 5
    out(f"{start_step}/{total_steps}. Populating the server with bed files...")
    overarching_path = populate_binding_sites(big_storage, rna_info, data_load_sources, main_rbp)
    # out("done!"
    # out(""

    out(f"{start_step + 1}/{total_steps}. Converting all the bed files to bb files...")
    convert_bed_to_bb(overarching_path, data_load_sources)
    # out("Done!"
    # out(""

    out(f"{start_step + 2}/{total_steps}. Generating density plot for RBP binding...")
    rbp_no_dict = density_plot(big_storage, rna_info, data_load_sources, overarching_path, return_rbp_no=True)
    # out("done!"
    # out(""

    out(f"{start_step + 3}/{total_steps}. Converting all the wig files to bw files...")
    convert_wig_to_bw(overarching_path, data_load_sources)
    # out("Done!"
    # out(""

    out(f"{start_step + 4}/{total_steps}. Uploading the bigBed and bigWig files to a local track hub server...")
    local_dir = "./website/output-data/ucsc/ucsc-genome-track-fake/"

    date_time_folder_name = overarching_path.split("/")[-2]
    local_stage = local_dir + date_time_folder_name + "/"
    rbp_peaks = {k: max(big_storage[k].sum_over_all().return_depth()) for k in big_storage}
    hub_name = populate_local_track_hub(overarching_path, main_rbp, rna_info, local_stage, rbp_no_dict, rbp_peaks)
    # out("done!"

    # out(""
    out(f"{start_step + 5}/{total_steps}. Copying the local track files to a global server...")
    online_dir = "s3://rnpfind-data/ucsc-trackhub/"
    upload_online(local_dir, online_dir)
    # out("done!"
    # out(""
    out("The link to your hub has been generated:")

    hub_url = "https://rnpfind-data.s3-us-west-1.amazonaws.com/ucsc-trackhub/" + \
              date_time_folder_name + "/" + hub_name.replace(" ", "+") + ".hub.txt"

    ucsc_url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" + genome_version + "&hubUrl=" + \
               hub_url + "&position=chr" + str(rna_info['chr_n']) + ":" + str(rna_info['start_coord']) + "-" + str(
                rna_info['end_coord'])

    out(ucsc_url)
    return ucsc_url
