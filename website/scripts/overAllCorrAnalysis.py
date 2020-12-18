# File dedicated to the correlation matrix generation analysis function


def generate_csv(symmetric_corr_table, rna_info, data_load_sources, stringency):
    RNA = rna_info['official_name']

    path_to_save = "./website/data/csv_files/" + RNA + "_" + "_".join(data_load_sources) + "_" + str(stringency) + ".csv"
    with open(path_to_save, 'w+') as csv_file:
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
    del rna_info
    print("HEAT MAP GENERATED")
    return


# As of now, takes two rna arguments only.
def overall_correlation_analysis(big_storage, rna_info):
    print("For this analysis method, we need the nucleotide base distance stringency threshold number.")
    print('''This number indicates the number of bases away that two binding sites can be before they are ''' +
          '''considered "too far"''')
    print("")
    threshold = ""
    while not threshold.isdigit():
        print("Could you please type in your preferred base distance stringency threshold number? (e.g. 30)")
        threshold = input(">")

    threshold = int(threshold)
    print("Thanks, this could take some time...")
    symmetric_corr_tables = {}
    data_load_sources = big_storage.keys()
    for data_load_source in data_load_sources:
        print("")
        print("Working on ", data_load_source, "data...")

        storage = big_storage[data_load_source]

        _, symmetric_corr_table, _, _ = storage.self_analysis(bp_threshold=threshold)
        symmetric_corr_tables[data_load_source] = symmetric_corr_table
        print("Done!")

    output_type = ""
    while output_type != "0" and output_type != "1" and output_type != "01" and output_type != "10":
        print("We can either generate a csv file or show you a heat map, what would you prefer?")
        print("[0]: Heat map")
        print("[1]: CSV file")
        print('''Write any combination [e.g. "01"]''')
        output_type = input(">")

    if "1" in output_type:
        print("Creating CSV files...")
        for data_load_source in data_load_sources:
            symmetric_corr_table = symmetric_corr_tables[data_load_source]
            path_saved = generate_csv(symmetric_corr_table, rna_info, [data_load_source], threshold)
            print("The file was saved at", path_saved + ".")
        print("Done!")

    if "0" in output_type:
        print("Generating heat map...")
        for data_load_source in data_load_sources:
            symmetric_corr_table = symmetric_corr_tables[data_load_source]
            generate_heat_map(symmetric_corr_table, rna_info, data_load_sources, threshold)
        print("Done!")
