from .config import data_load_sources_supported, data_load_sources_supported_short_form
from .analysis_functions import analysis_methods_supported_long, analysis_methods_supported_short
from .gene_coordinates import gene_to_coord


def get_rna_coord(rna_gene_name):
    """

    :param rna_gene_name: 

    """
    user_input_required = True
    deduction_attempt_made = False

    rna_info = gene_to_coord(rna_gene_name)

    is_success = rna_info['success']
    if is_success:
        RNA_chr_no = rna_info['chr_n']
        RNA_start_chr_coord = rna_info['start_coord']
        RNA_end_chr_coord = rna_info['end_coord']

        deduction_attempt_made = True
        print("We have automatically deduced that this gene lies on chromosome " + str(RNA_chr_no) + " from " +
              str(RNA_start_chr_coord) + " to " + str(RNA_end_chr_coord) + " (with length " +
              str(RNA_end_chr_coord - RNA_start_chr_coord) + " bases)")
        print("Are you okay with the above coordinates? [y/n]: ")
        s = input()
        while len(s) != 1 or s.lower() not in "yn":
            print("")
            print("Please type 'y' or 'n'")
            print("Are you okay with the above coordinates? [y/n]: ")
            s = input()
        print("")
        if s.lower() == "y":
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
                "Sorry, we are having trouble figuring out the location of this gene on the genome. Could you tell us?")
            print("")

        print("Chromosome number (1-23): > ")
        RNA_chr_no = input()
        print("")
        print("Thanks! What about the start coordinate of this gene on chromosome " + RNA_chr_no + "?:")
        print("")
        print("Start coordinate: > ")
        RNA_start_chr_coord = input()
        print("")
        print("Thanks! What about the end coordinate of this gene on chromosome " + RNA_chr_no + "?:")
        print("")
        print("End coordinate: > ")
        RNA_end_chr_coord = input()
    return RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord


def get_user_rna_preference() -> dict:
    """ """
    print("")
    print("")
    print("Welcome to RNPFind!")
    print("")
    print("")
    print("Which RNA transcript would you like to analyse today?")
    RNA = input()  # "Neat1"

    RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = get_rna_coord(RNA)
    rna_info = {}
    rna_info['official_name'] = RNA.upper()
    rna_info['chr_n'] = RNA_chr_no
    rna_info['start_coord'] = int(RNA_start_chr_coord)
    rna_info['end_coord'] = int(RNA_end_chr_coord)
    return rna_info


def get_user_data_source_preference():
    """ """
    print("")
    print("")
    print("Which sources of data would you like to collect RBP binding data from today?")
    print("")
    for i, source in enumerate(data_load_sources_supported):
        print("[" + str(i) + "]: " + source)
    print("")
    s = ""
    while not (s.isdigit() and all([int(c) < len(data_load_sources_supported_short_form) for c in s])):
        print("Please choose any combination from above as you like (e.g. 124)")
        print(">")
        s = input()
    print("Thank you")
    return [data_load_sources_supported_short_form[int(i)] for i in s]


def get_user_analysis_preference():
    """ """
    print("")
    print("")
    print("Which method of analysis on the data would you like to employ today?")
    print("")
    for i, source in enumerate(analysis_methods_supported_long):
        print("[" + str(i) + "]: " + source)
    print("")
    print("Please pick just one analysis method and write the number associated with it (e.g. 3)")
    print(">")
    s = input()
    while not (len(s) == 1 and s.isdigit() and int(s) < len(analysis_methods_supported_short)):
        print("Please pick just one analysis method and write the number associated with it (e.g. 3)")
        print(">")
        s = input()
    print("Thank you")
    return analysis_methods_supported_short[int(s)]
