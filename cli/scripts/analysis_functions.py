"""
This file is responsible for managing all the analysis functions that RNPFind
supports. If you want to develop and incorporate your own analysis function into
RNPFind, this is the file to edit!

The way any analysis function works in RNPFind is as follows. RNPFind will call
the analysis function like this:

analysis_method_function_name(big_storage, rna_info)

where big_storage is a dictionary that maps a data_source string to a Storage
variable (which in turn maps rbp names to a list of binding sites; for more
details see load_data.py), and rna_info is a dictionary that stores the RNA name
and its chromosomal location on the hg38 chromosome. If you preferred to have
the RNA sequence, get_human_seq() from pwm_scan should help you (see
attract_data_load.py for an example of that).

Upon getting the above input, an analysis function is free to do with the
information as it pleases. The information it outputs to the user could be in
any form, such as text on the CLI, csv files saved on the local directory, heat
maps produced on a figure generated, or a URL to web-generated content. Once
done with the work, the function should return so that RNPFind can continue with
its process.

If an analysis function is made to the above specification, the below lines
should be edited to incorporate the function on RNPFind!
# """

# First, import the function from wherever it is defined:
from .over_all_corr_analysis import overall_correlation_analysis
from .ucsc_visualize import ucsc_visualize

# Give your method a short name:
analysis_methods_supported_short = ["csv", "bed"]

# Give your method a long name:
analysis_methods_supported_long = [
    "Binding correlation csv",
    "Bed files with binding data",
]

# Map your short name to the variable imported above that corresponds to your
# function!
analysis_method_functions = {
    "csv": overall_correlation_analysis,
    "bed": ucsc_visualize,
}
