"""
This file is responsible for managing all the data loading functions that
RNPFind supports. If you want to develop and incorporate your own data loading
function into RNPFind, this is the file to edit!

The way any data loading function works in RNPFind is as follows. RNPFind will
call the data loading function like this: data_load_function_name(rna_info)

where rna_info is a list that stores the RNA name and its chromosomal location
on the hg38 chromosome. If you preferred to have the RNA sequence,
get_human_seq() from pwm_scan should help you (see attract_data_load.py for an
example of that).

Upon getting the above input, a data loading function must return an iterable
object that contains elements of the following form that each represent a
binding site of an rbp on the RNA molecule: (rbp, start, end, annotation). For
example, the function might return a list that looks like
[
    (rbp_1, start_1, end_1, annotation_1),
    (rbp_2, start_2, end_2, annotation_2), ...
    (rbp_n, star_n, end_n, annotation_n)
].

See attract_data_load.py for an example where "yield" was used to make an
iterator instead that serves the same purpose (more out of a laziness of
changing code than any demonstrable efficiency).

Above, for each binding site, "rbp" is the name of an RBP molecule (a string).

The "start" refers to the starting coordinate of a binding site relative to the
start of the RNA (where "start" is always the thing on the left side, doesn't
matter if the gene is on anti-sense strand) on the hg38 chromosome, and "end"
refers to the nucleotide number where the site ends (again relative to where the
RNA starts on the chromosome).

"Annotation" is any annotation you might want to add to the binding site. I have
struggled with standardizing this, switching back and forth between strings with
delimiters or just a list of annotations. For now, I have stuck with the
following standardization: the data loading function should give a list of
"features" that it wants to associate with the binding site. For example, a
binding site may be associated with the features "database id", "binding score",
and "experiment_type", so the function would return something like
(rbp_1, start_1, end_1, ["nm00042", 0.9, "clip seq"]). Internally, RNPFind does
stuff like converting it to a string separated by delimiters, but that's just a
relic from the past and my lack of firm will in changing the code.

Above, annotation to the binding sites are given just as data, but what about
the "column headings", so to speak? That's exactly a part of what this file
does!

If you make a data loading function to the above specification, the below lines
should be edited to incorporate the function on RNPFind!
"""

from .attract_data_load import (
    ATTRACT_DEFAULT_MOUSE_OVER_INDEX,
    attract_column_descriptions,
    attract_column_names,
    attract_data_load,
    attract_default_label_index,
)

# These imports are more obvious once you see below
from .colors import blue, dark_green, orange
from .custom_data_load import custom_data_load
from .postar_data_load import (
    POSTAR_DEFAULT_MOUSE_OVER_INDEX,
    postar_column_descriptions,
    postar_column_names,
    postar_data_load,
    postar_default_label_index,
)
from .rbpdb_data_load import (
    RBPDB_DEFAULT_MOUSE_OVER_INDEX,
    rbpdb_column_descriptions,
    rbpdb_column_names,
    rbpdb_data_load,
    rbpdb_default_label_index,
)

# First, add a short name for your data source method to this list
data_load_sources_supported_short = ["rbpdb", "attract", "postar"]

# Now give it a long name
data_load_sources_supported_long = [
    "RBPDB (computational)",
    "ATTRACT (computational)",
    "POSTAR (experimental)",
]

# Link your short name to the actual data loading function that complies with
# the above specifications here.

# This is where importing is useful!
data_load_sources_functions = {
    "rbpdb": rbpdb_data_load,
    "attract": attract_data_load,
    "postar": postar_data_load,
    "custom": custom_data_load,
}

# Finally, for the sake of UCSC visualization, give a little bit info about
# your data source function's annotations. Map the short name of your data load
# function to another dictionary, which has the following mappings:
#   "names": mapped to a list of column headings (strings) that the annotations
#       have
#   "default_label": mapped to the index of the column heading that would be
#       preferred to be displayed next to each binding site on the UCSC view.
#       This feature is NOT implemented yet.
#   "default_mouse_over": mapped to the index of the column heading that would
#       be preferred to be displayed when a mouse is hovered on a binding site
#       on the UCSC view. This feature is NOT implemented yet.
#   "descriptions: mapped to a list of descriptions for each of the column
#       headings. This could be the same as names, if the names are
#       self-explanatory.

column_data = {
    "postar": {
        "names": postar_column_names,
        "default_label": postar_default_label_index,
        "default_mouse_over": POSTAR_DEFAULT_MOUSE_OVER_INDEX,
        "descriptions": postar_column_descriptions,
    },
    "attract": {
        "names": attract_column_names,
        "default_label": attract_default_label_index,
        "default_mouse_over": ATTRACT_DEFAULT_MOUSE_OVER_INDEX,
        "descriptions": attract_column_descriptions,
    },
    "rbpdb": {
        "names": rbpdb_column_names,
        "default_label": rbpdb_default_label_index,
        "default_mouse_over": RBPDB_DEFAULT_MOUSE_OVER_INDEX,
        "descriptions": rbpdb_column_descriptions,
    },
}

# When the data source is viewed on UCSC, what color should the binding sites
# be? Give a tuple with three numbers, indicating the rgb values of the color
# (each ranging from 0 to 255, as is common).

data_load_source_colors = {
    "postar": orange,
    "attract": blue,
    "rbpdb": dark_green,
}

# To consider: right now, I have mainly filtered through the "columns of
# interest" for each of the data loading functions WITHIN the data loading
# functions, just to be selective about which columns are interesting to view
# on the UCSC visualization feature. This seems like bad practice; maybe in the
# future someone will develop more analysis functions that are interested in
# some of the other columns of the data; it might be wise to consider loading
# all columns of data at the data loading level and letting the analysis
# functions individually deal with which of the columns they are interested in.
