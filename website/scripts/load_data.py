# This file is responsible for loading RBP-RNA interaction data from the needed data sources

# Explanation of data structure used in this program:

# Top level: big_storage (dict)
# Branches (keys): data_load_sources:
# 'postar', 'attract', 'rbpdb', etc.

# Next level:
# big_storage[source]
# Value: Storage variable (RBP -> RNA intervals mapping)
# Branches (keys): RBPs
# 'HNRNPD', 'HNRNPC', etc.

# Next level:
# big_storage[source][RBP]
# Value: a BindingSites variable,
# a bunch of RBP binding sites (intervals)

# E.g. To get the binding sites of HNRNPC on
# the template RNA as determined experimentally by postar:

# big_storage['postar']['HNRNPC']

from .data_load_functions import data_load_sources_functions
from .merge_annotation_funcs import generate_merge_func
from .bind_analysis import BindingSites, Storage
from .config import experimental_binding_site_acceptable_coverage_ratio, annotation_row_delimiter, \
    annotation_column_delimiter
# from synonym_dict_build import deal_with_dictionary_building


def load_data(data_load_sources, rna_info: dict, out=None, total_steps=7):
    # This function creates a big_storage variable that maps data sources to storage variables that store binding data
    # retrieved from the data source
    # TODO: investigate if the synonym_func is still relevant.
    # synonym_func = deal_with_dictionary_building()
    if not out:
        out = lambda s:print(*s)
    
    big_storage = {}
    for i, data_load_source in enumerate(data_load_sources):
        out(f"{i+1}/{total_steps}. Loading binding sites from {data_load_source}")
        # TODO: is the merge func still relevant?
        merge_func = generate_merge_func(data_load_source)
        storageSpace = Storage(annotation_merge_func=merge_func)
        big_storage[data_load_source] = storageSpace
        collected_data = data_load_sources_functions[data_load_source](rna_info, out=out)

        for rbp, start, end, annotation in collected_data:
            if rbp not in storageSpace:
                storageSpace[rbp] = BindingSites(overlap_mode=True)
            storageSpace[rbp].add((start, end, annotation))
        # Now we merge all the binding sites that overlap.
        # TODO: fix the implementation of overlap_collapse so annotations are not lost

        # Get max coverage
        max_coverage = max([binding_site.base_cover() for rbp, binding_site in storageSpace.items()])
        # Filter the allowed amount
        allowed_coverage = experimental_binding_site_acceptable_coverage_ratio * max_coverage

        for binding_site in storageSpace.values():
            binding_site.overlap_collapse("baseCoverNumber", allowed_coverage, in_place=True,
                                          annotation_merger=lambda t: annotation_row_delimiter.join(t))

    return big_storage


def annotation_to_columns(annotation):
    rows = annotation.split(annotation_row_delimiter)
    array = [tuple(r.split(annotation_column_delimiter)) for r in rows]
    array = list(set(array))
    len_row = len(array[0])
    no_of_rows = len(array)
    array_of_strings = ['______'.join([array[i][j] for i in range(no_of_rows)]) for j in range(len_row)]

    return array_of_strings


keys = ["rbpdb", "attract", "rbpmap", "postar", "custom"]
data_source_annotation_to_columns = {k: annotation_to_columns for k in keys}

if __name__ == '__main__':
    pass
