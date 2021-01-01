"""
This file is responsible for loading RBP-RNA interaction data from the needed
data sources.

Explanation of data structure used in RNPFind  (and thus generated by the
functions of this module):

Top level: big_storage (dict)
Branches (keys): data_load_sources:
'postar', 'attract', 'rbpdb', etc.

Next level:
big_storage[source]
Value: Storage variable (RBP -> RNA intervals mapping)
Branches (keys): RBPs
'HNRNPD', 'HNRNPC', etc.

Next level:
big_storage[source][RBP]
Value: a BindingSites variable,
a bunch of RBP binding sites (intervals)

E.g. To get the binding sites of HNRNPC on
the template RNA as determined experimentally by postar:

big_storage['postar']['HNRNPC']

"""

from .data_load_functions import data_load_sources_functions
from .merge_annotation_funcs import generate_merge_func
from .bind_analysis import (
    BindingSites,
    Storage
)
from .config import (
    ANNOTATION_COLUMN_DELIMITER,
    ANNOTATION_ROW_DELIMITER,
    EXPERIMENTAL_BINDING_SITE_ACCEPTABLE_COVERAGE_RATIO,
)

def load_data(data_load_sources, rna_info: dict, out=None, total_steps=7):
    """
    Goes over a list of data sources of interest for a particular RNA and
    populates a Storage instance (for each of the data sources) with binding
    sites obtained from the data source. Returns a big_storage variable that
    maps data sources to storage variables that store binding data retrieved
    from the data source.

    :param data_load_sources: a list containing data sources, such as 'rbpdb',
        'postar', etc.
    :param rna_info: dict: a dictionary consisting of information about the RNA
        of interest, such as its name and genomic location.
    :param out: redirects progress status output if specified.
        (Default value = None)
    :param total_steps: will be deprecated soon, don't use. (Default value = 7)
    :returns: a dictionary mapping data load source to a Storage instance
        containing binding sites obtained from that data source.

    """

    if not out:
        out = print

    big_storage = {}
    for i, data_load_source in enumerate(data_load_sources):
        out(
            f"{i+1}/{total_steps}."
            f" Loading binding sites from {data_load_source}"
        )

        # TODO: is the merge func still relevant?
        merge_func = generate_merge_func(data_load_source)

        storage_space = Storage(annotation_merge_func=merge_func)
        big_storage[data_load_source] = storage_space
        collected_data = (
            data_load_sources_functions[data_load_source](rna_info, out=out)
        )

        for rbp, start, end, annotation in collected_data:
            if rbp not in storage_space:
                storage_space[rbp] = BindingSites(overlap_mode=True)
            storage_space[rbp].add((start, end, annotation))


        # Now we merge all the binding sites that overlap.

        # TODO: fix the implementation of overlap_collapse so annotations are
        # not lost

        # Get max coverage
        max_coverage = max(
            [
                binding_site.base_cover()
                for rbp, binding_site in storage_space.items()
            ]
        )

        # Filter the allowed amount
        allowed_coverage = (
            EXPERIMENTAL_BINDING_SITE_ACCEPTABLE_COVERAGE_RATIO * max_coverage
        )

        for binding_site in storage_space.values():
            binding_site.overlap_collapse(
                "baseCoverNumber", allowed_coverage, in_place=True,
                annotation_merger=ANNOTATION_ROW_DELIMITER.join
            )

    return big_storage


def annotation_to_columns(annotation):
    """
    Defines an annotation-to-columns conversion function for UCSC visualization
    data analysis function. This function takes an annotation as input and
    returns the annotation split into arrays, so that UCSC can display details
    of a binding site as a table.

    :param annotation: input annotation from a binding site.

    """

    rows = annotation.split(ANNOTATION_ROW_DELIMITER)
    array = [tuple(r.split(ANNOTATION_COLUMN_DELIMITER)) for r in rows]
    array = list(set(array))
    len_row = len(array[0])
    no_of_rows = len(array)

    array_of_strings = [
        '______'.join([array[i][j] for i in range(no_of_rows)])
        for j in range(len_row)
    ]

    return array_of_strings


keys = ["rbpdb", "attract", "rbpmap", "postar", "custom"]
data_source_annotation_to_columns = {k: annotation_to_columns for k in keys}

if __name__ == '__main__':
    pass
