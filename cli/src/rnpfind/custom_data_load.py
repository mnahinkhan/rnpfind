"""
This module defines a data load function. In particular, it loads a custom,
user-specified data that is available in custom_binding_data.py
"""
from .custom_binding_data import custom_data


def custom_data_load(rna_info):
    """
    Loads custom binding data specified in custom_binding_data.py and returns
    a Generator / Iterator object with the binding sites.

    :param rna_info: A dictionary containing the name of the RNA molecule,
        specified under its 'official_name' key.

    """
    rna = rna_info["official_name"]
    for rbp, binding_sites in custom_data[rna].items():
        for start, end in binding_sites:
            annotation = ""
            yield rbp, start, end, annotation
