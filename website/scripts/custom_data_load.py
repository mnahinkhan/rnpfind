from .custom_binding_data import custom_data


def custom_data_load(rna_info):
    RNA = rna_info['official_name']
    for rbp, binding_sites in custom_data[RNA].items():
        for start, end in binding_sites:
            annotation = ""
            yield rbp, start, end, annotation
