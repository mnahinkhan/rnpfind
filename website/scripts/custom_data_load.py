from .custom_binding_data import custom_data


def custom_data_load(rna_info):
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    for rbp, binding_sites in custom_data[RNA].items():
        for start, end in binding_sites:
            annotation = ""
            yield rbp, start, end, annotation
