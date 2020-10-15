from .config import annotation_column_delimiter
from .picklify import picklify
from .pwm_scan import get_human_seq, str_to_pwm, pwm_scan

attract_all_column_names = ["Gene_name", "Gene_id", "Mutated", "Organism", "Motif", "Len", "Experiment_description",
                            "Database", "Pubmed", "Experiment", "Family", "Matrix_id", "attractScore"]
attract_all_column_descriptions = attract_all_column_names
attract_columns_of_interest = [2, 4, 6, 7, 8, 9, 10, 11, 12]
attract_default_label_index = [8]
attract_default_mouse_over_index = 9

attract_column_names = [attract_all_column_names[i] for i in attract_columns_of_interest]
attract_column_descriptions = [attract_all_column_descriptions[i] for i in attract_columns_of_interest]


def generate_matrix_to_pwm_dict():
    attract_pwm_file_path = "./website/data/attract-pwm.txt"
    matrix_to_pwm_dict = {}
    with open(attract_pwm_file_path) as handle:
        s = handle.readline()
        while s:
            matrix_id = s.split()[0][1:]
            raw_pwm_str = ""
            s = handle.readline()
            while s and s[0] != ">":
                raw_pwm_str += s
                s = handle.readline()
            # Now we have the raw text, we convert it to pwm and add to dictionary
            matrix_to_pwm_dict[matrix_id] = str_to_pwm(raw_pwm_str)
    return matrix_to_pwm_dict


def attract_data_load(rna_info, out=None):
    attract_protein_file_path = "./website/data/ATtRACT_db.txt"
    RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord = rna_info
    rna_seq = get_human_seq(RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord)
    matrix_to_pwm_dict = picklify(generate_matrix_to_pwm_dict)
    with open(attract_protein_file_path) as handle:
        columns = handle.readline().strip().split('\t')
        assert (columns == ["Gene_name", "Gene_id", "Mutated", "Organism", "Motif", "Len", "Experiment_description",
                            "Database", "Pubmed", "Experiment_description", "Family", "Matrix_id", "Score"])
        protein_columns = handle.readline().replace("\n", "").split('\t')
        while protein_columns != ['']:
            # Warning: Score ends with \n here, maybe remove using strip or indexing. For now, we don't care about
            # score as it seems to be about literature

            # We only care about human RBPs for now.
            if protein_columns[3] != "Homo_sapiens":
                protein_columns = handle.readline().replace("\n", "").split('\t')
                continue
            annotation = annotation_column_delimiter.join([protein_columns[i] for i in attract_columns_of_interest])

            rbp = protein_columns[0]

            matrix_id = protein_columns[11]

            pwm = matrix_to_pwm_dict[matrix_id]
            sites = pwm_scan(rna_seq, pwm)
            if not sites:
                protein_columns = handle.readline().replace("\n", "").split('\t')
                continue

            for start, end in sites:
                yield rbp, start, end, annotation

            protein_columns = handle.readline().replace("\n", "").split('\t')
