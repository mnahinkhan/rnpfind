# # This file deals with making dictionaries that map gene names to a "default" gene name among its group of synonyms.
# # I do not know if this file is serving a purpose in the project anymore or not.

# # import pandas as pd

# def build(ncbi_gene_files):
#     # Declare a dictionary that maps gene synonyms to an official name:
#     synonym_dict = {}

#     syn_sheets = map(lambda f: pd.read_excel(f, keep_default_na=False),
#                      ncbi_gene_files)

#     for syn_sheet in syn_sheets:
#         for row in syn_sheet[["Symbol", "Synonyms"]].itertuples():
#             official_symbol = row[1].upper()

#             # All official symbols should map to themselves:
#             if official_symbol not in synonym_dict:
#                 synonym_dict[official_symbol] = official_symbol
#             elif official_symbol not in synonym_dict[official_symbol]:
#                 synonym_dict[official_symbol] = (synonym_dict[official_symbol]
#                                                  + "," + official_symbol)

#             # All synonyms represent the official symbol:
#             synonyms = row[2].upper().split("|")

#             for synonym in synonyms:

#                 # Blank
#                 if synonym == "-":
#                     continue

#                 # Otherwise synonym represents official symbol:
#                 if synonym not in synonym_dict:
#                     synonym_dict[synonym] = official_symbol
#                 elif official_symbol not in synonym_dict[synonym]:
#                     synonym_dict[synonym] = (synonym_dict[synonym]
#                                              + "," + official_symbol)

#     return synonym_dict


# # The function below uses the dictionary built above to convert any
# # given gene into the official representative gene name.

# def load_synonym_func(synonym_dict):
#     def synonym_func(gene):
#         if gene in synonym_dict:
#             return synonym_dict[gene]

#         # Some names contain hyphens and can be simplified like this:
#         if gene.replace("-", "") in synonym_dict:
#             return synonym_dict[gene.replace("-", "")]

#         # Some of the entries in the ATTRACT database end with the "*" symbol
#         # to denote mutations:
#         if gene[-1] == "*" and gene[:-1] in synonym_dict:
#             return synonym_func(gene[:-1]) + "*"

#         # This is one entry I looked up the simplification for:
#         if gene == "YBX2-A":
#             return synonym_dict["YBX2"]

#         # Just warn of genes that are not being dealt with:
#         print("WARNING:", gene, "is unaccounted for in the synonyms")

#         return gene

#     return synonym_func
