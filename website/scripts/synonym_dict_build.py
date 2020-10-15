# # Full confession: I'm not entirely sure if this piece of code and its related instances scattered throughout the
# # project are still relevant or not. The original purpose was to allow for data sources to add RBPs that are
# # seemingly different but are in fact the same, since in biology we have a lot of gene synonyms. Another use could be
# # at the querying level, where we could fetch the binding sites of an RBP from a Storage container variable but allow
# # for the query to be a synonym that is automatically standardized. I'm not sure if that was done, but it may have been.
# from config import dictionaryBuild, refreshSynonymDict, ncbi_gene_files, ncbi_gene_dir
# import pickle
# import gene_synonyms


# def deal_with_dictionary_building():
#     if dictionaryBuild and not ('synonym_dict' in vars() and not refreshSynonymDict):

#         print("Building the synonym dictionary...")

#         if refreshSynonymDict:

#             ncbi_gene_file_paths = [ncbi_gene_dir + s for s in ncbi_gene_files]

#             synonym_dict = gene_synonyms.build(ncbi_gene_file_paths)
#             # TODO: Consider rerouting this file to a better location
#             with open('gene_synonyms.pickle', 'wb') as handle:
#                 pickle.dump(synonym_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

#             print("Building complete!")

#         else:
#             with open('gene_synonyms.pickle', 'rb') as handle:
#                 synonym_dict = pickle.load(handle)

#         synonym_func = gene_synonyms.load_synonym_func(synonym_dict)

#     elif not dictionaryBuild:
#         print("Synonym dictionary building skipped as specified...")
#         synonym_func = -1
#     else:
#         print("Synonym dictionary already exists, skipping building..."
#               + "set refreshSynonymDict to True to change this!")
#         synonym_func = -1

#     return synonym_func
