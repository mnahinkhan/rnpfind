"""
An old file that I did not integrate into RNPFind. But maybe it has potential?
(probably not to anyone else)
"""

# from bind_analysis import BindingSites
from operator import itemgetter

# from merge_annotation_funcs import generate_merge_func
# import statistics  # Standard deviation mostly
# import operator

firstTwoItems = itemgetter(0, 1)
firstItem = itemgetter(0)
secondItem = itemgetter(1)
thirdItem = itemgetter(2)


def per_binding_site_analysis(big_storage, rna_info):
    """
    An analysis function that generates information on each binding site bound
    to on an RNA of interest by an RBP of interest. Information for each binding
    site includes the binding sites of other RBPs, competitive/cooperative
    relations, etc.

    :param big_storage: a dictionary indexed by data sources that correspond to
        Storage instances containing binding site data extracted from said data
        source.
    :param rna_info: A dictionary containing data on RNA of interest, such as
        its name and location on the human genome.

    """
    # This function previously required these inputs
    # bigStorage, AUF1Filter, analysis_per_binding_site_window,
    # analysis_per_binding_site_competitive_range, analysis_sources

    # Variable that comes useful later.
    print("unimplemented for now!")
    print(big_storage, rna_info, "", "")

    # TODO:This function needs to be reconsidered with respect to RNPFind and
    # generalized for one template

    # binding_info = {}
    #
    # # CHECK IF THIS SELECTION IS IMPORTANT OR NOT
    # # Is this selection important??
    # # storageSpace = select(bigStorage,analysis_sources)
    #
    # # For each lncRNA, we want to analyze the AUF1 binding sites
    # lncRNA = RNA
    # lncRNA_storage = storageSpace
    #
    # print("")
    # print("*****************")
    # print("analysing AUF1 sites on", lncRNA + "...\n")
    # print("*****************")
    #
    # # Get the RBPs that bind close to sites on AUF1
    # filtered_storage_dict = lncRNA_storage.sites_analysis(
    #     "AUF1", bp_threshold=analysis_per_binding_site_window
    # )

    # # Sort AUF1 sites by readCount
    # auf1_sites = sorted(filtered_storage_dict.keys(),
    #                     key=lambda k: int(k[2].split()[1]), reverse=True)
    #
    # # Variables for later analysis
    # distance_read_collection = {}
    #
    # # Analyse each site
    # for site in auf1_sites:
    #
    #     # Number of times the binding site was detected in RNASeq
    #     readCount = int(site[2].split()[1])
    #
    #     # Get collection of RBPs and their binding sites close to the current
    #     # AUF1 binding site being considered
    #     dB_storage = filtered_storage_dict[site]
    #
    # print(
    #     "\nanalysing binding site from",
    #     str(site[0]) + "bps to " + str(site[1]) + "bps; readCount =",
    #     readCount, 'from', site[2].split()[0]
    # )

    #     # We'll collect data in a list then sort them later.
    #     listOfThingsToPrint = []
    #
    #     # For each RBP and its binding sites:
    # for rbp, binding_sites in dB_storage.items():
    #     # Get the distances away from the current AUF1 site being considered
    #     distances = list(
    #         map(lambda s: BindingSites.distance(s, site), binding_sites)
    #     )

    #     # Separate binding sites and annotation:
    #     closest_sites = list(map(firstTwoItems, binding_sites))
    #     annotations = list(map(thirdItem, binding_sites))

    #     # Check if the RBP has experimental evidence for binding AUF1:
    #     isBindAUF1 = AUF1Filter(rbp)

    #     listOfThingsToPrint.append(
    #         (rbp, closest_sites, distances, annotations, isBindAUF1))

    #     # Variables to analyse later
    #     distance_read_collection[rbp] = (
    #         distance_read_collection.get(rbp, [])
    #         + list(map(lambda d: (d, readCount), distances))
    #     )

    # # Now sort and print!
    # for toPrint in sorted(listOfThingsToPrint,
    #                       key=lambda k: min(k[2])):

    #     # Unpack
    #     rbp, closest_sites, distances, annotations, isBindAUF1 = toPrint

    #     # Modify
    #     closest_sites = str(closest_sites)[1:-1]
    #     isBindAUF1 = 'yes' if isBindAUF1 else 'no'
    #     isAllCompetitive = (
    #         all(map(
    #         lambda k: k <= analysis_per_binding_site_competititive_range,
    #         distances
    #         ))
    #     )

    #     isAllCooperative = (
    #         all(map(
    #         lambda k: k >= analysis_per_binding_site_competititive_range,
    #         distances))
    #     )
    #     isAllCooperativeORCompetitive = (
    #         'competitive' if isAllCompetitive else
    #         ('cooperative' if isAllCooperative else 'neither')
    #     )
    #     distances = str(distances)[1:-1]

    #     # TODO: Consider methods of tracking datasources despite multiple
    #     # source merges

    #     data_load_source = analysis_sources[0]
    #     merge_func = generate_merge_func(data_load_source)

    #     annotations = merge_func(annotations)

    #     # Repack
    #     toPrint = (
    #         rbp, closest_sites, distances, annotations, isBindAUF1,
    #         isAllCooperativeORCompetitive
    #     )

    #         if False:
    #             # Print
    #             print(";".join(map(str, toPrint)))
    #
    #     print("*****************")
    #     print("*****************")
    #     print("*****************")
    #
    #     if True:
    #         info_collection = {}
    #
    #         def sorter(item):
    #             rbp = item[0]
    #             dist_read = item[1]
    #             flat_read_list = []
    #             for dist, read in dist_read:
    #                 flat_read_list += [dist] * read
    #
    #             reads = list(map(secondItem, dist_read))
    #             deviation = (
    #                 statistics.stdev(flat_read_list) if len(reads) > 1 else 100
    #             )
    #             return (
    #                 -sum(reads) / deviation if deviation != 0 else -sum(reads)
    #             )

    #             # return (-sum(reads), deviation)

    #         for rbp, dist_read in sorted(
    #             distance_read_collection.items(), key=sorter
    #         ):

    #             # Only because weighted sum standard deviation is difficult to
    #             # calculate:
    #             flat_read_list = []
    #             for dist, read in dist_read:
    #                 flat_read_list += [dist] * read

    #             reads = list(map(secondItem, dist_read))
    #             distance = list(map(firstItem, dist_read))
    #             deviation = (
    #                 statistics.stdev(flat_read_list) if len(reads) > 1 else 100
    #             )

    #             info_collection[rbp] = (
    #                 rbp, statistics.median(flat_read_list), deviation,
    #                 sum(reads), len(set(reads)), len(reads), min(distance),
    #                 max(distance), ['no', 'yes'][AUF1Filter(rbp)],
    #                 str(dist_read)[1:-1][:100]
    #             )

    #             if True:
    #                 print(";".join(map(str, info_collection[rbp])))

    #     binding_info[lncRNA] = (filtered_storage_dict, auf1_sites,
    #                             distance_read_collection, info_collection)

    # countofInvolvement = {}
    # for lncRNA in ["Neat1", "Malat1"]:
    #     distance_read_collection = binding_info[lncRNA][2]
    #     for rbp, dist_read in distance_read_collection.items():
    #         competitiveScore = (
    #             sum(map(
    #             lambda t: (
    #                 t[1] if t[0] < analysis_per_binding_site_competititive_range
    #                 else 0
    #             ),
    #             dist_read
    #             ))
    #         )

    #         cooperativeScore = (
    #             sum(map(
    #             lambda t: (
    #                 t[1] if t[0] > analysis_per_binding_site_competititive_range
    #                 else 0
    #             ),
    #             dist_read
    #             ))
    #         )

    #         if rbp not in countofInvolvement: countofInvolvement[rbp] = [
    #             0, 0, 0, 0
    #         ]

    #         if lncRNA == "Malat1":
    #             countofInvolvement[rbp][0] += competitiveScore
    #             countofInvolvement[rbp][1] += cooperativeScore
    #         elif lncRNA == "Neat1":
    #             countofInvolvement[rbp][2] += competitiveScore
    #             countofInvolvement[rbp][3] += cooperativeScore
    #         else:
    #             raise ValueError('Fatal Error, please debug')

    # # normalization = [0]*4
    # # for scores in countofInvolvement.values():
    # # 	for j in [0,1,2,3]:
    # # 		normalization[j] = max(scores[j],normalization[j])

    # total_neat1_reads = (
    #     sum(map(lambda k: int(k[2].split()[1]), storageSpace['Neat1']["AUF1"]))
    # )
    # total_malat1_reads = (
    #     sum(map(lambda k: int(k[2].split()[1]), storageSpace['Malat1']["AUF1"]))
    # )
    # normalization = [total_malat1_reads, total_malat1_reads,
    #                  total_neat1_reads, total_neat1_reads]
    # print(normalization, "-->normalization")

    # def sorter(item):
    #     rbp = item[0]
    #     scores = item[1]
    #     norm_scores = list(map(operator.truediv, scores, normalization))
    # What should I prioritize?
    # scores [a , b, c, d]
    # 		 ^   ^  ^  ^
    # 		 |   |  |  |
    # 		/   /    \  \
    # 	Malat1 		  Neat1
    # comp., cooper     comp., cooperative...

    # Presumably, Malat1 and Neat1 have a similar level
    # of binding to AUF1... In that case, perhaps there isn't
    #     # as much competitive inhibition, and the thing that
    #     # destabiilzes Neat1 but not Malat1 is a cooperative helper
    #     # that's with Neat1 but not on Malat1.
    #     # Alternatively, if this presumption is wrong, we can expect
    #     # a competitive RBP to be affecting on Malat1 but not on
    #     # Neat1.
    #     # We conclude two ways to evaluate interest:
    #     #	1. Maximize (Cooperation on Neat1 - Cooperation on Malat1)
    #     #	2. Maximize (Competitive on Malat1 - Competitive on Neat1)
    #
    #     # For now, i've chosen to maximize both at the same time:
    #
    #     comp_malat1 = norm_scores[0]
    #     coOP_malat1 = norm_scores[1]
    #     comp_neat1 = norm_scores[2]
    #     coOP_neat1 = norm_scores[3]
    # return max(abs(coOP_neat1 - coOP_malat1), abs(comp_malat1 - comp_neat1))
    #
    # for rbp, scores in sorted(countofInvolvement.items(),
    #                           key=sorter, reverse=True):
    #     norm_scores = list(map(operator.truediv, scores, normalization))
    #     item = ([rbp] + scores + ['yes' if AUF1Filter(rbp) else 'no'] +
    #             norm_scores)
    #     print(";".join(map(str, item)))
