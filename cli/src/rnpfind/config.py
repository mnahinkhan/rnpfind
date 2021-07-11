"""
Stores most of the configurable variables for RNPFind for ease of access.
"""

from pathlib import Path

# Building dictionary takes time. Set to false to use an identity function
# instead.
# dictionaryBuild = True
# ncbi_gene_dir = "./website/data/NCBI Gene Info/"

ncbi_gene_files = [
    "Human Genes and Synonyms.xlsx",
    "Saccharomyces_cerevisiae.xlsx",
    "Drosophila_melanogaster.xlsx",
]

# refreshSynonymDict = False

# Loading the data for RBP binding sites on lncRNA, as well as BIOGRID
# protein-protein interactions.
# dataLoad = True


# Analysis takes time. Set to true to analyse:
# analysis_overall_correlation = False
# What is the base pair stringency you want to analyse with?
# analysis_threshold_bps = [10, 15, 30, 50]

# The newer form of analysis: per binding site
# analysis_per_binding_site = False
# Window of bp range to look outside each AUF1 site
# analysis_per_binding_site_window = 56
# Threshold for how much is considered to be competitive
# analysis_per_binding_site_competititive_range = 15

# Some of the variables that take long to load will not be loaded again if this
# parameter is set to False. Namely, neat1_storage, malat1_storage, and
# synonym_dict.

# refreshStorages = True

# Add modified custom data?
# custom_data_add = True

EXPERIMENTAL_BINDING_SITE_ACCEPTABLE_COVERAGE_RATIO = (
    1  # Professor Ihab has decided
)

# Genome version being used
GENOME_VERSION = "hg38"

UCSC_TRACK_VISIBILITY = "dense"

# Used for csv output format
DEFAULT_BASE_STRINGENCY = 30

# This should be False normally. Only set to True if you want to make a
# dedicated directory for your RNA of interest.
# Note that this flag only affects UCSC browser visualization method.
DEDICATED_ANALYSIS = False

# Exactly how strong is a letter when written as a motif, compared with other
# letters? How many times more likely?

# I believe any value from 2-10 makes no difference in practice
RBPDB_MOTIF_PWM_LETTER_STRENGTH = 4
RBPDB_MOTIF_N_REPEAT_REQ = 6

PWM_SCAN_CUT_OFF_PERCENTAGE = 0.80

# This is data loading specific. If you are implementing a data load function,
# make sure your data does not have
# the following substrings in their annotations for binding sites. If they do,
# you should change the strings below!
ANNOTATION_COLUMN_DELIMITER = ",,,,,"
ANNOTATION_ROW_DELIMITER = ";;;;;"

RO_DATA_TAR_NAME = Path(__file__).parent / "all.tar.gz"
RO_DATA_URL = "https://rnpfind.com/ro-data/all.tar.gz"

# Path for all cached data
CACHE_PATH = Path(__file__).parent / "cache"
# Path for pickled data
PICKLE_PATH = f"{CACHE_PATH}/pickles"
# Path for all (input) read-only data
RO_DATA_PATH = Path(__file__).parent / "ro-data"
# Path for RBPDB data
RBPDB_PATH = f"{RO_DATA_PATH}/rbpdb"
# Path for ATTRACT data
ATTRACT_PATH = f"{RO_DATA_PATH}/attract"
# Path for POSTAR data
POSTAR_PATH = f"{RO_DATA_PATH}/postar"
# Path for AutoSQL files
AUTOSQL_PATH = f"{RO_DATA_PATH}/autosql"
# Path for UCSC tools
UCSCTOOL_PATH = Path(__file__).parent / "ucsc-tools"
