"""
A file dedicated to useful functions related to position weight matrices that
describe the binding profiles of RBPs, such as scanning an RNA sequence to
detect where the RBP might bind on it. Possibly the file I had the most pleasure
writing, despite its simplicity.

What is a PWM? A PWM is a position weight matrix that summarizes a collection of
strings. In the context of biology, an example could be a collection of
sequences that an RBP was found to bind to protein. Say, for example, that an
RBP called HNRNPK was found to bind to 100 RNA sequences, each of length 5
bases. We could then summarize, for each of those 5 bases, the frequency of
occurrence of a nucleotide base in that position. An example binding profile
could be:

Base / Position
  1   2   3   4   5
A 50  100 100 0   0
G 50  0   0   60  100
C 0   0   0   30  0
T 0   0   0   10  0

Notice how each column adds up to 100, since we had a collection of 100. We can
conclude that this RBP binds to a 5-letter sequence that looks like
(A/G)AA(G/C/T)G, with the distribution shown above (e.g. first letter is equally
likely to be A or G, and so on). Technically, what is shown above is a position
frequency matrix. It is common to normalize the above so each column adds up to
1, to generate a position probability matrix. For some purposes, one might take
the log values of the above too. I think all these things can be called PWMs.
Thus in this project, I have referred to any form of position frequency matrix
(normalized or not) as a PWM. We do not deal with logs here!

Clearly, pwms are useful for making predictions about where an RBP might bind to
an RNA sequence. I did not do extensive research on statistics here to make such
predictions; instead I used a simplest model described on a Wikipedia page,
where one assumes an equiprobable base distribution for prior belief and
independence between motif base positions in terms of nucleotide preference.

Specifically, we decide if an rbp binds to a specific location on an RNA
sequence by associating a "score" with the rbp binding profile on a specific
location, calculated by multiplying all the position weights that show up on the
RNA molecule. The RBP is considered to bind if the score calculated is more than
80% of the highest possible score attainable (cutoff percentage can be adjusted
in config.py).

We could have chosen a number of ways to represent PWM's. Somewhat arbitrarily,
I went with a pwm variable being represented by a dictionary that maps the
letters of the nucleotide bases ("A", "C", etc.) to a list of length equal to
the length of the motif. The list has numbers that represent frequency (not
assumed to be probabilities) of the base at that position. For example,
pwm["T"][6] is the frequency of the base T on the 7th nucleotide of the motif
represented by "pwm".

"""

# TODO: consider defining a class for PWM's for ease of interface


from functools import partial
from pathlib import Path

from .config import GENOME_VERSION, PWM_SCAN_CUT_OFF_PERCENTAGE, RO_DATA_PATH

bases = ["A", "G", "C", "T"]

seqs = {}


def get_human_seq(rna_info):
    """
    Returns the human genomic sequence of the specified chromosome from the
    start to the end coordinate.

    :param rna_info: a dictionary of information regarding the location on the
        human genome from which the nucleotide sequence is sought. Keys provided
        should include:
            'chr_n' : chromosome number
            'start_coord' : start coordinate on the chromosome
            'end_coord' : end coordinate on the chromosome

    """
    chr_no = rna_info["chr_n"]
    chr_start = rna_info["start_coord"]
    chr_end = rna_info["end_coord"]

    if (chr_no, chr_start, chr_end) in seqs:
        return seqs[(chr_no, chr_start, chr_end)]

    chr_files_dir = Path(RO_DATA_PATH) / f"{GENOME_VERSION}-human-genome"
    chr_file = chr_files_dir / f"chr{chr_no}.fa"

    with open(chr_file) as handle:
        # Need to skip the first line, and then a bunch more to reach the line
        # that has the chr_start numbered base. Note that the .fa file comes in
        # lines of 50 nucleotides with "\n" at the end. Admittedly, a part of
        # the below formula was trial and error until it was empirically correct
        # (maybe the -1 part).

        handle.seek(
            4
            + len(str(chr_no))
            + 1
            + 51 * (chr_start // 50)
            + chr_start % 50
            - 1
        )

        # Once at the required position, keep adding lines of bases as long as
        # too short. Then adjust in case we overshot the required length.
        gene = ""
        while len(gene) < chr_end - chr_start:
            gene += handle.readline().strip()
        gene = gene[: chr_end - chr_start]

    seq = gene.upper()
    seqs[(chr_no, chr_start, chr_end)] = seq
    return seq


def product(list_of_numbers):
    """
    Returns a product of numbers.

    :param list_of_numbers: the input list of numbers whose product is sought.

    """
    answer = 1
    for number in list_of_numbers:
        answer *= number
    return answer


def pwm_scan_naive_brute_force(gene, pwm):
    """
    Scans gene, a string of nucleotide bases, for positions where pwm could bind
    using a naive brute force approach. Slower than pwm_scan most of the times,
    but much faster when the degree of freedom (total number of strings
    represented) of pwm is very high (recommended to use this if degree of
    freedom of pwm is > 2000).

    :param gene: a string of nucleotide bases
    :param pwm: a position weight matrix (for format see (no clue))

    """
    len_gene = len(gene)
    len_pwm = len(pwm["A"])
    highest_scores = [max([pwm[b][i] for b in bases]) for i in range(len_pwm)]
    max_score = product(highest_scores)
    cut_off_percentage = 0.80

    cut_off_threshold = cut_off_percentage * max_score

    binding_sites = []
    for i in range(len_gene - len_pwm + 1):
        score = 1
        for j in range(len_pwm):
            score *= pwm[gene[i + j]][j]
            if score < cut_off_threshold:
                break

        if score >= cut_off_threshold:
            binding_sites += [(i, i + len_pwm)]
    return binding_sites


def pwm_degree_of_freedom(pwm):
    """
    Counts the number of binding strings described by a pwm. Can be useful in
    picking between pwm_scan_naive_brute_force and pwm_scan

    :param pwm: a position weight matrix (for format see (no clue))

    """
    # I think that this function can sometimes overestimate the
    # "degree of freedom", because it's just multiplying by the number of bases
    # possibly described at each position by checking if it's above the cutoff
    # percentage of the maximum score at that position, instead of looking
    # holistically. So this function should be considered an approximation.
    len_pwm = len(pwm["A"])
    freedom = 1
    for i in range(len_pwm):
        max_score = max([pwm[b][i] for b in bases])
        freedom *= len(
            [
                1
                for b in bases
                if pwm[b][i] >= PWM_SCAN_CUT_OFF_PERCENTAGE * max_score
            ]
        )
    return freedom


def findall(needle, haystack):
    """
    Yields all the positions of the pattern needle in the string haystack.
    :param needle: substring to be searched for
    :param haystack: the string to be searched in (for needle)

    """
    # Adopted from: https://stackoverflow.com/a/34445090/8551394
    i = haystack.find(needle)
    while i != -1:
        yield i
        i = haystack.find(needle, i + 1)


def pwm_scan(gene, pwm):
    """
    Scans gene, a string of nucleotide bases, for positions where pwm could
    bind. The approach used here is to generate all possible strings of bases
    that the pwm could represent, and use Python's inbuilt substring searching
    algorithm to get all places where the collection of strings could bind on
    the RNA string. Usually fast, but if the pwm represents too many strings
    (e.g. your pwm has a stretch of 15 bases with equi-possible bases
    (think "N"), this algorithm will likely not terminate ever). Consider using
    pwm_scan_naive_brute_force instead in such occasions.

    :param gene: a string of nucleotide bases
    :param pwm: a position weight matrix (for format see (no clue))

    """
    len_pwm = len(pwm["A"])
    cut_off_percentage = PWM_SCAN_CUT_OFF_PERCENTAGE

    possible_seqs = [("", 1)]

    for i in range(len_pwm):
        new_possible_seqs = []
        max_base_score = max([pwm[b][i] for b in bases])
        for seq, score in possible_seqs:
            for base in bases:
                if score * pwm[base][i] / max_base_score >= cut_off_percentage:
                    new_possible_seqs += [
                        (seq + base, score * pwm[base][i] / max_base_score)
                    ]
        possible_seqs = new_possible_seqs

    possible_seqs = [x for x, y in possible_seqs]

    binding_sites = []
    for substring in possible_seqs:
        binding_sites += [(i, i + len_pwm) for i in findall(substring, gene)]
    return binding_sites


def motif_to_pwm(motif, letter_strength=4):
    """
    Converts a motif (a string like "AGCNNNYWS") to a corresponding pwm. In
    making the pwm, letter_strength (defaulted to 4) represents the number of
    times more likely a base is considered to be when represented by the motif
    string at a particular position. For example, in "AGCNNNYWS", the first
    position is A but we don't represent it as a 100% A frequency, instead we
    give A a probability of 4/7 while G,C,T each get a probability of 1/7
    (when letter_strength = 4).

    :param motif: a string representing a nucleotide motif, using standard
        bioinformatics notation (insert ref. link here)
    :param letter_strength: the number of times more likely a base is considered
        to be when represented by the motif string at a particular position.
        (Default value = 4)

    """

    # In practice, the letter_strength might not matter much if it's 4, since
    # our cutoff threshold overall for binding is 80% (so the moment a letter is
    # only 25% it fails). But maybe that's good!

    # I think the above comment talks about use case in RNPFind

    motif = motif.strip()
    motif = motif.replace("U", "T")

    # TODO: Clean up the code below to avoid repeating!

    pwm = {b: [] for b in bases}

    for nucleotide in motif:
        for base in bases:
            if nucleotide == "N":
                pwm[base].append(0.25)
            elif nucleotide == "Y":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2)
                    if base in "TC"
                    else 1 / (2 * letter_strength + 2)
                )
            elif nucleotide == "W":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2)
                    if base in "TA"
                    else 1 / (2 * letter_strength + 2)
                )
            elif nucleotide == "S":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2)
                    if base in "GC"
                    else 1 / (2 * letter_strength + 2)
                )
            elif nucleotide == "K":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2)
                    if base in "GT"
                    else 1 / (2 * letter_strength + 2)
                )
            elif nucleotide == "M":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2)
                    if base in "AC"
                    else 1 / (2 * letter_strength + 2)
                )
            elif nucleotide == "R":
                pwm[base].append(
                    letter_strength / (2 * letter_strength + 2)
                    if base in "GA"
                    else 1 / (2 * letter_strength + 2)
                )
            elif nucleotide == "D":
                pwm[base].append(
                    letter_strength / (3 * letter_strength + 1)
                    if base != "C"
                    else 1 / (3 * letter_strength + 1)
                )
            elif nucleotide == "H":
                pwm[base].append(
                    letter_strength / (3 * letter_strength + 1)
                    if base != "G"
                    else 1 / (3 * letter_strength + 1)
                )
            elif nucleotide == "B":
                pwm[base].append(
                    letter_strength / (3 * letter_strength + 1)
                    if base != "A"
                    else 1 / (3 * letter_strength + 1)
                )
            else:
                assert nucleotide in bases
                pwm[base].append(
                    letter_strength / (letter_strength + 3)
                    if nucleotide == base
                    else 1 / (letter_strength + 3)
                )

    return pwm


def pwm_summary(pwm):
    """
    Returns a string representing the most likely string represented by the pwm

    :param pwm: a position weight matrix (for format see (no clue))

    """
    len_pwm = len(pwm["A"])

    def get_pwm_base_index(base, index):
        return pwm[base][index]

    return_new = "".join(
        [
            max(bases, key=partial(get_pwm_base_index, index=i))
            for i in range(len_pwm)
        ]
    )

    return return_new


def str_to_pwm(raw_motif_str, is_transpose=False):
    """
    Converts a string containing numbers representing motif frequencies into a
    pwm. is_transpose allows for working with strings oriented differently. Go
    with default if your string has many rows and 4 columns; otherwise set
    is_transpose=True

    :param raw_motif_str: a string containing numbers representing motif
        frequencies
    :param is_transpose: set False if string has many rows and 4 columns,
        otherwise set to True. (Default value = False)

    """
    motif = raw_motif_str.split()
    motif = [float(k) for k in motif]
    base_index = {"A": 0, "G": 1, "C": 2, "T": 3}
    pwm = {b: [] for b in bases}

    for base in bases:
        for i in range(len(motif) // 4):
            if is_transpose:
                pwm[base].append(
                    motif[i + (len(motif) // 4) * base_index[base]]
                )
            else:
                pwm[base].append(motif[i * 4 + base_index[base]])

    return pwm


if __name__ == "__main__":

    # Example usage
    # from user_input import get_user_rna_preference

    # [RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord] = user_input()

    # start = timer()
    # rna_info = RNA, RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord

    # RNA_sequence = get_human_seq(
    #     RNA_chr_no, RNA_start_chr_coord, RNA_end_chr_coord
    # )

    # print(RNA_sequence)

    # example_motif = (
    #  '''0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538
    #     0.971153846154	0.00961538461538	0.00961538461538	0.00961538461538'''
    # )

    # test_pwm = pwm_str_to_dict(example_motif)

    # sites = pwm_scan(RNA_sequence, test_pwm)
    # for s, t in sites:
    #     print("binding site found!")
    #     print("Position: ", s, "to", t - 1)
    #     print("Sequence: ", RNA_sequence[s:t])

    # print("Number of sites:", len(sites))

    # end = timer()
    # print(end - start, "seconds elapsed")

    # print(pwm_motif_to_dict("(U)(2)"))

    print(
        pwm_scan(
            "CCCCCCTTTTCGCGCGCGCTTTTCGCGCGCGCGCGTTTTTTT", motif_to_pwm("TTTT")
        )
    )
    print(
        pwm_scan_naive_brute_force(
            "CCCCCCTTTTCGCGCGCGCTTTTCGCGCGCGCGCGTTTTTTT", motif_to_pwm("TTTT")
        )
    )
    print(pwm_degree_of_freedom(motif_to_pwm("NNNNN")))
    print(pwm_degree_of_freedom(motif_to_pwm("ACGTGTGTCG")))
    print(pwm_degree_of_freedom(motif_to_pwm("ACGTGTGDDDH")))
