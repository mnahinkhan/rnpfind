"""
Defines functions that are involved in getting some data from the user.
Used for the command line interface variant of RNPFind.
"""

import sys

from hgfind import WrongGeneName, hgfind

from .config import GENOME_VERSION
from .gene_coordinates import Chromosome


def parse_genome_coord(transcript: str) -> dict:
    """
    Given a string of the form <chr>:<start>-<end> representing genomic
    coordinates, extracts each of the chromosome number, start, and end
    coordinates.

    The input may be malformed.

    :param transcript: input genomic coordinates
    :returns: a dictionary containing the following keys and associated values:
        'success': whether the conversion worked or not. If this is False, the
            other keys below are not guaranteed to exist.
        'chr_n': the Chromosome object on which the gene lies
        'start_coord': the start coordinate of the gene on the chromosome
        'end_coord': the end coordinate of the gene on the chromosome
    """

    fail_output = {"success": False}

    # import pdb
    # pdb.set_trace()

    # Check if structure seems correct
    if len([c for c in transcript if c == ":"]) != 1:
        return fail_output
    if len([c for c in transcript if c == "-"]) != 1:
        return fail_output
    if transcript.find(":") > transcript.find("-"):
        return fail_output

    chr_n = transcript.split(":")[0]
    try:
        chr_n = Chromosome(chr_n)
        start_coord = transcript.split(":")[1].split("-")[0]
        end_coord = transcript.split(":")[1].split("-")[1]
        start_coord = int(start_coord)
        end_coord = int(end_coord)
    except ValueError:
        return fail_output

    if start_coord > end_coord:
        return fail_output

    success_output = {
        "success": True,
        "chr_n": chr_n,
        "start_coord": start_coord,
        "end_coord": end_coord,
    }

    return success_output


def get_rna_coord(rna_gene_name):
    """
    Given a gene name, returns its hg38 coordinate.
    The input may also be of the form <chr>:<start>-<end>.

    :param rna_gene_name: input RNA gene name (string)

    """

    try:
        rna_info = hgfind(rna_gene_name)
    except WrongGeneName:
        # Transcript might have been given in <chr no>:<start>-<end> form
        rna_info = parse_genome_coord(rna_gene_name)
        if not rna_info["success"]:
            print(
                "Error parsing transcript. Please give a correct gene name"
                "or a genomic coordinate of the form <chr>:<start>-<end>",
                file=sys.stderr,
            )
            sys.exit(1)
        if not "strand" in rna_info:
            rna_info["strand"] = "+"

    rna_chr_no = rna_info["chr_n"]
    rna_start_chr_coord = rna_info["start_coord"]
    rna_end_chr_coord = rna_info["end_coord"]
    rna_strand = rna_info["strand"]

    print(
        f"Analyzing {GENOME_VERSION}"
        f" {rna_chr_no}:{rna_start_chr_coord}-{rna_end_chr_coord}"
        f" ({rna_strand})"
        f" (length = {rna_end_chr_coord - rna_start_chr_coord + 1} bases)",
        file=sys.stderr,
    )

    return rna_chr_no, rna_start_chr_coord, rna_end_chr_coord, rna_strand


def get_user_rna_preference(transcript: str) -> dict:
    """
    Asks the user for the RNA molecule they wish to analyze (the name of the
    gene is expected). Converts the gene name into coordinates.

    :returns: a dictionary containing the official name of the RNA, along with
              its chrosome number, start coordinate, and end coordinate on the
              human genome.

    Start and end coordinates are 1-based, fully closed.
    """
    (
        rna_chr_no,
        rna_start_chr_coord,
        rna_end_chr_coord,
        rna_strand,
    ) = get_rna_coord(transcript)
    rna_info = {}
    rna_info["official_name"] = (
        transcript.upper()
        if ":" not in transcript
        else f"chr{rna_chr_no}-gene"
    )
    rna_info["chr_n"] = rna_chr_no
    rna_info["start_coord"] = int(rna_start_chr_coord)
    rna_info["end_coord"] = int(rna_end_chr_coord)
    rna_info["strand"] = rna_strand
    return rna_info
