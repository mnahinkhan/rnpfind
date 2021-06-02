# pylint: disable=no-member
"""
This module defines functions used to generate and serve resources when clients
make requests to specific URLs (the mappings are defined in urls.py)

Important views needed for RNPFind are:
 1. A view used to serve the index page
 2. A view used for clients to request analysis of a gene
 3. A view used for clients to ask on the status of a previously requested
    analysis
"""
from django.shortcuts import render

# from time import sleep

from django.http import (
    # HttpResponse,
    # HttpResponseRedirect,
    JsonResponse,
)

# from django.urls import reverse

# RNPFind source script files
from website.scripts.gene_coordinates import gene_to_coord
from website.scripts.load_data import load_data
from website.scripts.analysis_functions import analysis_method_functions

from .models import AnalysisStatus, BindingSummaryInfo, Gene

# total number of steps it takes to analyse a gene - useful for outputting
# status messages.
TOTAL_STEPS = 13


def index(request, bad_gene=""):
    """
    Serve the index page to the client. Provide a list of analysed genes from
    the database. Optionally, provide a "bad gene" name in case the user
    recently requested analysis of a non-existent gene.
    """
    return render(
        request,
        "website/index.html",
        {"bad_gene": bad_gene, "analysed_genes": Gene.objects.all()},
    )


def gene_page_redirector(request):
    """
    For when a GET request is made by the client from the search bar, instead of
    the needed gene name being in the URL. This view extracts the requested
    gene name from the body of the GET rquest and redirects to the gene-page
    view.
    """
    gene_name = request.GET["gene_name"]
    return gene_page(request, gene_name)


def gene_page(request, gene_name):
    """
    Serves up a page dedicated to the requested gene name. If the requested gene
    is not a real gene, it serves the index page with a "bad gene" marker.
    """
    rna_info = gene_to_coord(gene_name)
    success = rna_info["success"]

    if not success:
        return index(request, bad_gene=gene_name)

    official_name = rna_info["official_name"]
    is_in_database = len(Gene.objects.filter(name=official_name)) != 0

    return render(
        request,
        "website/gene_page.html",
        {
            "cached": is_in_database,
            "gene_official_name": official_name,
            "gene_typed_name": gene_name,
            "gene_cached_data": (
                Gene.objects.get(name=official_name) if is_in_database else 0
            ),
        },
    )


def analysis_status(request, request_id):
    """
    Returns the current status of a request as a JSON object.
    """
    del request

    if len(AnalysisStatus.objects.filter(request_id=request_id)) > 0:
        status = AnalysisStatus.objects.get(request_id=request_id).status
    else:
        status = ""

    return JsonResponse(
        {"status": status if status else "0/13. Sending request to server for analysis"}
    )


def out_fn_gen(request_id):
    """
    Generates a function that allows for emitting status updates of an analysis,
    given the request id for the analysis. The "out" function returned may be
    called multiple times.

    Internally, the out function updates a database entry with every output
    (status update) message emitted. This allows for requests that check on
    the current status of an analysis.
    """

    def out(msg):
        if len(AnalysisStatus.objects.filter(request_id=request_id)) > 0:
            analysis_status_obj = AnalysisStatus.objects.get(request_id=request_id)
            analysis_status_obj.status = msg
            analysis_status_obj.save()
        else:
            AnalysisStatus(request_id=request_id, status=msg).save()

    return out


def analysis_request(request, gene):
    """
    Given a gene name, analyses its transcript.

    This function is adopted from analysis_script() in main.py (the command
    line tool). It specifically populates the binding sites for the gene from
    various (three) data sources, and then returns a JSON object to confirm
    completion. The client may then refresh their page to see the newly
    populated data for their transcript of interest.
    """
    del request

    rna_info = gene_to_coord(gene)
    assert rna_info["success"]
    assert rna_info["official_name"] == gene
    assert len(Gene.objects.filter(name=rna_info["official_name"])) == 0

    if len(AnalysisStatus.objects.filter(request_id=gene)) > 0:
        return JsonResponse({"process_complete": False})
    gene_entry = Gene(name=gene)

    request_id = gene
    out = out_fn_gen(request_id)

    # collect binding sites
    data_load_sources = ["rbpdb", "attract", "postar"]
    big_storage = load_data(
        data_load_sources, rna_info, out=out, total_steps=TOTAL_STEPS
    )

    out(f"4/{TOTAL_STEPS}. Creating visualizations on UCSC")

    analysis_method = "ucsc"
    analysis_method_function = analysis_method_functions[analysis_method]
    ucsc_url = analysis_method_function(
        big_storage, rna_info, out=out, total_steps=TOTAL_STEPS
    )
    gene_entry.ucsc_url = ucsc_url

    out(f"11/{TOTAL_STEPS}. Caching computed values")
    gene_entry.chromosome_number = rna_info["chr_n"]
    gene_entry.start_coord = rna_info["start_coord"]
    gene_entry.end_coord = rna_info["end_coord"]

    gene_entry.save()

    out(f"12/{TOTAL_STEPS}. Generating summary data")
    # Save binding site summary info
    for data_source in big_storage:
        no_unique_rbps, no_unqiue_sites = big_storage[data_source].summary()
        BindingSummaryInfo(
            data_source_type=data_source,
            gene=gene_entry,
            number_of_sites=no_unqiue_sites,
            number_of_rbps=no_unique_rbps,
        ).save()

    assert len(Gene.objects.filter(name=rna_info["official_name"])) == 1
    out(
        f"13/{TOTAL_STEPS}. Complete! Refreshing your page now..."
        " If it does not automatically, please refresh this page!"
    )
    return JsonResponse({"process_complete": True})
