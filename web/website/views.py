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
import sys

from django.http import JsonResponse  # HttpResponse,; HttpResponseRedirect,
from django.shortcuts import redirect, render
from hgfind import WrongGeneName, hgfind

from .models import AnalysisStatus, BindingSummaryInfo, Gene
from .tasks import analyze_gene, analyze_gene_done


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
    if not gene_name.strip():
        return redirect(index)
    return redirect(gene_page, gene_name=gene_name)


def gene_page(request, gene_name):
    """
    Serves up a page dedicated to the requested gene name. If the requested gene
    is not a real gene, it serves the index page with a "bad gene" marker.
    """
    try:
        rna_info = hgfind(gene_name)
    except WrongGeneName:
        return index(request, bad_gene=gene_name)

    official_name = rna_info["official_name"]

    if gene_name != official_name.lower():
        return redirect(gene_page, gene_name=official_name.lower())

    is_in_database = len(Gene.objects.filter(name=official_name)) != 0
    dbg_print(official_name)
    dbg_print(is_in_database)

    return render(
        request,
        "website/gene_page.html",
        {
            "cached": is_in_database,
            "gene_official_name": official_name,
            "gene_cached_data": (
                Gene.objects.get(name=official_name) if is_in_database else 0
            ),
            "binding_site_data": (
                BindingSummaryInfo.objects.filter(gene=official_name).exclude(
                    data_source_type=BindingSummaryInfo.TOTAL
                )
                if is_in_database
                else 0
            ),
        },
    )


def analysis_status(request, request_id):
    """
    Returns the current status of a request as a JSON object.
    """
    del request

    task = analyze_gene.AsyncResult(request_id.lower())

    if (
        task.state == "SUCCESS"
        and len(Gene.objects.filter(name=request_id.upper())) == 1
    ):
        state = "SUCCESS"
    else:
        state = "PENDING"

    if len(AnalysisStatus.objects.filter(request_id=request_id.upper())) > 0:
        status = AnalysisStatus.objects.get(
            request_id=request_id.upper()
        ).output
    else:
        status = "No analysis record"
        state = "UNREQUESTED"

    return JsonResponse({"message": status, "status": state})


def dbg_print(string: str):
    """
    Prints string to stderr
    """
    print(string, file=sys.stderr)


def analysis_result(request, gene):
    """
    Returns the status of the analysis as a JSON object
    """
    del request

    if len(Gene.objects.filter(name=gene.upper())) == 0:
        if len(AnalysisStatus.objects.filter(request_id=gene.upper())) == 0:
            return reject("unrequested")

        return reject("analysis ongoing")

    gene_obj = Gene.objects.get(name=gene.upper())
    return JsonResponse(
        {
            "request_accepted": True,
            "chr_no": gene_obj.chromosome_number,
            "start_coord": gene_obj.start_coord,
            "end_coord": gene_obj.end_coord,
            "ucsc_url": gene_obj.ucsc_url,
        }
    )


def reject(reason):
    """
    Returns a query rejection JSON object with customizable message
    """
    return JsonResponse({"request_accepted": False, "reason": reason})


def analysis_request(request, gene):
    """
    API callpoint for requesting analysis of genes
    """
    del request

    try:
        rna_info = hgfind(gene)
    except WrongGeneName:
        return reject("unrecognized gene")

    if rna_info["official_name"].lower() != gene.lower():
        return reject(f"Try {rna_info['official_name'].lower()} instead")

    if len(Gene.objects.filter(name=gene.upper())) != 0:
        return reject(f"{gene.upper()} has already been analyzed")

    if len(AnalysisStatus.objects.filter(request_id=gene.upper())) > 0:
        return reject(f"{gene.upper()} is currently being analyzed")

    task = analyze_gene.apply_async(
        args=[gene.upper()], task_id=gene.lower(), link=analyze_gene_done.s()
    )
    assert task.id == gene.lower()
    AnalysisStatus(
        request_id=gene.upper(), output="analysis request queued"
    ).save()
    assert len(AnalysisStatus.objects.filter(request_id=gene.upper())) > 0
    return JsonResponse({"request_accepted": True, "task_id": task.id})
