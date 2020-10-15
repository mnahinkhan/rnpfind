from django.shortcuts import render
from .models import Gene, BindingSummaryInfo, AnalysisStatus
from time import sleep

# Create your views here.
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse
from django.shortcuts import render
from django.urls import reverse

# RNPFind source script files
from website.scripts.gene_coordinates import gene_to_coord
from website.scripts.load_data import load_data
from website.scripts.analysis_functions import analysis_method_functions


def index(request, bad_gene=""):
    # bad_gene = request.GET['gene_name'] if 'gene_name' in request.GET else ""
    return render(request, "website/index.html", {
        'bad_gene': bad_gene,
        'analysed_genes': Gene.objects.all()
    })

def gene_page_redirector(request):
    gene_name = request.GET['gene_name']
    return gene_page(request, gene_name)

def gene_page(request, gene_name):
    # Gene.objects.all().delete();
    # gene_name = request.GET['gene_name']
    success, _,_,_, official_name = gene_to_coord(gene_name)
    if not success:
        return index(request, bad_gene=gene_name)

    is_in_database = len(Gene.objects.filter(name=official_name)) != 0

    return render(request, "website/gene_page.html", {
        'cached': is_in_database,
        'gene_official_name': official_name,
        'gene_typed_name': gene_name,
        'gene_cached_data': Gene.objects.get(name=official_name) if is_in_database else 0
    })

def process(gene_name):
    return gene_name[::-1]


def analysis_status(request, request_id):
    status = AnalysisStatus.objects.get(request_id=request_id).status if len(AnalysisStatus.objects.filter(request_id=request_id)) > 0 else ""
    return JsonResponse({"status": status if status else "0/13. Sending request to server for analysis"})

def out_fn_gen(request_id):
    def out(msg):
        if len(AnalysisStatus.objects.filter(request_id=request_id)) > 0:
            analysis_status_obj = AnalysisStatus.objects.get(request_id=request_id)
            analysis_status_obj.status = msg
            analysis_status_obj.save()
        else:
            AnalysisStatus(request_id=request_id, status=msg).save()
    return out

def analysis_request(request, gene):
    success, chr_n, start_coord, end_coord, official_name = gene_to_coord(gene)
    assert success
    assert official_name == gene
    assert len(Gene.objects.filter(name=official_name)) == 0

    if len(AnalysisStatus.objects.filter(request_id=gene)) > 0:
        return JsonResponse({"process_complete": False})
    gene_entry = Gene(name=gene)

    # # get reversed
    # reversed_gene = process(gene)
    # gene_entry.reverse_name = reversed_gene

    total_steps = 13

    request_id = gene
    out = out_fn_gen(request_id)
    # out("Loading binding sites")

    # collect binding sites
    data_load_sources = ['rbpdb', 'attract', 'postar']
    rna_info = [gene, chr_n, start_coord, end_coord]
    big_storage = load_data(data_load_sources, rna_info, out=out, total_steps=total_steps)
    # out("Done loading binding sites!")

    
    out(f"4/{total_steps}. Creating visualizations on UCSC")

    analysis_method = 'ucsc'
    analysis_method_function = analysis_method_functions[analysis_method]
    ucsc_url = analysis_method_function(big_storage, rna_info, out=out, total_steps=total_steps)
    gene_entry.ucsc_url = ucsc_url

    

    out(f"11/{total_steps}. Caching computed values")
    gene_entry.chromosome_number = chr_n
    gene_entry.start_coord = start_coord
    gene_entry.end_coord = end_coord

    gene_entry.save()

    out(f"12/{total_steps}. Generating summary data")
    # Save binding site summary info
    for data_source in big_storage:
        no_unique_rbps, no_unqiue_sites = big_storage[data_source].summary(is_return=True) 
        BindingSummaryInfo(data_source_type=data_source, gene=gene_entry, number_of_sites=no_unqiue_sites, number_of_rbps=no_unique_rbps).save()

    assert len(Gene.objects.filter(name=official_name)) == 1
    out(f"13/{total_steps}. Complete! Refreshing your page now")
    return JsonResponse({"process_complete": True})
