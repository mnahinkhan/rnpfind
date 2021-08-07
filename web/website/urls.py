"""
Django files used to define URLs supported and corresponding view-functions
used to serve requests made to URLs
"""
from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path(
        "gene-search", views.gene_page_redirector, name="gene-page-redirector"
    ),
    path("gene/<str:gene_name>", views.gene_page, name="gene-page"),
    path(
        "analysis-request/<str:gene>",
        views.analysis_request,
        name="analysis-request",
    ),
    path(
        "analysis-status/<str:request_id>",
        views.analysis_status,
        name="analysis-status",
    ),
    path(
        "analysis-result/<str:gene>",
        views.analysis_result,
        name="analysis-result",
    ),
    path("about", views.about, name="about"),
    path("cli/README.md", views.cli_docs, name="cli_doc"),
]
