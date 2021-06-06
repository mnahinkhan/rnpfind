"""
Django files used to define URLs supported and corresponding view-functions
used to serve requests made to URLs
"""
from django.urls import path

from . import views

urlpatterns = [
    path("", views.index, name="index"),
    path("gene", views.gene_page_redirector, name="gene-page-redirector"),
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
]
