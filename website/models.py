from django.db import models

# Create your models here.


from django.utils.translation import gettext_lazy as _
from django.db import models


class Gene(models.Model):
    name = models.CharField(max_length=20)
    chromosome_number = models.CharField(max_length=3)
    start_coord = models.IntegerField()
    end_coord = models.IntegerField()
    ucsc_url = models.URLField()


class BindingSummaryInfo(models.Model):
    RBPDB = 'RB'
    ATTRACT = 'AT'
    POSTAR = 'PO'
    DATABASE_CHOICES = [
        (RBPDB, 'RBPDB'),
        (ATTRACT, 'ATTRACT'),
        (POSTAR, 'POSTAR'),
    ]
    data_source_type = models.CharField(
        max_length=2,
        choices=DATABASE_CHOICES,
    )
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE, related_name="binding_summaries")
    number_of_sites = models.IntegerField()
    number_of_rbps = models.IntegerField()

class AnalysisStatus(models.Model):
    request_id = models.CharField(max_length=20)
    status = models.CharField(max_length=200)