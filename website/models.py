from django.db import models

# Create your models here.


from django.utils.translation import gettext_lazy as _
from django.db import models


class Gene(models.Model):
    name = models.CharField(max_length=200)
    chromosome_number = models.CharField(max_length=200)
    start_coord = models.IntegerField()
    end_coord = models.IntegerField()
    ucsc_url = models.URLField()

    def total_sites(self):
        return sum(summary.number_of_sites for summary in self.binding_summaries.all())

    def max_unique_rbps(self):
        return max(summary.number_of_rbps for summary in self.binding_summaries.all())
    
    def size(self):
        return self.end_coord  - self.start_coord



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
        max_length=200,
        choices=DATABASE_CHOICES,
    )
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE, related_name="binding_summaries")
    number_of_sites = models.IntegerField()
    number_of_rbps = models.IntegerField()

class AnalysisStatus(models.Model):
    request_id = models.CharField(max_length=200)
    status = models.CharField(max_length=200)