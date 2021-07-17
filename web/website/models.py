# pylint: disable=no-member
"""
Defines the database models used by RNPFind.
Classes defined are:
    Gene - stores collected data for a gene like name, coordinate, sites, etc.
    BindingSummaryInfo - stores summary data for a gene (no. of binding sites,
                                                         etc.)
    AnalysisStatus - stores real time analysis status data (which step is being
                                                            done, etc.)
"""
from django.db import models

# Create your models here.


class Gene(models.Model):
    """
    stores collected data for a gene like name, coordinate, sites, etc.
    """

    name = models.CharField(max_length=200)
    chromosome_number = models.CharField(max_length=200)
    start_coord = models.IntegerField()
    end_coord = models.IntegerField()
    ucsc_url = models.URLField()

    def size(self):
        """
        size of gene in bases
        """
        return self.end_coord - self.start_coord


class AnalysisStatus(models.Model):
    """
    model for helping store real-time data on computation status of a gene
    analysis request
    """

    request_id = models.CharField(max_length=200)
    output = models.CharField(max_length=200)
