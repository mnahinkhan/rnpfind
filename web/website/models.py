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

    def total_sites(self):
        """
        computes total number of binding sites stored for the gene
        """
        return BindingSummaryInfo.objects.get(
            gene=self.name, data_source_type=BindingSummaryInfo.TOTAL
        ).number_of_sites

        # return sum(summary.number_of_sites
        #                             for summary in self.binding_summaries.all())

    def num_unique_rbps(self):
        """
        returns the number of unique RBPs that bound to this gene, amongst
        the various data sources that were used
        """
        return BindingSummaryInfo.objects.get(
            gene=self.name, data_source_type=BindingSummaryInfo.TOTAL
        ).number_of_rbps

    def size(self):
        """
        size of gene in bases
        """
        return self.end_coord - self.start_coord + 1


class BindingSummaryInfo(models.Model):
    """
    class for defining possible data source types.
    """

    TOTAL = "TOTAL"
    # RBPDB = 'RB'
    # ATTRACT = 'AT'
    # POSTAR = 'PO'
    # TOTAL = 'TO'
    # DATABASE_CHOICES = [
    #     (RBPDB, 'RBPDB'),
    #     (ATTRACT, 'ATTRACT'),
    #     (POSTAR, 'POSTAR'),
    #     (TOTAL, 'TOTAL'),
    # ]
    data_source_type = models.CharField(
        max_length=200,
        # choices=DATABASE_CHOICES,
    )
    gene = models.CharField(max_length=80)
    number_of_sites = models.IntegerField()
    number_of_rbps = models.IntegerField()


class AnalysisStatus(models.Model):
    """
    model for helping store real-time data on computation status of a gene
    analysis request
    """

    request_id = models.CharField(max_length=200)
    output = models.CharField(max_length=200)
