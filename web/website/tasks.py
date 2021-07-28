"""
Module for defining celery tasks
"""

# pylint: disable=no-member
# Create your tasks here

import os
import sys

from celery import shared_task
from hgfind import hgfind

from rnpfind import GENOME_VERSION, rnpfind

from .models import AnalysisStatus, BindingSummaryInfo, Gene


class DbFileObject:
    """
    A file-like object to redirect output from gene analysis tasks.
    The object writes outputs as records on a database model.

    Inspiration for this idea:
    https://stackoverflow.com/a/27253257/8551394
    """

    def __init__(self, _key: str):
        self.key = _key
        assert len(AnalysisStatus.objects.filter(request_id=self.key)) > 0
        # AnalysisStatus(request_id=self.key, output="starting analysis").save()
        assert len(AnalysisStatus.objects.filter(request_id=self.key)) > 0
        analysis_status_obj = AnalysisStatus.objects.get(request_id=self.key)
        analysis_status_obj.output = "starting analysis"
        analysis_status_obj.save()

    def write(self, data):
        """
        Writes data to a database
        """
        if len(data) < 4:
            print("too short not saving to db", file=sys.__stderr__)
            return

        print(
            f"Writing to db: {data}",
            file=sys.__stderr__,
        )

        if "Collected" in data:
            assert len([c for c in data if c == "("]) == 1
            assert len([c for c in data if c == ")"]) == 1
            assert data.find("(") < data.find(")")
            nums = [int(s) for s in data.split() if s.isdigit()]
            source_text = data[data.find("(") + 1 : data.find(")")]

            num_rbps = nums[0]
            num_sites = nums[1]

            if "total" in source_text:
                data_source = BindingSummaryInfo.TOTAL
            else:
                data_source = source_text.split()[1]

            BindingSummaryInfo(
                data_source_type=data_source,
                gene=self.key,
                number_of_sites=num_sites,
                number_of_rbps=num_rbps,
            ).save()

        assert len(AnalysisStatus.objects.filter(request_id=self.key)) > 0

        analysis_status_obj = AnalysisStatus.objects.get(request_id=self.key)

        analysis_status_obj.output = data
        analysis_status_obj.save()

        assert AnalysisStatus.objects.get(request_id=self.key).output == data

    def close(self):
        """
        Needed to simulate files
        """


@shared_task(rate_limit="4/m")
def analyze_gene(gene):
    """
    Analyzes a given gene; redirects output to write to a database
    """
    prev_stderr = sys.stderr
    sys.stderr = DbFileObject(gene)
    rnpfind(
        gene,
        methods=["bed"],
        is_trackhub_only=True,
        out_dir=f"/app/staticfiles/{gene.lower()}",
    )

    sys.stderr = prev_stderr
    return gene


@shared_task
def analyze_gene_done(gene):
    """
    Writes the results of gene analysis (to mark its completion)
    """
    # Save the result in db
    rna_info = hgfind(gene)
    assert rna_info["official_name"] == gene

    host_url = os.environ.get("HOST_URL")
    hub_url = (
        f"{host_url}/static/{gene.lower()}/trackhub"
        f"/rbps-on-{gene.lower()}.hub.txt"
    )

    ucsc_url = (
        f"https://genome.ucsc.edu/cgi-bin/hgTracks?db={GENOME_VERSION}&hubUrl="
        f"{hub_url}&position=chr{rna_info['chr_n']}:{rna_info['start_coord']}"
        f"-{rna_info['end_coord']}"
    )

    Gene(
        name=gene,
        chromosome_number=rna_info["chr_n"],
        start_coord=rna_info["start_coord"],
        end_coord=rna_info["end_coord"],
        ucsc_url=ucsc_url,
    ).save()

    print("Gene record saved!", file=sys.stderr)
    print(Gene.objects.get(name=gene).total_sites())
    print(Gene.objects.get(name=gene).num_unique_rbps())
