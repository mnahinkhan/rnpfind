# pylint: disable=no-member
# Create your tasks here

import sys

from celery import shared_task
from hgfind import hgfind

from rnpfind import GENOME_VERSION, rnpfind

from .models import AnalysisStatus, Gene


# Credit for this idea:
# https://stackoverflow.com/a/27253257/8551394
class DbFileObject:
    def __init__(self, _key: str):
        self.key = _key
        assert len(AnalysisStatus.objects.filter(request_id=self.key)) > 0
        # AnalysisStatus(request_id=self.key, output="starting analysis").save()
        assert len(AnalysisStatus.objects.filter(request_id=self.key)) > 0
        analysis_status_obj = AnalysisStatus.objects.get(request_id=self.key)
        analysis_status_obj.output = "starting analysis"
        analysis_status_obj.save()

    def write(self, data):
        if len(data) < 4:
            print("too short not saving to db", file=sys.__stderr__)
            return

        print(
            f"Writing to db: {data}", file=sys.__stderr__,
        )

        assert len(AnalysisStatus.objects.filter(request_id=self.key)) > 0

        analysis_status_obj = AnalysisStatus.objects.get(request_id=self.key)

        analysis_status_obj.output = data
        analysis_status_obj.save()

        assert AnalysisStatus.objects.get(request_id=self.key).output == data

    def close(self):
        pass


@shared_task
def analyze_gene(gene):
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
    # Save the result in db
    rna_info = hgfind(gene)
    assert rna_info["official_name"] == gene

    hub_url = f"https://rnpfind.com/static/{gene.lower()}/trackhub/rbps-on-{gene.lower()}.hub.txt"

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
