from website.models import AnalysisStatus, Gene

bad_records = [
    x
    for x in AnalysisStatus.objects.all()
    if not Gene.objects.filter(name=x.request_id).exists()
]

for bad_record in bad_records:
    bad_record.delete()
