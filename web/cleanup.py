"""
Removes 'bad analysis status records' from the database.

Such bad records are those which exist despite their counterpart
records not existing among the 'Gene records', which store the results of
the analysis after it is complete.

These records are bad as their exitence from the outset might fool the
system into thinking that the gene is currently being analyzed, when in
fact it is not.
"""
import sys

from website.models import AnalysisStatus, Gene

bad_records = [
    x
    for x in AnalysisStatus.objects.all()
    if not Gene.objects.filter(name=x.request_id).exists()
]

for bad_record in bad_records:
    print(
        f"{bad_record.request_id} deleted as it is a bad record",
        file=sys.stderr,
    )
    bad_record.delete()
