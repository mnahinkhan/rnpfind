"""
Removes malat1 entry from the database.

This allows testing the analysis of malat1 both in production and development
"""
import sys

from website.models import Gene

if Gene.objects.filter(name="MALAT1").exists():
    malat1_obj = Gene.objects.get(name="MALAT1")
    print(
        f"{malat1_obj.name} deleted for debugging...",
        file=sys.stderr,
    )
    malat1_obj.delete()
