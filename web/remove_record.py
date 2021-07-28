"""
Removes a gene entry from the database.

This allows testing the analysis of that gene both in production and
development
"""
import sys

from website.models import Gene

TO_REMOVE = "SPRR4"

if Gene.objects.filter(name=TO_REMOVE).exists():
    gene_obj = Gene.objects.get(name=TO_REMOVE)
    print(
        f"{gene_obj.name} deleted for debugging...",
        file=sys.stderr,
    )
    gene_obj.delete()
