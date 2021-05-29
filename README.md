# RNPFind

## What
RNPFind is a bioinformatics tool. It takes a gene name as input and gives a
list of RNA binding proteins (RBPs) that bind to the input gene's RNA
transcript and their binding sites on the transcript as an output. In addition
to collecting the data needed to produce this output, RNPFind provides options
for analysing the data collected (for example, suggesting associations among
the RBPs).

## When
This tool was produced at iYLab at CMUQ to help both my molecular biology
research project as well as that of others in the lab. Originally a CLI tool
that needed to be installed on the user's computer, it was made into a webapp
later for ease of use and access.

## Where
The webapp is avaiable for use at https://rnp-find.herokuapp.com

## How
RNPFind collects its data from three main databases: ATTRACT, RBPDB, and
POSTAR. Each database stores binding sites of RBPs deduced through various
methods, such as various experimental methods as well as computational methods.
RBPDB and ATTRACT provide binding site patterns for RBPs, so we scan the input
gene's sequence to deduce the binding sites. POSTAR stores experimentally
deduced binding sites for the transcriptome for all RBPs, so we simply query.

The various analysis methods are detailed elsewhere.



