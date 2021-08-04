#!/usr/bin/env bash
cat common_genes.txt | xargs -I _ curl -4 https://rnpfind.com/analysis-request/_
