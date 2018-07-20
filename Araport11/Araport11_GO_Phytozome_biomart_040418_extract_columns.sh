#!/bin/bash

# Extract columns 1 and 3
tail -n +2 Araport11_GO_Phytozome_biomart_040418.txt | \
awk 'BEGIN {OFS="\t"}; {print $1, $3}' | \
grep "GO:" > Araport11_GO_Phytozome_biomart_040418_geneID_GOann_tmp.txt

# Sort file and retain only unique lines
sort Araport11_GO_Phytozome_biomart_040418_geneID_GOann_tmp.txt | uniq > Araport11_GO_Phytozome_biomart_040418_geneID_GOann.txt

