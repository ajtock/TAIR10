#!/bin/bash

# Extract column 1 values containing text matching the following pattern
awk '$1 ~ "^AT[1-5|M|C]G"' TAIR10_ATH_GO_GOSLIM_030418.txt | \
awk 'BEGIN {OFS="\t"}; {print $1}' | \
# Remove ecotype info
sed 's/-\w\+$//g' > TAIR10_ATH_GO_GOSLIM_030418_column1.txt

# Extract GO annotations as separate column
awk '$1 ~ "^AT[1-5|M|C]G"' TAIR10_ATH_GO_GOSLIM_030418.txt | \
awk '{
  for(i = 1; i <= NF; i++) {
    if ($i ~ /^GO:/) {
      printf "%s ", $i, $(i + 1)
    }
  }
  print ""
}' | \
# Where multiple GO annotations are present, separate with comma
sed -e 's/ GO:/,GO:/g' | \
sed -e 's/[|]/,/g' > TAIR10_ATH_GO_GOSLIM_030418_columnGO.txt

# Create file with gene ID (column 1) and GO annotations (column 2)
paste -d "\t" TAIR10_ATH_GO_GOSLIM_030418_column1.txt TAIR10_ATH_GO_GOSLIM_030418_columnGO.txt > TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann_tmp.txt
# Sort file and retain only unique lines
sort TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann_tmp.txt | uniq > TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann.txt


# Remove intermediate files
#rm TAIR10_ATH_GO_GOSLIM_030418_column1.txt TAIR10_ATH_GO_GOSLIM_030418_columnGO.txt TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann_tmp.txt
