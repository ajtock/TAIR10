library(rtracklayer)
ann <- import.gff3("/projects/ajt200/TAIR10/Araport11/Araport11_GFF3_genes_transposons.201606.gff")
ann <- ann[seqnames(ann) != "ChrC" & seqnames(ann) != "ChrM"]
ann_mRNA <- ann[ann$type == "mRNA"]

