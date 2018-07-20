# Reshape reformatted TAIR10 GOSLIM annotation file for use with R package topGO

library(data.table)

GOann <- read.table(file = "/projects/ajt200/TAIR10/TAIR10_GO_Phytozome_biomart_040418_geneID_GOann.txt", header = F)
colnames(GOann) <- c("gene_ID", "GO_ann")
DT <- data.table(GOann, key = "gene_ID")
GOannReshaped <- DT[ , list("GO_ann" = list(GO_ann)), by = "gene_ID"]
GOannReshapedDF <- data.frame(cbind(as.character(GOannReshaped$gene_ID),
                                    vapply(GOannReshaped$GO_ann, paste,
                                           collapse = ",", character(1L))))
write.table(GOannReshapedDF, "/projects/ajt200/TAIR10/TAIR10_GO_Phytozome_biomart_040418_geneID_GOann_reshaped.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

