library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
#library(BSgenome.Athaliana.TAIR.TAIR9)

## Chromosome sequence definitions
#chr1 <- Athaliana$Chr1
#chr2 <- Athaliana$Chr2
#chr3 <- Athaliana$Chr3
#chr4 <- Athaliana$Chr4
#chr5 <- Athaliana$Chr5

GTFfile <- "/projects/ajt200/TAIR10/Arabidopsis_thaliana.TAIR10.38.gtf"
FASTAfile <- "/projects/ajt200/TAIR10/TAIR10_chr_all.fa"
outDir <- "/projects/ajt200/TAIR10/union_gene_model_210318/"

# Load the annotation and reduce it
GTF <- import.gff3(GTFfile, format = "gtf", feature.type = c("exon", "five_prime_utr", "three_prime_utr"))
seqlevels(GTF) <- sub("Mt", "mitochondria", seqlevels(GTF))
seqlevels(GTF) <- sub("Pt", "chloroplast", seqlevels(GTF))
GTF <- GTF[GTF$gene_biotype == "protein_coding" | GTF$gene_biotype == "nontranslating_CDS"]
GTFreduceSplitGeneID <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(GTFreduceSplitGeneID, use.names = T)

elementMetadata(reducedGTF)$gene_id <- rep(names(GTFreduceSplitGeneID), elementNROWS(GTFreduceSplitGeneID))

# Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

elementMetadata(reducedGTF)$seq <- getSeq(FASTA, reducedGTF)
splitReducedGTF <- split(reducedGTF, elementMetadata(reducedGTF)$gene_id)
geneFasta <- NULL
geneFastaDNAStringSet <- NULL
for(i in 1:length(splitReducedGTF)) {
  geneSeq <- NULL
  geneSeqDNAStringSet <- NULL
  for(j in 1:length(splitReducedGTF[[i]])) {
    print(j)
    if(j == 1 & j+1 > length(splitReducedGTF)) {
      geneSeq <- splitReducedGTF[[i]]$seq[j]
      geneSeqDNAStringSet <- splitReducedGTF[[i]]$seq[j]
    }
    if(j == 1 & j+1 <= length(splitReducedGTF)) {
      geneSeq <- paste0(splitReducedGTF[[i]]$seq[j], splitReducedGTF[[i]]$seq[j+1])
      geneSeqDNAStringSet <- xscat(splitReducedGTF[[i]]$seq[j], splitReducedGTF[[i]]$seq[j+1])
    }
    if(j > 2) {
      geneSeq <- paste0(geneSeq, splitReducedGTF[[i]]$seq[j])
      geneSeqDNAStringSet <- xscat(geneSeqDNAStringSet, splitReducedGTF[[i]]$seq[j])
    }
  }
  geneFasta <- rbind(geneFasta, geneSeq)
  #geneFastaDNAStringSet <- rbind(geneFastaDNAStringSet, geneSeqDNAStringSet)
}

  geneSeq <- NULL
  for(j in 1:length(splitReducedGTF[[i]])) {
    print(j)
    if(j == 1) { 
      geneSeq <- xscat(splitReducedGTF[[i]]$seq[j], splitReducedGTF[[i]]$seq[j+1])
    }
    if(j > 2) {
      geneSeq <- xscat(geneSeq, splitReducedGTF[[i]]$seq[j])
    }
  }


  geneFasta <- rbind(geneFasta, paste0(">", splitReducedGTF[[i]]$gene_id))
  geneFasta <- rbind(geneFasta, paste0(splitReducedGTF[[i]]$seq[]))
}


getSeq(FASTA, reducedGTF)
elementMetadata(reducedGTF)$seq <- getSeq(FASTA, reducedGTF)

length(reducedGTF[grepl("AT[1-9]G[0-9+]", reducedGTF$gene_id)]$gene_id)

reducedGTF[!grepl("AT", reducedGTF$gene_id)]$gene_id

## Seb
# I assume you are mapping against the genome rather the transcriptome, since for the later the length would be trivial.
#
# Assuming the first, I think not only the coding sections should be included but also the UTR, since reads can map against them which is what we ultimately care about.
#
#In general, I found gene annotation files (e.g. gff or gtf) can be inconsistent in terms of naming, so it's good practice to inspect and double check. Below is some R code to import the annotation and calculate isoform lengths:
#
library(rtracklayer)
# Reading into a GRanges object
anno <- import.gff3("annotation.gff")
# Filtering exons and UTRs
exons <- anno[anno@elementMetadata$type %in% c("exon","five_prime_UTR","three_prime_UTR"),]

# Depending on the annotation at hand, the most sensible is probably best to count the length of each isoform which are often contained in the "Parent" column of the annotation file:
# splitting up isoforms as preparation for the next step
tmp <- split(exons,as.character(exons$Parent))
# for each isoform, calculate the sum of all reduced exons
Gene_length <- sum(width(reduce(tmp)))

# Note, reduce merges overlapping intervals together, since UTRs can "contain" bits of exons which would be otherwise double counted.

# This code can of course be adapted mainly by changing the "Parent", "exon" etc.


## Devon Ryan
# If you're filtering for exons then you needn't include the UTRs. You're not hurting anything since you reduce() them out anyway, but you could have freed up some memory. You can also do the filtering directly in import.gff3 (see the first link in my answer).

