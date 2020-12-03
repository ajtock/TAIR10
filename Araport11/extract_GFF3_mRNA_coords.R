#!/applications/R/R-3.5.0/bin/Rscript

# Create table of representative genes from GFF3 file
# and generate random loci of the same number and width distribution
# Write as GFF3 and BED files 

# Usage:
# ./extract_GFF3_mRNA_coords.R Araport11_GFF3_genes_transposons.201606.gff

args <- commandArgs(trailingOnly = T)
inFile <- args[1]

library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(data.table)

# Genomic definitions
fai <- read.table("TAIR10_chr_all.fa.fai", header = F)
chrs <- paste0("Chr", fai[,1])[1:5]
chrLens <- fai[,2][1:5]
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genomeGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = rep(1, length(chrs)),
                                     end = chrLens),
                    strand = "*")

genes <- readGFF(inFile)
genes <- genes[!(genes$seqid %in% c("ChrC", "ChrM")),]
mRNA <- genes[genes$type == "mRNA",]
print(dim(mRNA))
#[1] 52060    25

# Obtain frequency of occurrence of each gene parent ID
n_occur <- data.frame(table(unlist(mRNA$Parent)))

# Obtain mRNA records for which the gene parent ID occurs only once
mRNA_unique <-  as.data.frame(
  mRNA[ unlist(mRNA$Parent)
    %in% n_occur$Var1[n_occur$Freq == 1],
  ]
)
# Obtain mRNA records for which the gene parent ID occurs more than once
mRNA_multi <- as.data.frame(
  mRNA[ unlist(mRNA$Parent)
    %in% n_occur$Var1[n_occur$Freq > 1],
  ]
)

# For each gene parent ID in mRNA_multi, obtain the mRNA record with the
# longest transcript
# If multiple mRNA records have the longest transcript,
# keep the first reported one only
mRNA_multi_list <- mclapply(seq_along(mRNA_multi[,1]), function(h) {
  mRNA_multi_ID_all <- mRNA_multi[ unlist(mRNA_multi$Parent)
                         == unlist(mRNA_multi[h,]$Parent),
                       ]
  mRNA_multi_ID_all[ mRNA_multi_ID_all$end-mRNA_multi_ID_all$start
    == max(mRNA_multi_ID_all$end-mRNA_multi_ID_all$start),
  ][1,]
}, mc.cores = detectCores())

# Collapse mRNA_multi_list into single data.frame and remove duplicates
mRNA_multi_dup <- rbindlist(mRNA_multi_list)
mRNA_multi_rep <- unique(as.data.frame(mRNA_multi_dup))


# Combine into one representative set of mRNA entries, order,
# and output in GFF3 and BED formats
mRNA_rep <- rbind(mRNA_unique, mRNA_multi_rep)
stopifnot(dim(mRNA_rep)[1] == dim(n_occur)[1])
mRNA_rep <- mRNA_rep[ order(mRNA_rep$seqid,
                            mRNA_rep$start,
                            mRNA_rep$end), ]
write.table(mRNA_rep[,1:14],
            file = "Araport11_representative_mRNA.gff3",
            quote = F, sep = "\t", row.names = F, col.names = F)
mRNA_rep_bed <- data.frame(chr = as.character(mRNA_rep[,1]),
                           start = as.integer(mRNA_rep[,4]-1),
                           end = as.integer(mRNA_rep[,5]),
                           name = as.integer(1:length(mRNA_rep[,1])),
                           score = as.numeric(mRNA_rep[,6]),
                           strand = as.character(mRNA_rep[,7]))
write.table(mRNA_rep_bed,
            file = "Araport11_representative_mRNA.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)

mRNA_repGR <- GRanges(seqnames = mRNA_rep$seqid,
                      ranges = IRanges(start = mRNA_rep$start,
                                       end = mRNA_rep$end),
                      strand = mRNA_rep$strand)

# Define function to select randomly positioned loci of the same
# width distribution as mRNA_repGR
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Define seed so that random selections are reproducible
set.seed(93750174)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as mRNA_repGR
ranLocGR <- GRanges()
for(i in 1:length(chrs)) {
  mRNA_repChrGR <- mRNA_repGR[seqnames(mRNA_repGR) == chrs[i]]
  genomeChrGR <- genomeGR[seqnames(genomeGR) == chrs[i]]
  # Contract genomeChrGR so that random loci end coordinates do not extend beyond chromosome ends
  end(genomeChrGR) <- end(genomeChrGR)-max(width(mRNA_repChrGR))
  ranLocChrStart <- ranLocStartSelect(coordinates = (start(genomeChrGR):end(genomeChrGR)),
                                      n = length(mRNA_repChrGR))
  ranLocChrGR <- GRanges(seqnames = chrs[i],
                         ranges = IRanges(start = ranLocChrStart,
                                          width = width(mRNA_repChrGR)),
                         strand = strand(mRNA_repChrGR))
  ranLocGR <- append(ranLocGR, ranLocChrGR)
}
ranLoc_bed <- data.frame(chr = as.character(seqnames(ranLocGR)),
                         start = as.integer(start(ranLocGR)-1),
                         end = as.integer(end(ranLocGR)),
                         name = as.integer(1:length(ranLocGR)),
                         score = rep("NA", length(ranLocGR)),
                         strand = as.character(strand(ranLocGR)))
write.table(ranLoc_bed,
            file = "Araport11_representative_mRNA_randomLoci.bed",
            quote = F, sep = "\t", row.names = F, col.names = F)
