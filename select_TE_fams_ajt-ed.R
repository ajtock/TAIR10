#################################################
# select and annotate analysis classes/families #
#################################################
#########################################################################################################
# The following element families were selected from Buisine transposon annotation                       #
# RNA: LTR/Gypsy, LTR/Copia, LINE/L1, [SINE,RathE1_cons,RathE2_cons,RathE3_cons]                        #
# DNA: Helitron, MuDR, En-Spm, hAT, PIF/Harbinger, [Pogo,Tc1,Mariner]                                   #
# Unassigned were dropped (n=129), LINE? dropped	, DNA dropped as caused spiking in AT-TSS and TTS plots	#
#########################################################################################################

outDirDNA <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
outDirRNA <- "/projects/ajt200/TAIR10/TE_classes/RNA/"
outDirAll <- "/projects/ajt200/TAIR10/TE_classes/"

all.tes <- read.table(file="/projects/ajt200/TAIR10/TAIR10_Buisine_TEs_strand_tab_ann.txt",header=T)
print(dim(all.tes))
all.tes <- all.tes[-which(all.tes[,7]=="Unassigned"),]
print(dim(all.tes))
all.tes <- all.tes[-which(all.tes[,7]=="LINE?"),]
print(dim(all.tes))
all.tes <- all.tes[-which(all.tes[,7]=="DNA"),]
print(dim(all.tes))
#[1] 29150     7
super.fams <- as.character(unique(all.tes[,7]))

################
# DNA elements #
################

heli <- all.tes[which(all.tes[,7]=="RC/Helitron"),]
heli <- cbind(heli,rep("heli",length(heli[,1])))
colnames(heli) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(heli))
#[1] 12945     8
sum(heli[,3]-heli[,2])
#[1] 7535358
write.table(heli, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_heli.txt"), sep = "\t", row.names = F, quote = F)

pogo <- all.tes[which(all.tes[,7]=="DNA/Pogo"),]
tc1 <- all.tes[which(all.tes[,7]=="DNA/Tc1"),]
mari <- all.tes[which(all.tes[,7]=="DNA/Mariner"),]
ptmari <- rbind(pogo,tc1,mari)
ptmari <- cbind(ptmari,rep("ptmari",length(ptmari[,1])))
colnames(ptmari) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(ptmari))
#[1] 590   8
sum(ptmari[,3]-ptmari[,2])
#[1] 188850
write.table(ptmari, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_ptmari.txt"), sep = "\t", row.names = F, quote = F)

mudr <- all.tes[which(all.tes[,7]=="DNA/MuDR"),]
mudr <- cbind(mudr,rep("mudr",length(mudr[,1])))
colnames(mudr) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(mudr))
#[1] 5410    8
sum(mudr[,3]-mudr[,2])
#[1] 4137203
write.table(mudr, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_mudr.txt"), sep = "\t", row.names = F, quote = F)

enspm <- all.tes[which(all.tes[,7]=="DNA/En-Spm"),]
enspm <- cbind(enspm,rep("enspm",length(enspm[,1])))
colnames(enspm) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(enspm))
#[1] 941   8
sum(enspm[,3]-enspm[,2])
#[1] 1176532
write.table(enspm, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_enspm.txt"), sep = "\t", row.names = F, quote = F)

hat <- all.tes[which(all.tes[,7]=="DNA/HAT"),]
hat <- cbind(hat,rep("hat",length(hat[,1])))
colnames(hat) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(hat))
#[1] 1035    8
sum(hat[,3]-hat[,2])
#[1] 460061
write.table(hat, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_hat.txt"), sep = "\t", row.names = F, quote = F)

harbinger <- all.tes[which(all.tes[,7]=="DNA/Harbinger"),]
harbinger <- cbind(harbinger,rep("harbinger",length(harbinger[,1])))
colnames(harbinger) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(harbinger))
#[1] 379   8
sum(harbinger[,3]-harbinger[,2])
#[1] 248007
write.table(harbinger, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_harbinger.txt"), sep = "\t", row.names = F, quote = F)

dna <- rbind(heli,ptmari,mudr,enspm,hat,harbinger)
dna <- cbind(dna,rep("dna",length(dna[,1])))
colnames(dna) <- c("chr","start","end","strand","name","family","superfamily","class", "dna.rna")
print(dim(dna))
#[1] 21300     9
sum(dna[,3]-dna[,2])
#[1] 13746011
write.table(dna, file = paste0(outDirDNA, "TAIR10_Buisine_TEs_strand_tab_ann_dna.txt"), sep = "\t", row.names = F, quote = F)


################
# RNA elements #
################

gypsy <- all.tes[which(all.tes[,7]=="LTR/Gypsy"),]
gypsy <- cbind(gypsy,rep("gypsy",length(gypsy[,1])))
colnames(gypsy) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(gypsy))
#[1] 4181    8
sum(gypsy[,3]-gypsy[,2])
#[1] 6970284
write.table(gypsy, file = paste0(outDirRNA, "TAIR10_Buisine_TEs_strand_tab_ann_gypsy.txt"), sep = "\t", row.names = F, quote = F)

copia <- all.tes[which(all.tes[,7]=="LTR/Copia"),]
copia <- cbind(copia,rep("copia",length(copia[,1])))
colnames(copia) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(copia))
#[1] 1781    8
sum(copia[,3]-copia[,2])
#[1] 1917157
write.table(copia, file = paste0(outDirRNA, "TAIR10_Buisine_TEs_strand_tab_ann_copia.txt"), sep = "\t", row.names = F, quote = F)

linel1 <- all.tes[which(all.tes[,7]=="LINE/L1"),]
linel1 <- cbind(linel1,rep("linel1",length(linel1[,1])))
colnames(linel1) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(linel1))
#[1] 1366    8
sum(linel1[,3]-linel1[,2])
#[1] 1247521
write.table(linel1, file = paste0(outDirRNA, "TAIR10_Buisine_TEs_strand_tab_ann_linel1.txt"), sep = "\t", row.names = F, quote = F)

sine <- all.tes[which(all.tes[,7]=="SINE"),]
sine <- cbind(sine,rep("sine",length(sine[,1])))
colnames(sine) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(sine))
# 131
rathe1 <- all.tes[which(all.tes[,7]=="RathE1_cons"),]
rathe1 <- cbind(rathe1,rep("sine",length(rathe1[,1])))
colnames(rathe1) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(rathe1))
# 213
rathe2 <- all.tes[which(all.tes[,7]=="RathE2_cons"),]
rathe2 <- cbind(rathe2,rep("sine",length(rathe2[,1])))
colnames(rathe2) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(rathe2))
# 74
rathe3 <- all.tes[which(all.tes[,7]=="RathE3_cons"),]
rathe3 <- cbind(rathe3,rep("sine",length(rathe3[,1])))
colnames(rathe3) <- c("chr","start","end","strand","name","family","superfamily","class")
print(dim(rathe3))
# 104
sine <- rbind(sine,rathe1,rathe2,rathe3)
print(dim(sine))
# 552
sum(sine[,3]-sine[,2])
# 86402
write.table(sine, file = paste0(outDirRNA, "TAIR10_Buisine_TEs_strand_tab_ann_sine.txt"), sep = "\t", row.names = F, quote = F)

rna <- rbind(gypsy,copia,linel1,sine)
rna <- cbind(rna,rep("rna",length(rna[,1])))
colnames(rna) <- c("chr","start","end","strand","name","family","superfamily","class", "dna.rna")
print(dim(rna))
# 7850
sum(rna[,3]-rna[,2])
# 10221364
write.table(rna, file = paste0(outDirRNA, "TAIR10_Buisine_TEs_strand_tab_ann_rna.txt"), sep = "\t", row.names = F, quote = F)

all.tes <- rbind(dna,rna)
print(dim(all.tes))
# 29150
sum(all.tes[,3]-all.tes[,2])
# 23967375
print(unique(all.tes[,8]))
# heli ptmari mudr enspm hat harbinger gypsy copia linel1 sine
write.table(all.tes, file = paste0(outDir, "TAIR10_Buisine_TEs_strand_tab_ann_all.txt"), sep = "\t", row.names = F, quote = F)

