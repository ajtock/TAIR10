#################################################
# select and annotate analysis classes/families #
#################################################
#########################################################################################################
# The following element families were selected from Buisine transposon annotation                       #
# RNA: LTR/Gypsy, LTR/Copia, LINE/L1, [SINE,RathE1_cons,RathE2_cons,RathE3_cons]                        #
# DNA: Helitron, MuDR, En-Spm, hAT, PIF/Harbinger, [Pogo,Tc1,Mariner]                                   #
# Unassigned were dropped (n=129), LINE? dropped	, DNA dropped as caused spiking in AT-TSS and TTS plots	#
#########################################################################################################

all.tes <- read.table(file="TAIR10_Buisine_TEs_strand.txt",header=T)
print(dim(all.tes))
all.tes <- all.tes[-which(all.tes[,7]=="Unassigned"),]
print(dim(all.tes))
all.tes <- all.tes[-which(all.tes[,7]=="LINE?"),]
print(dim(all.tes))
all.tes <- all.tes[-which(all.tes[,7]=="DNA"),]
print(dim(all.tes))
# 29150
super.fams <- as.character(unique(all.tes[,7]))

################
# DNA elements #
################

heli <- all.tes[which(all.tes[,7]=="RC/Helitron"),]
heli <- cbind(heli,rep("heli",length(heli[,1])))
colnames(heli) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 12945
sum(heli[,3]-heli[,2])
# 7535358

pogo <- all.tes[which(all.tes[,7]=="DNA/Pogo"),]
tc1 <- all.tes[which(all.tes[,7]=="DNA/Tc1"),]
mari <- all.tes[which(all.tes[,7]=="DNA/Mariner"),]
ptmari <- rbind(pogo,tc1,mari)
ptmari <- cbind(ptmari,rep("ptmari",length(ptmari[,1])))
colnames(ptmari) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 590
sum(ptmari[,3]-ptmari[,2])
# 188850

mudr <- all.tes[which(all.tes[,7]=="DNA/MuDR"),]
mudr <- cbind(mudr,rep("mudr",length(mudr[,1])))
colnames(mudr) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 5410
sum(mudr[,3]-mudr[,2])
# 4137203

enspm <- all.tes[which(all.tes[,7]=="DNA/En-Spm"),]
enspm <- cbind(enspm,rep("enspm",length(enspm[,1])))
colnames(enspm) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 941 
sum(enspm[,3]-enspm[,2])
# 1176532

hat <- all.tes[which(all.tes[,7]=="DNA/HAT"),]
hat <- cbind(hat,rep("hat",length(hat[,1])))
colnames(hat) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 1035
sum(hat[,3]-hat[,2])
# 460061

harbinger <- all.tes[which(all.tes[,7]=="DNA/Harbinger"),]
harbinger <- cbind(harbinger,rep("harbinger",length(harbinger[,1])))
colnames(harbinger) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 379
sum(harbinger[,3]-harbinger[,2])
# 248007

dna <- rbind(heli,ptmari,mudr,enspm,hat,harbinger)
dna <- cbind(dna,rep("dna",length(dna[,1])))
colnames(dna) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class","dna.rna")
# 21300
sum(dna[,3]-dna[,2])
# 13746011

################
# RNA elements #
################

gypsy <- all.tes[which(all.tes[,7]=="LTR/Gypsy"),]
gypsy <- cbind(gypsy,rep("gypsy",length(gypsy[,1])))
colnames(gypsy) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 4181
sum(gypsy[,3]-gypsy[,2])
# 6970284

copia <- all.tes[which(all.tes[,7]=="LTR/Copia"),]
copia <- cbind(copia,rep("copia",length(copia[,1])))
colnames(copia) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 1781
sum(copia[,3]-copia[,2])
# 1917157

linel1 <- all.tes[which(all.tes[,7]=="LINE/L1"),]
linel1 <- cbind(linel1,rep("linel1",length(linel1[,1])))
colnames(linel1) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 1366
sum(linel1[,3]-linel1[,2])
# 1247521

sine <- all.tes[which(all.tes[,7]=="SINE"),]
sine <- cbind(sine,rep("sine",length(sine[,1])))
colnames(sine) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 131
rathe1 <- all.tes[which(all.tes[,7]=="RathE1_cons"),]
rathe1 <- cbind(rathe1,rep("sine",length(rathe1[,1])))
colnames(rathe1) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 213
rathe2 <- all.tes[which(all.tes[,7]=="RathE2_cons"),]
rathe2 <- cbind(rathe2,rep("sine",length(rathe2[,1])))
colnames(rathe2) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 74
rathe3 <- all.tes[which(all.tes[,7]=="RathE3_cons"),]
rathe3 <- cbind(rathe3,rep("sine",length(rathe3[,1])))
colnames(rathe3) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class")
# 104
sine <- rbind(sine,rathe1,rathe2,rathe3)
# 552
sum(sine[,3]-sine[,2])
# 86402

rna <- rbind(gypsy,copia,linel1,sine)
rna <- cbind(rna,rep("rna",length(rna[,1])))
colnames(rna) <- c("Chr","start","end","strand","Transposon_Name","Transposon_Family","Transposon_Super_Family","class","dna.rna")
# 7850
sum(rna[,3]-rna[,2])
# 10221364

all.tes <- rbind(dna,rna)
# 29150
sum(all.tes[,3]-all.tes[,2])
# 23967375
print(unique(all.tes[,8]))
# heli ptmari mudr enspm hat harbinger gypsy copia linel1 sine


