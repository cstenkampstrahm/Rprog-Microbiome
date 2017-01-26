# Enrique Doster
# Script for Chloe's data
# 9/18/2016


library(metagenomeSeq)
library("biom")
library(vegan)

## Read biome file, convert to MRexperiment and add sample metadata ############
setwd("C:\\Users\\enrique\\Dropbox\\Projects\\MISC_projects\\Chloe\\")
# import with default parameters, specify a biome file in my set wd
biom_file <- "otu_table_mc2_w_tax_no_pynast_failuresfiltered.biom"
x = read_biom(biom_file)

#check properties of biome file
show(x)
print(x)
header(x)
biom_data(x)
biom_shape(x)
nrow(x) #159354 rows
ncol(x) #301 samples
observation_metadata(x) # this is all the taxonomic metadata
sample_metadata(x)

#convert biom to MRexperiment 
x1=biom2MRexperiment(x)
#View(fData(x1))
#View(pData(x1))
#View(MRcounts(x1))

#Looks like the current MRexperiment does not have the full array of sample metadata, just removed the "#" from the first column
pheno = load_phenoData("QIIME_map_all_corrected.txt", sep = "\t", tran=TRUE)
#this removes the rows that don't match our count table
ord=match(colnames(MRcounts(x1)), rownames(pheno))
pheno=pheno[ord, ]
head(pheno)

pheno=AnnotatedDataFrame(pheno)
features=AnnotatedDataFrame(fData(x1))
# This makes our main object for metagenomeseq with the updated sample metadata
microbiome = newMRexperiment(MRcounts(x1),phenoData=pheno, featureData = features)

#removal of sparse features (rows)
sparseFeatures = which(rowSums(MRcounts(microbiome) > 0)< 3)
length(sparseFeatures) # this seems like way too many sparse features 97,779
microbiometrim = microbiome[-sparseFeatures, ]

# normalization step
p = cumNormStatFast(microbiometrim, pFlag=TRUE, rel=0.1)
microbiometrimnorm = cumNorm(microbiometrim, p = p)


# This is the code for aggregating to each taxanomic level ##########
MicroPhylum = as.character(fData(microbiometrim)[, 2])
Phylummicrobiome= aggTax(MRcounts(microbiometrim, norm=TRUE, log=FALSE),norm=TRUE, log=FALSE, lvl = MicroPhylum, out = "MRexperiment")
ordp=match(colnames(MRcounts(Phylummicrobiome)), rownames(pheno))
phenop=pheno[ordp, ]
Phylummicrobiome = newMRexperiment(MRcounts(Phylummicrobiome),phenoData=phenop, featureData=NULL)

MicroClass = as.character(fData(microbiometrim)[, 3])
Classmicrobiome = aggTax(MRcounts(microbiometrim, norm=TRUE, log=FALSE),norm=TRUE, log=FALSE, lvl = MicroClass, out = "MRexperiment")
ordc=match(colnames(MRcounts(Classmicrobiome)), rownames(pheno))
phenoc=pheno[ordc, ]
Classmicrobiome = newMRexperiment(MRcounts(Classmicrobiome), phenoData=phenoc, featureData=NULL)

MicroOrder = as.character(fData(microbiometrim)[, 4])
Ordermicrobiome = aggTax(MRcounts(microbiometrim, norm=TRUE, log=FALSE),norm=TRUE, log=FALSE, lvl = MicroOrder, out = "MRexperiment")
ordo=match(colnames(MRcounts(Ordermicrobiome)), rownames(pheno))
phenoo=pheno[ordo, ]
Ordermicrobiome = newMRexperiment(MRcounts(Ordermicrobiome), phenoData=phenoo, featureData=NULL)

MicroFamily = as.character(fData(microbiometrim)[, 5])
Familymicrobiome = aggTax(MRcounts(microbiometrim, norm=TRUE, log=FALSE),norm=TRUE, log=FALSE, lvl = MicroFamily, out = "MRexperiment")
ordf=match(colnames(MRcounts(Familymicrobiome)), rownames(pheno))
phenof=pheno[ordf, ]
Familymicrobiome = newMRexperiment(MRcounts(Familymicrobiome), phenoData=phenof, featureData=NULL)

MicroGenus = as.character(fData(microbiometrim)[, 6])
Genusmicrobiome = aggTax(MRcounts(microbiometrim, norm=TRUE, log=FALSE),norm=TRUE, log=FALSE, lvl = MicroGenus, out = "MRexperiment")
ordg=match(colnames(MRcounts(Genusmicrobiome)), rownames(pheno))
phenog=pheno[ordg, ]
Genusmicrobiome = newMRexperiment(MRcounts(Genusmicrobiome), phenoData=phenog, featureData=NULL)

MicroSpecies = as.character(fData(microbiometrim)[, 7])
Speciesmicrobiome = aggTax(MRcounts(microbiometrim, norm=TRUE, log=FALSE),norm=TRUE, log=FALSE, lvl = MicroSpecies, out = "MRexperiment")
ords=match(colnames(MRcounts(Speciesmicrobiome)), rownames(pheno))
phenos=pheno[ords, ]
Speciesmicrobiome = newMRexperiment(MRcounts(Speciesmicrobiome), phenoData=phenos, featureData=NULL)
 
# Make sample groups  ####################
envirosamples = which(is.na(pData(microbiometrimnorm)$Individual_animal))

#Cow_samples
Cow_samples= microbiometrimnorm[, - envirosamples]
Cow_samples_Phylum= Phylummicrobiome[, - envirosamples]
Cow_samples_Class= Classmicrobiome[, - envirosamples]

#Descriptive categories##############
cow_enviro=pData(Cow_samples)$Cow_enviro
DIM=pData(Cow_samples)$DIM
Farm=pData(Cow_samples)$Farm
Pathotype=pData(Cow_samples)$Pathotype
Pattern=pData(Cow_samples)$Pattern
Treatment=pData(Cow_samples)$Treatment
Disease=pData(Cow_samples)$Disease
Cworker=pData(Cow_samples)$Cow_worker

Indanimal=pData(Cow_samples)$Individual_animal
sample_ID=pData(Cow_samples)$sample_ID


# Ordination code #####################
pdf("ordi_micro_cows_Phylum.pdf") # his is to make one PDF per taxa level (you can remove the comments from each part if you want particular pdfs)
#pdf("ordi_micro_cows_Phylum_byEnviro.pdf")
dist_micro_cow_byenviro <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"),"bray")
metaMDS_micro_cow_byenviro <- metaMDS(dist_micro_cow_byenviro, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byenviro, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byenviro, display="sites", pch=20, col=as.numeric(cow_enviro))
text(metaMDS_micro_cow_byenviro, lab=row.names(metaMDS_micro_cow_byenviro))
groupz <- sort(unique(cow_enviro))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byenviro, cow_enviro,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(cow_enviro),pch = 20, title = "cow_enviro",bty = "n")
#dev.off()
anosim(dist_micro_cow_byenviro,cow_enviro)

#pdf("ordi_micro_cows_Phylum_byDIM.pdf")
dist_micro_cow_byDIM <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byDIM <- metaMDS(dist_micro_cow_byDIM, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byDIM, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byDIM, display="sites", pch=20, col=as.numeric(DIM))
text(metaMDS_micro_cow_byDIM, lab=row.names(metaMDS_micro_cow_byDIM))
groupz <- sort(unique(DIM))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byDIM, DIM,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(DIM),pch = 20, title = "DIM",bty = "n")
#dev.off()
anosim(dist_micro_cow_byDIM,DIM)

#pdf("ordi_micro_cows_Phylum_byFarm.pdf")
dist_micro_cow_byFarm <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byFarm <- metaMDS(dist_micro_cow_byFarm, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byFarm, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byFarm, display="sites", pch=20, col=as.numeric(Farm))
text(metaMDS_micro_cow_byFarm, lab=row.names(metaMDS_micro_cow_byFarm))
groupz <- sort(unique(Farm))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byFarm, Farm,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Farm),pch = 20, title = "Farm",bty = "n")
#dev.off()
anosim(dist_micro_cow_byFarm,Farm)

#pdf("ordi_micro_cows_Phylum_byPathotype.pdf")
dist_micro_cow_byPathotype <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byPathotype <- metaMDS(dist_micro_cow_byPathotype, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byPathotype, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byPathotype, display="sites", pch=20, col=as.numeric(Pathotype))
text(metaMDS_micro_cow_byPathotype, lab=row.names(metaMDS_micro_cow_byPathotype))
groupz <- sort(unique(Pathotype))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byPathotype, Pathotype,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Pathotype),pch = 20, title = "Pathotype",bty = "n")
#dev.off()
anosim(dist_micro_cow_byPathotype,Pathotype)

#pdf("ordi_micro_cows_Phylum_byTreatment.pdf")
dist_micro_cow_byTreatment <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byTreatment <- metaMDS(dist_micro_cow_byTreatment, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byTreatment, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byTreatment, display="sites", pch=20, col=as.numeric(Treatment))
text(metaMDS_micro_cow_byTreatment, lab=row.names(metaMDS_micro_cow_byTreatment))
groupz <- sort(unique(Treatment))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byTreatment, Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Treatment),pch = 20, title = "Treatment",bty = "n")
#dev.off()
anosim(dist_micro_cow_byTreatment,Treatment)

#pdf("ordi_micro_cows_Phylum_byDisease.pdf")
dist_micro_cow_byDisease <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byDisease <- metaMDS(dist_micro_cow_byDisease, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byDisease, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byDisease, display="sites", pch=20, col=as.numeric(Disease))
text(metaMDS_micro_cow_byDisease, lab=row.names(metaMDS_micro_cow_byDisease))
groupz <- sort(unique(Disease))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byDisease, Disease,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Disease),pch = 20, title = "Disease",bty = "n")
#dev.off()
anosim(dist_micro_cow_byDisease,Disease)

#pdf("ordi_micro_cows_Phylum_byCworker.pdf")
dist_micro_cow_byCworker <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byCworker <- metaMDS(dist_micro_cow_byCworker, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byCworker, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byCworker, display="sites", pch=20, col=as.numeric(Cworker))
text(metaMDS_micro_cow_byCworker, lab=row.names(metaMDS_micro_cow_byCworker))
groupz <- sort(unique(Cworker))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byCworker, Cworker,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Cworker),pch = 20, title = "Cworker",bty = "n")
#dev.off()
anosim(dist_micro_cow_byCworker,Cworker)

#pdf("ordi_micro_cows_Phylum_byIndanimal.pdf")
dist_micro_cow_byIndanimal <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_byIndanimal <- metaMDS(dist_micro_cow_byIndanimal, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byIndanimal, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byIndanimal, display="sites", pch=20, col=as.numeric(Indanimal))
text(metaMDS_micro_cow_byIndanimal, lab=row.names(metaMDS_micro_cow_byIndanimal))
groupz <- sort(unique(Indanimal))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byIndanimal, Indanimal,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Indanimal),pch = 20, title = "Indanimal",bty = "n")
#dev.off()
anosim(dist_micro_cow_byIndanimal,Indanimal)

#pdf("ordi_micro_cows_Phylum_bysample_ID.pdf")
dist_micro_cow_bysample_ID <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, norm=TRUE)), "hell"))
metaMDS_micro_cow_bysample_ID <- metaMDS(dist_micro_cow_bysample_ID, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_bysample_ID, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_bysample_ID, display="sites", pch=20, col=as.numeric(sample_ID))
text(metaMDS_micro_cow_bysample_ID, lab=row.names(metaMDS_micro_cow_bysample_ID))
groupz <- sort(unique(sample_ID))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_bysample_ID, sample_ID,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(sample_ID),pch = 20, title = "sample_ID",bty = "n")
#dev.off()
anosim(dist_micro_cow_bysample_ID,sample_ID)

dev.off() # this corresponds to the first pdf label that covers all Phylum graphs

## Ordination by Class ##########
pdf("ordi_micro_cows_Class.pdf")
#pdf("ordi_micro_cows_Class_byEnviro.pdf")
dist_micro_cow_byenviro <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"),"bray")
metaMDS_micro_cow_byenviro <- metaMDS(dist_micro_cow_byenviro, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byenviro, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byenviro, display="sites", pch=20, col=as.numeric(cow_enviro))
text(metaMDS_micro_cow_byenviro, lab=row.names(metaMDS_micro_cow_byenviro))
groupz <- sort(unique(cow_enviro))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byenviro, cow_enviro,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(cow_enviro),pch = 20, title = "cow_enviro",bty = "n")
#dev.off()
anosim(dist_micro_cow_byenviro,cow_enviro)

#pdf("ordi_micro_cows_Class_byDIM.pdf")
dist_micro_cow_byDIM <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byDIM <- metaMDS(dist_micro_cow_byDIM, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byDIM, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byDIM, display="sites", pch=20, col=as.numeric(DIM))
text(metaMDS_micro_cow_byDIM, lab=row.names(metaMDS_micro_cow_byDIM))
groupz <- sort(unique(DIM))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byDIM, DIM,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(DIM),pch = 20, title = "DIM",bty = "n")
#dev.off()
anosim(dist_micro_cow_byDIM,DIM)

#pdf("ordi_micro_cows_Class_byFarm.pdf")
dist_micro_cow_byFarm <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byFarm <- metaMDS(dist_micro_cow_byFarm, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byFarm, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byFarm, display="sites", pch=20, col=as.numeric(Farm))
text(metaMDS_micro_cow_byFarm, lab=row.names(metaMDS_micro_cow_byFarm))
groupz <- sort(unique(Farm))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byFarm, Farm,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Farm),pch = 20, title = "Farm",bty = "n")
#dev.off()
anosim(dist_micro_cow_byFarm,Farm)

#pdf("ordi_micro_cows_Class_byPathotype.pdf")
dist_micro_cow_byPathotype <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byPathotype <- metaMDS(dist_micro_cow_byPathotype, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byPathotype, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byPathotype, display="sites", pch=20, col=as.numeric(Pathotype))
text(metaMDS_micro_cow_byPathotype, lab=row.names(metaMDS_micro_cow_byPathotype))
groupz <- sort(unique(Pathotype))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byPathotype, Pathotype,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Pathotype),pch = 20, title = "Pathotype",bty = "n")
#dev.off()
anosim(dist_micro_cow_byPathotype,Pathotype)

#pdf("ordi_micro_cows_Class_byTreatment.pdf")
dist_micro_cow_byTreatment <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byTreatment <- metaMDS(dist_micro_cow_byTreatment, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byTreatment, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byTreatment, display="sites", pch=20, col=as.numeric(Treatment))
text(metaMDS_micro_cow_byTreatment, lab=row.names(metaMDS_micro_cow_byTreatment))
groupz <- sort(unique(Treatment))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byTreatment, Treatment,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Treatment),pch = 20, title = "Treatment",bty = "n")
#dev.off()
anosim(dist_micro_cow_byTreatment,Treatment)

#pdf("ordi_micro_cows_Class_byDisease.pdf")
dist_micro_cow_byDisease <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byDisease <- metaMDS(dist_micro_cow_byDisease, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byDisease, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byDisease, display="sites", pch=20, col=as.numeric(Disease))
text(metaMDS_micro_cow_byDisease, lab=row.names(metaMDS_micro_cow_byDisease))
groupz <- sort(unique(Disease))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byDisease, Disease,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Disease),pch = 20, title = "Disease",bty = "n")
#dev.off()
anosim(dist_micro_cow_byDisease,Disease)

#pdf("ordi_micro_cows_Class_byCworker.pdf")
dist_micro_cow_byCworker <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byCworker <- metaMDS(dist_micro_cow_byCworker, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byCworker, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byCworker, display="sites", pch=20, col=as.numeric(Cworker))
text(metaMDS_micro_cow_byCworker, lab=row.names(metaMDS_micro_cow_byCworker))
groupz <- sort(unique(Cworker))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byCworker, Cworker,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Cworker),pch = 20, title = "Cworker",bty = "n")
#dev.off()
anosim(dist_micro_cow_byCworker,Cworker)

#pdf("ordi_micro_cows_Class_byIndanimal.pdf")
dist_micro_cow_byIndanimal <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_byIndanimal <- metaMDS(dist_micro_cow_byIndanimal, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_byIndanimal, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byIndanimal, display="sites", pch=20, col=as.numeric(Indanimal))
text(metaMDS_micro_cow_byIndanimal, lab=row.names(metaMDS_micro_cow_byIndanimal))
groupz <- sort(unique(Indanimal))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byIndanimal, Indanimal,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(Indanimal),pch = 20, title = "Indanimal",bty = "n")
#dev.off()
anosim(dist_micro_cow_byIndanimal,Indanimal)

#pdf("ordi_micro_cows_Class_bysample_ID.pdf")
dist_micro_cow_bysample_ID <- vegdist(decostand(t(MRcounts(Cow_samples_Class, norm=TRUE)), "hell"))
metaMDS_micro_cow_bysample_ID <- metaMDS(dist_micro_cow_bysample_ID, distance="none",symmetric=TRUE,)
plot(metaMDS_micro_cow_bysample_ID, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_bysample_ID, display="sites", pch=20, col=as.numeric(sample_ID))
text(metaMDS_micro_cow_bysample_ID, lab=row.names(metaMDS_micro_cow_bysample_ID))
groupz <- sort(unique(sample_ID))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_bysample_ID, sample_ID,font=2, cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(sample_ID),pch = 20, title = "sample_ID",bty = "n")
#dev.off()
anosim(dist_micro_cow_bysample_ID,sample_ID)

dev.off()

## Example of using "zero-inflated gaussian model" #####
settings = zigControl(maxit=20, verbose=TRUE)

designMicroTime = model.matrix(~0+ Time + Treatment)
dup_resMicroTime <- duplicateCorrelation(MRcounts(Microbiometrimnorm),block=Animal, design= designMicroTime)
resMicroTime = fitZig(obj= Microbiometrimnorm, mod = designMicroTime, control = settings, useCSSoffset=FALSE, useMixedModel=dup_resMicroTime$consensus)
zigFitMicroTime = resMicroTime$fit
finalModMicroTime = resMicroTime$fit$design
contrastMicroTime = makeContrasts(TimeDay11-TimeArrival, levels=finalModMicroTime)
resMicro2Time = contrasts.fit(zigFitMicroTime, contrastMicroTime)
resMicro2EBTime = eBayes(resMicro2Time)
write.table(topTable(resMicro2EBTime, coef=1, adjust.method="BH",number = 1000), file = "Day11vsArrival_ALL_Micro_block.tsv", sep="\t")

## Comparisons of Richness and Shannons, with glm and wilcoxon #####
Arrival_Micro_S_norm_Order = specnumber(t(MRcounts(MicroArrival_Order, norm=TRUE)))
Arrival_Micro_S_glm_norm_Order = glm(Arrival_Micro_S_norm_Order ~ TreatmentArrival)
summary(Arrival_Micro_S_glm_norm_Order)
Arrival_Micro_H_norm_Order = diversity(t(MRcounts(MicroArrival_Order, norm=TRUE)))
Arrival_Micro_H_glm_norm_Order = glm(Arrival_Micro_H_norm_Order ~ TreatmentArrival)
summary(Arrival_Micro_H_glm_norm_Order)

Arrival_Micro_S = specnumber(t(MRcounts(MicroArrival)))
Arrival_Micro_S_wilcox = wilcox.test(Arrival_Micro_S ~ TreatmentArrival)
summary(Arrival_Micro_S_wilcox)
Arrival_Micro_H = diversity(t(MRcounts(MicroArrival)))
Arrival_Micro_H_wilcox = wilcox.test(Arrival_Micro_H ~ TreatmentArrival)
summary(Arrival_Micro_H_wilcox)


# Some figures for exploring your data #######
require(interactiveDisplay) 
display(Cow_samples_Phylum)

heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(Indanimal))]
heatmapCols = colorRampPalette(brewer.pal(9,"RdBu"),bias=1)(50)

nrow(fData(Cow_samples_Phylum)) 
plotMRheatmap(obj = Cow_samples_Phylum, n =15, cexRow = 0.4, cexCol = 0.4,trace = "none", col = heatmapCols, ColSideColors = heatmapColColors, norm=FALSE)

heatmap.2(MRcounts(Cow_samples_Phylum),
          cellnote = MRcounts(Cow_samples_Phylum),  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=heatmapCols,       # use on color palette defined earlier 
          dendrogram="none",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering

col_distance = dist(t(MRcounts(Cow_samples_Phylum)), method ="manhattan")
col_cluster = hclust(col_distance, method = "ward.D")
heatmap.2(MRcounts(Cow_samples_Phylum),
          main= "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=heatmapCols,       # use on color palette defined earlier
          dendrogram="column",     # only draw a row dendrogram
          Colv=as.dendrogram(col_cluster))            # turn off column clustering

library(igraph)
library(gplot)
library(ggplot2)
Cow_samples_phylum_net2 <- graph.incidence(t(MRcounts(Cow_samples_Phylum,norm=TRUE)),weighted = TRUE) #this is a weighted network
V(Cow_samples_phylum_net2)$name
V(Cow_samples_phylum_net2)$type
E(Cow_samples_phylum_net2)

plot(Cow_samples_phylum_net2, layout= layout.reingold.tilford(Cow_samples_phylum_net2,circular=TRUE), vertex.size=20, vertex.color="yellow", edge.width=E(Cow_samples_phylum_net2)$weight)
plot(Cow_samples_phylum_net2, layout= layout.reingold.tilford(Cow_samples_phylum_net2,circular=TRUE), vertex.size=20, vertex.color="yellow", edge.width=E(Cow_samples_phylum_net2)$weight)

plot(Cow_samples_phylum_net2, layout=-layout.bipartite(Cow_samples_phylum_net2)[,2:1], vertex.size=30, vertex.shape=ifelse(V(Cow_samples_phylum_net2)$type, "rectangle", "circle"), vertex.color=ifelse(V(Cow_samples_phylum_net2)$type, "red", "cyan"))

A <- get.adjacency(Cow_samples_phylum_net2, sparse=FALSE)
library(network)
g <- network::as.network.matrix(A)
library(sna)
pdf("network_Cow_samples_Phylum.pdf")
sna::gplot.target(g, degree(g), main="Degree", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Cow_samples_phylum_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
sna::gplot.target(g, evcent(g)$vector, main="Eigen Vector", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Cow_samples_phylum_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
sna::gplot.target(g, betweenness(g), main="Eigen Vector", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Cow_samples_phylum_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
sna::gplot.target(g, closeness(g), main="Eigen Vector", circ.col="skyblue",usearrows = FALSE,
                  vertex.col=ifelse(V(Cow_samples_phylum_net2)$type, "red", "cyan"),
                  edge.col="darkgray")
dev.off()
