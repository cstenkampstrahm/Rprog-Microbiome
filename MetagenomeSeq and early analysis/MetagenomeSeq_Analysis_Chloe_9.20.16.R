# Chloe Stenkamp
# Script for Microbiome data analysis in biome, metagenomeSeq, vegan
# Following source code Enrique produced 9/20/16

library(metagenomeSeq)
library("biom")
library(vegan)

# to find out where the library is being stored:
.libPaths()

# if the library is not saved in the appropriate working directory location, 
# can indicate the package to call for an R script by listing script and then
# function. Example:
# y2 <- metagenomeSeq::biom2MRexperiment(x)

# Read in a biom file produced by QIIME, convert to sample MRexperiment 
# R-object, add sample metadata
setwd("/Users/Chloe/Dropbox/Chloe")
biom_file <- "otu_table_mc2_w_tax_no_pynast_failuresfiltered.biom"
biom <- read_biom(biom_file)

# Look at different properties of the biom file read in to R
show(biom)
print(biom)
header(biom)
biom_data(biom)
biom_shape(biom) # shows a truncated version of a biom matrix?
nrow(biom)
ncol(biom)
observation_metadata(biom) # gives you the taxonomy lists found in the biom 
# file, and total counts

# convert biom file to an MRexperiment (biomformat library may automatically 
# load on its own?
library(biomformat)
biom <- biom2MRexperiment(biom)
#View(fData(biom))
#View(pData(biom))
#View(MRcounts(biom))
# the above are three different objects inherent to the MRexperiment files. 
# you can call on View to print a matrix of these items to the console
# fData is 'feature data'or the taxonomy associated with the counts
# pData is 'phenotype data' or the metadata associated with the samples 
# (wont be any until you load the phenotypic data)
# MRCounts are the counts of each of the taxa classes by sample number 
# (row, column)

# load phenotypic (aka sample metadata) data to group with biom file counts 
# and taxonomy:
# new file (made 9/19/16) has a column specified as sample_type 
# (cow, worker_button, air, soil)
# need to specify 'tran' in the following function because covariates are 
# in columns and samples along rows
pheno_metadata <- load_phenoData("new_QIIME_map_Rprog.txt", sep = '\t', 
                                 tran= TRUE)

# remove rows that don't match the count table (a row in the metadata is the 
# sample number, a column in the biom file called MRcounts is the sample number)
# then pull out the row 'ord' in pheno_metadata
ord<-match(colnames(MRcounts(biom)), rownames(pheno_metadata))
pheno_metadata <- pheno_metadata[ord, ]
dim(pheno_metadata)

# need to tie the biom count table (fData of the biom MR experiment) to the 
# sample metadata (pheno_metadata) by first converting with AnnotatedDataFrame:
pheno_metadata=AnnotatedDataFrame(pheno_metadata)
features=AnnotatedDataFrame(fData(biom))

# Now that all are annotaded data frames, tie MRcounts, fData, and 
# pheno_metadata together in a new MR experiment:
cow_microbiome = newMRexperiment(MRcounts(biom),phenoData=pheno_metadata, 
                                 featureData = features)

# omit sample taxa/OTU that have lower than 3 total counts when summed across 
# all samples
# sparseFeatures = which(rowSums(MRcounts(cow_microbiome) > 0)< 3) 
# not sure how this annotation works, >0 and less than 3??
# maybe another option : 
newsparsefeatures = which(rowSums(MRcounts(cow_microbiome))<= 3)  : 
# this gives us 73,593 OTUs as sparse features
length(newsparsefeatures) # this seems like way too many sparse features 99,042
cow_microbiome_trimmed = cow_microbiome[-newsparsefeatures, ]
cow_microbiome_trimmed = cow_microbiome[-sparseFeatures, ]
View(MRcounts(cow_microbiome_trimmed)) # shows that now there are only 60,312 OTUs

# normalization step. metagenomeSeq uses cumulative sum scaling; calculates 
# scaling factors equal to the sum of counts up to a particular quantile. 
# This aims to account for sparsity due to undersampling
# first step is to calculate proper percentile (p) by which to normalize counts: 
p = cumNormStatFast(cow_microbiome_trimmed)
# then use that percentile to calculate scaling factors:
microbiometrimnorm <- cumNorm(cow_microbiome_trimmed, p = p)

# to export normalized count matrices, to the working directory:
mat <- MRcounts(microbiometrimnorm, norm = TRUE, log = TRUE)[1:5, 1:5] # specify 
# the row and columns for only the first few lines of matrix. Now can export the 
# normalized count matrices:
exportMat(mat, file = "microbiometrim.tsv")

# to look at sample statistics like scaling factor, quantile value, number of 
# id'd features and library size:
exportStats(microbiometrimnorm[ , 1:5], file = "microbiomenormstats.tsv")

# to view the first few features in the saved normalization file:
head(read.csv(file = "microbiomenormstats.tsv"), sep = "\t")


# First thing to do before this is to aggregate the counts to the different 
# taxonomic levels! This is not annotated in the MetagenomeSeq file, but we want
# to detect things that are differentially abundant by their level and not be 
# effected by whatever level the sequencing actually hit to. If we have things 
# that aren't classified to that level though (say if we want to look by genus) 
# those will alot to a zero count and not be included in the analysis: 

# This is the code for aggregating to each taxanomic level, will just do phylum
# level analysis to start:
MicroPhylum <- as.character(fData(microbiometrimnorm)[, 2])
Phylummicrobiome <- aggTax(MRcounts(microbiometrimnorm, norm=TRUE, log=FALSE),
                         norm=TRUE, log=FALSE, lvl = MicroPhylum, out = 
                           "MRexperiment")
ordp <- match(colnames(MRcounts(Phylummicrobiome)), rownames(pheno_metadata))
phenop <- pheno_metadata[ordp, ]
Phylum_cow_microbiome <- newMRexperiment(MRcounts(Phylummicrobiome),
                                         phenoData=phenop, featureData=NULL)

# Now will make categories- need to remove the environmental samples, do it 
# based on the presence of an NA in the 'fecal' column:
envirosamples <- which(is.na(pData(microbiometrimnorm)$Fecal))
Cow_samples= microbiometrimnorm[, - envirosamples]
Cow_samples_Phylum= Phylum_cow_microbiome[, - envirosamples]
# do we need to do cumulative sum scaling AFTER we remove the duplicate samples??
# seems like the CSS would be different with just cows and their taxa vs all
# make a data set with only laluna cows by milk weight
non_milk_samples <- which(is.na(pData(microbiometrimnorm)$Week_avglbs))
laluna_samples <- microbiometrimnorm[, - envirosamples]
laluna_samples <- laluna_samples[, - non_milk_samples]
Cow_samples_Phylum_LL <- Cow_samples_Phylum_LL[, - non_milk_samples]
levels(pData(Cow_samples_Phylum_LL)$Week_avglbs)
levels(pData(laluna_samples)$Week_avglbs)
is.na(pData(Cow_samples_Phylum_LL)$Week_avglbs) # ugh still shows NA in here!
non_milk_samples <- which(is.na(pData(Cow_samples_Phylum_LL)$Week_avglbs))
#new_milk_samples <- filter(pData(microbiometrimnorm), !is.na(Week_avglbs)) # this 
# script works better as a single script compared to all the others
Cow_samples_Phylum_LL <- Cow_samples_Phylum_LL[, -non_milk_samples]
# these last two levels should have the same amount of categories in them
# need to pull out LL samples only prior to CSS??
# Make descriptive categories based on the pData:
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
milk_wt=pData(Cow_samples_Phylum_LL)$Week_avglbs
parity=pData(Cow_samples)$Parity



# Ordination code #####################
# his is to make one PDF per taxa level 
# (you can remove the comments from each part if you want particular pdfs)
# per Enrique looking at ANOSIM or analysis of similarities- the R statistic 
# can either be high or low, and both may be significant. One means that there 
# is more within sample differences and one means there is more between group
# differences?? Need to look it up
# first look at ordi by phylum grouped by cow
pdf("ordi_micro_cows_Phylum.pdf")
dist_micro_cow_byenviro <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                        norm=TRUE)), "hell"),
                                   "bray")
metaMDS_micro_cow_byenviro <- metaMDS(dist_micro_cow_byenviro, 
                                      distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_byenviro, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byenviro, display="sites", pch=20, 
       col=as.numeric(cow_enviro))
text(metaMDS_micro_cow_byenviro, lab=row.names(metaMDS_micro_cow_byenviro))
groupz <- sort(unique(cow_enviro))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byenviro, 
                                  cow_enviro,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(cow_enviro), col=factor(cow_enviro),
       pch = 20, title = "cow_enviro", bty = "n")
dev.off() #turn this off if you want to read through and make multiple pdfs at 
# one time
anosim(dist_micro_cow_byenviro,cow_enviro) #statistically assesses difference in 
# two or more sampling units

#now look at ordi with phylum grouped by DIM
pdf("ordi_micro_cows1_Phylum_byDIM.pdf")
dist_micro_cow_byDIM <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                     norm=TRUE)), "hell"))
metaMDS_micro_cow_byDIM <- metaMDS(dist_micro_cow_byDIM, distance="none",
                                   symmetric=TRUE)
plot(metaMDS_micro_cow_byDIM, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byDIM, display="sites", pch=20, col=as.numeric(DIM))
text(metaMDS_micro_cow_byDIM, lab=row.names(metaMDS_micro_cow_byDIM))
groupz <- sort(unique(DIM))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byDIM, DIM,font=2, 
                                  cex=1.5, col=i, show.groups=groupz[i])}
legend("topright", legend=levels(DIM), fill=levels(DIM), pch = 20, 
       title = "DIM",bty = "n")
dev.off()
anosim(dist_micro_cow_byDIM,DIM)

# now look at ordi by phylum grouped by Farm
pdf("ordi_micro_cows_Phylum_byFarm.pdf")
dist_micro_cow_byFarm <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                      norm=TRUE)), "hell"))
metaMDS_micro_cow_byFarm <- metaMDS(dist_micro_cow_byFarm, distance="none",
                                    symmetric=TRUE)
plot(metaMDS_micro_cow_byFarm, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byFarm, display="sites", pch=20, col=as.numeric(Farm))
text(metaMDS_micro_cow_byFarm, lab=row.names(metaMDS_micro_cow_byFarm))
groupz <- sort(unique(Farm))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byFarm, 
                                  Farm,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(Farm),fill=levels(Farm), pch = 20, 
       title = "Farm",bty = "n")
dev.off()
anosim(dist_micro_cow_byFarm,Farm)

# now by pathotype
pdf("ordi_micro_cows_Phylum_byPathotype.pdf")
dist_micro_cow_byPathotype <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                        norm=TRUE)), "hell"))
metaMDS_micro_cow_byPathotype <- metaMDS(dist_micro_cow_byPathotype, 
                                         distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_byPathotype, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byPathotype, display="sites", pch=20, 
                col=as.numeric(Pathotype))
legend("topright", legend=levels(Pathotype), col=as.numeric(Pathotype), 
       pch = 20, title = "Pathotype", bty = "n")
text(metaMDS_micro_cow_byPathotype, lab=row.names(metaMDS_micro_cow_byPathotype))
groupz <- sort(unique(Pathotype))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byPathotype, 
                                    Pathotype,font=2, cex=1.5, col=i, 
                                    show.groups=groupz[i])}

dev.off()
anosim(dist_micro_cow_byPathotype,Pathotype)

# now look at pattern:
pdf("ordi_micro_cows_Phylum_byPattern.pdf")
dist_micro_cow_byPattern <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                         norm=TRUE)), "hell"))
metaMDS_micro_cow_byPattern <- metaMDS(dist_micro_cow_byPattern, 
                                  distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_byPattern, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byPattern, display="sites", pch=20, 
         col=as.numeric(Pattern))
text(metaMDS_micro_cow_byPattern, lab=row.names(metaMDS_micro_cow_byPattern))
groupz <- sort(unique(Pattern))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byPattern, 
                                    Pattern,font=2, cex=1.5, col=i, 
                                    show.groups=groupz[i])}
legend("topright", legend=levels(Pattern), pch = 20, col=as.numeric(Pattern),
         title = "Pattern", bty = "n")
dev.off()
anosim(dist_micro_cow_byPattern,Pattern)

# now by treatment:
pdf("ordi_micro_cows_Phylum_byTreatment.pdf")
dist_micro_cow_byTreatment <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                           norm=TRUE)), "hell"))
metaMDS_micro_cow_byTreatment <- metaMDS(dist_micro_cow_byTreatment, 
                                         distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_byTreatment, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byTreatment, display="sites", pch=20, 
       col=as.numeric(Treatment))
text(metaMDS_micro_cow_byTreatment, lab=row.names(metaMDS_micro_cow_byTreatment))
groupz <- sort(unique(Treatment))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byTreatment, 
                                  Treatment,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(Treatment),pch = 20, title = "Treatment",
       bty = "n")
dev.off()
anosim(dist_micro_cow_byTreatment,Treatment)

# now by disease: 
pdf("ordi_micro_cows_Phylum_byDisease.pdf")
dist_micro_cow_byDisease <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                         norm=TRUE)), "hell"))
metaMDS_micro_cow_byDisease <- metaMDS(dist_micro_cow_byDisease, 
                                       distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_byDisease, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byDisease, display="sites", pch=20, 
       col=as.numeric(Disease))
text(metaMDS_micro_cow_byDisease, lab=row.names(metaMDS_micro_cow_byDisease))
groupz <- sort(unique(Disease))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byDisease, 
                                  Disease,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(Disease),pch = 20, title = "Disease",bty = "n")
dev.off()
anosim(dist_micro_cow_byDisease,Disease)

# in enriques code also have cworker, ind_animal listed. These would make the 
# same output as cow_enviro did (just classified based on the cow)
# now will look by sample_ID (this is A, B, C, D, E for the diff days)
pdf("ordi_micro_cows_Phylum_bysample_ID.pdf")
dist_micro_cow_bysample_ID <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                           norm=TRUE)), "hell"))
metaMDS_micro_cow_bysample_ID <- metaMDS(dist_micro_cow_bysample_ID, 
                                         distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_bysample_ID, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_bysample_ID, display="sites", pch=20, 
       col=as.numeric(sample_ID))
text(metaMDS_micro_cow_bysample_ID, lab=row.names(metaMDS_micro_cow_bysample_ID))
groupz <- sort(unique(sample_ID))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_bysample_ID, 
                                  sample_ID,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(sample_ID), fill=levels(as.numeric(sample_ID)), 
       pch = 20, title = "sample_ID", bty = "n")
dev.off()
anosim(dist_micro_cow_bysample_ID,sample_ID)

# now will look by parity:
pdf("ordi_micro_cows_Phylum_byparity.pdf")
dist_micro_cow_byparity <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum, 
                                                           norm=TRUE)), "hell"))
metaMDS_micro_cow_byparity <- metaMDS(dist_micro_cow_byparity, 
                                         distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_byparity, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_byparity, display="sites", pch=20, 
       col=as.numeric(parity))
text(metaMDS_micro_cow_byparity, lab=row.names(metaMDS_micro_cow_byparity))
groupz <- sort(unique(parity))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_byparity, 
                                  parity,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(parity), fill=levels(parity), pch = 20, 
       title = "Parity", bty = "n")
dev.off()
anosim(dist_micro_cow_byparity,parity)

# now phylum by milk weight of la luna only cows:
pdf("ordi_micro_cows_Phylum_bymilkwt.pdf")
dist_micro_cow_bymilkwt <- vegdist(decostand(t(MRcounts(Cow_samples_Phylum_LL, 
                                                           norm=TRUE)), "hell"))
metaMDS_micro_cow_bymilkwt <- metaMDS(dist_micro_cow_bymilkwt, 
                                         distance="none",symmetric=TRUE)
plot(metaMDS_micro_cow_bymilkwt, type="none", display=c("sites"),axes=TRUE)
points(metaMDS_micro_cow_bymilkwt, display="sites", pch=20, 
       col=as.numeric(milk_wt))
text(metaMDS_micro_cow_bymilkwt, lab=row.names(metaMDS_micro_cow_bymilkwt))
groupz <- sort(unique(milk_wt))
for(i in seq(groupz)) {ordispider(metaMDS_micro_cow_bymilkwt, 
                                  milk_wt,font=2, cex=1.5, col=i, 
                                  show.groups=groupz[i])}
legend("topright", legend=levels(milk_wt), fill =levels(milk_wt), pch = 20, 
       title = "Milk Avg Lbs",
       bty = "n")
dev.off()
anosim(dist_micro_cow_bymilkwt,milk_wt)

#dev.off # could put this code right here if wanted all Pdfs to write to the 
# same file (like a file with all the phylum comparisons in it)
# now can go back through and do NMDS for the different categories, at each of
# the taxonomic levels
# can look at differential abundance with zero inflated gaussian. One of the 
# first steps based on the metagenomeseq file seem to be calculate normalization
# factors using the trimmed and cumulative sum scaled data:
normfactor = normFactors(Cow_samples)
normfactor <- log2(normfactor/median(normfactor)+1)
# the notation also discusses remove features (taxa) based on the number of 
# estimated effective samples (see calculateEffectiveSamples- this estimates the effective samples per feature.
# Have not done this
settings = zigControl(maxit=20, verbose=TRUE)
designPathoCow = model.matrix(~Pathotype + cow_enviro + normfactor)

#change the above to the variables of interest in modeling.....
dup_resMicroTime <- duplicateCorrelation(MRcounts(Microbiometrimnorm),
                    block=sample_ID, design= designMicroTime)res
resMicroTime = fitZig(obj= Microbiometrimnorm, mod = designMicroTime, 
                      control = settings, useCSSoffset=FALSE, 
                      useMixedModel=dup_resMicroTime$consensus)
# the above is in ED code not sure what it does? Cant find function duplicate 
# correlation anywhere... not sure how this changes things? if cow_enviro or 
# samp_ID or whatever we want to block for is already in the model.......
# per ED on 10/3: this is another way to account for blocking by animal, it 
# results in more conservative p-values for differential abundance estimates. 
# ED and Morley group has emailed Paulsen about this, is it best to include the 
# variable as a covariate in the model or to create a mixed model separately and 
# incorporate it into the fitZig adjustment?

zigFitPathoCow = resPathoCow$fit
finalModPathoCow = resPathoCow$fit$design
# now can specify contrasts and F tests to compare multiple groups
contrastPathoCow = makeContrasts(Pattern, parity,
      levels=finalModPathoCow)
resPatho2Cow = contrasts.fit(zigFitPathoCow, contrastPathoCow)
resPatho2EBCow = eBayes(resPatho2Cow)
topTable(resPatho2EBCow) # I think this prints the contrast F stats and p values?
write.table(topTable(resPatho2EBCow, coef=1, adjust.method="BH",number = 1000), 
            file = "Pathotype-Pattern", sep="\t")
