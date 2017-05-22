# Per committee meeting, better to do negative binomial model with DeSeq2 than to
# keep going with the zero inflated gaussian. Log-normal would prob be the best
# but not working with MetagenomeSeq. 
# Looking at DeSeq2 package. Looks like you can import straight from phyloseq.
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
library(phyloseq)
load("Phyloseq files/Cowonly")
Cowonly

# set the reference sample class of Pathotype_1 to "0", make factors for design matrix:
sample_data(Cowonly)$Pathotype_1 <- factor(sample_data(Cowonly)$Pathotype_1) #makes a factor
sample_data(Cowonly)$Individual_animal <- factor(sample_data(Cowonly)$Individual_animal)#makes a factor
sample_data(Cowonly)$Pathotype_1 <- relevel(sample_data(Cowonly)$Pathotype_1, "0")

# convert from phyloseq to DESeq2 with a specified design matrix
# this line of code estimates dispersions based on the design matrix
deseqCowonly <- phyloseq_to_deseq2(Cowonly, ~Individual_animal + Pathotype_1)
# following the McMurdie script on 
#https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
# hers uses a local fit, which is regression of log dispersions over log base mean
# versus just a dispersion mean in parametric
# below changed to parametric and used betaPrior = TRUE to get the zero mean normal prior on 
# the estimated coefficients
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseqCowonly), 1, gm_mean)
deseqCowonly = estimateSizeFactors(deseqCowonly, geoMeans = geoMeans)
deseqCowonly = DESeq(deseqCowonly, fitType="local")
res = results(deseqCowonly)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
# no things that have adjusted p values <0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Cowonly)[rownames(sigtab), ],
                                            "matrix"))
head(sigtab)
## since did not get anything significant, would be good to aggregate the taxa 
## as done prior with metagenomeseq. how best to do this? did we figure out how to 
## do it in phyloseq with Zaid's code?
library(pheatmap)
library(vegan)
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
library(cluster)
library(NbClust)
source(file="files from zaid/functions_Chloe.R")
OTUS <- otu_table(Cowonly)
sampledata <- sample_data(Cowonly)
taxa <- tax_table(Cowonly)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Cowonly)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Cowonly)

taxanew <- data.frame(taxa)
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)

taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]
library(tidyverse)
taxa.df <-mutate(taxa.df, OTUs = row.names(taxa.df))
# from saved previously
load("splitcowtaxa")
head(taxa1.df)
nm.vt = colnames(OTUS)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(OTUS)=nm.vt
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
phylum.df = combined.ls[[1]]
head(phylum.df)
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=3)
class.df = combined.ls[[1]]
head(class.df)
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=4)
order.df = combined.ls[[1]]
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=5)
family.df = combined.ls[[1]]
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=6)
genus.df = combined.ls[[1]]
combined.ls = taxa.split.combine.qimme.ftn(taxa.df=taxa1.df,d.df=OTUS,1=7)
# for some reason can't get the species level aggregation to work. need to look
# into the code from zaid (1=7 in code above does not work)
# now need to make this the count data for the phyloseq object
# original phyloseq OTU table is sample in the columns, and bacteria in the rows
# need to switch rows and column to get back, use transpose function
class.df <- as.data.frame(t(class.df))
genus.df <- as.data.frame(t(genus.df))
family.df <- as.data.frame(t(family.df))
order.df <- as.data.frame(t(order.df))
phylum.df <- as.data.frame(t(phylum.df))

## GENUS LEVEL AGGREGATION  ##
sampledata <- sample_data(Cowonly)
taxatable <- tax_table(Cowonly)
genus.df <- as.matrix(genus.df)
otutable <- otu_table(genus.df, taxa_are_rows=TRUE)
cowgenus <- merge_phyloseq(otutable, sampledata)
sample_data(cowgenus)$Pathotype_1 <- factor(sample_data(cowgenus)$Pathotype_1)
sample_data(cowgenus)$Individual_animal <- factor(sample_data(cowgenus)$Individual_animal)
sample_data(cowgenus)$Pathotype_1 <- relevel(sample_data(cowgenus)$Pathotype_1, "0")
deseqcowgenus <- phyloseq_to_deseq2(cowgenus, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowgenus), 1, gm_mean)
deseqcowgenus = estimateSizeFactors(deseqcowgenus, geoMeans = geoMeans)
deseqcowgenus = DESeq(deseqcowgenus, fitType="local")
res = results(deseqcowgenus)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
# nothing is significant again, after the multiple testing correction


## ORDER LEVEL AGGREGATION ##
order.df <- as.matrix(order.df)
otutable <- otu_table(order.df, taxa_are_rows=TRUE)
coworder <- merge_phyloseq(otutable, sampledata)
sample_data(coworder)$Pathotype_1 <- factor(sample_data(coworder)$Pathotype_1)
sample_data(coworder)$Individual_animal <- factor(sample_data(coworder)$Individual_animal)
sample_data(coworder)$Pathotype_1 <- relevel(sample_data(coworder)$Pathotype_1, "0")
deseqcoworder <- phyloseq_to_deseq2(coworder, ~Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcoworder), 1, gm_mean)
deseqcoworder = estimateSizeFactors(deseqcoworder, geoMeans = geoMeans)
deseqcoworder = DESeq(deseqcoworder, fitType="local")
res = results(deseqcoworder)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
# once again, nothing is significant after correcting for multiple testing
# maybe need to filter to OTUs that are only present in 50 (25%) of the cows?
# in metagenomeseq files from committee_meeting_5.4.17.R file,
# d.famfilt, d.genfilt, and d.sppfilt are those count tables filtered to taxa present
# in 50 samples. Metagenomseq assayData is not that easily converted to different 
# formats than the eSet ones. But can convert these tables to biom ones....
# lots of trouble doing this
# just using MRcounts call to get the tables to matrices
# THESE HAVE BEEN NORMALIZED WITH CUMULATIVE SUM SCALING! #
genus <- MRcounts(d.genfilt)
# Non-normalized
d.gen1 = aggTax(MRexp_cowonly,lvl="Rank6", norm = FALSE)
d.genfilt1 = filterData(d.gen1, present = 50, depth = 1)
genus1 <- MRcounts(d.genfilt1)
genus2 <- MRcounts(d.gen1)

d.spp1 = aggTax(MRexp_cowonly,lvl="Rank7", norm = FALSE) 
d.sppfilt1 = filterData(d.spp1, present = 50, depth = 1)
spp1 <- MRcounts(d.sppfilt1)
spp2 <- MRcounts(d.spp1)

d.fam1 = aggTax(MRexp_cowonly,lvl="Rank5", norm = FALSE)
d.famfilt1 = filterData(d.fam1, present = 50, depth = 1)
family1 <- MRcounts(d.famfilt1)
family2 <- MRcounts(d.fam1)

d.otu1 = filterData(MRexp_cowonly, present = 50, depth = 1)
otu1 <- MRcounts(d.otu1)
otu2 <- MRcounts(MRexp_cowonly)

## GENUS LEVEL WITH AGGREG THROUGH METAGENOMESEQ AND PRESENT IN 50 COWS ##
genusotu <- otu_table(genus1, taxa_are_rows=TRUE)
cowgenus1 <- merge_phyloseq(genusotu, sampledata)
sample_data(cowgenus1)$Pathotype_1 <- factor(sample_data(cowgenus1)$Pathotype_1)
sample_data(cowgenus1)$Individual_animal <- factor(sample_data(cowgenus1)$Individual_animal)
sample_data(cowgenus1)$Pathotype_1 <- relevel(sample_data(cowgenus1)$Pathotype_1, "0")
deseqcowgenus1 <- phyloseq_to_deseq2(cowgenus1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowgenus1), 1, gm_mean)
deseqcowgenus1 = estimateSizeFactors(deseqcowgenus1, geoMeans = geoMeans)
deseqcowgenus1 = DESeq(deseqcowgenus1, fitType="parametric", betaPrior = TRUE)
res = results(deseqcowgenus1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)
# nothing sig diff with multiple testing
alpha = 0.2
sigtab = res[(res$pvalue < alpha), ]
sigtab
unique(res$pvalue)
# pvalues not corrected for multiple testing show a couple <0.06

## SPECIES LEVEL WITH AGGREG THROUGH METAGENOMESEQ AND PRESENT IN 50 COWS ##
speciesotu <- otu_table(spp1, taxa_are_rows=TRUE)
cowspecies1 <- merge_phyloseq(speciesotu, sampledata)
sample_data(cowspecies1)$Pathotype_1 <- factor(sample_data(cowspecies1)$Pathotype_1)
sample_data(cowspecies1)$Individual_animal <- factor(sample_data(cowspecies1)$Individual_animal)
sample_data(cowspecies1)$Pathotype_1 <- relevel(sample_data(cowspecies1)$Pathotype_1, "0")
deseqcowspecies1 <- phyloseq_to_deseq2(cowspecies1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowspecies1), 1, gm_mean)
deseqcowspecies1 = estimateSizeFactors(deseqcowspecies1, geoMeans = geoMeans)
deseqcowspecies1 = DESeq(deseqcowspecies1, fitType="parametric", betaPrior = TRUE)
res = results(deseqcowspecies1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)
# nothing significant across the adjusted p values
alpha = 0.2
sigtab = res[(res$pvalue < alpha), ]
sigtab
unique(res$pvalue)


## FAMILY LEVEL WITH AGGREG THROUGH METAGENOMESEQ AND PRESENT IN 50 COWS ##
familyotu <- otu_table(family1, taxa_are_rows=TRUE)
cowfamily1 <- merge_phyloseq(familyotu, sampledata)
sample_data(cowfamily1)$Pathotype_1 <- factor(sample_data(cowfamily1)$Pathotype_1)
sample_data(cowfamily1)$Individual_animal <- factor(sample_data(cowfamily1)$Individual_animal)
sample_data(cowfamily1)$Pathotype_1 <- relevel(sample_data(cowfamily1)$Pathotype_1, "0")
deseqcowfamily1 <- phyloseq_to_deseq2(cowfamily1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowfamily1), 1, gm_mean)
deseqcowfamily1 = estimateSizeFactors(deseqcowfamily1, geoMeans = geoMeans)
deseqcowfamily1 = DESeq(deseqcowfamily1, fitType="parametric", betaPrior = TRUE)
res = results(deseqcowfamily1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)
# nothing significant across the adjusted p values
alpha = 0.2
sigtab = res[(res$pvalue < alpha), ]
sigtab
unique(res$pvalue)

# Going to try again with OTUs filtered to those present in 50 samples
otuonly <- otu_table(otu1, taxa_are_rows=TRUE)
cowotu1 <- merge_phyloseq(otuonly, sampledata)
sample_data(cowotu1)$Pathotype_1 <- factor(sample_data(cowotu1)$Pathotype_1)
sample_data(cowotu1)$Individual_animal <- factor(sample_data(cowotu1)$Individual_animal)
sample_data(cowotu1)$Pathotype_1 <- relevel(sample_data(cowotu1)$Pathotype_1, "0")
deseqcowotu1 <- phyloseq_to_deseq2(cowotu1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowotu1), 1, gm_mean)
deseqcowotu1 = estimateSizeFactors(deseqcowotu1, geoMeans = geoMeans)
deseqcowotu1 = DESeq(deseqcowotu1, fitType="parametric", betaPrior = TRUE)
res = results(deseqcowotu1)
restest = results(deseqcowotu1, contrast=c("Pathotype_1", "0", "1")) # get the same thing when specifying 
# contrasts this way
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)
# nothing significant across the adjusted p values
alpha = 0.2
sigtab = res[(res$pvalue < alpha), ]
sigtab
unique(res$pvalue)

##### Going to try the above code, but use cumulative sum scaled data to begin 
# with:
d.gen1norm = aggTax(MRexp_cowonly,lvl="Rank6", norm = TRUE)
d.genfilt1norm = filterData(d.gen1norm, present = 50, depth = 1)
genus1norm <- MRcounts(d.genfilt1norm)

genusotu <- otu_table(genus1norm, taxa_are_rows=TRUE)
cowgenus1 <- merge_phyloseq(genusotu, sampledata)
sample_data(cowgenus1)$Pathotype_1 <- factor(sample_data(cowgenus1)$Pathotype_1)
sample_data(cowgenus1)$Individual_animal <- factor(sample_data(cowgenus1)$Individual_animal)
sample_data(cowgenus1)$Pathotype_1 <- relevel(sample_data(cowgenus1)$Pathotype_1, "0")
deseqcowgenus1 <- phyloseq_to_deseq2(cowgenus1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowgenus1), 1, gm_mean)
deseqcowgenus1 = estimateSizeFactors(deseqcowgenus1, geoMeans = geoMeans)
deseqcowgenus1 = DESeq(deseqcowgenus1, fitType="parametric", betaPrior = TRUE)
#deseqcowgenus12 = DESeq(deseqcowgenus1, fitType="parametric")
res = results(deseqcowgenus1)
#res12 = results(deseqcowgenus12) # the results from this, with no betaPrior are the same
# as for betaPrior = TRUE
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)



d.spp1norm = aggTax(MRexp_cowonly,lvl="Rank7", norm = TRUE) 
d.sppfilt1norm = filterData(d.spp1norm, present = 50, depth = 1)
spp1norm <- MRcounts(d.sppfilt1norm)

speciesotu <- otu_table(spp1norm, taxa_are_rows=TRUE)
cowspecies1 <- merge_phyloseq(speciesotu, sampledata)
sample_data(cowspecies1)$Pathotype_1 <- factor(sample_data(cowspecies1)$Pathotype_1)
sample_data(cowspecies1)$Individual_animal <- factor(sample_data(cowspecies1)$Individual_animal)
sample_data(cowspecies1)$Pathotype_1 <- relevel(sample_data(cowspecies1)$Pathotype_1, "0")
deseqcowspecies1 <- phyloseq_to_deseq2(cowspecies1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowspecies1), 1, gm_mean)
deseqcowspecies1 = estimateSizeFactors(deseqcowspecies1, geoMeans = geoMeans)
deseqcowspecies1 = DESeq(deseqcowspecies1, fitType="parametric", betaPrior = TRUE)
res = results(deseqcowspecies1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)

d.fam1norm = aggTax(MRexp_cowonly,lvl="Rank5", norm = TRUE)
d.famfilt1norm = filterData(d.fam1norm, present = 50, depth = 1)
family1norm <- MRcounts(d.famfilt1norm)

familyotu <- otu_table(family1norm, taxa_are_rows=TRUE)
cowfamily1 <- merge_phyloseq(familyotu, sampledata)
sample_data(cowfamily1)$Pathotype_1 <- factor(sample_data(cowfamily1)$Pathotype_1)
sample_data(cowfamily1)$Individual_animal <- factor(sample_data(cowfamily1)$Individual_animal)
sample_data(cowfamily1)$Pathotype_1 <- relevel(sample_data(cowfamily1)$Pathotype_1, "0")
deseqcowfamily1 <- phyloseq_to_deseq2(cowfamily1, ~ Individual_animal + Pathotype_1)
geoMeans = apply(counts(deseqcowfamily1), 1, gm_mean)
deseqcowfamily1 = estimateSizeFactors(deseqcowfamily1, geoMeans = geoMeans)
deseqcowfamily1 = DESeq(deseqcowfamily1, fitType="parametric", betaPrior = TRUE)
res = results(deseqcowfamily1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab
res
unique(res$padj)

# Looks like betaPrior = TRUE gets us the same results for log2fold change...
# want to make sure that there is a 0 in the intercept when the forumal created.
# when attempting to do this, the call for DESeq forces you to put betaprior=false
deseqcowfamily2 <- phyloseq_to_deseq2(cowfamily1, ~Individual_animal + Pathotype_1 - 1)
deseqcowfamily2 <- DESeq(deseqcowfamily2, fitType = "parametric", betaPrior = FALSE)
results2 <- results(deseqcowfamily2, contrast=c("Pathotype_1", "1", "0"))
# these results appear to be different than what was in res12 and res from above family test,
# need to go through and redo with the appropriate contrasts using a model matrix that is 
# having a -1 in the right column
Data sets





## FIXED WHAT WAS BEING TRIED BELOW- NEED BETAPRIOR=TRUE and OKAY ##
## Going to try to get a -1 in the model call.....,specify design~1 in phylo_to_des
# then DESeq put in design
familytest <- phyloseq_to_deseq2(cowfamily1, ~1)
familytest = estimateSizeFactors(familytest, geoMeans = geoMeans)
familytest = DESeq(familytest, fitType="local")



















# following the joey bioconductor script on 
# http://joey711.github.io/phyloseq-extensions/DESeq2.html
deseqCowonly <- DESeq(deseqCowonly, test="Wald", fitType="parametric")
res <- results(deseqCowonly, cooksCutoff = FALSE)
alpha <- 0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Cowonly)[rownames(sigtab), ], "matrix"))
head(sigtab)
# do I need to variance stabilize prior to running the negative binomial? It looks like no,
# based on the notation of what variance stabilization is/is used for, and the fact that
# the McMurdie script does not use it and shes the guru of doing this stuff with 
# microbiome specific data sets
