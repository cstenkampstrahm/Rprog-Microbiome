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
combined.ls = taxa.split.combine.qimme.ftn(taxa.df=taxa1.df,d.df=OTUS,1=7
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


