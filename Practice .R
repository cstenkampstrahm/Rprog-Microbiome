# 3/12/17- going to practice aggregating taxa with the file zaid created for me
# using the phyloseq workflow, can generate the cowdatarich (nonnormed) phyloseq object
# and also the Normcowdatanew (normalized with CSS) phyloseq object
# additionally cowdata is phyloseq with everything (enviro and non) and cowonly and environly
# have only those samples, respectively (there is just different data in their metadata). saved
# files to load easily later

#save(Cowonly, file="Cowonly")
#save(Enviroonly, file="Enviroonly")
#save(Cowdatarich, file="Cowdata")
#save(Normcowonly, file="Normcowdatanew")
#source(file="Phyloseq_workflow.R")
load("Phyloseq files/Cowdata")
Cowdatarich
OTUS <- otu_table(Cowdatarich)
sampledata <- sample_data(Cowdatarich)
taxa <- tax_table(Cowdatarich)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Cowdatarich)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Cowdatarich)

taxanew <- data.frame(taxa)
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)

taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]
library(tidyverse)
taxa.df <-mutate(taxa.df, OTUs = row.names(taxa.df))

# from Zaid's analysis file:
library(pheatmap)
library(vegan)
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
library(cluster)
library(NbClust)
source(file="files from zaid/functions_Chloe.R")

taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df)
#save(taxa1.df, file="splitcowtaxa")
nm.vt = colnames(OTUS)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(OTUS)=nm.vt
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
combined.df = combined.ls[[1]]
#will get a taxa table with aggregates by sample up to the Phylum level