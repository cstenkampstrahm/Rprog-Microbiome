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

taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df) # added a note in here so can
# save this generated taxa table with taxa parsed out by level
#save(taxa1.df, file="splitcowtaxa")
#giving the OTUS in the split taxa names based on those in the OTU table
nm.vt = colnames(OTUS)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(OTUS)=nm.vt
# now combo the OTU table and the taxa table based on the level desired
## PHYLUM
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
phylum.df = combined.ls[[1]]
#will get a taxa table with aggregates by sample up to the Phylum level. doing this
#for all different levels now to generate data files to work with
## CLASS
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=3)
class.df = combined.ls[[1]]
## ORDER
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=4)
order.df = combined.ls[[1]]
## FAMILY
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=5)
family.df = combined.ls[[1]]
## GENUS
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=6)
genus.df = combined.ls[[1]]


# want to look at spread of data for Table 1:
library("xlsx")
table_1 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormed.xlsx", 1)
pattern <- table_1 %>% group_by(Pattern_1, Farm) %>% summarise(pattnval = n())
pathotype <- table_1 %>% group_by(Pathotype_1, Farm) %>% summarise(pathnval = n())
evernever <- table_1 %>% group_by(EvNev_1, Farm) %>% summarise(evnevnval = n())
