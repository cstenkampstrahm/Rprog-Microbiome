# Want to make a basic overview of microbiome data for manuscript
# will attempt to make a stacked bar chart with all day 1 samples, a subset of 
# samples by pathotype, and then include ordination plots as well
library("tidyverse")
library("ggplot2")
library(phyloseq)
source(file="files from zaid/functions_Chloe.R")
load("Phyloseq files/Normcowdata")
Normcowdata
### double checked that the above table has metagenomeSeq normalized counts in it
### pull out day 1 samples only
Day1 <- c("10A", "14A", "15A", "18A", "1A", "20A", "24A", "29A", "34A", "35A", "38A", "39A", "42A", 
            "44A", "45A", "48A", "50A", "52A", "53A", "55A", "58A", "61A", "62A", "63A", "64A", 
            "65A", "66A", "67A", "68A", "69A", "6A", "70A", "71A", "73A", "74A", "8A")
Normcowdataday1 <- prune_samples(Day1, Normcowdata)
Normcowdataday1
## pull out only non shedding or shedding samples, only one per cow (max # poss for shed = 14)
Sheddingsubset <- c("10A", "18E", "1D", "38A", "45D", "48B", "50E", "57B", "58C", "63E", "67A", 
                      "6D", "74C", "8C",
                      "14E", "15A", "20B", "24C", "29D", "34E", "35A", "39B", "42C", "44D", 
                      "52E", "53A", "55B", "60C")
Normcowdatashedsub <- prune_samples(Sheddingsubset, Normcowdata)
Normcowdatashedsub

#### DAY 1 FIRST ####
## Generation of vegan data frames (OTU file, taxa file, treatment file; reducing
## by OTUs that have nothing in them and also reducing taxa file to only taxa present)
OTUS <- otu_table(Normcowdataday1)
sampledata <- sample_data(Normcowdataday1)
taxa <- tax_table(Normcowdataday1)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Normcowdataday1)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Normcowdataday1)

taxanew <- data.frame(taxa)
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)

taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]
taxa.df <-mutate(taxa.df, OTUs = row.names(taxa.df))
taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df)

## phylum level
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
phylumday1.df = combined.ls[[1]]
bar.taxa.sample.ftn(phylumday1.df,cutoff=.01)
bar.taxa.ftn(phylumday1.df,cutoff=.01)

## family level
combined.ls3 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=5)
familyday1.df = combined.ls3[[1]]
bar.taxa.sample.ftn(familyday1.df,cutoff=.01)
bar.taxa.ftn(familyday1.df,cutoff=.01)


#### NOW SHEDDING SUBSET ####
## Generation of vegan data frames (OTU file, taxa file, treatment file; reducing
## by OTUs that have nothing in them and also reducing taxa file to only taxa present)
OTUS2 <- otu_table(Normcowdatashedsub)
sampledata2 <- sample_data(Normcowdatashedsub)
taxa2 <- tax_table(Normcowdatashedsub)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Normcowdatashedsub)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS2 <- vegan_otu(Normcowdatashedsub)

taxanew2 <- data.frame(taxa2)
d.cs2 = apply(OTUS2,2,sum)
OTUS2 = OTUS2[,d.cs2>0]
otu.nm2 = colnames(OTUS2)

taxa.df2 = taxanew2[row.names(taxanew2)%in%otu.nm2,]
taxa.df2 <-mutate(taxa.df2, OTUs = row.names(taxa.df2))
taxa1.df2 = taxa.cleanup.qiime.ftn(taxa.df=taxa.df2)

## phylum level
combined.ls2 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df2,d.df=OTUS2,l=2)
phylumshedsub.df = combined.ls2[[1]]
bar.taxa.sample.ftn(phylumshedsub.df,cutoff=.01)
bar.taxa.ftn(phylumshedsub.df,cutoff=.01)

## family level
combined.ls4 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df2,d.df=OTUS2,l=5)
familyshedsub.df = combined.ls4[[1]]
bar.taxa.sample.ftn(familyshedsub.df,cutoff=.01)
bar.taxa.ftn(familyshedsub.df,cutoff=.01)

## How to get the shedding and non-shedding samples in a row?
## Did the below before for diff abundance plot... not working with this plot
xaxistitles <- c("10A", "18E", "1D", "38A", "45D", "48B", "50E", "57B", "58C", "63E", "67A", 
                 "6D", "74C", "8C",
                 "14E", "15A", "20B", "24C", "29D", "34E", "35A", "39B", "42C", "44D", 
                 "52E", "53A", "55B", "60C")
plot1 + scale_x_discrete(limits = xaxistitles)

# will try to do them as two separate data frames and then cut and merge the two

shedding <- c("10A", "18E", "1D", "38A", "45D", "48B", "50E", "57B", "58C", "63E", "67A", 
              "6D", "74C", "8C")
nonshedding <- c("14E", "15A", "20B", "24C", "29D", "34E", "35A", "39B", "42C", "44D", 
                 "52E", "53A", "55B", "60C")

## shedding
Normcowdatashed <- prune_samples(shedding, Normcowdata)
Normcowdatashed
OTUS3 <- otu_table(Normcowdatashed)
sampledata3 <- sample_data(Normcowdatashed)
taxa3 <- tax_table(Normcowdatashed)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Normcowdatashed)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS3 <- vegan_otu(Normcowdatashed)

taxanew3 <- data.frame(taxa3)
d.cs3 = apply(OTUS3,2,sum)
OTUS3 = OTUS3[,d.cs3>0]
otu.nm3 = colnames(OTUS3)

taxa.df3 = taxanew3[row.names(taxanew3)%in%otu.nm3,]
taxa.df3 <-mutate(taxa.df3, OTUs = row.names(taxa.df3))
taxa1.df3 = taxa.cleanup.qiime.ftn(taxa.df=taxa.df3)

## phylum level
combined.ls5 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df3,d.df=OTUS3,l=2)
phylumshed.df = combined.ls5[[1]]
bar.taxa.sample.ftn(phylumshed.df,cutoff=.01)
bar.taxa.ftn(phylumshed.df,cutoff=.01)
## family level
combined.ls6 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df3,d.df=OTUS3,l=5)
familyshed.df = combined.ls6[[1]]
bar.taxa.sample.ftn(familyshed.df,cutoff=.01)
bar.taxa.ftn(familyshed.df,cutoff=.01)

## nonshedding
Normcowdatanon <- prune_samples(nonshedding, Normcowdata)
Normcowdatanon
OTUS4 <- otu_table(Normcowdatanon)
sampledata4 <- sample_data(Normcowdatanon)
taxa4 <- tax_table(Normcowdatanon)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Normcowdatanon)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS4 <- vegan_otu(Normcowdatanon)

taxanew4 <- data.frame(taxa4)
d.cs4 = apply(OTUS4,2,sum)
OTUS4 = OTUS4[,d.cs4>0]
otu.nm4 = colnames(OTUS4)

taxa.df4 = taxanew4[row.names(taxanew4)%in%otu.nm4,]
taxa.df4 <-mutate(taxa.df4, OTUs = row.names(taxa.df4))
taxa1.df4 = taxa.cleanup.qiime.ftn(taxa.df=taxa.df4)

## phylum level
combined.ls7 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df4,d.df=OTUS4,l=2)
phylumnon.df = combined.ls7[[1]]
bar.taxa.sample.ftn(phylumnon.df,cutoff=.01)
bar.taxa.ftn(phylumnon.df,cutoff=.01)
## family level
combined.ls7 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df4,d.df=OTUS4,l=5)
familynon.df = combined.ls7[[1]]
bar.taxa.sample.ftn(familynon.df,cutoff=.01)
bar.taxa.ftn(familynon.df,cutoff=.01)

plot1 + bar.taxa.ftn(familynon.df,cutoff=.01)

## the above will not suffice because slightly different abundances from total
## in shedding and non shedding. will splice out original stacked bars




# Getting relative abundance (after the CSS normalization) to compare totals:
# want to use the whole normed data set; Normcowdata

abundnormdata  = transform_sample_counts(Normcowdata, function(x) x / sum(x) )
abundnormfilt = filter_taxa(abundnormdata, function(x) mean(x) > 1e-5, TRUE)
# now in relative abundance by sample
rank_names(abundnormfilt)

famglob <- tax_glom(abundnormdata, taxrank="Rank5")

famtaxglob <- tax_table(famglob)[,"Rank5"]
famotu <- otu_table(famglob)
famglomTable = merge(famotu,famtaxglob,by=0,all=TRUE)
library(xlsx)
write.xlsx(famglomTable, "FamTable.xlsx")

phyloglob <- tax_glom(abundnormdata, taxrank="Rank2")
phylotaxglob <- tax_table(phyloglob)[,"Rank2"]
phylootu <- otu_table(phyloglob)
phyloglomTable = merge(phylootu,phylotaxglob,by=0,all=TRUE)
write.xlsx(phyloglomTable, "PhyloTable.xlsx")

## took to excel and divided shedders and nonshedders, 
## 176 shedders and 22 non shedders
## divided each row sum by this to get the average % by shedding or
## non-shedding
