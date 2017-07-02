# Want to make a basic overview of microbiome data for manuscript
# will attempt to make a stacked bar chart with all day 1 samples, a subset of 
# samples by pathotype, and then include ordination plots as well
save(sampledata, file = "sampledata")
load("sampledata")
sampledata
library("tidyverse")
library("ggplot2")
library(phyloseq)
load("Phyloseq files/Cowonly")
Cowonly

library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
MRexp_cowonly

d.gen = aggTax(MRexp_cowonly,lvl="Rank6", norm = TRUE) 
d.spp = aggTax(MRexp_cowonly,lvl="Rank7", norm = TRUE) 
d.fam = aggTax(MRexp_cowonly,lvl="Rank5", norm = TRUE)
d.phyl = aggTax(MRexp_cowonly,lvl="Rank2", norm = TRUE)
d.class = aggTax(MRexp_cowonly,lvl="Rank3", norm = TRUE)
d.ord = aggTax(MRexp_cowonly,lvl="Rank4", norm = TRUE)


genus_cts <- MRcounts(d.gen)
genus_table <- otu_table(genus_cts, taxa_are_rows = TRUE)
genus <- merge_phyloseq(genus_table, sampledata)
View(otu_table(genus))

phylum_cts <- MRcounts(d.phyl)
phylum_table <- otu_table(phylum_cts, taxa_are_rows = TRUE)
phylum <- merge_phyloseq(phylum_table, sampledata)
View(otu_table(phylum))
View(sample_data(phylum))

# check to see what the example data set looks like 
plot_bar(phylum, x = "SampleID", y = "Abundance", fill = row.names(otu_table(phylum)))
data("GlobalPatterns")
Cowonly
# normalized phyloseq object:
Normcowdatanew
plot_bar(Normcowdatanew, x = "Pathotype_1", y = "Abundance", fill = "Genus")

###################################
## Generation of vegan data frames (OTU file, taxa file, treatment file; reducing
## by OTUs that have nothing in them and also reducing taxa file to only taxa present)
OTUS <- otu_table(Normcowdatanew)
sampledata <- sample_data(Normcowdatanew)
taxa <- tax_table(Normcowdatanew)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Normcowdatanew)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Normcowdatanew)

taxanew <- data.frame(taxa)
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]

taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]
library(tidyverse)
taxa.df <-mutate(taxa.df, OTUs = row.names(taxa.df))

source(file="files from zaid/functions_Chloe.R")
#taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df)
load("splitcowtaxa")
taxa1.df
nm.vt = colnames(OTUS)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(OTUS)=nm.vt

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
phylum.df = combined.ls[[1]]
bar.taxa.sample.ftn(phylum.df,cutoff=.05)
bar.taxa.ftn(phylum.df,cutoff=.01)

combined.ls1 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=6)
genus.df = combined.ls1[[1]]
bar.taxa.sample.ftn(genus.df,cutoff=.09)
bar.taxa.ftn(genus.df,cutoff=.01)

combined.ls2 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=3)
class.df = combined.ls2[[1]]
bar.taxa.sample.ftn(class.df,cutoff=.01)
bar.taxa.ftn(class.df,cutoff=.01)

combined.ls3 = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=5)
family.df = combined.ls3[[1]]
bar.taxa.sample.ftn(family.df,cutoff=.01)
bar.taxa.ftn(family.df,cutoff=.01)
