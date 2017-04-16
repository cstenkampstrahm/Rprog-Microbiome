# starting with phyloseq object generated in phyloseq workflow, called 'Normcowdata'
# this file was generated and has counts normalized via CSS in the OTU table. the 
# normalized counts were made without the 18A and 7A samples
library(pheatmap)
library(vegan)
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
library(cluster)
library(NbClust)
library(tidyverse)
source(file="files from zaid/functions_Chloe.R")
load("Phyloseq files/Normcowdata") 
Normcowdata
OTUS <- otu_table(Normcowdata)
sampledata <- sample_data(Normcowdata)
taxa <- tax_table(Normcowdata)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Normcowdata)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Normcowdata)

taxanew <- data.frame(taxa)
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)

taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]
taxa.df <-mutate(taxa.df, OTUs = row.names(taxa.df))
taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df) 

nm.vt = colnames(OTUS)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(OTUS)=nm.vt

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
phylum.df = combined.ls[[1]]

#bar.taxa.sample.ftn(phylum.df,llab="Phylum")
#NbClust(data=phylum.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
#pheatmap.1.ftn(phylum.df,cutoff=3,cellh=3,cellw=10,fonts=7)

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=3)
class.df = combined.ls[[1]]

#bar.taxa.sample.ftn(class.df,llab="Class")
#NbClust(data=class.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
#pheatmap.1.ftn(class.df,cutoff=3,cellh=3,cellw=5,fonts=7)

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=4)
order.df = combined.ls[[1]]

#bar.taxa.sample.ftn(order.df,llab="Order")
#NbClust(data=order.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
#pheatmap.1.ftn(order.df,cutoff=3,cellh=3,cellw=5,fonts=7)

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=5)
family.df = combined.ls[[1]]

#bar.taxa.sample.ftn(family.df,llab="Family")
#NbClust(data=family.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
#pheatmap.1.ftn(family.df,cutoff=3,cellh=1,cellw=1,fonts=5)

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=6)
genus.df = combined.ls[[1]]

#bar.taxa.sample.ftn(genus.df,llab="Genus")
#NbClust(data=genus.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
#pheatmap.1.ftn(genus.df,cutoff=3,cellh=2,cellw=2,fonts=7)

# need to get the treatment file in and lined up like Zaids

View(sampledata)
sample_data <- data.frame(sampledata)
sample_data <- sample_data %>% mutate(Disease = as.factor(Disease),
                                      DIM = as.numeric(DIM),
                                  Parity = as.factor(Parity),
                                  Farm = as.factor(Farm),
                                  Pattern_1 = as.factor(Pattern_1),
                                  EvNev_1 = as.factor(EvNev_1),
                                  Pathotype_1 = as.factor(Pathotype_1),
                                  Individual_animal = as.factor(Individual_animal))
row.names(sample_data) <- sample_data[,1]
View(sample_data)


##Trying some code from Zaid's class:######
#bar.ftn(class.df,sample_data[,23]) # doesn't work well
#pheatmap.2.ftn(class.df,sample_data[,23:24],cutoff=3,cellh=3,cellw=3,fonts=7) # kind of
# works to see a heatmap of OTUs by a sample category
#ord.plot.ftn(class.df,sample_data[,23],x=c(-0.5,0.5),y=c(-.3,.3),type=6,cex=0.5) # doesn't work

########### NOW ORDINATING BASED on CSS NORMALIZED CLASS LEVEL ##################
class.otu = phyloseq(otu_table(class.df,taxa_are_rows = FALSE),sample_data(sample_data))
sample_data(sample_data)
class.ord = ordinate(class.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(class.otu,class.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue"))

# need to make a phyloseq tree using the 'ape' package according to phyloseq 
# github page, and save as the phytree object to be used in Unifrac creation.
library(ape)
random_tree = rtree(ncol(class.df), rooted=TRUE, tip.label= colnames(class.df))
plot(random_tree)
class.otu <- merge_phyloseq(otu_table(class.df,taxa_are_rows = FALSE),
                            sample_data(sample_data), random_tree)
class.ordU = ordinate(class.otu, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(class.otu, class.ordU, color="Pathotype_1", shape="Parity")
# getting errory that phy_tree slot is empty
