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

combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=7)
all.df = combined.ls[[1]]

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
#### Ordination with NMDS first ####
######## CLASS ########
class.otu = phyloseq(otu_table(class.df,taxa_are_rows = FALSE),sample_data(sample_data))
sample_data(sample_data)
class.ord = ordinate(class.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(class.otu,class.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Class Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(class.otu,class.ord,"samples",color="Pattern_1") + 
  scale_color_manual(values = c("red", "blue", "green")) + ggtitle("Class Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(class.otu,class.ord,"samples",color="EvNev_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Class Level NMDS with Bray-Curtis
                                                                   (CSS normalized)")

##### ORDER ######
order.otu = phyloseq(otu_table(order.df,taxa_are_rows = FALSE),sample_data(sample_data))
sample_data(sample_data)
order.ord = ordinate(order.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(order.otu,order.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Order Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(order.otu,order.ord,"samples",color="Pattern_1") + 
  scale_color_manual(values = c("red", "blue", "green")) + ggtitle("Order Level NMDS with Bray-Curtis
                                                                   (CSS normalized)")
plot_ordination(order.otu,order.ord,"samples",color="EvNev_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Order Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
####### FAMILY ########
family.otu = phyloseq(otu_table(family.df,taxa_are_rows = FALSE),sample_data(sample_data))
sample_data(sample_data)
family.ord = ordinate(family.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(family.otu,family.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Family Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(family.otu,family.ord,"samples",color="Pattern_1") + 
  scale_color_manual(values = c("red", "blue", "green")) + ggtitle("Family Level NMDS with Bray-Curtis
                                                                   (CSS normalized)")
plot_ordination(family.otu,family.ord,"samples",color="EvNev_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Family Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
###### GENUS ######
genus.otu = phyloseq(otu_table(genus.df,taxa_are_rows = FALSE),sample_data(sample_data))
sample_data(sample_data)
genus.ord = ordinate(genus.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(genus.otu,genus.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Genus Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(genus.otu,genus.ord,"samples",color="Pattern_1") + 
  scale_color_manual(values = c("red", "blue", "green")) + ggtitle("Genus Level NMDS with Bray-Curtis
                                                                   (CSS normalized)")
plot_ordination(genus.otu,genus.ord,"samples",color="EvNev_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Genus Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
##### CLUSTER UP TO SPECIES #####
all.otu = phyloseq(otu_table(all.df,taxa_are_rows = FALSE),sample_data(sample_data))
sample_data(sample_data)
all.ord <- ordinate(all.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(all.otu,all.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Species Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(all.otu,all.ord,"samples",color="Pattern_1") + 
  scale_color_manual(values = c("red", "blue", "green")) + ggtitle("Species Level NMDS with Bray-Curtis
                                                                   (CSS normalized)")
plot_ordination(all.otu,all.ord,"samples",color="EvNev_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("Species Level NMDS with Bray-Curtis
                                                          (CSS normalized)")

#### NON AGGREGATED #####
orig.otu = phyloseq(otu_table(Normcowdata),sample_data(sample_data))
sample_data(sample_data)
orig.ord <- ordinate(orig.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(orig.otu,orig.ord,"samples",color="Pathotype_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("OTU Level NMDS with Bray-Curtis
                                                          (CSS normalized)")
plot_ordination(orig.otu,orig.ord,"samples",color="Pattern_1") + 
  scale_color_manual(values = c("red", "blue", "green")) + ggtitle("OTU Level NMDS with Bray-Curtis
                                                                   (CSS normalized)")
plot_ordination(orig.otu,orig.ord,"samples",color="EvNev_1") + 
  scale_color_manual(values = c("red", "blue")) + ggtitle("OTU Level NMDS with Bray-Curtis
                                                          (CSS normalized)")


#### PCoA with Unifrac now @@@@
# need to make a phyloseq tree using the 'ape' package according to phyloseq 
# github page, and save as the phytree object to be used in Unifrac creation.
#library(ape)
#random_tree = rtree(ncol(class.df), rooted=TRUE, tip.label= colnames(class.df))
#plot(random_tree)
#treefile <- read.tree("rep_set.tre")
# will bring in original OTU table and perform this analysis on the non-aggregated
# data, because the differences between tree tips are what is measured (closer things
# will be less different). Added this code to the phyloseq workflow and saved the 
# phyloseq object as Normcowdatatree
# need to bring in that phyloseq object, but parse it out to mutate the sampledata
# to make factors for plots. then will put back together:

library(ape)
treefile <- read.tree("rep_set.tre")
load("Phyloseq files/Normcowdatatree")
Normcowdatatree
OTUS.1 <- otu_table(Normcowdatatree)
sampledata.1 <- sample_data(Normcowdatatree)
taxa.1 <- tax_table(Normcowdatatree)
sampledata.1 <- data.frame(sampledata.1)
sampledata.1 <- sampledata.1 %>% mutate(Disease = as.factor(Disease),
                                      DIM = as.numeric(DIM),
                                      Parity = as.factor(Parity),
                                      Farm = as.factor(Farm),
                                      Pattern_1 = as.factor(Pattern_1),
                                      EvNev_1 = as.factor(EvNev_1),
                                      Pathotype_1 = as.factor(Pathotype_1),
                                      Individual_animal = as.factor(Individual_animal))
row.names(sampledata.1) <- sampledata.1[,1]
sampledata.1 <- sample_data(sampledata.1)
Normcowdatatree <- merge_phyloseq(OTUS.1, taxa.1, sampledata.1, treefile)

# need the treefile to only contain the colnames of the otu_table. (could have 
# been some in the enviro samples we don't want....)
#OTUS.1 <- data.frame(OTUS.1)
#OTUskeep <- rownames(OTUS.1)
#tips.tree <- treefile$tip.label
#tips.tree <- data.frame(tips.tree)
#OTUskeep <- data.frame(OTUskeep)
#list <- tips.tree[tips.tree%in%OTUskeep,]
#tree <- drop.tip(treefile, list, trim.internal = TRUE)
# after running the above code, there are no OTUs that are being omitted


#class.otu <- merge_phyloseq(otu_table(class.df,taxa_are_rows = FALSE),
#                           sample_data(sample_data), treefile)

# chose to do unweighted because data has been normalized already
ordU = ordinate(Normcowdatatree, "PCoA", "unifrac", weighted=FALSE)
plot_ordination(Normcowdatatree, ordU, color="Pathotype_1") +
  ggtitle("Un-Weighted Unifrac PCoA")
plot_ordination(Normcowdatatree, ordU, color="Pattern_1") +
  ggtitle("Un-Weighted Unifrac PCoA")
plot_ordination(Normcowdatatree, ordU, color="EvNev_1") +
  ggtitle("Un-Weighted Unifrac PCoA")
# what does weighted look like?
ordUW = ordinate(Normcowdatatree, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(Normcowdatatree, ordUW, color="Pathotype_1") +
  ggtitle("Weighted Unifrac PCoA")
plot_ordination(Normcowdatatree, ordUW, color="Pattern_1") +
  ggtitle("Weighted Unifrac PCoA")
plot_ordination(Normcowdatatree, ordUW, color="EvNev_1") +
  ggtitle("Weighted Unifrac PCoA")

# making cluster dendrograms with hellinger transformed data, euclidian distances
# do I need to use non-normalized data prior to hellinger? methinks no (still need
# to account for different sequencing depths) but maybe should check....
clustotu <- vegan_otu(Normcowdata)
clustotu.h <- decostand(clustotu, "hellinger")
clustotu.h.d <- vegdist(clustotu.h, "euclidean")
plot(hclust(clustotu.h.d, method="average"), hang = 0.4, cex = 0.6, 
  main = "Cluster Diagram Using Hellinger-transfrormed Distance")

# How much does it matter that we look at different levels? Class, Order, Family
# Genus etc? Can sum the column counts for the all.df data frame?
dim(genus.df)
#[1] 196 490
dim(family.df)
#[1] 196 235
dim(order.df)
#[1] 196 109
dim(class.df)
#[1] 196  57
dim(phylum.df)
#[1] 196  28
dim(all.df)
#[1] 196 550

# many of the original OTUs (159,000) aggregate to the same level in the taxa 
# table (to the order Clostridiales, for instance)
