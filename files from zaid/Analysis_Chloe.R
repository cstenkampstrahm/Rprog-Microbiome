### Prepared by Zaid Abdo/Metagenomics Class-2017
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq','metagenomeSeq')
library(pheatmap)
library(vegan)
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
library(cluster)
library(NbClust)
source(file="functions.R")
load("Chloe_taxa1.RData")
taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df)
nm.vt = colnames(d.df)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(d.df)=nm.vt
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=d.df,l=2)
combined.df = combined.ls[[1]]










### Reading the data and checking the samples
taxa.df = read.table(file="results/results.0.03.cons.taxonomy",header=TRUE)
data.df = read.table(file="results/results.shared",header=TRUE)
data.df = data.frame(data.df[,-c(1:3)],row.names=data.df[,2])

### Mock check (if you have any) and if you called it "Mock" (use whatever else you called it)
mock.nm = "Mock"
data.df[mock.nm,data.df["Mock",]>0] ## Going to assume 1 and 2 counts more towards error
data.df = data.df - 2
data.df[data.df<0] = 0
data.df = data.df[row.names(data.df)!=mock.nm,]

### Remove otus that have nothing in them
d.cs = apply(data.df,2,sum)
data.df = data.df[,d.cs>0]
otu.nm = names(data.df)
taxa.df = taxa.df[taxa.df$OTU%in%otu.nm,]

####### OTU Analysis
### (7 = otu/species) otu analysis
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=7)
d.df = d.taxa.ls[[1]]
taxa.vt = d.taxa.ls[[2]]

### check distribution of depth of coverage
d.rs = apply(data.df,1,sum)
hist(d.rs)

### did we sample these communities well? (Rarifaction curves)
rarecurve(x=d.df,col=1:length(d.rs))

### A look at proportions before normalizing
bar.taxa.sample.ftn(d.df,cutoff=.01)

### Diversity analysis (old fashioned ecological approach [richness and diversity])
d.otu = otu_table(d.df,taxa_are_rows = FALSE)
d.rich = estimate_richness(d.otu)
rich.mt = t(rarefy(d.df,min(d.rs),se=TRUE))
rich.df = data.frame(cbind(d.rich[,c(1,6,8)],rich.mt[,1]))
names(rich.df) = c("Observed","Shannon","InvSimpson","Richness")
df.plot.ftn(rich.df,type="point")

### Using the smallest depth with proportions
s.depth = min(d.rs)
d.sd.df = round((d.df/d.rs)*s.depth,0)
d.cs = apply(d.sd.df,2,sum)
d.sd.df=d.sd.df[,d.cs>0]
rarecurve(x=d.sd.df,col=1:length(d.rs))
d.otu = otu_table(d.sd.df,taxa_are_rows = FALSE)
d.rich = estimate_richness(d.otu)

### CSS normalization (using this all way through)
norm.vt = cumNormStatFast(as.matrix(t(d.df)))
norm.mt = round(t(cumNormMat(as.matrix(t(d.df)),norm.vt)))
d.otu = otu_table(norm.mt,taxa_are_rows = FALSE)
d.rich = estimate_richness(d.otu)

### Clustering
## Identifying the number of clusters
NbClust(data=norm.mt,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
## plotting the clusters
pheatmap.1.ftn(norm.mt,cutoff=3,cellh=10,cellw=3,fonts=5)

###### Higher level plotting
### Genus (6)
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=norm.mt,l=6)
d.df = d.taxa.ls[[1]]

# A look at proportions before normalizing
bar.taxa.sample.ftn(d.df,llab="Genus")

# Clustering
NbClust(data=d.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
pheatmap.1.ftn(d.df,cutoff=3,cellh=50,cellw=10,fonts=7)

### Family (5)
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=norm.mt,l=5)
d.df = d.taxa.ls[[1]]

# A look at proportions before normalizing
bar.taxa.sample.ftn(d.df,llab="Family")

# Clustering
NbClust(data=d.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
pheatmap.1.ftn(d.df,cutoff=3,cellh=50,cellw=10,fonts=7)

### Phylum (2)
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=norm.mt,l=2)
d.df = d.taxa.ls[[1]]

# A look at proportions before normalizing
bar.taxa.sample.ftn(d.df,llab="Phylum")

# Clustering
NbClust(data=d.df,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
pheatmap.1.ftn(d.df,cutoff=3,cellh=10,cellw=10,fonts=7)