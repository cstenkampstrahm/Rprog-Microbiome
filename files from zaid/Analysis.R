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

### Reading the data and checking the samples
taxa.df = read.table(file="files from zaid/results/results.0.03.cons.taxonomy",header=TRUE)
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
### (6 = otu/species) otu analysis
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
#plot_richness(d.otu)
d.rich = estimate_richness(d.otu)
rich.mt = t(rarefy(d.df,min(d.rs),se=TRUE))
rich.mt
rich.df = data.frame(cbind(d.rich[,c(1,6,8)],rich.mt[,1]))
names(rich.df) = c("Observed","Shannon","InvSimpson","Richness")
df.plot.ftn(rich.df)

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
#plot_richness(d.otu)

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
pheatmap.1.ftn(d.df,cutoff=3,cellh=50,cellw=10,fonts=7)



################################## Experiment
###### Looking at an experiment
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=norm.mt,l=5)
d.df = d.taxa.ls[[1]]
trt.df = read.table(file="metadata/samples.file",header=FALSE)
trt.df = data.frame(trt.df[,c(1,2)],row.names=trt.df[,3])
names(trt.df) = c("Trt","rep")

bar.ftn(d.df,trt.df)
pheatmap.2.ftn(d.df,trt.df,cutoff=3,cellh=30,cellw=10,fonts=7)
ord.plot.ftn(d.df,trt.df[,1],x=c(-0.5,0.5),y=c(-.3,.3),type=6,cex=0.5)
d.otu = phyloseq(otu_table(d.df,taxa_are_rows = FALSE),sample_data(trt.df))
sample_data(trt.df)
d.ord = ordinate(d.otu, "CCA", "bray")
theme_set(theme_bw())
plot_ordination(d.otu,d.ord,"samples",color="Trt")

##### Simple analysis
## Diversity analysis (old fashioned ecological approach [richness and diversity])
d.otu = otu_table(d.df,taxa_are_rows = FALSE)
d.rich = estimate_richness(d.otu)
rich.mt = t(rarefy(d.df,min(d.rs),se=TRUE))
rich.df = data.frame(trt.df[,1],d.rich[,c(1,6,8)],rich.mt[,1])
names(rich.df) = c("trt","Observed","Shannon","InvSimpson","Richness")

## simple anova
rich.lm = lm(Shannon~trt,data=rich.df)
anova(rich.lm)

## redundancy analysis
data.hellinger <- decostand(d.df, "hellinger") #Hellinger transformation to linearize relationship
mod1 <- rda(data.hellinger~ trt.df[,1])
mod0 <- rda(data.hellinger ~ 1)
mod <- step(mod0, scope = formula(mod1), test = "perm")
#Backward elimination
modb <- step(mod1, scope = formula(mod0), test = "perm")


rda<- rda(data.hellinger~ trt.df[,1], distance= "euclidean")
rda
plot(rda, type="text", xlim=c(-0.5,0.5), ylim=c(-0.5,0.5)) #2D plot based on best model including woman as factor