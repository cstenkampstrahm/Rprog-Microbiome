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

### Reading the data and checking the samples 9 data points
## reading the taxonomy file
taxa.df = read.table(file="results2/results.0.03.cons.taxonomy",header=TRUE)
## reading the data table
data.df = read.table(file="results2/results.shared",header=TRUE)
## removing the first 3 columns and using the second of those as row names
data.df = data.frame(data.df[,-c(1:3)],row.names=data.df[,2])

### Mock check (if you have any) and if you called it "Mock" (use whatever else you called it)
## name of row with mock for this data set is Mock, we create a variable that stores this name
mock.nm = "Mock"
## number of mock bacteria in the data is about 20 so we expect to have 20 OTUs
## with large numbers (here we look at anything > 0)
data.df[mock.nm,data.df["Mock",]>0] ## Going to assume 1 and 2 counts more towards error
## we can see 19 that meets the large size expected and some OTUs with 1 and 2
## we use 2 as the minimum number of reads allowed. We remove it from the data
data.df = data.df - 2
## we set any value in the data frame that is less than 0 to zero
data.df[data.df<0] = 0
## we remove the mock community row
data.df = data.df[row.names(data.df)!=mock.nm,]

### Remove otus that have nothing in them
## the above would have created OTUs with all zero in them in all samples 
## after removing the 1's and 2's
## the following sums over all columns
d.cs = apply(data.df,2,sum)
## then uses this sum (which has zeros corresponding to columns with nothing
## in them) to remove the empty columns (number of OTUS is down from 189 to 160)
data.df = data.df[,d.cs>0]
## it next extracts the names of these otus
otu.nm = names(data.df)
## and matches the taxa to these names to reduce the taxa.df data frame to match
## the taxa remaining in the collection
taxa.df = taxa.df[taxa.df$OTU%in%otu.nm,]

### Reading the treatment file for a simple experiment (one factor)
## file samples.file includes 5 columns: 1) treatment levels, 2) replicate number,
## 3) sample name, 4) name of R1 files and 5) name of R2 files
trt.df = read.table(file="results2/samples.file",header=FALSE)
## we are only interested in the treatment levels and the replicates
## we create a dataframe that takes only the frist two columns and 
## uses the sample names as row names
trt.df = data.frame(trt.df[,c(1,2)],row.names=trt.df[,4])
## we name the columns as "Trt" and "rep", for treatment and replicate 
## respectively
names(trt.df) = c("Trt","rep")

### Making sure that the treatments match the data
data.df = data.df[row.names(trt.df),]

### Choosing level for analysis 
##level = (1=Domain, 2=Phylum,3=Class, 4=Order, 5=Family, 6=Genus, 7=otu/species)
level = 5
## function taxa.split.combine.ftn does this splitting and puts the data in
## a list that has two parts: 1) otu table and 2) taxa table
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=level) # Phulum
## we extract and use the otu table only
d.df = d.taxa.ls[[1]]

### CSS normalization (using this all way through)
## we use cumulative sum scaling to normalize utilizing metagenomeSeq
norm.vt = cumNormStatFast(as.matrix(t(d.df)))
norm.mt = round(t(cumNormMat(as.matrix(t(d.df)),norm.vt)))

### Looking at an experiment
## we use function bar.ftn to create the bar graphs
bar.ftn(norm.mt,trt.df)
## we use NbClust to identify the number of clusters
NbClust(data=norm.mt,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
## we use the identified number of clusters to set the cutoff and create the
## the dendrogram using the pheatmap.2.ftn, cellh=cell height, cellw=cell width
## fonts = font size to allow for control of the size of the dendrogram
pheatmap.2.ftn(norm.mt,trt.df,cutoff=6,cellh=20,cellw=10,fonts=7) # Using the normalized counts
## we use ord.plot.ftn to create an ordination in vegan using bray distance
## the first entry is the normalized data, the second is only the treatment levels
## which is the first column of the treatment table, the x an y entries 
## identify the size of the plot, type (between 1 and 7) allows for manipulation
## of the looks of the graph (try it) and cex controls the font size.
## It generates and NMDS ordination
ord.plot.ftn(norm.mt,trt.df[,1],x=c(-0.5,0.5),y=c(-.3,.3),type=6,cex=0.5,distance = "bray")
## we de the same next using phyloseq (first create a phyloseq object)
## that includes the normalized data=otu_table and the metadata=sample_data
d.otu = phyloseq(otu_table(norm.mt,taxa_are_rows = FALSE),sample_data(trt.df))
## and we create and ordination object 
d.ord = ordinate(d.otu, "NMDS", "bray")
## change the plot's back groun
theme_set(theme_bw())
## and plot the ordination
plot_ordination(d.otu,d.ord,"samples",color="Trt")+geom_polygon(aes(alpha=0.5,fill=Trt))

####### OTU Analysis
#### First analysis on diversity measures and that is done on the otu level
## recovering the otus levels = 7 (7 = otu/species) otu analysis
level = 7
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=level)
d.df = d.taxa.ls[[1]]
## calculating the row sums
d.rs = apply(d.df,1,sum)

### generating the diversity measures 
## we use phyloseq first
d.otu = otu_table(d.df,taxa_are_rows = FALSE)
d.rich = estimate_richness(d.otu)
## look at them
d.rich

## estimating richness using rarefaction and the minimum sample coverage
## rarefy takes the data = d.df here, and the min sample coverage min(d.rs)
## and can calculate the standard error if you tell it to se=TRUE
rich.mt = t(rarefy(d.df,min(d.rs),se=TRUE))
## look at it
rich.mt
## we are only interested in observed, shannon, invers simpson and the estimated 
## richness from vegan so we extract those and combine them to make a table
## and we add the treatment levels as the first column
rich.df = data.frame(trt.df[,1],d.rich[,c(1,6,8)],rich.mt[,1])
## then we name these columns
names(rich.df) = c("trt","Observed","Shannon","InvSimpson","Richness")
## and plot the data
rich.plot.ftn(d.df=rich.df)

### simple anova
## to run an analysis of variance we use the function lm in R
## the equation we use is Shannon ~ trt or InvSimpson ~ trt or Richness ~ trt
## with data = rich.df
rich.lm = lm(Shannon~trt,data=rich.df)
## summary gives estimates of the parameters of the model
summary(rich.lm)
## anova gives the anoalysis of variance table
## last column of the table provides significance
anova(rich.lm)

### redundancy analysis
##level = (1=Domain, 2=Phylum,3=Class, 4=Order, 5=Family, 6=Genus, 7=otu/species)
level = 5
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=level) # Phulum
d.df = d.taxa.ls[[1]]
## CSS normalize
norm.vt = cumNormStatFast(as.matrix(t(d.df)))
norm.mt = round(t(cumNormMat(as.matrix(t(d.df)),norm.vt)))
## create full model using the normalized matrix and the against the treatment levels
mod1 = rda(norm.mt ~ trt.df[,1],distance = "bray")
## create thre reduced model with nothing in it
mod0 = rda(norm.mt ~ 1, distance = "bray")
## Backward elimination for model selection and look at AIC (small is good)
modb = step(mod1, scope = formula(mod0), test = "perm")

### we want to plot the ordination resulting from the rda analysis
## first we extract the treatment levels
trt = trt.df[,1]
## we run the rda analysis using the best model from the previous step
rda = rda(norm.mt ~ trt, distance= "bray")
## choosing colors of the treatment lvels
c = 1:length(levels(trt))
## we plot c is fo collors here and cex is for font size
plot(rda, type="n")
ordiellipse(rda, trt,col=c,lwd=1,label=TRUE)
ordispider(rda, trt,col=c,label = FALSE)
text(rda,display="site",cex=0.5) 
text(rda,display="species",cex=0.5) 


##################################### Two factors anslysis
### Reading the data and checking the samples from the full results in directory
## results2
taxa.df = read.table(file="results2/results.0.03.cons.taxonomy",header=TRUE)
data.df = read.table(file="results2/results.shared",header=TRUE)
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

### Reading the treatment file for a simple experiment (two factors)
trt.df = read.table(file="metadata2/samples.file",header=FALSE)
## we want to extract the first three columns: 1) treatment 1, 2) treatment 2 and 3)
## replicates per treatment 1 x treatment 2 levels (open file samples.file in results2 and look at it)
## we use sample name as row name
trt.df = data.frame(trt.df[,1:3],row.names=trt.df[,4])
## add a header to the new data frame
names(trt.df) = c("Trt1","Trt2","rep")
## look at it
trt.df

### Making sure that the treatments match the data
data.df = data.df[row.names(trt.df),]

## level = (1=Domain, 2=Phylum,3=Class, 4=Order, 5=Family, 6=Genus, 7=otu/species)
level = 5
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=level) # Phulum
d.df = d.taxa.ls[[1]]

### CSS normalization (using this all way through)
norm.vt = cumNormStatFast(as.matrix(t(d.df)))
norm.mt = round(t(cumNormMat(as.matrix(t(d.df)),norm.vt)))

### Looking at an experiment
## bar plot
bar.ftn(norm.mt,trt.df)
## identify number of clusters
NbClust(data=norm.mt,distance="euclidean",method="ward.D",index="silhouette",max.nc = 8)
# use number of clusters identified as cutoff for the dendrogram
pheatmap.2.ftn(norm.mt,trt.df,cutoff=2,cellh=5,cellw=10,fonts=7) # Using the normalized counts
## ordination per each treatment 1 (trt.df[,1])
ord.plot.ftn(norm.mt,trt.df[,1],x=c(-0.5,0.5),y=c(-.3,.3),type=6,cex=0.5)
## ordination per each treatment 2 (trt.df[,2])
ord.plot.ftn(norm.mt,trt.df[,2],x=c(-0.5,0.5),y=c(-.3,.3),type=6,cex=0.5)
### ordination for both treatments 
## first we combine the treatments in trt.cm
trt.cm = factor(paste(trt.df[,1],trt.df[,2],sep=""))
## then we plot
ord.plot.ftn(norm.mt,trt.cm,x=c(-.6,.6),y=c(-.6,.6),type=6,cex=0.5)
## using phyloseq to ordinate and plot
d.otu = phyloseq(otu_table(norm.mt,taxa_are_rows = FALSE),sample_data(trt.df))
d.ord = ordinate(d.otu, "NMDS", "bray")
theme_set(theme_bw())
plot_ordination(d.otu,d.ord,"samples",color="Trt1",shape = "Trt2")+geom_polygon(aes(alpha=0.5,fill=Trt1))

####### OTU Analysis using diversity measures
### (7 = otu/species) otu analysis
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=7)
d.df = d.taxa.ls[[1]]
d.rs = apply(d.df,1,sum)

### Diversity analysis (old fashioned ecological approach [richness and diversity])
## Simple analysis
d.otu = otu_table(d.df,taxa_are_rows = FALSE)
d.rich = estimate_richness(d.otu)
rich.mt = t(rarefy(d.df,min(d.rs),se=TRUE))
## this is like above, except that we add treatments 1 and 2 to the data frame
## instead of only treatment 1
rich.df = data.frame(trt.df[,1:2],d.rich[,c(1,6,8)],rich.mt[,1])
names(rich.df) = c("trt1","trt2","Observed","Shannon","InvSimpson","Richness")
## look at it
rich.df

## plot
rich.plot.ftn(d.df=rich.df)

## simple anova but now we have two treatments
## * indicates the full model y = treatment 1 + treatment 2 + the interaction
rich.lm = lm(Shannon~trt1*trt2,data=rich.df)
anova(rich.lm)

#### redundancy analysis like above
## level = (1=Domain, 2=Phylum,3=Class, 4=Order, 5=Family, 6=Genus, 7=otu/species)
level = 5
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa.df,d.df=data.df,l=level) # Phulum
d.df = d.taxa.ls[[1]]

### CSS normalization (using this all way through)
norm.vt = cumNormStatFast(as.matrix(t(d.df)))
norm.mt = round(t(cumNormMat(as.matrix(t(d.df)),norm.vt)))

mod1 <- rda(norm.mt~ trt.df[,1]*trt.df[,2],distance="bray")
mod0 <- rda(norm.mt ~ 1, distance="bray")
#Backward elimination
modb <- step(mod1, scope = formula(mod0), test = "perm")

## using best model to plot (here I use treatement 1 to demonstrate; best model
## is the empty model)
trt = trt.df[,1]
c=1:length(levels(trt))
rda<- rda(norm.mt ~ trt, distance= "bray")
## plot
plot(rda,type="n")
text(rda,display="site",cex=0.5) 
text(rda,display="species",cex=0.5) 
ordiellipse(rda, trt,col=c,lwd=1,label=TRUE)
ordispider(rda, trt,col=c,label = FALSE)



