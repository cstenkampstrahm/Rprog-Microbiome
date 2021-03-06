library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(metagenomeSeq)
library(cluster)
library(NbClust)
library(tidyverse)
source(file="files from zaid/functions_1.R")

### using non normalized phyloseq experiment object from phyloseq_workflow.R
OTUS <- otu_table(Cowdatarich)
sampledata <- sample_data(Cowdatarich)
taxa <- tax_table(Cowdatarich)

## need to coerce the OTU table into a data frame for vegan (??)
## found on https://github.com/joey711/phyloseq/issues/190
vegan_otu <- function(physeq) {
  OTU <- otu_table(Cowdatarich)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Cowdatarich)
# looks like the function made a matrix.....
# following Zaids code from class 2/17, data.df is the samps and OTU cts, taxa.df
# is taxa table, trt.df is metadata file. For me, OTUS is OTU cts, taxa is tax table
# and sampledata is metadata file.
taxanew <- data.frame(taxa)

### Remove otus that have nothing in them
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)
# made the OTU values a new column. Now try to select only the ones left in my 
# OTU table after pruning:
taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]

### Before continuing on with this analysis, need to generate the taxa.df similar
### to the version that zaid has, where there is numbered rows, a column of OTUs,
### their size in the next column, and then the taxa full name in the next column

# going to try to convert to an annotated data frame and back again to use agg tax
names <- sampledata[,1]
OTUSpruned <- OTUS[names$SampleID,] # this didn't get rid of any taxa but checks 
# to see that they are the same samp names between OTU table and sampledata
OTUSpruned = t(OTUSpruned) # what does this do?\
# looks like I need to get the sample names in the column position?
sampledata = t(sampledata)
sampledata = data.frame(sampledata)
trt.an = AnnotatedDataFrame(sampledata)
taxa.an = AnnotatedDataFrame(taxa.df)
d.ms = newMRexperiment(OTUSpruned,phenoData = trt.an,featureData = taxa.an) #NOT WORKING
pData(d.ms) #phenotype data that describes the metadata
fData(d.ms) #this is the otudata 
head(MRcounts(d.ms))#the count data



taxa1.df <- unite(taxa.df, col = bacteria, Rank1:Rank7, sep = " ;")
#changed the 1s to 2s for column spec, changed the 3 to 1 for column spec, also
# changed the strspl to "\\;" to match space and semi colon in 
#taxa split combine function in functions_1 script. Having trouble with 'temp'
#line, unable to decipher what this means exactly in the function....

####### OTU Analysis
### (6 = otu/species) otu analysis
d.taxa.ls = taxa.split.combine.ftn(taxa=taxa1.df,d.df=OTUS,l=level)
d.df = d.taxa.ls[[1]]
taxa.vt = d.taxa.ls[[2]]

### check distribution of depth of coverage
d.rs = apply(OTUS,1,sum)
hist(d.rs)

### did we sample these communities well? (Rarifaction curves)
rarecurve(x=d.df,col=1:length(d.rs)) 
## this is bogging down computer, just 
## ends up running forever. Will try to specify the first day samples
rarecurve(x=d.df[c(2,14,23,25,33,38,39,43,51,56,57,58,68,71,83,90,94,
                    95,106,114,118,122,129,135,140,141,143,153,156,158,
                    163,166,173,176,178,180,193,197),],col=1:38)
## Finally generated rarefaction curves, after -8 or -9 hours. Going to remove
## samples 18A and 7A that have an extremely low depth of sampling (197 and 193)
day1rarecurve <- rarecurve(x=d.df[c(2,14,23,25,33,38,39,43,51,56,57,58,68,71,83,90,94,
                   95,106,114,118,122,129,135,140,141,143,153,156,158,
                   163,166,173,176,178,180),],col=1:36)
### A look at proportions before normalizing
bar.taxa.sample.ftn(d.df,cutoff=.01)
### Try to make a bar taxa plot just with the day 1 cows, without the 18 & 7.....
day1bartaxa <- bar.taxa.sample.ftn(d.df[c(2,14,23,25,33,38,39,43,51,56,57,58,68,71,83,90,94,
                           95,106,114,118,122,129,135,140,141,143,153,156,158,
                           163,166,173,176,178,180),],cutoff=.01)
### What does it mean for a cutoff of 0.01? Tried to understand by looking at function,
### and still not sure. For the day 1 cows have 68 in the legend, for all cows have 119...
### Diversity analysis (old fashioned ecological approach [richness and diversity])
d.otu = otu_table(d.df,taxa_are_rows = FALSE)
#plot_richness(d.otu)
d.rich = estimate_richness(d.otu)
rich.mt = t(rarefy(d.df,min(d.rs),se=TRUE))
rich.df = data.frame(cbind(d.rich[,c(1,6,8)],rich.mt[,1]))
names(rich.df) = c("Observed","Shannon","InvSimpson","Richness")
df.plot.ftn(rich.df)

### now do the same diversity analysis with only the day 1 cows:
d.otu = otu_table(d.df[c(2,14,23,25,33,38,39,43,51,56,57,58,68,71,83,90,94,
                         95,106,114,118,122,129,135,140,141,143,153,156,158,
                         163,166,173,176,178,180),],taxa_are_rows = FALSE)
#plot_richness(d.otu)
d.rich = estimate_richness(d.otu)
rich.mt = t(rarefy(d.df[c(2,14,23,25,33,38,39,43,51,56,57,58,68,71,83,90,94,
                          95,106,114,118,122,129,135,140,141,143,153,156,158,
                          163,166,173,176,178,180),],min(d.rs),se=TRUE))
rich.df = data.frame(cbind(d.rich[,c(1,6,8)],rich.mt[,1]))
names(rich.df) = c("Observed","Shannon","InvSimpson","Richness")
df.plot.ftn(rich.df)

