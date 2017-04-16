### Prepared by Zaid Abdo/Metagenomics Class-2017
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq','metagenomeSeq')
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
source(file="functions.R")

##################################### Two factors anslysis
### Reading the data and checking the samples from the full results in directory
## results2
taxa.df = read.table(file="results2/results.0.03.cons.taxonomy",header=TRUE,as.is=TRUE)
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

### Making sure that the taxonomic table matches the data table
otu.nm = names(data.df)
taxa.df = taxa.df[taxa.df$OTU%in%otu.nm,]

### Making sure the taxonomic table matches the required format for metagenomSeq
taxa.df = taxa.ftn(taxa.df)

### Reading the treatment file for a simple experiment (two factors)
## and making sure it matches the metagenome seq format
trt.df = read.table(file="results2/samples.file",header=FALSE)
r.nm = trt.df[,4]
trt.df = data.frame(trt.df[,2:3],row.names=r.nm)
names(trt.df) = c("Trt1","Trt2")

### Getting the data in a metagenomeSeq friendly format
data.df = data.df[r.nm,]
data.df = t(data.df)

### Creating the MRexperiment object
trt.an = AnnotatedDataFrame(trt.df)
taxa.an = AnnotatedDataFrame(taxa.df)
d.ms = newMRexperiment(data.df,phenoData = trt.an,featureData = taxa.an)
pData(d.ms) #phenotype data that describes the metadata
fData(d.ms) #this is the otudata 
head(MRcounts(d.ms))#the count data

### Choosing the OTUs that are found in at least half the samples
# d.ms = filterData(d.ms, present = 18, depth = 1) # present in 18 samples, at least one 'count'
d.ms = filterData(d.ms, present = 9, depth = 1) # present in 9 samples, at least one 'count'
### CSS normalization (using this all way through)
p = cumNormStatFast(d.ms)
d.ms = cumNorm(d.ms,p)
nf = normFactors(d.ms)

#### metagenomeSeq analysis, looking for differental OTU abundance, testing a linear model fit by 
## treatment of each OTU one at a time. Because there are many OTUs, there are corrections for mult
## testing. in MetagenomeSeq they use FDR (false discovery rate)- less conservative (assumes some are
## and some are not significant, and finds somewhere in between- bonferroni and shafayes (?) and family-wise
## and others assume that there is no significance and must prove otehrwise!)
t.df = pData(d.ms)
mod = model.matrix(~Trt1-1, t.df) # minus one is taking out the mean, is going to contrast the variable
# you specify with all of the other ones in the data frame that you have
# mod = model.matrix(~Trt1+Trt2-1, t.df) # another model
d.fit = fitZig(d.ms, mod)
head(MRcoefs(d.fit))

#### metagenomeSeq analysis on family level
d.fam = aggTax(d.ms,lvl="Family",norm = TRUE)
d.fam = cumNorm(d.fam,p=1)
mod = model.matrix(~Trt1-1, t.df)
colnames(mod) = levels(t.df$Trt1)
d.fit = fitZig(d.fam, mod,useCSSoffset = FALSE)
head(MRcoefs(d.fit)) # these are the coefficients
head(MRtable(d.fit))
zigFit = d.fit$fit
finalMod = d.fit$fit$design
contrast.matrix = makeContrasts(Trt1 - Trt3, Trt1 - Trt2, levels = finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 # gives you all sorts of info, p vals assoc for each contrast, stddev, sigmas, residuals, coeffs
fit2 = eBayes(fit2) # looking at the error of the residuals here
topTable(fit2)
