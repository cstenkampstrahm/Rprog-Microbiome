# MetagenomeSeq workflow: will use this page to build models of differential abundance 

# can generate the cowdatarich (nonnormed- although renamed simply Cowdata in files) phyloseq object
# and also the Normcowdatanew (normalized with CSS) phyloseq object
# additionally cowdata is phyloseq with everything (enviro and non) and cowonly and environly
# have only those samples, respectively (there is just different data in their metadata). saved
# files to load easily later
library("xlsx")
library("lme4")
library("nnet") # multinomial modeling
library("gee") # if using generalized estimating equations
library("tidyverse")
#save(Cowonly, file="Cowonly")
#save(Enviroonly, file="Enviroonly")
#save(Cowdatarich, file="Cowdata")
#save(Normcowonly, file="Normcowdatanew")
#source(file="Phyloseq_workflow.R")
library(phyloseq)
load("Phyloseq files/Cowdata")
Cowdatarich
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowdatarich)
MRexp_cowonly

fData(MRexp_cowonly)
pData(MRexp_cowonly)
head(MRcounts(MRexp_cowonly))
colSums(MRcounts(MRexp_cowonly))

# filter OTUs to present in at least 150 (75%) cow samples, depth of 1:
d.ms = filterData(MRexp_cowonly, present = 150, depth = 1)
# went from 159,000 to 815 features
# do CSS normalization
# p = cumNormStatFast(d.ms)
# d.ms = cumNorm(d.ms,p)
# got error warning sample with one or zero features
# need to remove sample that has hardly any counts; which sample? 7A
d.ms = cumNorm(d.ms,p = cumNormStatFast(d.ms))
nf = normFactors(d.ms)
t.df = pData(d.ms)
# specify model to look at differential abundance of things based on pathotype
mod = model.matrix(~ Pathotype_1, t.df)
d.fit = fitZig(d.ms, mod)
head(MRcoefs(d.fit))
head(MRtable(d.fit))
View(MRcoefs(d.fit))
View(MRtable(d.fit))

# need to control for animal though. how to set up the matrix to do this?

mod1.0 = model.matrix(~ Pathotype_1 + (1|Individual_animal), t.df)
d.fit1.0 = fitZig(d.ms, mod1.0)
head(MRcoefs(d.fit1.0))
head(MRtable(d.fit1.0))
View(MRcoefs(d.fit1.0))
View(MRtable(d.fit1.0))

mod1 = model.matrix(~ Pathotype_1 - 1, t.df)
d.fit1 = fitZig(d.ms, mod1)
head(MRcoefs(d.fit1))
head(MRtable(d.fit1))
View(MRcoefs(d.fit1))
View(MRtable(d.fit1))

# what does subtracting 1 do?

mod = model.matrix(~ EvNev_1, t.df)
d.fit = fitZig(d.ms, mod)
head(MRcoefs(d.fit))
head(MRtable(d.fit))
View(MRcoefs(d.fit))
View(MRtable(d.fit))

test.mod <- lmer(~ Pathotype_1 + (1|Individual_animal), 
                  data = t.df) 
                 
