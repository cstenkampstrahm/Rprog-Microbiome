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
load("Phyloseq files/Cowonly")
Cowonly
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
MRexp_cowonly

fData(MRexp_cowonly)
pData(MRexp_cowonly)
head(MRcounts(MRexp_cowonly))
colSums(MRcounts(MRexp_cowonly))
# met with Zaid on 4.21.17. For modeling best to use fitFeatureModel (as this is
# zero inflated log normal model of the count data). Mine is currently not working,
# and the Gaussian distribution of the method using FitZig is not ideal, because our
# data is not normal and does not go below zero. In theory if we log transformed our data
# we would have the correct distribution, but log(0) = 1 and log(1)=0. So log transforming
# and running FitZig does not seem to be a doable option. Best idea to figure out how
# to get the fitFeatureModel to work and run the OTU data. The CSS offset should always
# be a model factor, so don't put it as FALSE (we are always going to use CSS scaling
# before analyzing the samples). Probably best to start with the aggregation at the 
# lowest level possible (spp or genus) and go from there. In models should always specify
# in the model.matrix statement "-1" because it gets rid of the equation general mean.
# this way we are looking at the deviation of the mean of each cell that meets our 
# category criteria, versus a total mean for the whole population. For the Pathotype, we 
# want to run a model as Pathotype_1 + Individual_animal - 1 to take into account the individual
# animal and also get rid of the population mean. For Day prior versus day of shedding event
# we will extract the rows of the appropriate cows and make two separate OTU matrices (matched up. We will
# substract the Day prior matrix - Day of matrix and then do a nonparametric Wilcoxon Rank Sum test
# on the output matrix asking whether or not the values are different from zero. Will use FDRtools to 
# correct this output for multiple testing.  For day prior vs never shedding cows, will go through the
# data and select one 'day prior' sample from all available cows with a day prior option. Then of non-shedding 
# cows will select a random day's sample (any of the days) to use in the analysis (so we only ever have 
# a single value from any one cow). Then in the fitFeatureModel will run DPvNS - 1 and look at the model 
# output for significant OTU changes between these categories. 
d.ms = cumNorm(MRexp_cowonly, p = cumNormStatFast(MRexp_cowonly))
nf = normFactors(d.ms)
t.df = pData(d.ms)
# now need to aggregate the taxa
# not sure why norm = TRUE here, because already normed. Following Zaid's code from class:
d.gen = aggTax(d.ms,lvl="Rank6",norm = TRUE) # get 590 to the genus
d.spp = aggTax(d.ms,lvl="Rank7", norm = TRUE) # only get 195 to the spp
d.fam = aggTax(d.ms,lvl="Rank5", norm = TRUE) # get 290 that agg to the family
#d.fam = cumNorm(d.fam,p=1)
# filter OTUs to present in at least 150 (75%) cow samples, depth of 1:
d.genfilt = filterData(d.gen, present = 100, depth = 1)
# went from 159,000 to 51 features at genus level
d.genfiltnorm = cumNorm(d.genfilt, p = cumNormStatFast(d.genfilt)) # need to normalize again per the
# fitFeatureModel below
# specify model to look at differential abundance of things based on pathotype

mod = model.matrix(~ Pathotype_1 + Individual_animal - 1, t.df) 
d.fit = fitFeatureModel(d.genfiltnorm, mod)
head(MRcoefs(d.fit))
head(MRtable(d.fit))
View(MRcoefs(d.fit))
View(MRtable(d.fit))

# try to make levels of pathotype 1 so can contrast them. isn't working!:
library(tidyverse)
t.df <- mutate(t.df, Pathotype_name = factor(Pattern_1, levels = c(0, 1),
                                             labels = c("None", "O157")))
mod = model.matrix(~ Pathotype_name + Individual_animal - 1, t.df) 
d.fit = fitFeatureModel(d.genfiltnorm, mod)
# says Error in fitFeatureModel(d.genfiltnorm, mod) : Can't analyze currently.
#In addition: Warning message:
  #In cbind(mod, log(normFactors(obj)/median(normFactors(obj)))) :
  #number of rows of result is not a multiple of vector length (arg 2)
#head(MRcoefs(d.fit))
#Fit1 = d.fit$fit
#finalMod = d.fit$fit$design
#contrast.matrix = makeContrasts(0 - 1, levels = finalMod)
#fit2 = contrasts.fit(zigFit, contrast.matrix)

library(tidyverse)
t.df <- mutate(t.df, Pattern_name = factor(Pattern_1, levels = c(0, 1, 2),
                                           labels = c("No_shed", "Intermitt", "Multi")))

mod1.0 = model.matrix(~ Pattern_name, t.df)
colnames(mod1.0) = levels(t.df$Pattern_name) # if I add in animal, I can't name the 
#columns here because they don't match
d.fit1.0 = fitZig(d.genfiltnorm, mod1.0)
head(MRcoefs(d.fit1.0))
head(MRtable(d.fit1.0))
#View(MRcoefs(d.fit1.0))
View(MRtable(d.fit1.0))
zigFit = d.fit1.0$fit
finalMod = d.fit1.0$fit$design
contrast.matrix = makeContrasts(Multi - No_shed, Intermitt - No_shed, levels = finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 # gives you all sorts of info, p vals assoc for each contrast, stddev, sigmas, residuals, coeffs
fit2 = eBayes(fit2) # looking at the error of the residuals here
topTable(fit2)

# try to export the taxa names
taxa = sapply(strsplit(as.character(fData(MRexp_cowonly)$taxa), split = ";"),
              function(i) {
                i[length(i)]
          })  ## this code is not working....
head(MRcoefs(d.fit1.0, taxa=taxa))

              



                 
