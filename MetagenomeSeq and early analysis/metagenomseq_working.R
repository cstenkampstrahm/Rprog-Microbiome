library("xlsx")
library("lme4")
library("nnet") # multinomial modeling
library("gee") # if using generalized estimating equations
library("tidyverse")
library(phyloseq)
load("Phyloseq files/Cowonly")
Cowonly
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
MRexp_cowonly

d.gen = aggTax(MRexp_cowonly,lvl="Rank6", norm = TRUE) 
d.spp = aggTax(MRexp_cowonly,lvl="Rank7", norm = TRUE) 
d.fam = aggTax(MRexp_cowonly,lvl="Rank5", norm = TRUE)
d.class = aggTax(MRexp_cowonly, lvl="Rank3", norm = TRUE)
d.ord = aggTax(MRexp_cowonly, lvl="Rank4", norm = TRUE)
t.df = pData(MRexp_cowonly)
dim(d.gen) # features = 590
dim(d.spp) # features = 195
dim(d.fam) # features = 290
dim(d.class) # features = 111
dim(d.ord) # features = 191

d.genfilt = filterData(d.gen, present = 50, depth = 1)
dim(d.genfilt) # 79 features in a total of 50 cows
d.genfiltnorm = cumNorm(d.genfilt, p = cumNormStatFast(d.genfilt))
modgen = model.matrix(~ Pathotype_1 + Individual_animal - 1, t.df) 
d.fitgen = fitFeatureModel(d.genfiltnorm, modgen)
View(d.fitgen)
head(d.fitgen$fit)
head(MRcoefs(d.fitgen), counts = 1)
head(MRtable(d.fitgen), counts = 1)
head(MRcoefs(d.fitgen))
head(MRtable(d.fitgen))
View(MRcoefs(d.fitgen))
View(MRtable(d.fitgen))

d.famfilt = filterData(d.fam, present = 50, depth = 1)
dim(d.famfilt) # 68 features in a total of 50 cows
d.famfiltnorm = cumNorm(d.famfilt, p = cumNormStatFast(d.famfilt))
modfam = model.matrix(~ Pathotype_1 + Individual_animal - 1, t.df) 
modfamnew= model.matrix (~ Pathotype_1, t.df) # says incorrect number of dimensions when -1 added, 
# can't analyze currently when individual animal added, for fitfeature model. runs with fitzig
# but get zero counts with -1 or individual animal
d.fitfam = fitFeatureModel(d.famfiltnorm, modfam)
d.fitfamnew = fitFeatureModel(d.famfiltnorm, modfamnew)
d.fitfamzig <- fitZig(d.famfiltnorm, modfam)
d.fitfamnewzig <- fitZig(d.famfiltnorm, modfamnew)

head(MRcoefs(d.fitfam))
head(MRtable(d.fitfam))
View(MRcoefs(d.fitfam))
View(MRtable(d.fitfam))
head(MRcoefs(d.fitfamnew))
head(MRtable(d.fitfamnew))
View(MRcoefs(d.fitfamnew))
View(MRtable(d.fitfamnew))
head(MRcoefs(d.fitfamzig))
head(MRtable(d.fitfamzig))
View(MRcoefs(d.fitfamzig))
View(MRtable(d.fitfamzig))
head(MRcoefs(d.fitfamnewzig))
head(MRtable(d.fitfamnewzig))
View(MRcoefs(d.fitfamnewzig))
View(MRtable(d.fitfamnewzig))

modother <- model.matrix(~Pathotype_1+Individual_animal, t.df)
famother <- fitZig(d.famfiltnorm, modother)
famother1 <- fitFeatureModel(d.famfiltnorm, modother) # cant analyze currently
View(MRtable(famother))

Pathotype_new <- factor(t.df$Pathotype_1) #makes it a factor
trymod <- model.matrix(~ Pathotype_new + Individual_animal - 1, t.df) # makes actual contrasts for pathotype
d.fitfamtry = fitFeatureModel(d.famfiltnorm, trymod) #cant analyze currently (with or without -1!)
tryfit = d.fitfamtry$fit
finalMod = d.fitfam$design
contrast.matrix = makeContrasts(Pathotype_1 - Pathotype_1, levels = finalMod)
fit2 = contrasts.fit(tryfit, contrast.matrix)
model.matrix(~Pathotype_1, t.df)

d.sppfilt = filterData(d.spp, present = 50, depth = 1)
dim(d.sppfilt) # only 17 species to compare
d.sppfiltnorm = cumNorm(d.sppfilt, p = cumNormStatFast(d.sppfilt))
modspp = model.matrix(~ Pathotype_1 + Individual_animal - 1, t.df) 
d.fitspp = fitFeatureModel(d.sppfiltnorm, modspp)
head(MRcoefs(d.fitspp))
head(MRtable(d.fitspp))
View(MRcoefs(d.fitspp))
View(MRtable(d.fitspp))

# Need to parse down cowonly data set to the samples that are DP vs NS for 
# comparisons. Went through data set and wrote samples that want to keep
# only 13 total samples that were day prior to shedding event. One cow that had two
# samps prior to shedding were available, 8b and 8d. Kept one for subset data set. 
# For subset data set picked one sample each from the day matching the DP samps 
# of the available NS samps. Have a subset data set and a total one, with the total
# one can put individual animal in the model

DPNSsubset <- c("18D", "1C", "30C", "45C", "48A", "50D", "58B", "63D", "64D", "6C", "74B", "8B",
               "14D", "15C", "20C", "24C", "29A", "34D", "35B", "39D", "42D", "44C", "52B", "53B", "55A", 
               "60B", "61C", "62D", "65E", "66A", "68B", "69C", "80D", "73E", "81A", "7B")
DPNSsamps <- c("18D", "1C", "30C", "45C", "48A", "50D", "58B", "63D", "64D", "6C", "74B", 
               "8B", "14A", "14B", "14C", "14D", "14E", "15A", "15B", "15C", "15D", "15E", "20A", "20B", "20C", "20D", 
               "20E", "24A", "24B", "24C", "24D", "24E", "29A", "29B", "29C", "29D", "29E", "34A", "34B", "34C", "34D", 
               "34E", "35A", "35B", "35C", "35D", "35E", "39A", "39B", "39C", "39D", "39E", "42A", "42B", "42C", "42D", 
               "42E", "44A", "44B", "44C", "44D", "44E", "52A", "52B", "52C", "52D", "52E", "53A", "53B", "53C", 
               "53D", "53E", "55A", "55B", "55C", "55D", "55E", "60B", "60C", "60D", "60E", "61A", "61B", "61C", 
               "61D", "61E", "62A", "62B", "62C", "62D", "62E", "65A", "65B", "65C", "65D", "65E", "66A", "66B", 
               "66C", "66D", "66E", "68A", "68B", "68C", "68D", "68E", "69A", "69B", "69C", "69D", "69E", "70A", 
               "70B", "70C", "70D", "70E", "71A", "71B", "71C", "71D", "73A", "73B", "73C", "73D", "73E", 
               "7B", "7C", "7D", "7E")

DPNSsubset <- prune_samples(DPNSsubset, Cowonly)
DPNSsamps <- prune_samples(DPNSsamps, Cowonly)
MRexp_DPNSsubset <- phyloseq_to_metagenomeSeq(DPNSsubset)
MRexp_DPNSsamps <- phyloseq_to_metagenomeSeq(DPNSsamps)










###### TRYING OTHER THINGS WITH fitFEATUREMODEL ######
# try farm
modgenfarm = model.matrix(~ Farm + Individual_animal - 1, t.df) 
d.fitgenfarm = fitFeatureModel(d.genfiltnorm, modgenfarm)
View(MRcoefs(d.fitgenfarm))
View(MRtable(d.fitgenfarm))
# no change by farm at the genus level
# re-ran the above with filtering to 100 cows or to 150 cows. Same genus come
# up in the table, in some cases their p values are slightly higher though (ones
# that are 0.09 are now 0.1 for instance)


# try to mutate the variable to get the different heading for pathotype in the output
library(tidyverse)
t.df <- mutate(t.df, Pathotype_name = factor(Pathotype_1, levels = c(0, 1),
                                             labels = c("None", "O157")))
t.df <- mutate(t.df, Pathotype_name = as.integer(Pathotype_name))
mod = model.matrix(~ Pathotype_name + Individual_animal - 1, t.df) 
d.fit = fitFeatureModel(d.genfiltnorm, mod)
# how come when I run this code I get different output for 1 and 2 than I got for 
# 0 and 1 when I simply ran pathotype?


### trying other ways to specify the variables like from 2013 metagenomeseq annotation
# for fitzig model
d.genfilt = filterData(d.gen, present = 50, depth = 1)
dim(d.genfilt) # 79 features in a total of 50 cows
d.genfiltnorm = cumNorm(d.genfilt, p = cumNormStatFast(d.genfilt))
pathotype <- pData(MRexp_cowonly)$Pathotype_1
individual <- pData(MRexp_cowonly)$Individual_animal
modgen = model.matrix(~ pathotype + individual -1) 
d.fitgen = fitFeatureModel(d.genfiltnorm, modgen)
d.fitgen
head(MRcoefs(d.fitgen))
head(MRtable(d.fitgen))
View(MRcoefs(d.fitgen))
View(MRtable(d.fitgen))

### making a heatmap  

trials <- pData(d.genfilt)$Pathotype_1
trials <- sort(trials, decreasing = TRUE)

pData <-pData(d.genfilt)
pData <- pData[order(pData[,33],decreasing=T),]
trials <- pData$Pathotype_1
heatmapColColors=brewer.pal(12,"Set3")[trials]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj=d.genfilt, n = 10, fun = mad, cexRow = 0.4, cexCol = 0.4,
              trace = "none", col = heatmapCols, ColSideColors = heatmapColColors,
              Colv=TRUE, dendrogram = "row")


save(MRexp_cowonly, file="cowexperiment")
