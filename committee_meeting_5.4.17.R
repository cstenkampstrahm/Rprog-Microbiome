### Code used for differential abundance outputs for committee meeting on 5/4/17

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
Pathotype_new <- factor(t.df$Pathotype_1) #makes it a factor
Individual_animalnew <- factor(t.df$Individual_animal)
d.genfiltnorm = cumNorm(d.genfilt, p = cumNormStatFast(d.genfilt))
trymod <- model.matrix(~ Pathotype_new + Individual_animalnew - 1, t.df) # makes actual contrasts for pathotype
d.genfit = fitZig(d.genfiltnorm, trymod)
View(MRcoefs(d.genfit))
View(MRtable(d.genfit))
zigFitgen = d.genfit$fit
finalMod = d.genfit$fit$design
contrast.matrix = makeContrasts(Pathotype_new0-Pathotype_new1, levels = finalMod)
fit2 = contrasts.fit(zigFitgen, contrast.matrix)
fit2 # gives you all sorts of info, p vals assoc for each contrast, stddev, sigmas, residuals, coeffs
fit2 = eBayes(fit2) # looking at the error of the residuals here
topTable(fit2)
View(topTable(fit2))



d.famfilt = filterData(d.fam, present = 50, depth = 1)
Pathotype_new <- factor(t.df$Pathotype_1) #makes it a factor
Individual_animalnew <- factor(t.df$Individual_animal)
d.famfiltnorm = cumNorm(d.famfilt, p = cumNormStatFast(d.famfilt))
trymod <- model.matrix(~ Pathotype_new + Individual_animalnew - 1, t.df) # makes actual contrasts for pathotype
d.famfit = fitZig(d.famfiltnorm, trymod)
View(MRcoefs(d.famfit))
View(MRtable(d.famfit))
zigFitfam = d.famfit$fit
finalMod = d.famfit$fit$design
contrast.matrix = makeContrasts(Pathotype_new0-Pathotype_new1, levels = finalMod)
fit3 = contrasts.fit(zigFitfam, contrast.matrix)
fit3 # gives you all sorts of info, p vals assoc for each contrast, stddev, sigmas, residuals, coeffs
fit3 = eBayes(fit3) # looking at the error of the residuals here
topTable(fit3)
View(topTable(fit3))

d.sppfilt = filterData(d.spp, present = 50, depth = 1)
Pathotype_new <- factor(t.df$Pathotype_1) #makes it a factor
Individual_animalnew <- factor(t.df$Individual_animal)
d.sppfiltnorm = cumNorm(d.sppfilt, p = cumNormStatFast(d.sppfilt))
trymod <- model.matrix(~ Pathotype_new + Individual_animalnew - 1, t.df) # makes actual contrasts for pathotype
d.sppfit = fitZig(d.sppfiltnorm, trymod)
View(MRcoefs(d.sppfit))
View(MRtable(d.sppfit))
zigFitspp = d.sppfit$fit
finalMod = d.sppfit$fit$design
contrast.matrix = makeContrasts(Pathotype_new0-Pathotype_new1, levels = finalMod)
fit4 = contrasts.fit(zigFitspp, contrast.matrix)
fit4 # gives you all sorts of info, p vals assoc for each contrast, stddev, sigmas, residuals, coeffs
fit4 = eBayes(fit4) # looking at the error of the residuals here
topTable(fit4)
View(topTable(fit4))

