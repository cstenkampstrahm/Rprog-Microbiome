### placing all used code to measure differential abundance here (taken from 
## other files during working analysis). Code includes metagenomeSeq and 
## DeSeq2 models for Pathotype variable, and Wilcoxon rank testing of DPvsDO and
## DOvsDA variables


#### FROM final_metagenomseq.R
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
t.df = pData(MRexp_cowonly)

dim(d.gen) # features = 590
dim(d.spp) # features = 195
dim(d.fam) # features = 290

d.otufilt = filterData(MRexp_cowonly, present = 50, depth = 1)
dim(d.otufilt) # features = 3531
Pathotype_new <- factor(t.df$Pathotype_1) #makes it a factor
Individual_animalnew <- factor(t.df$Individual_animal)
d.otufiltnorm = cumNorm(d.otufilt, p = cumNormStatFast(d.otufilt))
trymod <- model.matrix(~ Pathotype_new + Individual_animalnew - 1, t.df) # makes actual contrasts for pathotype
d.otufit = fitZig(d.otufiltnorm, trymod)
View(MRcoefs(d.otufit))
View(MRtable(d.otufit))
zigFitotu = d.otufit$fit
finalMod = d.otufit$fit$design
contrast.matrix = makeContrasts(Pathotype_new0-Pathotype_new1, levels = finalMod)
fit5 = contrasts.fit(zigFitotu, contrast.matrix)
fit5 # gives you all sorts of info, p vals assoc for each contrast, stddev, sigmas, residuals, coeffs
fit5 = eBayes(fit5) # looking at the error of the residuals here
topTable(fit5)
View(topTable(fit5))

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
topTable(fit2, confint = TRUE)
View(topTable(fit2, confint = TRUE))



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
topTable(fit3, confint = TRUE)
View(topTable(fit3, confint = TRUE))

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
topTable(fit4, confint = TRUE)
View(topTable(fit4, confint = TRUE))


##### FROM final_deseq2.R
library(DESeq2)
library(phyloseq)
library(metagenomeSeq)
load("Phyloseq files/Cowonly")
Cowonly
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
MRexp_cowonly

## TO REMAKE THESE COUNT TABLES ##
d.gen1 = aggTax(MRexp_cowonly,lvl="Rank6", norm = FALSE)
d.genfilt1 = filterData(d.gen1, present = 50, depth = 1)
genus1 <- MRcounts(d.genfilt1)
dim(genus1) # 79 (same as filtered w metagenomeseq fitzig)

d.spp1 = aggTax(MRexp_cowonly,lvl="Rank7", norm = FALSE) 
d.sppfilt1 = filterData(d.spp1, present = 50, depth = 1)
spp1 <- MRcounts(d.sppfilt1)
dim(spp1) # 17 (same as filtered w metagenomeseq fitzig)


d.fam1 = aggTax(MRexp_cowonly,lvl="Rank5", norm = FALSE)
d.famfilt1 = filterData(d.fam1, present = 50, depth = 1)
family1 <- MRcounts(d.famfilt1)
dim(family1) # 68 (same as filtered w metagenomeseq fitzig)

d.otu1 = filterData(MRexp_cowonly, present = 50, depth = 1)
otu1 <- MRcounts(d.otu1)
dim(otu1) # 3531 (same as filtered w metagenomeseq fitzig)


#McMurdie's calculation of the geometric mean 

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## GENUS FILTERED
genusfilt <- otu_table(genus1, taxa_are_rows=TRUE)
genusfilt <- merge_phyloseq(genusfilt, sampledata)
sample_data(genusfilt)$Pathotype_1 <- factor(sample_data(genusfilt)$Pathotype_1)
sample_data(genusfilt)$Individual_animal <- factor(sample_data(genusfilt)$Individual_animal)
sample_data(genusfilt)$Pathotype_1 <- relevel(sample_data(genusfilt)$Pathotype_1, "0")
genusfiltDESEQ <- phyloseq_to_deseq2(genusfilt, ~Individual_animal + Pathotype_1 - 1)
geoMeans = apply(counts(genusfiltDESEQ), 1, gm_mean)
genusfiltDESEQg = estimateSizeFactors(genusfiltDESEQ, geoMeans = geoMeans)
genusfiltDESEQg = DESeq(genusfiltDESEQg, fitType="local", betaPrior = FALSE)
genusfiltresultsG <- results(genusfiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res5.1 = genusfiltresultsG
res5.1 = res5.1[order(res5.1$padj, na.last=NA), ]
alpha = 0.2
sigtab5.1 = res5.1[(res5.1$padj < alpha), ]
sigtab6.1 = res5.1[(res5.1$pvalue < alpha), ]
sigtab5.1
sigtab6.1

## SPECIES FILTERED
sppfilt <- otu_table(spp1, taxa_are_rows=TRUE)
sppfilt <- merge_phyloseq(sppfilt, sampledata)
sample_data(sppfilt)$Pathotype_1 <- factor(sample_data(sppfilt)$Pathotype_1)
sample_data(sppfilt)$Individual_animal <- factor(sample_data(sppfilt)$Individual_animal)
sample_data(sppfilt)$Pathotype_1 <- relevel(sample_data(sppfilt)$Pathotype_1, "0")
sppfiltDESEQ <- phyloseq_to_deseq2(sppfilt, ~Individual_animal + Pathotype_1 - 1)
geoMeans = apply(counts(sppfiltDESEQ), 1, gm_mean)
sppfiltDESEQg = estimateSizeFactors(sppfiltDESEQ, geoMeans = geoMeans)
sppfiltDESEQg = DESeq(sppfiltDESEQg, fitType="local", betaPrior = FALSE)
sppfiltresultsG <- results(sppfiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res8.1 = sppfiltresultsG
res8.1 = res8.1[order(res8.1$padj, na.last=NA), ]
alpha = 0.2
sigtab11.1 = res8.1[(res8.1$padj < alpha), ]
sigtab12.1 = res8.1[(res8.1$pvalue < alpha), ]
sigtab11.1
sigtab12.1

## FAMILY FILTERED
famfilt <- otu_table(family1, taxa_are_rows=TRUE)
famfilt <- merge_phyloseq(famfilt, sampledata)
sample_data(famfilt)$Pathotype_1 <- factor(sample_data(famfilt)$Pathotype_1)
sample_data(famfilt)$Individual_animal <- factor(sample_data(famfilt)$Individual_animal)
sample_data(famfilt)$Pathotype_1 <- relevel(sample_data(famfilt)$Pathotype_1, "0")
famfiltDESEQ <- phyloseq_to_deseq2(famfilt, ~Individual_animal + Pathotype_1 - 1)
geoMeans = apply(counts(famfiltDESEQ), 1, gm_mean)
famfiltDESEQg = estimateSizeFactors(famfiltDESEQ, geoMeans = geoMeans)
famfiltDESEQg = DESeq(famfiltDESEQg, fitType="local", betaPrior = FALSE)
famfiltresultsG <- results(famfiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res11.1 = famfiltresultsG
res11.1 = res11.1[order(res11.1$padj, na.last=NA), ]
alpha = 0.2
sigtab17.1 = res11.1[(res11.1$padj < alpha), ]
sigtab18.1 = res11.1[(res11.1$pvalue < alpha), ]
sigtab17.1
sigtab18.1

##### FROM final_DPvDO and DOvsDA.R
library(tidyverse)
library(phyloseq)
load("Phyloseq files/Cowonly")
Cowonly
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
MRexp_cowonly

# aggregate taxa and normalize in metagenomeseq
d.gen = aggTax(MRexp_cowonly,lvl="Rank6", norm = TRUE) 
d.spp = aggTax(MRexp_cowonly,lvl="Rank7", norm = TRUE) 
d.fam = aggTax(MRexp_cowonly,lvl="Rank5", norm = TRUE)

t.df = pData(MRexp_cowonly)

# convert metagenomeseq to phyloseq objects
genus1 <- MRcounts(d.gen)
species1 <- MRcounts(d.spp)
family1 <- MRcounts(d.fam)

genusphylo <- otu_table(genus1, taxa_are_rows=TRUE)
genusphylo <- merge_phyloseq(genusphylo, sampledata)

speciesphylo <- otu_table(species1, taxa_are_rows=TRUE)
speciesphylo <- merge_phyloseq(speciesphylo, sampledata)

familyphylo <- otu_table(family1, taxa_are_rows=TRUE)
familyphylo <- merge_phyloseq(familyphylo, sampledata)

# prune phyloseq experiments to have DPDO and DODA data
DPDO <- c("18D", "18E", "1C", "1D", "30C", "30D", "45C", "45D", "48A",
          "48B", "50D", "50E", "58B", "58C", "63D", "63E", "64D", "64E",
          "6C", "6D", "74B", "74C", "8B", "8C")
DODA <- c("10C", "10D", "30D", "30E", "38A", "38B", "45D", "45E", "57B",
          "57C", "58C", "58D", "67A", "67B", "6D", "6E", "74C", "74D", 
          "8C", "8D")

DPDOfamphylo <- prune_samples(DPDO, familyphylo)
DPDOsppphylo <- prune_samples(DPDO, speciesphylo)
DPDOgenphylo <- prune_samples(DPDO, genusphylo)

DODAfamphylo <- prune_samples(DODA, familyphylo)
DODAsppphylo <- prune_samples(DODA, speciesphylo)
DODAgenphylo <- prune_samples(DODA, genusphylo)

# convert to metagenomeseq de novo and use metagenomeseq to filter
# for some reason at this step the CSS normalization goes away?
MR_DPDOfamphylo <- phyloseq_to_metagenomeSeq(DPDOfamphylo)
MR_DPDOsppphylo <- phyloseq_to_metagenomeSeq(DPDOsppphylo)
MR_DPDOgenphylo <- phyloseq_to_metagenomeSeq(DPDOgenphylo)

MR_DODAfamphylo <- phyloseq_to_metagenomeSeq(DODAfamphylo)
MR_DODAgenphylo <- phyloseq_to_metagenomeSeq(DODAgenphylo)
MR_DODAsppphylo <- phyloseq_to_metagenomeSeq(DODAsppphylo)

MR_DPDOfamphylofilt = filterData(MR_DPDOfamphylo, present = 8, depth = 1)
MR_DPDOsppphylofilt = filterData(MR_DPDOsppphylo, present = 8, depth = 1)
MR_DPDOgenphylofilt = filterData(MR_DPDOgenphylo, present = 8, depth = 1)

MR_DODAfamphylofilt = filterData(MR_DODAfamphylo, present = 6, depth = 1)
MR_DODAsppphylofilt = filterData(MR_DODAsppphylo, present = 6, depth = 1)
MR_DODAgenphylofilt = filterData(MR_DODAgenphylo, present = 6, depth = 1)

# convert back to phyloseq
# put norm = true to CSS again...

DPDOgenus1 <- MRcounts(MR_DPDOgenphylofilt, norm = TRUE)
DPDOfamily1 <- MRcounts(MR_DPDOfamphylofilt, norm = TRUE)
DPDOspecies1 <- MRcounts(MR_DPDOsppphylofilt, norm = TRUE) #some seem nonnormed?

DODAgenus1 <- MRcounts(MR_DODAgenphylofilt, norm = TRUE)
DODAfamily1 <- MRcounts(MR_DODAfamphylofilt, norm = TRUE)
DODAspecies1 <- MRcounts(MR_DODAsppphylofilt, norm = TRUE) #some seem nonnormed?


DPDOgenusphylo <- otu_table(DPDOgenus1, taxa_are_rows=TRUE)
DPDOgenusfilt <- merge_phyloseq(DPDOgenusphylo, sampledata)

DPDOfamilyphylo <- otu_table(DPDOfamily1, taxa_are_rows=TRUE)
DPDOfamilyfilt <- merge_phyloseq(DPDOfamilyphylo, sampledata)

DPDOspeciesphylo <- otu_table(DPDOspecies1, taxa_are_rows=TRUE)
DPDOspeciesfilt <- merge_phyloseq(DPDOspeciesphylo, sampledata)

DODAgenusphylo <- otu_table(DODAgenus1, taxa_are_rows=TRUE)
DODAgenusfilt <- merge_phyloseq(DODAgenusphylo, sampledata)

DODAfamilyphylo <- otu_table(DODAfamily1, taxa_are_rows=TRUE)
DODAfamilyfilt <- merge_phyloseq(DODAfamilyphylo, sampledata)

DODAspeciesphylo <- otu_table(DODAspecies1, taxa_are_rows=TRUE)
DODAspeciesfilt <- merge_phyloseq(DODAspeciesphylo, sampledata)

# prune samples to DP1, DO2, DO3, DA4 in phyloseq
DP1 <- c("18D", "1C", "30C", "45C", "48A", "50D", "58B", "63D", "64D",
         "6C", "74B", "8B")
DO2 <- c("18E", "1D", "30D", "45D", "48B", "50E", "58C", "63E", "64E", 
         "6D", "74C", "8C")
DO3 <- c("10C", "30D", "38A", "45D", "57B", "58C", "67A", "6D", "74C", "8C")
DA4 <- c("10D", "30E", "38B", "45E", "57C", "58D", "67B", "6E", "74D", "8D")

DP1famphylo <- prune_samples(DP1, DPDOfamilyfilt)
DP1sppphylo <- prune_samples(DP1, DPDOspeciesfilt)
DP1genphylo <- prune_samples(DP1, DPDOgenusfilt)

DO2famphylo <- prune_samples(DO2, DPDOfamilyfilt)
DO2sppphylo <- prune_samples(DO2, DPDOspeciesfilt)
DO2genphylo <- prune_samples(DO2, DPDOgenusfilt)

DO3famphylo <- prune_samples(DO3, DODAfamilyfilt)
DO3sppphylo <- prune_samples(DO3, DODAspeciesfilt)
DO3genphylo <- prune_samples(DO3, DODAgenusfilt)

DA4famphylo <- prune_samples(DA4, DODAfamilyfilt)
DA4sppphylo <- prune_samples(DA4, DODAspeciesfilt)
DA4genphylo <- prune_samples(DA4, DODAgenusfilt)

##############WILCOX TESTING##################################################
#### FAMILY DPDO ####
# pull out OTU table, make sure order matches bw DP and DO
DP1fam <- otu_table(DP1famphylo, taxa_are_rows = TRUE)
DP1fam <- DP1fam[,c("30C", "6C", "63D", "45C", "58B", "50D", "48A", "8B", 
                    "1C", "64D", "74B", "18D")]
DO2fam <- otu_table(DO2famphylo, taxa_are_rows = TRUE)
head(DP1fam)
head(DO2fam)

# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DP1fam, DO2fam, paired = TRUE, mu = 0) 
# pvalue = 1.97e-8
DPDOfampvals <- sapply(1:nrow(DP1fam), function(i){
  wilcox.test(as.numeric(DP1fam[i,]), as.numeric(DO2fam[i,]), 
              paired = TRUE, mu = 0, exact = FALSE)$p.value
}) 
DPDOfampvals <- as.data.frame(DPDOfampvals, row.names=row.names(DP1fam))
pval <- DPDOfampvals[(DPDOfampvals$DPDOfampvals<0.1),]
rownames <- rownames(DPDOfampvals)[DPDOfampvals<0.1]
DPDOfamresults <- data.frame(pval, row.names = rownames)
DPDOfamresults
#pval
#f__[Mogibacteriaceae]  0.09168053
#f__Desulfovibrionaceae 0.05191296
#f__Methanobacteriaceae 0.04545875
#f__Pirellulaceae       0.02249427


####FDR CORRECTION####
DPDOfamresults <- mutate(DPDOfamresults, "padj" = p.adjust(DPDOfamresults$pval, 
                                                           method = "fdr", n = 12))
DPDOfamresults
#pval      padj
#1 0.09168053 0.2750416
#2 0.05191296 0.2076519
#3 0.04545875 0.2076519
#4 0.02249427 0.2076519


#### GENUS DPDO ####
DP1gen <- otu_table(DP1genphylo, taxa_are_rows = TRUE)
DP1gen <- DP1gen[,c("30C", "6C", "63D", "45C", "58B", "50D", "48A", "8B", 
                    "1C", "64D", "74B", "18D")]
DO2gen <- otu_table(DO2genphylo, taxa_are_rows = TRUE)
head(DP1gen)
head(DO2gen)
# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DP1gen, DO2gen, paired = TRUE, mu = 0) 
# pvalue = 0.2597
DPDOgenpvals <- sapply(1:nrow(DP1gen), function(i){
  wilcox.test(as.numeric(DP1gen[i,]), as.numeric(DO2gen[i,]), 
              paired = TRUE, mu = 0, exact = FALSE)$p.value
}) 
DPDOgenpvals <- as.data.frame(DPDOgenpvals, row.names=row.names(DP1gen))
pval <- DPDOgenpvals[(DPDOgenpvals$DPDOgenpvals<0.1),]
rownames <- rownames(DPDOgenpvals)[DPDOgenpvals<0.1]
DPDOgenresults <- data.frame(pval, row.names = rownames)
DPDOgenresults
#pval
#g__Anaeroplasma  0.06525726
#g__Bulleidia     0.04544695
#g__Butyrivibrio  0.07755617
#g__Paenibacillus 0.06654572

####FDR CORRECTION####
DPDOgenresults <- mutate(DPDOgenresults, "padj" = p.adjust(DPDOgenresults$pval, 
                                                           method = "fdr", n = 12))
DPDOgenresults
#pval      padj
#1 0.06525726 0.2326685
#2 0.04544695 0.2326685
#3 0.07755617 0.2326685
#4 0.06654572 0.2326685



#### SPECIES DPDO ####
DP1spp <- otu_table(DP1sppphylo, taxa_are_rows = TRUE)
DP1spp <- DP1spp[,c("30C", "6C", "63D", "45C", "58B", "50D", "48A", "8B", 
                    "1C", "64D", "74B", "18D")]
DO2spp <- otu_table(DO2sppphylo, taxa_are_rows = TRUE)
head(DP1spp)
head(DO2spp)
# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DP1spp, DO2spp, paired = TRUE, mu = 0) 
# pvalue = 0.4859
DPDOspppvals <- sapply(1:nrow(DP1spp), function(i){
  wilcox.test(as.numeric(DP1spp[i,]), as.numeric(DO2spp[i,]), 
              paired = TRUE, mu = 0, exact = FALSE)$p.value
}) 
DPDOspppvals <- as.data.frame(DPDOspppvals, row.names=row.names(DP1spp))
pval <- DPDOspppvals[(DPDOspppvals$DPDOspppvals<0.1),]
rownames <- rownames(DPDOspppvals)[DPDOspppvals<0.1]
DPDOsppresults <- data.frame(pval, row.names = rownames)
DPDOsppresults
#pval
#<0 rows> (or 0-length row.names)


### SPECIES DODA ###
DO3spp <- otu_table(DO3sppphylo, taxa_are_rows = TRUE)
DO3spp <- DO3spp[,c("10C", "30D", "74C", "8C", "38A", "67A", "58C", "45D", 
                    "57B", "6D")]
DA4spp <- otu_table(DA4sppphylo, taxa_are_rows = TRUE)
head(DO3spp)
head(DA4spp)
# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DO3spp, DA4spp, paired = TRUE, mu = 0) 
# pvalue = 0.9098
DODAspppvals <- sapply(1:nrow(DO3spp), function(i){
  wilcox.test(as.numeric(DO3spp[i,]), as.numeric(DA4spp[i,]), 
              paired = TRUE, mu = 0, exact = FALSE)$p.value
}) 
DODAspppvals <- as.data.frame(DODAspppvals, row.names=row.names(DO3spp))
pval <- DODAspppvals[(DODAspppvals$DODAspppvals<0.1),]
rownames <- rownames(DODAspppvals)[DODAspppvals<0.1]
DODAsppresults <- data.frame(pval, row.names = rownames)
DODAsppresults
#[1] pval
#<0 rows> (or 0-length row.names)

### GENUS DODA ###
DO3gen <- otu_table(DO3genphylo, taxa_are_rows = TRUE)
DO3gen <- DO3gen[,c("10C", "30D", "74C", "8C", "38A", "67A", "58C", "45D", 
                    "57B", "6D")]
DA4gen <- otu_table(DA4genphylo, taxa_are_rows = TRUE)
head(DO3gen)
head(DA4gen)
# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DO3gen, DA4gen, paired = TRUE, mu = 0) 
# pvalue = 0.07743
DODAgenpvals <- sapply(1:nrow(DO3gen), function(i){
  wilcox.test(as.numeric(DO3gen[i,]), as.numeric(DA4gen[i,]), 
              paired = TRUE, mu = 0, exact = FALSE)$p.value
}) 
DODAgenpvals <- as.data.frame(DODAgenpvals, row.names=row.names(DO3gen))
pval <- DODAgenpvals[(DODAgenpvals$DODAgenpvals<0.1),]
rownames <- rownames(DODAgenpvals)[DODAgenpvals<0.1]
DODAgenresults <- data.frame(pval, row.names = rownames)
DODAgenresults
#pval
#g__Corynebacterium 0.04231527
#g__rc4-4           0.06654572
DODAgenresults <- mutate(DODAgenresults, "padj" = p.adjust(DODAgenresults$pval, 
                                                           method = "fdr", n = 10))
DODAgenresults
#pval      padj
#1 0.04231527 0.3327286
#2 0.06654572 0.3327286


### FAMILY DODA ###
DO3fam <- otu_table(DO3famphylo, taxa_are_rows = TRUE)
DO3fam <- DO3fam[,c("10C", "30D", "74C", "8C", "38A", "67A", "58C", "45D", 
                    "57B", "6D")]
DA4fam <- otu_table(DA4famphylo, taxa_are_rows = TRUE)
head(DO3fam)
head(DA4fam)
# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DO3fam, DA4fam, paired = TRUE, mu = 0) 
# pvalue = 0.004145
DODAfampvals <- sapply(1:nrow(DO3fam), function(i){
  wilcox.test(as.numeric(DO3fam[i,]), as.numeric(DA4fam[i,]), 
              paired = TRUE, mu = 0, exact = FALSE)$p.value
}) 
DODAfampvals <- as.data.frame(DODAfampvals, row.names=row.names(DO3fam))
pval <- DODAfampvals[(DODAfampvals$DODAfampvals<0.1),]
rownames <- rownames(DODAfampvals)[DODAfampvals<0.1]
DODAfamresults <- data.frame(pval, row.names = rownames)
DODAfamresults
#pval
#f__Corynebacteriaceae 0.02997397
DODAfamresults <- mutate(DODAfamresults, "padj" = p.adjust(DODAfamresults$pval, 
                                                           method = "fdr", n = 10))
DODAfamresults
#pval      padj
#1 0.02997397 0.2997397
