load("Phyloseq files/Cowonly")
Cowonly
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
MRexp_cowonly

############ Using metagenomSeq and FitZig to look at diff abundance 
####between O157- and O157+ samples
d.gen = aggTax(MRexp_cowonly,lvl="Rank6", norm = TRUE) 
dim(d.gen) # features = 590
d.genfilt = filterData(d.gen, present = 50, depth = 1)
Pathotype_new <- factor(t.df$Pathotype_1) #makes it a factor
Individual_animalnew <- factor(t.df$Individual_animal)
d.genfiltnorm = cumNorm(d.genfilt, p = cumNormStatFast(d.genfilt))
trymod <- model.matrix(~ Pathotype_new + Individual_animalnew - 1, t.df) 
d.genfit = fitZig(d.genfiltnorm, trymod)
View(MRcoefs(d.genfit))
View(MRtable(d.genfit))
zigFitgen = d.genfit$fit
finalMod = d.genfit$fit$design
contrast.matrix = makeContrasts(Pathotype_new0-Pathotype_new1, levels = finalMod)
fit2 = contrasts.fit(zigFitgen, contrast.matrix)
fit2 
fit2 = eBayes(fit2) 
topTable(fit2)
View(topTable(fit2))

########## Using DeSeq2 to look at diff abundance between O157- and O157+ samps
# using McMurdie's code to calc geomeans first, then use a 'local' fit
d.gen1 = aggTax(MRexp_cowonly,lvl="Rank6", norm = FALSE)
d.genfilt1 = filterData(d.gen1, present = 50, depth = 1)
genus1 <- MRcounts(d.genfilt1)
genus2 <- MRcounts(d.gen1)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
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

########### Using wilcox test to look at paired differences in abundance between
#day prior to shedding versus day of shedding samples (DPDO)
d.gen = aggTax(MRexp_cowonly,lvl="Rank6", norm = TRUE) 
genus1 <- MRcounts(d.gen)
genusphylo <- otu_table(genus1, taxa_are_rows=TRUE)
genusphylo <- merge_phyloseq(genusphylo, sampledata)
DPDO <- c("18D", "18E", "1C", "1D", "30C", "30D", "45C", "45D", "48A",
          "48B", "50D", "50E", "58B", "58C", "63D", "63E", "64D", "64E",
          "6C", "6D", "74B", "74C", "8B", "8C", "8D", "8E")
DPDOgenphylo <- prune_samples(DPDO, genusphylo)
MR_DPDOgenphylo <- phyloseq_to_metagenomeSeq(DPDOgenphylo)
MR_DPDOgenphylofilt = filterData(MR_DPDOgenphylo, present = 8, depth = 1)
DPDOgenus1 <- MRcounts(MR_DPDOgenphylofilt, norm = TRUE)
DPDOgenusphylo <- otu_table(DPDOgenus1, taxa_are_rows=TRUE)
DPDOgenusfilt <- merge_phyloseq(DPDOgenusphylo, sampledata)
DP1 <- c("18D", "1C", "30C", "45C", "48A", "50D", "58B", "63D", "64D",
         "6C", "74B", "8B", "8D")
DO2 <- c("18E", "1D", "30D", "45D", "48B", "50E", "58C", "63E", "64E", 
         "6D", "74C", "8C", "8E")
DP1genphylo <- prune_samples(DP1, DPDOgenusfilt)
DO2genphylo <- prune_samples(DO2, DPDOgenusfilt)
#### GENUS DPDO ####
DP1gen <- otu_table(DP1genphylo, taxa_are_rows = TRUE)
DP1gen <- DP1gen[,c("30C", "6C", "63D", "45C", "58B", "50D", "48A", "8B", 
                    "1C", "64D", "74B", "18D", "8D")]
DO2gen <- otu_table(DO2genphylo, taxa_are_rows = TRUE)
DP1gen
DO2gen
# wilcox test overall, then wilcox test giving p values by row and pulling out <0.1
wilcox.test(DP1gen, DO2gen, paired = TRUE, mu = 0) 
# pvalue = 0.5814
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
#g__Atopobium      0.07556057
#g__Bacteroides    0.06921297
#g__Bulleidia      0.03763288
#g__Methanosphaera 0.09349248

####FDR CORRECTION####
DPDOgenresults <- mutate(DPDOgenresults, "padj" = p.adjust(DPDOgenresults$pval, 
                                                           method = "fdr", n = 13))
DPDOgenresults
#pval      padj
#1 0.07556057 0.3038506
#2 0.06921297 0.3038506
#3 0.03763288 0.3038506
#4 0.09349248 0.3038506


