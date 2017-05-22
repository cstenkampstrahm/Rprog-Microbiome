# Notes about files from working deseq 2 r script:
# Cowonly = phyloseq experiment with only cows in it, non-normalized
# genus1, family1, species1, otu1 are MRexp count tables filtered to those taxa
# present in 50 animals at least, non normalized
# genus2, species2, family2, otu2 are MRexp count tables non filtered and non normalized
## TO REMAKE THESE COUNT TABLES ##
d.gen1 = aggTax(MRexp_cowonly,lvl="Rank6", norm = FALSE)
d.genfilt1 = filterData(d.gen1, present = 50, depth = 1)
genus1 <- MRcounts(d.genfilt1)
genus2 <- MRcounts(d.gen1)

d.spp1 = aggTax(MRexp_cowonly,lvl="Rank7", norm = FALSE) 
d.sppfilt1 = filterData(d.spp1, present = 50, depth = 1)
spp1 <- MRcounts(d.sppfilt1)
spp2 <- MRcounts(d.spp1)

d.fam1 = aggTax(MRexp_cowonly,lvl="Rank5", norm = FALSE)
d.famfilt1 = filterData(d.fam1, present = 50, depth = 1)
family1 <- MRcounts(d.famfilt1)
family2 <- MRcounts(d.fam1)

d.otu1 = filterData(MRexp_cowonly, present = 50, depth = 1)
otu1 <- MRcounts(d.otu1)
otu2 <- MRcounts(MRexp_cowonly)

#McMurdie's calculation of the geometric mean 

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# options with DESeq2:
# non filtered vs filtered to those in 50 cows, (these do come up with values slightly different)
# local fit vs parametric fit (these do come up with values slightly different)
# estimating geometric means as mcmurdie prior to DESeq, vs not estimating geometric means prior
# (estimating gm_means does come up with values slightly different)


#### GENUS LEVEL ######
## NON-FILTERED
genusnofilt <- otu_table(genus2, taxa_are_rows=TRUE)
genusnofilt <- merge_phyloseq(genusnofilt, sampledata)
sample_data(genusnofilt)$Pathotype_1 <- factor(sample_data(genusnofilt)$Pathotype_1)
sample_data(genusnofilt)$Individual_animal <- factor(sample_data(genusnofilt)$Individual_animal)
sample_data(genusnofilt)$Pathotype_1 <- relevel(sample_data(genusnofilt)$Pathotype_1, "0")
genusnofiltDESEQ <- phyloseq_to_deseq2(genusnofilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
genusnofiltDESEQP <- DESeq(genusnofiltDESEQ, fitType = "parametric", betaPrior = FALSE)
genusnofiltresultsP <- results(genusnofiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res1 = genusnofiltresultsP
res1 = res1[order(res1$padj, na.last=NA), ]
alpha = 0.2
sigtab1 = res1[(res1$padj < alpha), ]
sigtab2 = res1[(res1$pvalue < alpha), ]
sigtab1
sigtab2
# now with local
genusnofiltDESEQL <- DESeq(genusnofiltDESEQ, fitType = "local", betaPrior = FALSE)
genusnofiltresultsL <- results(genusnofiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res2 = genusnofiltresultsL
res2 = res2[order(res2$padj, na.last=NA), ]
alpha = 0.2
sigtab3 = res2[(res2$padj < alpha), ]
sigtab4 = res2[(res2$pvalue < alpha), ]
sigtab3
sigtab4
# now with geomeans first
geoMeans = apply(counts(genusnofiltDESEQ), 1, gm_mean) 
genusnofiltDESEQg = estimateSizeFactors(genusnofiltDESEQ, geoMeans = geoMeans)
genusnofiltDESEQg = DESeq(genusnofiltDESEQg, fitType="local", betaPrior = FALSE)
genusnofiltresultsG <- results(genusnofiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res5 = genusnofiltresultsG
res5 = res5[order(res5$padj, na.last=NA), ]
alpha = 0.2
sigtab5 = res5[(res5$padj < alpha), ]
sigtab6 = res5[(res5$pvalue < alpha), ]
sigtab5
sigtab6

## Filtered
genusfilt <- otu_table(genus1, taxa_are_rows=TRUE)
genusfilt <- merge_phyloseq(genusfilt, sampledata)
sample_data(genusfilt)$Pathotype_1 <- factor(sample_data(genusfilt)$Pathotype_1)
sample_data(genusfilt)$Individual_animal <- factor(sample_data(genusfilt)$Individual_animal)
sample_data(genusfilt)$Pathotype_1 <- relevel(sample_data(genusfilt)$Pathotype_1, "0")
genusfiltDESEQ <- phyloseq_to_deseq2(genusfilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
genusfiltDESEQP <- DESeq(genusfiltDESEQ, fitType = "parametric", betaPrior = FALSE)
genusfiltresultsP <- results(genusfiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res1.1 = genusfiltresultsP
res1.1 = res1.1[order(res1.1$padj, na.last=NA), ]
alpha = 0.2
sigtab1.1 = res1.1[(res1.1$padj < alpha), ]
sigtab2.1 = res1.1[(res1.1$pvalue < alpha), ]
sigtab1.1
sigtab2.1
# now with local
genusfiltDESEQL <- DESeq(genusfiltDESEQ, fitType = "local", betaPrior = FALSE)
genusfiltresultsL <- results(genusfiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res2.1 = genusfiltresultsL
res2.1 = res2.1[order(res2.1$padj, na.last=NA), ]
alpha = 0.2
sigtab3.1 = res2.1[(res2.1$padj < alpha), ]
sigtab4.1 = res2.1[(res2.1$pvalue < alpha), ]
sigtab3.1
sigtab4.1
# now with geomeans first
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
