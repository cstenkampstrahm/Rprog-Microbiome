# Notes about files from working deseq 2 r script:
# Cowonly = phyloseq experiment with only cows in it, non-normalized
# genus1, family1, species1, otu1 are MRexp count tables filtered to those taxa
# present in 50 animals at least, non normalized
# genus2, species2, family2, otu2 are MRexp count tables non filtered and non normalized
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
# look at genus, species, otu and family level
# non filtered vs filtered to those in 50 cows, (these do come up with values slightly different)
# local fit vs parametric fit (these do come up with values slightly different)
# estimating geometric means as mcmurdie prior to DESeq w local fit, vs not estimating 
# geometric means prior (estimating gm_means does come up with values slightly different)


#### GENUS LEVEL ######
## NON-FILTERED
genusnofilt <- otu_table(genus2, taxa_are_rows=TRUE)
genusnofilt <- merge_phyloseq(genusnofilt, sampledata) # somehow can acquire the sample data from 
# original metagenomeseq object
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

## FILTERED
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





#### SPECIES LEVEL ######
## NON-FILTERED
sppnofilt <- otu_table(spp2, taxa_are_rows=TRUE)
sppnofilt <- merge_phyloseq(sppnofilt, sampledata)
sample_data(sppnofilt)$Pathotype_1 <- factor(sample_data(sppnofilt)$Pathotype_1)
sample_data(sppnofilt)$Individual_animal <- factor(sample_data(sppnofilt)$Individual_animal)
sample_data(sppnofilt)$Pathotype_1 <- relevel(sample_data(sppnofilt)$Pathotype_1, "0")
sppnofiltDESEQ <- phyloseq_to_deseq2(sppnofilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
sppnofiltDESEQP <- DESeq(sppnofiltDESEQ, fitType = "parametric", betaPrior = FALSE)
sppnofiltresultsP <- results(sppnofiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res6 = sppnofiltresultsP
res6 = res6[order(res6$padj, na.last=NA), ]
alpha = 0.2
sigtab7 = res6[(res6$padj < alpha), ]
sigtab8 = res6[(res6$pvalue < alpha), ]
sigtab7
sigtab8
# now with local
sppnofiltDESEQL <- DESeq(sppnofiltDESEQ, fitType = "local", betaPrior = FALSE)
sppnofiltresultsL <- results(sppnofiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res7 = sppnofiltresultsL
res7 = res7[order(res7$padj, na.last=NA), ]
alpha = 0.2
sigtab9 = res7[(res7$padj < alpha), ]
sigtab10 = res7[(res7$pvalue < alpha), ]
sigtab9
sigtab10
# now with geomeans first
geoMeans = apply(counts(sppnofiltDESEQ), 1, gm_mean) 
sppnofiltDESEQg = estimateSizeFactors(sppnofiltDESEQ, geoMeans = geoMeans)
sppnofiltDESEQg = DESeq(sppnofiltDESEQg, fitType="local", betaPrior = FALSE)
sppnofiltresultsG <- results(sppnofiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res8 = sppnofiltresultsG
res8 = res8[order(res8$padj, na.last=NA), ]
alpha = 0.2
sigtab11 = res8[(res8$padj < alpha), ]
sigtab12 = res8[(res8$pvalue < alpha), ]
sigtab11
sigtab12

## FILTERED
sppfilt <- otu_table(spp1, taxa_are_rows=TRUE)
sppfilt <- merge_phyloseq(sppfilt, sampledata)
sample_data(sppfilt)$Pathotype_1 <- factor(sample_data(sppfilt)$Pathotype_1)
sample_data(sppfilt)$Individual_animal <- factor(sample_data(sppfilt)$Individual_animal)
sample_data(sppfilt)$Pathotype_1 <- relevel(sample_data(sppfilt)$Pathotype_1, "0")
sppfiltDESEQ <- phyloseq_to_deseq2(sppfilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
sppfiltDESEQP <- DESeq(sppfiltDESEQ, fitType = "parametric", betaPrior = FALSE)
sppfiltresultsP <- results(sppfiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res6.1 = sppfiltresultsP
res6.1 = res6.1[order(res6.1$padj, na.last=NA), ]
alpha = 0.2
sigtab7.1 = res6.1[(res6.1$padj < alpha), ]
sigtab8.1 = res6.1[(res6.1$pvalue < alpha), ]
sigtab7.1
sigtab8.1
# now with local
sppfiltDESEQL <- DESeq(sppfiltDESEQ, fitType = "local", betaPrior = FALSE)
sppfiltresultsL <- results(sppfiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res7.1 = sppfiltresultsL
res7.1 = res7.1[order(res7.1$padj, na.last=NA), ]
alpha = 0.2
sigtab9.1 = res7.1[(res7.1$padj < alpha), ]
sigtab10.1 = res7.1[(res7.1$pvalue < alpha), ]
sigtab9.1
sigtab10.1
# now with geomeans first
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







#### FAMILY LEVEL ######
## NON-FILTERED
famnofilt <- otu_table(family2, taxa_are_rows=TRUE)
famnofilt <- merge_phyloseq(famnofilt, sampledata)
sample_data(famnofilt)$Pathotype_1 <- factor(sample_data(famnofilt)$Pathotype_1)
sample_data(famnofilt)$Individual_animal <- factor(sample_data(famnofilt)$Individual_animal)
sample_data(famnofilt)$Pathotype_1 <- relevel(sample_data(famnofilt)$Pathotype_1, "0")
famnofiltDESEQ <- phyloseq_to_deseq2(famnofilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
famnofiltDESEQP <- DESeq(famnofiltDESEQ, fitType = "parametric", betaPrior = FALSE)
famnofiltresultsP <- results(famnofiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res9 = famnofiltresultsP
res9 = res9[order(res9$padj, na.last=NA), ]
alpha = 0.2
sigtab13 = res9[(res9$padj < alpha), ]
sigtab14 = res9[(res9$pvalue < alpha), ]
sigtab13
sigtab14
# now with local
famnofiltDESEQL <- DESeq(famnofiltDESEQ, fitType = "local", betaPrior = FALSE)
famnofiltresultsL <- results(famnofiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res10 = famnofiltresultsL
res10 = res10[order(res10$padj, na.last=NA), ]
alpha = 0.2
sigtab15 = res10[(res10$padj < alpha), ]
sigtab16 = res10[(res10$pvalue < alpha), ]
sigtab15
sigtab16
# now with geomeans first
geoMeans = apply(counts(famnofiltDESEQ), 1, gm_mean) 
famnofiltDESEQg = estimateSizeFactors(famnofiltDESEQ, geoMeans = geoMeans)
famnofiltDESEQg = DESeq(famnofiltDESEQg, fitType="local", betaPrior = FALSE)
famnofiltresultsG <- results(famnofiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res11 = famnofiltresultsG
res11 = res11[order(res11$padj, na.last=NA), ]
alpha = 0.2
sigtab17 = res11[(res11$padj < alpha), ]
sigtab18 = res11[(res11$pvalue < alpha), ]
sigtab17
sigtab18

## FILTERED
famfilt <- otu_table(family1, taxa_are_rows=TRUE)
famfilt <- merge_phyloseq(famfilt, sampledata)
sample_data(famfilt)$Pathotype_1 <- factor(sample_data(famfilt)$Pathotype_1)
sample_data(famfilt)$Individual_animal <- factor(sample_data(famfilt)$Individual_animal)
sample_data(famfilt)$Pathotype_1 <- relevel(sample_data(famfilt)$Pathotype_1, "0")
famfiltDESEQ <- phyloseq_to_deseq2(famfilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
famfiltDESEQP <- DESeq(famfiltDESEQ, fitType = "parametric", betaPrior = FALSE)
famfiltresultsP <- results(famfiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res9.1 = famfiltresultsP
res9.1 = res9.1[order(res9.1$padj, na.last=NA), ]
alpha = 0.2
sigtab13.1 = res9.1[(res9.1$padj < alpha), ]
sigtab14.1 = res9.1[(res9.1$pvalue < alpha), ]
sigtab13.1
sigtab14.1
# now with local
famfiltDESEQL <- DESeq(famfiltDESEQ, fitType = "local", betaPrior = FALSE)
famfiltresultsL <- results(famfiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res10.1 = famfiltresultsL
res10.1 = res10.1[order(res10.1$padj, na.last=NA), ]
alpha = 0.2
sigtab15.1 = res10.1[(res10.1$padj < alpha), ]
sigtab16.1 = res10.1[(res10.1$pvalue < alpha), ]
sigtab15.1
sigtab16.1
# now with geomeans first
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






##### OTU LEVEL #####
## NON-FILTERED
otunofilt <- otu_table(otu2, taxa_are_rows=TRUE)
otunofilt <- merge_phyloseq(otunofilt, sampledata)
sample_data(otunofilt)$Pathotype_1 <- factor(sample_data(otunofilt)$Pathotype_1)
sample_data(otunofilt)$Individual_animal <- factor(sample_data(otunofilt)$Individual_animal)
sample_data(otunofilt)$Pathotype_1 <- relevel(sample_data(otunofilt)$Pathotype_1, "0")
otunofiltDESEQ <- phyloseq_to_deseq2(otunofilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
otunofiltDESEQP <- DESeq(otunofiltDESEQ, fitType = "parametric", betaPrior = FALSE)
otunofiltresultsP <- results(otunofiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res12 = otunofiltresultsP
res12 = res12[order(res12$padj, na.last=NA), ]
alpha = 0.2
sigtab19 = res12[(res12$padj < alpha), ]
sigtab20 = res12[(res12$pvalue < alpha), ]
sigtab19
sigtab20
# now with local
otunofiltDESEQL <- DESeq(otunofiltDESEQ, fitType = "local", betaPrior = FALSE)
otunofiltresultsL <- results(otunofiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res13 = otunofiltresultsL
res13 = res12[order(res12$padj, na.last=NA), ]
alpha = 0.2
sigtab21 = res13[(res13$padj < alpha), ]
sigtab22 = res13[(res13$pvalue < alpha), ]
sigtab21
sigtab22
# now with geomeans first
geoMeans = apply(counts(otunofiltDESEQ), 1, gm_mean) 
otunofiltDESEQg = estimateSizeFactors(otunofiltDESEQ, geoMeans = geoMeans)
otunofiltDESEQg = DESeq(otunofiltDESEQg, fitType="local", betaPrior = FALSE)
otunofiltresultsG <- results(otunofiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res14 = otunofiltresultsG
res14 = res14[order(res14$padj, na.last=NA), ]
alpha = 0.2
sigtab23 = res14[(res14$padj < alpha), ]
sigtab24 = res14[(res14$pvalue < alpha), ]
sigtab23
sigtab24

## FILTERED
otufilt <- otu_table(otu1, taxa_are_rows=TRUE)
otufilt <- merge_phyloseq(otufilt, sampledata)
sample_data(otufilt)$Pathotype_1 <- factor(sample_data(otufilt)$Pathotype_1)
sample_data(otufilt)$Individual_animal <- factor(sample_data(otufilt)$Individual_animal)
sample_data(otufilt)$Pathotype_1 <- relevel(sample_data(otufilt)$Pathotype_1, "0")
otufiltDESEQ <- phyloseq_to_deseq2(otufilt, ~Individual_animal + Pathotype_1 - 1)
# with parametric
otufiltDESEQP <- DESeq(otufiltDESEQ, fitType = "parametric", betaPrior = FALSE)
otufiltresultsP <- results(otufiltDESEQP, contrast=c("Pathotype_1", "1", "0"))
res12.1 = otufiltresultsP
res12.1 = res12.1[order(res12.1$padj, na.last=NA), ]
alpha = 0.2
sigtab19.1 = res12.1[(res12.1$padj < alpha), ]
sigtab20.1 = res12.1[(res12.1$pvalue < alpha), ]
sigtab19.1
sigtab20.1
# now with local
otufiltDESEQL <- DESeq(otufiltDESEQ, fitType = "local", betaPrior = FALSE)
otufiltresultsL <- results(otufiltDESEQL, contrast=c("Pathotype_1", "1", "0"))
res13.1 = otufiltresultsL
res13.1 = res13.1[order(res13.1$padj, na.last=NA), ]
alpha = 0.2
sigtab21.1 = res13.1[(res13.1$padj < alpha), ]
sigtab22.1 = res13.1[(res13.1$pvalue < alpha), ]
sigtab21.1
sigtab22.1
# now with geomeans first
geoMeans = apply(counts(otufiltDESEQ), 1, gm_mean)
otufiltDESEQg = estimateSizeFactors(otufiltDESEQ, geoMeans = geoMeans)
otufiltDESEQg = DESeq(otufiltDESEQg, fitType="local", betaPrior = FALSE)
otufiltresultsG <- results(otufiltDESEQg, contrast=c("Pathotype_1", "1", "0"))
res14.1 = otufiltresultsG
res14.1 = res14.1[order(res14.1$padj, na.last=NA), ]
alpha = 0.2
sigtab23.1 = res14.1[(res14.1$padj < alpha), ]
sigtab24.1 = res14.1[(res14.1$pvalue < alpha), ]
sigtab23.1
sigtab24.1














###### Now want to look at DP vs NS #####
### Only going to make models with filtered or non filtered and the geomeans option
## Make phyloseq object with cows of interest
# DP samples and single random sample from NS cows (no duplicate cows in either,
# n = 12 and n = 24)
dpvsnsSubset <- c("18D", "1C", "30C", "45C", "48A", "50D", "58B", "63D", "64D", "6C", 
                  "74B", "8B", "14D", "15C", "20C", "24C", "29A", "34D", "35B", "39D", 
                  "42D", "44C", "52B", "53B", "55A", "60B", "61C", "62D", "65E", "66A", 
                  "68B", "69C", "70D", "73E", "71A", "7B")
dpvsnsSubset <- prune_samples(dpvsnsSubset, Cowonly)

# All DP samples, all NS samples (many duplicate cows, n = 13 and n = 117)
dpvsns <- c("18D", "1C", "30C", "45C", "48A", "50D", "58B", "63D", "64D", "6C", 
            "74B", "8B", "8D", "14A", "14B", "14C", "14D", "14E", "15A", "15B",
            "15C", "15D", "15E", "20A", "20B", "20C", "20D", "20E", "24A", "24B", 
            "24C", "24D", "24E", "29A", "29B", "29C", "29D", "29E", "34A", "34B",
            "34C", "34D", "34E", "35A", "35B", "35C", "35D", "35E", "39A", "39B", 
            "39C", "39D", "39E", "42A", "42B", "42C", "42D", "42E", "44A", "44B",
            "44C", "44D", "44E", "52A", "52B", "52C", "52D", "52E", "53A", "53B",
            "53C", "53D", "53E", "55A", "55B", "55C", "55D", "55E", "60B", "60C",
            "60D", "60E", "61A", "61B", "61C", "61D", "61E", "62A", "62B", "62C",
            "62D", "62E", "65A", "65B", "65C", "65D", "65E", "66A", "66B", "66C",
            "66D", "66E", "68A", "68B", "68C", "68D", "68E", "69A", "69B", "69C",
            "69D", "69E", "70A", "70B", "70C", "70D", "70E", "71A", "71B", "71C",
            "71D", "73A", "73B", "73C", "73D", "73E", "7B", "7C", "7D", "7E")
dpvsns <- prune_samples(dpvsns, Cowonly)

#### DPvsNSSubset ####

MRexp_dpvsnsSubset <- phyloseq_to_metagenomeSeq(dpvsnsSubset)
MRexp_dpvsnsSubset

gen1dnsub = aggTax(MRexp_dpvsnsSubset,lvl="Rank6", norm = FALSE)
genfiltdnsub = filterData(gen1dnsub, present = 12, depth = 1)
genus1dnsub <- MRcounts(genfiltdnsub)
genus2dnsub <- MRcounts(gen1dnsub)

spp1dnsub = aggTax(MRexp_dpvsnsSubset,lvl="Rank7", norm = FALSE) 
sppfiltdnsub = filterData(spp1dnsub, present = 12, depth = 1)
spp.1dnsub <- MRcounts(sppfiltdnsub)
spp.2dnsub <- MRcounts(spp1dnsub)

fam1dnsub = aggTax(MRexp_dpvsnsSubset,lvl="Rank5", norm = FALSE)
famfiltdnsub = filterData(fam1dnsub, present = 12, depth = 1)
family1dnsub <- MRcounts(famfiltdnsub)
family2dnsub <- MRcounts(fam1dnsub)

otu1dnsub = filterData(MRexp_dpvsnsSubset, present = 12, depth = 1)
otu1.dnsub <- MRcounts(otu1dnsub)
otu2.dnsub <- MRcounts(MRexp_dpvsnsSubset)

#### GENUS LEVEL ######
## NON-FILTERED
genusnofiltdnsub <- otu_table(genus2dnsub, taxa_are_rows=TRUE)
genusnofiltdnsub <- merge_phyloseq(genusnofiltdnsub, sampledata) # somehow can acquire the sample data from 
# original metagenomeseq object
sample_data(genusnofiltdnsub)$DPvNS <- factor(sample_data(genusnofiltdnsub)$DPvNS)
sample_data(genusnofiltdnsub)$Individual_animal <- factor(sample_data(genusnofiltdnsub)$Individual_animal) # don't need this line of code really...
sample_data(genusnofiltdnsub)$DPvNS <- relevel(sample_data(genusnofiltdnsub)$DPvNS, "0")
genusnofiltdnsubDESEQ <- phyloseq_to_deseq2(genusnofiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(genusnofiltdnsubDESEQ), 1, gm_mean) 
genusnofiltdnsubDESEQg = estimateSizeFactors(genusnofiltdnsubDESEQ, geoMeans = geoMeans)
genusnofiltdnsubDESEQg = DESeq(genusnofiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
genusnofiltdnsubresultsG <- results(genusnofiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub1 = genusnofiltdnsubresultsG
resdnsub1 = resdnsub1[order(resdnsub1$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub1 = resdnsub1[(resdnsub1$padj < alpha), ]
sigtabdnsub2 = resdnsub1[(resdnsub1$pvalue < alpha), ]
sigtabdnsub1
sigtabdnsub2

## FILTERED
genusfiltdnsub <- otu_table(genus1dnsub, taxa_are_rows=TRUE)
genusfiltdnsub <- merge_phyloseq(genusfiltdnsub, sampledata)
sample_data(genusfiltdnsub)$DPvNS <- factor(sample_data(genusfiltdnsub)$DPvNS)
sample_data(genusfiltdnsub)$Individual_animal <- factor(sample_data(genusfiltdnsub)$Individual_animal)
sample_data(genusfiltdnsub)$DPvNS <- relevel(sample_data(genusfiltdnsub)$DPvNS, "0")
genusfiltdnsubDESEQ <- phyloseq_to_deseq2(genusfiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(genusfiltdnsubDESEQ), 1, gm_mean)
genusfiltdnsubDESEQg = estimateSizeFactors(genusfiltdnsubDESEQ, geoMeans = geoMeans)
genusfiltdnsubDESEQg = DESeq(genusfiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
genusfiltdnsubresultsG <- results(genusfiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub3 = genusfiltdnsubresultsG
resdnsub3 = resdnsub3[order(resdnsub3$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub3= resdnsub3[(resdnsub3$padj < alpha), ]
sigtabdnsub4 = resdnsub3[(resdnsub3$pvalue < alpha), ]
sigtabdnsub3
sigtabdnsub4

#### SPECIES LEVEL ######
## NON-FILTERED
sppnofiltdnsub <- otu_table(spp.2dnsub, taxa_are_rows=TRUE)
sppnofiltdnsub <- merge_phyloseq(sppnofiltdnsub, sampledata)
sample_data(sppnofiltdnsub)$DPvNS <- factor(sample_data(sppnofiltdnsub)$DPvNS)
sample_data(sppnofiltdnsub)$Individual_animal <- factor(sample_data(sppnofiltdnsub)$Individual_animal)
sample_data(sppnofiltdnsub)$DPvNS <- relevel(sample_data(sppnofiltdnsub)$DPvNS, "0")
sppnofiltdnsubDESEQ <- phyloseq_to_deseq2(sppnofiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(sppnofiltdnsubDESEQ), 1, gm_mean) 
sppnofiltdnsubDESEQg = estimateSizeFactors(sppnofiltdnsubDESEQ, geoMeans = geoMeans)
sppnofiltdnsubDESEQg = DESeq(sppnofiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
sppnofiltdnsubresultsG <- results(sppnofiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub4 = sppnofiltdnsubresultsG
resdnsub4 = resdnsub4[order(resdnsub4$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub5 = resdnsub4[(resdnsub4$padj < alpha), ]
sigtabdnsub6 = resdnsub4[(resdnsub4$pvalue < alpha), ]
sigtabdnsub5
sigtabdnsub6

## FILTERED
sppfiltdnsub <- otu_table(spp.1dnsub, taxa_are_rows=TRUE)
sppfiltdnsub <- merge_phyloseq(sppfiltdnsub, sampledata)
sample_data(sppfiltdnsub)$DPvNS <- factor(sample_data(sppfiltdnsub)$DPvNS)
sample_data(sppfiltdnsub)$Individual_animal <- factor(sample_data(sppfiltdnsub)$Individual_animal)
sample_data(sppfiltdnsub)$DPvNS <- relevel(sample_data(sppfiltdnsub)$DPvNS, "0")
sppfiltdnsubDESEQ <- phyloseq_to_deseq2(sppfiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(sppfiltdnsubDESEQ), 1, gm_mean)
sppfiltdnsubDESEQg = estimateSizeFactors(sppfiltdnsubDESEQ, geoMeans = geoMeans)
sppfiltdnsubDESEQg = DESeq(sppfiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
sppfiltdnsubresultsG <- results(sppfiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub5 = sppfiltdnsubresultsG
resdnsub5 = resdnsub5[order(resdnsub5$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub7 = resdnsub5[(resdnsub5$padj < alpha), ]
sigtabdnsub8 = resdnsub5[(resdnsub5$pvalue < alpha), ]
sigtabdnsub7
sigtabdnsub8

#### FAMILY LEVEL ######
## NON-FILTERED
famnofiltdnsub <- otu_table(family2dnsub, taxa_are_rows=TRUE)
famnofiltdnsub <- merge_phyloseq(famnofiltdnsub, sampledata)
sample_data(famnofiltdnsub)$DPvNS <- factor(sample_data(famnofiltdnsub)$DPvNS)
sample_data(famnofiltdnsub)$Individual_animal <- factor(sample_data(famnofiltdnsub)$Individual_animal)
sample_data(famnofiltdnsub)$DPvNS <- relevel(sample_data(famnofiltdnsub)$DPvNS, "0")
famnofiltdnsubDESEQ <- phyloseq_to_deseq2(famnofiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(famnofiltdnsubDESEQ), 1, gm_mean) 
famnofiltdnsubDESEQg = estimateSizeFactors(famnofiltdnsubDESEQ, geoMeans = geoMeans)
famnofiltdnsubDESEQg = DESeq(famnofiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
famnofiltdnsubresultsG <- results(famnofiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub6 = famnofiltdnsubresultsG
resdnsub6 = resdnsub6[order(resdnsub6$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub9 = resdnsub6[(resdnsub6$padj < alpha), ]
sigtabdnsub10 = resdnsub6[(resdnsub6$pvalue < alpha), ]
sigtabdnsub9
sigtabdnsub10

## FILTERED
famfiltdnsub <- otu_table(family1dnsub, taxa_are_rows=TRUE)
famfiltdnsub <- merge_phyloseq(famfiltdnsub, sampledata)
sample_data(famfiltdnsub)$DPvNS <- factor(sample_data(famfiltdnsub)$DPvNS)
sample_data(famfiltdnsub)$Individual_animal <- factor(sample_data(famfiltdnsub)$Individual_animal)
sample_data(famfiltdnsub)$DPvNS <- relevel(sample_data(famfiltdnsub)$DPvNS, "0")
famfiltdnsubDESEQ <- phyloseq_to_deseq2(famfiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(famfiltdnsubDESEQ), 1, gm_mean)
famfiltdnsubDESEQg = estimateSizeFactors(famfiltdnsubDESEQ, geoMeans = geoMeans)
famfiltdnsubDESEQg = DESeq(famfiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
famfiltdnsubresultsG <- results(famfiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub7 = famfiltdnsubresultsG
resdnsub7 = resdnsub7[order(resdnsub7$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub11 = resdnsub7[(resdnsub7$padj < alpha), ]
sigtabdnsub12 = resdnsub7[(resdnsub7$pvalue < alpha), ]
sigtabdnsub11
sigtabdnsub12

##### OTU LEVEL #####
## NON-FILTERED
otunofiltdnsub <- otu_table(otu2.dnsub, taxa_are_rows=TRUE)
otunofiltdnsub <- merge_phyloseq(otunofiltdnsub, sampledata)
sample_data(otunofiltdnsub)$DPvNS <- factor(sample_data(otunofiltdnsub)$DPvNS)
sample_data(otunofiltdnsub)$Individual_animal <- factor(sample_data(otunofiltdnsub)$Individual_animal)
sample_data(otunofiltdnsub)$DPvNS <- relevel(sample_data(otunofiltdnsub)$DPvNS, "0")
otunofiltdnsubDESEQ <- phyloseq_to_deseq2(otunofiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(otunofiltdnsubDESEQ), 1, gm_mean) 
otunofiltdnsubDESEQg = estimateSizeFactors(otunofiltdnsubDESEQ, geoMeans = geoMeans)
otunofiltdnsubDESEQg = DESeq(otunofiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
otunofiltdnsubresultsG <- results(otunofiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub8 = otunofiltdnsubresultsG
resdnsub8 = resdnsub8[order(resdnsub8$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub13 = resdnsub8[(resdnsub8$padj < alpha), ]
sigtabdnsub14 = resdnsub8[(resdnsub8$pvalue < alpha), ]
sigtabdnsub13
sigtabdnsub14

## FILTERED
otufiltdnsub <- otu_table(otu1.dnsub, taxa_are_rows=TRUE)
otufiltdnsub <- merge_phyloseq(otufiltdnsub, sampledata)
sample_data(otufiltdnsub)$DPvNS <- factor(sample_data(otufiltdnsub)$DPvNS)
sample_data(otufiltdnsub)$Individual_animal <- factor(sample_data(otufiltdnsub)$Individual_animal)
sample_data(otufiltdnsub)$DPvNS <- relevel(sample_data(otufiltdnsub)$DPvNS, "0")
otufiltdnsubDESEQ <- phyloseq_to_deseq2(otufiltdnsub, ~DPvNS - 1)
# now with geomeans first
geoMeans = apply(counts(otufiltdnsubDESEQ), 1, gm_mean)
otufiltdnsubDESEQg = estimateSizeFactors(otufiltdnsubDESEQ, geoMeans = geoMeans)
otufiltdnsubDESEQg = DESeq(otufiltdnsubDESEQg, fitType="local", betaPrior = FALSE)
otufiltdnsubresultsG <- results(otufiltdnsubDESEQg, contrast=c("DPvNS", "1", "0"))
resdnsub9 = otufiltdnsubresultsG
resdnsub9 = resdnsub9[order(resdnsub9$padj, na.last=NA), ]
alpha = 0.2
sigtabdnsub15 = resdnsub9[(resdnsub9$padj < alpha), ]
sigtabdnsub16 = resdnsub9[(resdnsub9$pvalue < alpha), ]
sigtabdnsub15
sigtabdnsub16



















#### DPvsNS  ####
###### CAN'T LOOK AT THE DATA THIS WAY BECAUSE HAVE LINEAR COMBOS OF COLUMNS! HAVE TO DO IT
##### WITHOUT CONTROLLING FOR ANIMAL, BECAUSE THE ANIMAL WILL PREDICT THE DPVSNS CLASSIFICATION (LINEAR
MRexp_dpvsns <- phyloseq_to_metagenomeSeq(dpvsns)
MRexp_dpvsns


gen1dn = aggTax(MRexp_dpvsns,lvl="Rank6", norm = FALSE)
genfiltdn = filterData(gen1dn, present = 33, depth = 1)
genus1dn <- MRcounts(genfiltdn)
genus2dn <- MRcounts(gen1dn)

spp1dn = aggTax(MRexp_dpvsns,lvl="Rank7", norm = FALSE) 
sppfiltdn = filterData(spp1dn, present = 33, depth = 1)
spp.1dn <- MRcounts(sppfiltdn)
spp.2dn <- MRcounts(spp1dn)

fam1dn = aggTax(MRexp_dpvsns,lvl="Rank5", norm = FALSE)
famfiltdn = filterData(fam1dn, present = 33, depth = 1)
family1dn <- MRcounts(famfiltdn)
family2dn <- MRcounts(fam1dn)

otu1dn = filterData(MRexp_dpvsns, present = 33, depth = 1)
otu1.dn <- MRcounts(otu1dn)
otu2.dn <- MRcounts(MRexp_dpvsns)

#### GENUS LEVEL ######
## NON-FILTERED
genusnofiltdn <- otu_table(genus2dn, taxa_are_rows=TRUE)
genusnofiltdn <- merge_phyloseq(genusnofiltdn, sampledata) # somehow can acquire the sample data from 
# original metagenomeseq object
sample_data(genusnofiltdn)$DPvNS <- factor(sample_data(genusnofiltdn)$DPvNS)
sample_data(genusnofiltdn)$Individual_animal <- factor(sample_data(genusnofiltdn)$Individual_animal) # not full rank, bec one variable linear
