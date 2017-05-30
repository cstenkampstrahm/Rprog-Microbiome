## want to look at differential abundance of matched samples:
## Day prior vs Day of AND Day of vs Day after
## To do this need to take cowonly phyloseq obj with all the data, agg to family, spp, genus
# as done previously WITH NORMALIZATION (in metagenomeseq), convert to phyloseq and prune samples
# to day prior & day of (DPDO) and day of & day after (DODA). Put in metagenomeseq denovo 
# and filter (when together this ensures the same taxa are filtered for each set of data) 
# for DPDO n = 13 for each group, n = 26 total. Will filter to present in 8 (30%) samples
# for DODA n = 10 for each group, n = 20 total. WIll filter to present in 6 (30%) samples
# then put in phyloseq de novo and prune to DP1, DO2, DO3, DA4 data sets. These matrices
# can then be used for wilcox.test in stats package
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
