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
          "6C", "6D", "74B", "74C", "8B", "8C", "8D", "8E")
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

DPDOgenus1 <- MRcounts(MR_DPDOgenphylofilt)
DPDOfamily1 <- MRcounts(MR_DPDOfamphylofilt)
DPDOspecies1 <- MRcounts(MR_DPDOsppphylofilt)

DODAgenus1 <- MRcounts(MR_DODAgenphylofilt)
DODAfamily1 <- MRcounts(MR_DODAfamphylofilt)
DODAspecies1 <- MRcounts(MR_DODAsppphylofilt)


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
         "6C", "74B", "8B", "8D")
DO2 <- c("18E", "1D", "30D", "45D", "48B", "50E", "58C", "63E", "64E", 
         "6D", "74C", "8C", "8E")
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






