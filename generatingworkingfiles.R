### Cowdatarich phyloseq object (non-normalized cow only experiment, 
### with richness estimates in metadata). Can also generate enviro only here

library("phyloseq")
biom_file <- "otu_table_mc2_w_tax_no_pynast_failures.biom"
map_file <- "QIIME_map_all_corrected_new.txt"
biomot = import_biom(biom_file, parseFunction = parse_taxonomy_default)
bmsd = import_qiime_sample_data(map_file)
Cowdata <- merge_phyloseq(biomot, bmsd)

enviro_samps <- c("75A", "75B", "75C", "75D", "75E", "76A", "76B", "76C", 
                  "76D", "76E", "77A", "77B", "77C", "77D", "77E", "78A", 
                  "78B", "78C", "78D", "78E", "79A", "79B", "79C", "79D", 
                  "79E", "80A", "80B", "80C", "80D", "80E", "81A", "81B", 
                  "81C", "81D", "81E", "82B", "82C", "82D", "82E", "86A", 
                  "86B", "86C", "86D", "86E", "87E", "87C", "87D", "87E", 
                  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
                  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 32, 33, 34, 
                  35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 49, 50, 51, 52, 53, 
                  54, 55, 56, 57, 58, 59, 60, 61, 65, 63, 64, 66, 67, 68, 69, 
                  70, 71, 72, 73, 74, 75, 76)
cow_samps <- c("50A", "50B", "50C", "50D", "50E", "71A", "71B", "71C", "71D", "67A",
               "67B", "67C", "67D", "67E", "63A", "63B", "63C", "63D", "63E", "58A",
               "58B", "58C", "58D", "58E", "55A", "55B", "55C", "55D", "55E", "62A", "62B", "62C", "62D", 
               "62E", "65A", "65B", "65C", "65D", "65E", "66A", "66B", "66C", "66D", "66E", "70A", "70B", "70C",
               "70D", "70E", "69A", "69B", "69C", "69D", "69E", "68A", "68B", "68C", "68D", "68E", "61A", "61B",
               "61C", "61D", "61E", "74A", "74B", "74C", "74D", "74E", "73A", "73B", "73C", "73D", "73E", "53A",
               "53B", "53C", "53D", "53E", "60B", "60C", "60D", "60E", "52A", "52B", "52C", "52D", "52E", "64A",
               "64B", "64C", "64D", "64E", "57B", "57C", "57D", "57E", "42A", "42B", "42C", "42D", "42E", "45A",
               "45B", "45C", "45D", "45E", "20A", "20B", "20C", "20D", "20E", "34A", "34B", "34C", "34D",
               "34E", "38A", "38B", "38C", "38D", "38E", "44A", "44B", "44C", "44D", "44E", "18A", "18B", "18C",
               "18D", "18E", "35A", "35B", "35C", "35D", "35E", "1A", "1B", "1C", "1D", "1E", "10A", "10B",
               "10C", "10D", "10E", "6A", "6B", "6C", "6D", "6E", "24A", "24B", "24C", "24D", "24E", "7A",
               "7B", "7C", "7D", "7E", "39A", "39B", "39C", "39D", "39E", "14A", "14B", "14C", "14D",
               "14E", "48A", "48B", "48C", "48D", "48E", "15A", "15B", "15C", "15D", "15E", "8A", "8B",
               "8C", "8D", "8E", "29A", "29B", "29C", "29D", "29E", "30A", "30B", "30C", "30D", "30E")

Cowonly <- prune_samples(cow_samps, Cowdata)
Enviroonly <- prune_samples(enviro_samps, Cowdata)

## Have added the richness estimates
new_map_file <- "excel sheets/Cow_map_wrichness.txt"
bmsd = import_qiime_sample_data(new_map_file)
Cowdatarich <- merge_phyloseq(biomot, bmsd)

## Normcounts biom table. Normalize the biom file in MetagenomeSeq, create data set with cows only,
## exports the normalized biom as a .tsv
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowdatarich)
NormCounts <- MRcounts(MRexp_cowonly, norm=TRUE, log=FALSE)
exportMat(NormCounts, file = file.path("NormalizedBiom.tsv"))

## Normcowdata. Making a phyloseq object with a normalized biom table
norm_otu_table <- otu_table(NormCounts, taxa_are_rows = TRUE)
new_tax_table <- tax_table(Cowonly)
Normcowdata <- merge_phyloseq(norm_otu_table, bmsd, new_tax_table)
Normcowdata 

## Normcowdatanew. Making a phyloseq object with normalized biom table and 
## additional estimates of alpha diversity (note this is diff than other workflow)
Updated_map <- import_qiime_sample_data("excel sheets/Cow_map_wrichnshansnevennormed.txt")
Normcowdatanew <- merge_phyloseq(norm_otu_table, Updated_map, new_tax_table)


## Generation of vegan data frames (OTU file, taxa file, treatment file; reducing
## by OTUs that have nothing in them and also reducing taxa file to only taxa present)
library(tidyverse)
OTUS <- otu_table(Cowdatarich)
sampledata <- sample_data(Cowdatarich)
taxa <- tax_table(Cowdatarich)
vegan_otu <- function(physeq) {
  OTU <- otu_table(Cowdatarich)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

OTUS <- vegan_otu(Cowdatarich)
taxanew <- data.frame(taxa)
taxanew <- mutate(taxanew, OTUs = rownames(taxa))
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)
taxa.df = taxanew[taxanew$OTU%in%otu.nm,]
