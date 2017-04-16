# Workflow using Phyloseq package. Goal is to import biom, tie to mapping file 
# sample metadata, measure richness of samples

# import the OTU table, the mapping file, then merge them into a phyloseq file
library("phyloseq")
# based on phyloseq vignette says to also load the ggplot2 library, and set a
# black and white theme
library("ggplot2")
library("xlsx")
theme_set(theme_bw())
biom_file <- "otu_table_mc2_w_tax_no_pynast_failures.biom"
map_file <- "QIIME_map_all_corrected_new.txt"
biomot = import_biom(biom_file, parseFunction = parse_taxonomy_default)
bmsd = import_qiime_sample_data(map_file)
Cowdata <- merge_phyloseq(biomot, bmsd)
Cowdata

# want to remove samples from this data set that are not cows....
# making two RData files that have SampleID for enviro or cow samples to use in 
# future analyses
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
               "10C", "10D", "10E", "6A", "6B", "6C", "6D", "6E", "24A", "24B", "24C", "24D", "24E",
               "7B", "7C", "7D", "7E", "39A", "39B", "39C", "39D", "39E", "14A", "14B", "14C", "14D",
               "14E", "48A", "48B", "48C", "48D", "48E", "15A", "15B", "15C", "15D", "15E", "8A", "8B",
               "8C", "8D", "8E", "29A", "29B", "29C", "29D", "29E", "30A", "30B", "30C", "30D", "30E")
# 7A was omitted from above because no sequence data
Cowonly <- prune_samples(cow_samps, Cowdata)
Enviroonly <- prune_samples(enviro_samps, Cowdata)

# now going to estimate richness and shannons (observed only when not normalized) for cows:
obsrich_values <- estimate_richness(Cowonly, measures="Observed")
write.xlsx(obsrich_values, "obsrich_values.xlsx")
shannon_values <- estimate_richness(Cowonly, measures="Shannon")
write.xlsx(shannon_values, "shannon_values.xlsx")


# turn to metagenomseq to normalize
library(metagenomeSeq)
MRexp_cowonly <- phyloseq_to_metagenomeSeq(Cowonly)
head(MRcounts(MRexp_cowonly))
NormCounts <- cumNormMat(MRexp_cowonly, p = cumNormStat(MRexp_cowonly))

######## OPTIONAL TO VIEW INFO ABOUT NORMALIZATION
# calculate the percentile
#p = cumNormStat(MRexp_cowonly)
# use that percentile to normalize the counts
#MRexp_cowonlynormed = cumNorm(MRexp_cowonly,p)
#show normalization factors for each sample??
#nf = normFactors(MRexp_cowonlynormed)


## Turn back to phyloseq. Making a phyloseq object with a normalized biom table
norm_otu_table <- otu_table(NormCounts, taxa_are_rows = TRUE)
new_tax_table <- tax_table(Cowonly)
Normcowdata <- merge_phyloseq(norm_otu_table, bmsd, new_tax_table)
Normcowdata

save(Normcowdata, file="Phyloseq files/Normcowdata")  ## this is the file used in 
## later code to aggregate taxa, make heat maps, ordination plots, etc ##



# Measuring Shannon's in the normalized samples:

shannon_values <- estimate_richness(Normcowdata, measures="Shannon")
write.xlsx(shannon_values, "shannon_normed_values.xlsx")

# Got this error:
# Warning message:
#In estimate_richness(Normcowdata, measures = "Shannon") :
#  The data you have provided does not have
#any singletons. This is highly suspicious. Results of richness
# estimates (for example) are probably unreliable, or wrong, if you have already
# trimmed low-abundance taxa from the data.

# After reviewing Paulson et al., seems that it is innapropriate to normalize 
# the data prior to alpha diversity metrics. Added these measures to the cow 
# only mapping file from before, but not normalized.
# How to calculate evenness? 
# http://sciencing.com/calculate-species-evenness-2851.html
# Divide Shannon's diversity index H by natural logarithm of species richness ln(S) 
# to calculate the species evenness. In the example, 0.707 divided by 1.099 equals 
# 0.64. Note that species evenness ranges from zero to one, with zero signifying no 
# evenness and one, a complete evenness.
# calculated this in excel under a column called evenness. Saved as txt and exported 
# back with the now normalized OTU table

# After talking with Zaid on 2.6.17, decided that because we are modeling the alpha
# diversity measures, we must normalize them. Otherwise some may be sequenced 
# deeper than others, and we have a biased estimate in who is there as it can be
# due to the level at which is was sequenced.
# Because richness estimation only allows integers, need to put the OTU table
# with normalized counts to a single digit before using the estimate richness 
# command. rounding won't work because values <0.4 will be zero. Will use ceiling
# command:
integeronlyOTU <- ceiling(NormCounts)
rich_otu_table <- otu_table(integeronlyOTU, taxa_are_rows = TRUE)
norm_phyloseq_for_richness <- merge_phyloseq(rich_otu_table, bmsd, new_tax_table)
normrich_values <- estimate_richness(norm_phyloseq_for_richness, measures="Observed")
write.xlsx(normrich_values, "normrich_values.xlsx")

# added values to new xlsx and imported back as .txt with the same name as prior
# but with normed on the end; Cow_map_wrichnshansnevennormed.txt
Updated_map <- import_qiime_sample_data("excel sheets/Cow_map_wrichnshansnevennormed.txt")
Normcowdatanew <- merge_phyloseq(norm_otu_table, Updated_map, new_tax_table)

# another phyloseq file that is non-normed for use in the basic stats R script, basically the 
# first generated phyloseq file but with a better map:
new_map_file <- "excel sheets/Cow_map_wrichness.txt"
bmsd = import_qiime_sample_data(new_map_file)
Cowdatarich <- merge_phyloseq(biomot, bmsd) # although this is original OTU table, only
# cows end up being kept here

# want to make sure that evenness has been calculated correctly (as defined in the vegan
# ecological diversity vignette). In excell calc the shannons/natlog(richness). 
# however in vegan vignette says to do J <- H/log(species number).... so log and not the 
# natural log (these are different!) but also says that 'H', which is shannons, is commonly
# the natural log. Tried these different transformations in the xcel sheet and get 
# various values. Decided to reput in the excel sheet and with vegan loaded do the calc:
library("vegan")
practice <- read.xlsx(file="excel sheets/Cow_map_wrichnshansnevennormed.xlsx", 1)
J <- practice$Normshann/log(practice$Normrich)
J

# the values that I get are the ones that I have used all along (many around 0.7 and 
# bounded between 0 and 1). Going to stick with what I have for now, as that is 
# the calculation used in vegan. Will cite vegan in text of paper
citation("vegan")


## MESSING AROUND WITH PACKAGE A BIT:
# look at richness based on pathotype metric:
#p <- plot_richness(Cowonly, x = "Pathotype", color = "Parity", measures = "Observed")

# pq <- p + geom_boxplot(data=sample_data(Cowonly), aes(x=Pathotype, y=value, color=Null), alpha=0.1)
#plot_richness(Cowonly, x = "Pathotype", color = "Parity")
#plot_richness(Cowonly, x = "Pathotype_1", color = "Parity")
#plot_richness(Cowonly, x = "Pattern_1", color = "Parity")
#plot_richness(Cowonly, x = "EvNev_1", color = "Parity")

# Looking at abundance bar plots:
#TopOTUs <- names(sort(taxa_sums(Cowonly), TRUE)[1:15])
#path15 = prune_taxa(TopOTUs, Cowonly)


#plot_bar(path15, "Pathotype_1", fill = "Parity")

#plot_heatmap(Normcowdatanew, "NMDS", "bray", "SampleID", "Family")


