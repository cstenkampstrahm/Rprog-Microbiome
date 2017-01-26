# Workflow using Phyloseq package. Goal is to import biom, tie to mapping file 
# sample metadata, measure richness of samples

# import the OTU table, the mapping file, then merge them into a phyloseq file
library("phyloseq")
biom_file <- "otu_table_mc2_w_tax_no_pynast_failures.biom"
map_file <- "QIIME_map_all_corrected_new.txt"
biomot = import_biom(biom_file, parseFunction = parse_taxonomy_default)
bmsd = import_qiime_sample_data(map_file)
Cowdata <- merge_phyloseq(biomot, bmsd)
Cowdata

# look at richness based on pathotype metric:
plot_richness(Cowdata, x = "Pathotype", color = "Parity")

# without parity (greys because have NAs)
plot_richness(Cowdata, x = "Pathotype")

# want to remove samples from this data set that are not cows....
