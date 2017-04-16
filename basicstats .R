# 3/12/17- going to practice aggregating taxa with the file zaid created for me
# using the phyloseq workflow, can generate the cowdatarich (nonnormed) phyloseq object
# and also the Normcowdatanew (normalized with CSS) phyloseq object
# additionally cowdata is phyloseq with everything (enviro and non) and cowonly and environly
# have only those samples, respectively (there is just different data in their metadata). saved
# files to load easily later

#save(Cowonly, file="Cowonly")
#save(Enviroonly, file="Enviroonly")
#save(Cowdatarich, file="Cowdata")
#save(Normcowonly, file="Normcowdatanew")
#source(file="Phyloseq_workflow.R")
load("Phyloseq files/Cowdata") # load this for non normed data
load("Phyloseq files/Normcowonly") # load this for normed data
Cowdatarich
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
d.cs = apply(OTUS,2,sum)
OTUS = OTUS[,d.cs>0]
otu.nm = colnames(OTUS)

taxa.df = taxanew[row.names(taxanew)%in%otu.nm,]
library(tidyverse)
taxa.df <-mutate(taxa.df, OTUs = row.names(taxa.df))

# from Zaid's analysis file:
library(pheatmap)
library(vegan)
library(ggplot2)
library(phyloseq)
library(metagenomeSeq)
library(cluster)
library(NbClust)
source(file="files from zaid/functions_Chloe.R")

taxa1.df = taxa.cleanup.qiime.ftn(taxa.df=taxa.df) # added a note in here so can
# save this generated taxa table with taxa parsed out by level
#save(taxa1.df, file="splitcowtaxa")
#giving the OTUS in the split taxa names based on those in the OTU table
nm.vt = colnames(OTUS)
nm.vt = paste("OTU",nm.vt,sep="")
colnames(OTUS)=nm.vt
# now combo the OTU table and the taxa table based on the level desired
## PHYLUM
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=2)
phylum.df = combined.ls[[1]]
#will get a taxa table with aggregates by sample up to the Phylum level. doing this
#for all different levels now to generate data files to work with
## CLASS
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=3)
class.df = combined.ls[[1]]
## ORDER
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=4)
order.df = combined.ls[[1]]
## FAMILY
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=5)
family.df = combined.ls[[1]]
## GENUS
combined.ls = taxa.split.combine.qiime.ftn(taxa.df=taxa1.df,d.df=OTUS,l=6)
genus.df = combined.ls[[1]]



##### OLD CODE #######
# want to look at spread of data for Table 1:
library("xlsx")
table_1 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
pathotype_farm <- table_1 %>% group_by(Pathotype_1, Farm) %>% summarise(pathnval = n())
evernever_farm <- table_1 %>% group_by(EvNev_1, Farm) %>% summarise(evnevnval = n())
pattern_farm <- table_1 %>% group_by(Pattern_1, Farm) %>% summarise(pattnval = n())
pattern_DIM1 <- table_1 %>% group_by(Pattern_1, DIM) %>% summarise(pattnval = n())
# DIM day 1
table_2 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)
shapiro.test(table_2$DIM)
pattern_DIM <- table_2 %>% group_by(Pattern_1, DIM) %>% summarise(pattnval = n())
pattern_DIM_mean <- table_2 %>% group_by(Pattern_1) %>% summarise(mean = mean(DIM), sd = sd(DIM))
table_2 <- mutate(table_2, Pattern_1 = as.factor(Pattern_1))
anova(table_2$DIM, table_2$Pattern_1)
fisher.test(table_2$DIM, table_2$Pattern_1)
pathotype_DIM <- table_2 %>% group_by(Pathotype_1, DIM) %>% summarise(pathnval = n())
pathotype_DIM_mean <- table_2 %>% group_by(Pathotype_1) %>% summarise(mean = mean(DIM), sd = sd(DIM))
table_2 <- mutate(table_2, Pathotype_1 = as.factor(Pathotype_1))
t.test(table_2$DIM ~ table_2$Pathotype_1)
evernever_DIM <- table_2 %>% group_by(EvNev_1, DIM) %>% summarise(evnevnval = n())
evernevr_DIM_mean <- table_2 %>% group_by(EvNev_1) %>% summarise(mean = mean(DIM), sd = sd(DIM))
table_2 <- mutate(table_2, EvNev_1 = as.factor(EvNev_1))
t.test(table_2$DIM ~ table_2$EvNev_1)
pathotype_parity <- table_1 %>% group_by(Pathotype_1, Parity) %>% summarise(paritynval = n())
pathotype_parity1 <- table_1 %>% group_by(Pathotype_1, Parity_1) %>% summarise(paritynval = n())
pattern_parity <- table_1 %>% group_by(Pattern_1, Parity) %>% summarise(paritynval = n())
pattern_parity1 <- table_2 %>% group_by(Pattern_1, Parity_1) %>% summarise(paritynval = n())
evnev_parity <- table_1 %>% group_by(EvNev_1, Parity) %>% summarise(paritynval = n())
evnev_parity1 <- table_2 %>% group_by(EvNev_1, Parity_1) %>% summarise(paritynval = n())
pathotype_dx <- table_1 %>% group_by(Pathotype_1, Disease) %>% summarise(dxnval = n())
pattern_dx <- table_1 %>% group_by(Pattern_1, Disease) %>% summarise(dxnval = n())
evnev_dx <- table_1 %>% group_by(EvNev_1, Disease) %>% summarise(dxnval = n())
#### END OF OLD CODE ####



###### looking to redo table 1 on 4.10.17
table_1 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
table_2 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)

#Farm

pathotype_farm <- table_1 %>% group_by(Pathotype_1, Farm) %>% summarise(pathnval = n())
evernever_farm <- table_2 %>% group_by(EvNev_1, Farm) %>% summarise(evnevnval = n())
pattern_farm <- table_2 %>% group_by(Pattern_1, Farm) %>% summarise(pattnval = n())

chisq.test(table_1$Farm, table_1$Pathotype_1)
chisq.test(table_2$Farm, table_2$EvNev_1)
chisq.test(table_2$Farm, table_2$Pattern_1) # one cell has 5 values, so fisher
fisher.test(table_2$Farm, table_2$Pattern_1)

# DIM

shapiro.test(table_1$DIM) # no not normal
shapiro.test(table_2$DIM) # yes normal

table_1 <- mutate(table_1, Pathotype_1 = as.factor(Pathotype_1))
pathotype_DIM_IQR <- table_1 %>% group_by(Pathotype_1) %>% summarise(median = median(DIM), 
                          IQR = IQR(DIM)) 
wilcox.test(table_1$DIM ~ table_1$Pathotype_1)

evernevr_DIM_mean <- table_2 %>% group_by(EvNev_1) %>% summarise(mean = mean(DIM), sd = sd(DIM))
table_2 <- mutate(table_2, EvNev_1 = as.factor(EvNev_1))
t.test(table_2$DIM ~ table_2$EvNev_1)

pattern_DIM_mean <- table_2 %>% group_by(Pattern_1) %>% summarise(mean = mean(DIM), sd = sd(DIM))
table_2 <- mutate(table_2, Pattern_1 = as.factor(Pattern_1))
aov1 <- lm(DIM ~ Pattern_1, data = table_2)

# Disease

pathotype_dx <- table_1 %>% group_by(Pathotype_1, Disease) %>% summarise(dxnval = n())
chisq.test(table_1$Disease, table_1$Pathotype_1)

pattern_dx <- table_2 %>% group_by(Pattern_1, Disease) %>% summarise(dxnval = n())
fisher.test(table_2$Disease, table_2$Pattern_1)

evnev_dx <- table_2 %>% group_by(EvNev_1, Disease) %>% summarise(dxnval = n())
chisq.test(table_2$Disease, table_2$EvNev_1)

# Parity
pattern_parity1 <- table_2 %>% group_by(Pattern_1, Parity_1) %>% summarise(paritynval = n())
fisher.test(table_2$Pattern_1, table_2$Parity_1)

evnev_parity1 <- table_2 %>% group_by(EvNev_1, Parity_1) %>% summarise(paritynval = n())
fisher.test(table_2$EvNev_1, table_2$Parity_1)

pathotype_parity1 <- table_1 %>% group_by(Pathotype_1, Parity_1) %>% summarise(paritynval = n())
fisher.test(table_2$Pathotype_1, table_2$Parity_1)  
fisher.test(table_2$Parity_1, table_2$Pathotype_1)

                                                                            
# want to look at the counts for OTUs by sample. Can't figure out how to 
# group by metadata variables in phyloseq. Will output the sums and add to 
# metadata file (cow_map_wrichnshansnevennormednscaled.xlsx)
OTUcts <- sample_sums(Normcowdatanew)
OTUcts
write.xlsx(OTUcts, "OTUctsnormed.xlsx")
OTUcts <- sample_sums(Cowdatarich)
write.xlsx(OTUcts, "OTUcts.xlsx")

# non normalized counts
alpha_div <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
path_OTU <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathOTUavg = mean(OTUcts), pathOTUsd = sd(OTUcts),
            pathOTUmin = min(OTUcts),pathOTUmax = max(OTUcts), 
            pathnval = n()) %>% 
  ungroup()

evnev_OTU <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevOTUavg = mean(OTUcts), evnevOTUsd = sd(OTUcts), 
            evnevOTUmin = min(OTUcts),evnevOTUmax = max(OTUcts), 
            evnevnval = n()) %>% 
  ungroup()

patt_OTUs <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattOTUavg = mean(OTUcts), pattOTUsd = sd(OTUcts), 
            pattOTUmin = min(OTUcts),pattOTUmax = max(OTUcts), 
            pattnval = n())

