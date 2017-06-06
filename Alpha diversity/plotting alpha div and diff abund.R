### variance = avg of squared differences from mean, std deviation is square
### root of the variance
library("xlsx")
alpha_div <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
library(ggplot2)
library(tidyverse)

######## Old code, skip down for newer code for box and whisker per Sheryl recomm.######
#Richness
path_richness <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathrichavg = mean(Normrich), pathrichsd = sd(Normrich),
            pathrichmin = min(Normrich),pathrichmax = max(Normrich), 
            pathnval = n()) %>% 
  ungroup()

evnev_richness <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevrichavg = mean(Normrich), evnevrichsd = sd(Normrich), 
            evnevrichmin = min(Normrich),evnevrichmax = max(Normrich), 
            evnevnval = n()) %>% 
  ungroup()
  
patt_richness <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattrichavg = mean(Normrich), pattrichsd = sd(Normrich), 
            pattrichmin = min(Normrich),pattrichmax = max(Normrich), 
            pattnval = n())

#Shannons
path_shann <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathshanavg = mean(Normshann), pathshannsd = sd(Normshann), 
            pathshannmin = min(Normshann),pathshannmax = max(Normshann), 
            pathnval = n()) %>% 
  ungroup()

evnev_shann <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevshannavg = mean(Normshann), evnevshannsd = sd(Normshann), 
            evnevshannmin = min(Normshann),evnevshannmax = max(Normshann), 
            evnevnval = n()) %>% 
  ungroup()

patt_shann <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattshannavg = mean(Normshann), pattshannsd = sd(Normshann), 
            pattshannmin = min(Normshann),pattshannmax = max(Normshann), 
            pattnval = n())

# Evenness
path_even <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathevenavg = mean(Normeven), pathevensd = sd(Normeven), 
            pathevennmin = min(Normeven),pathevenmax = max(Normeven), 
            pathnval = n()) %>% 
  ungroup()

evnev_even <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevevenavg = mean(Normeven), evnevevensd = sd(Normeven), 
            evnevevenmin = min(Normeven),evnevevenmax = max(Normeven), 
            evnevnval = n()) %>% 
  ungroup()

patt_even <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattevenavg = mean(Normeven), pattevensd = sd(Normeven), 
            pattevenmin = min(Normeven),pattevenmax = max(Normeven), 
            pattnval = n())
###placed value output in the excel sheet that had the other test output for
#also added a column with n values and calculated the std error (sd/sqrt(n))
#per zaid recommendation, need se's on the graph 
###modeling in it
library(ggthemes)
library(xlsx)
plotting <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 4)
plotting <- mutate(plotting, errorbarSE = as.numeric(0.5*Error))

plot_1 <- ggplot(plotting, aes(x=Data, y=Value, fill = Type)) + 
                   geom_bar(stat="identity", colour="white",position="dodge", 
                            show.legend = FALSE, width=0.5) +
                    geom_errorbar(aes(ymin = Value-errorbarSE, 
                    ymax = Value+errorbarSE),
                    width=.1, color = "red") +
                   facet_grid(alpha_measure~Type, scales = "free") + 
                  xlab('') + ylab("Alpha Diversity Value") +
                  theme_bw() +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        panel.background = element_rect(fill = "white"))

 ### getting some stats for Table 1 based on shedding status of cows
library(tidyverse)
EvNevsummary <- alpha_div %>% group_by(Farm, Individual_animal) %>% 
  summarise(EvNev_1 = sum(EvNev_1))

patternsummary <- alpha_div %>% group_by(Farm, Individual_animal) %>%
  summarise(Pattern = sum(Pattern_1))

pathotypesummary <- alpha_div %>% group_by(Farm, Day) %>%
  summarise(Pathotype_1 = sum(Pathotype_1), n = n())
### double check the values
dgaf <- alpha_div %>% select(Farm, Day, Pathotype_1) %>% group_by(Farm, Day)  

##### New code for box and whisker plot on 4.11.17
table_1 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
table_2 <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)
#Evnev used averaged values
#Pattern used averaged values
#Pathotype used values not averaged
path_richness <- table_1 %>% group_by(Pathotype_1) %>% 
  summarise(ymin = min(Normrich), ymax = max(Normrich),
            median = median(Normrich), 
            lower = quantile(Normrich, probs=0.25),
            upper = quantile(Normrich, probs=0.75),
            pathnval = n()) %>%
  ungroup()

path_shannons <- table_1 %>% group_by(Pathotype_1) %>% 
  summarise(ymin = min(Normshann), ymax = max(Normshann),
            median = median(Normshann), 
            lower = quantile(Normshann, probs=0.25),
            upper = quantile(Normshann, probs=0.75),
            pathnval = n()) %>%
  ungroup()

path_even <- table_1 %>% group_by(Pathotype_1) %>% 
  summarise(ymin = min(Normeven), ymax = max(Normeven),
            median = median(Normeven), 
            lower = quantile(Normeven, probs=0.25),
            upper = quantile(Normeven, probs=0.75),
            pathnval = n()) %>%
  ungroup()

patt_richness <- table_2 %>% group_by(Pattern_1) %>% 
  summarise(ymin = min(Avgrich), ymax = max(Avgrich),
            median = median(Avgrich), 
            lower = quantile(Avgrich, probs=0.25),
            upper = quantile(Avgrich, probs=0.75),
            pattnval = n()) %>%
  ungroup()

patt_shann <- table_2 %>% group_by(Pattern_1) %>% 
  summarise(ymin = min(Avgshann), ymax = max(Avgshann),
            median = median(Avgshann), 
            lower = quantile(Avgshann, probs=0.25),
            upper = quantile(Avgshann, probs=0.75),
            pattnval = n()) %>%
  ungroup()

patt_even <- table_2 %>% group_by(Pattern_1) %>% 
  summarise(ymin = min(Avgeven), ymax = max(Avgeven),
            median = median(Avgeven), 
            lower = quantile(Avgeven, probs=0.25),
            upper = quantile(Avgeven, probs=0.75),
            pattnval = n()) %>%
  ungroup()

evnev_richness <- table_2 %>% group_by(EvNev_1) %>% 
  summarise(ymin = min(Avgrich), ymax = max(Avgrich),
            median = median(Avgrich), 
            lower = quantile(Avgrich, probs=0.25),
            upper = quantile(Avgrich, probs=0.75),
            pattnval = n()) %>%
  ungroup()

evnev_shann <- table_2 %>% group_by(EvNev_1) %>% 
  summarise(ymin = min(Avgshann), ymax = max(Avgshann),
            median = median(Avgshann), 
            lower = quantile(Avgshann, probs=0.25),
            upper = quantile(Avgshann, probs=0.75),
            pattnval = n()) %>%
  ungroup()

evnev_even <- table_2 %>% group_by(EvNev_1) %>% 
  summarise(ymin = min(Avgeven), ymax = max(Avgeven),
            median = median(Avgeven), 
            lower = quantile(Avgeven, probs=0.25),
            upper = quantile(Avgeven, probs=0.75),
            mid = quantile(Avgeven, probs=0.5),
            pattnval = n()) %>%
  ungroup()
# placed values into excel sheet
## still trying to change the order of the third facet- pattern- to be
## Never shed, intermittent, multi day rather than intermittent, multi-day, never shed
## (this would make each facet uniform with never shedding being first)

table_3 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 7)

ggplot(aes(y = median, x = Level, color= Alphadiv), data = table_3) +
  facet_grid(Alphadiv ~ O157metric, scales = "free") +
  geom_boxplot(aes(ymin = ymin, ymax = ymax, middle = median, 
                   lower = lower, upper = upper), stat = "identity", fill='#A4A4A4', color="black") +
  theme_classic() +
  #reorder(Pattern, c("Never shed", "Intermittent", "Multi-day")) +
  #reorder(Level, c("No O157", "O157", "Never shed", "Shed > 1 time", 
  #                 "Never shed", "Intermittent", "Multi-day")) +
  #scale_x_discrete(labels = table_3$Level) +
  #scale_x_discrete(labels = c("No O157", "O157", "Never shed", "Shed > 1 time", 
   #                                 "Never shed", "Intermittent", "Multi-day")) +
  #reorder(table_3, O157metric$Pattern, c("Never shed", "Intermittent", "Multi-day")) +
  labs(y = "IQR of Alpha Diversity", x = "O157 Metric") +
  theme(legend.position = "none")

### Plotting for Differential Abundance results
table_4 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 8)
limits <- aes(ymax = table_4$CI.R, ymin = table_4$CI.L)
titles <- c("Bacillus coagulans", "Blautia producta", "Clostridium neonatale",
            "Moryella spp", "Faecalibacterium spp", "Methanosphaera spp", 
            "Family Acetobacteraceae", "Family Corynebacteriaceae", "Family Pirellulaceae")
ggplot(table_4, aes(x=Title, y = LogFC, fill = level)) +
  geom_bar(stat= "identity", position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
  labs(y = "Log2 Fold Change from Non-O157 Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8, face = "italic")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = titles) +
  scale_fill_grey(start = 0.2, end = 0.9) +
  annotate("text",x=1,y=0.6,label="p=0.09") +
  annotate("text",x=4,y=-.65,label="p=0.04")
             
## what are the family/genus/species full names??
Cowonly
cowonlytaxa <- as.data.frame(tax_table(Cowonly))
View(cowonlytaxa)
coagulans <- subset(cowonlytaxa, Rank7 %in% "s__coagulans")
producta <- subset(cowonlytaxa, Rank7 %in% "s__producta")
neonatale <- subset(cowonlytaxa, Rank7 %in% "s__neonatale")
moryella <- subset(cowonlytaxa, Rank6 %in% "g__Moryella")
faecalibacterium <- subset(cowonlytaxa, Rank6 %in% "g__Faecalibacterium") # looks like only prausnitzii!
methanosphaera <- subset(cowonlytaxa, Rank6 %in% "g__Methanosphaera")
acetobacteraceae <- subset(cowonlytaxa, Rank5 %in% "f__Acetobacteraceae")
corynebacteriaceae <- subset(cowonlytaxa, Rank5 %in% "f__Corynebacteriaceae")
pirellulaceae <- subset(cowonlytaxa, Rank5 %in% "f__Pirellulaceae")

