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
## METAGENOMESEQ
table_4 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 8)
limits <- aes(ymax = table_4$CI.R, ymin = table_4$CI.L)
titles <- c("Bacillus coagulans", "Blautia producta", "Clostridium neonatale",
            "Moryella spp", "Faecalibacterium spp", "Methanosphaera spp", 
            "Family Acetobacteraceae", "Family Corynebacteriaceae", "Family Pirellulaceae")
ggplot(table_4, aes(x=Title, y = LogFC, fill = level)) +
  geom_bar(stat= "identity", position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9), width = 0.25) +
  labs(y = "Log2 Fold Change in a Non-O157 Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "italic")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = titles) +
  scale_fill_grey(start = 0.2, end = 0.9) +
  annotate("text",x=1,y=0.6,label="p=0.09") +
  annotate("text",x=4,y=-.65,label="p=0.04")

##DESEQ2
table_5 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 9)
limits_1 <- aes(ymax = table_5$CI.R, ymin = table_5$CI.L)
titles_1 <- c("Faecalibacterium prausnitzii", "Atopobium spp", "Corynebacterium spp",
            "Family Acetobacteraceae", "Family Pirellulaceae")
ggplot(table_5, aes(x=Title, y = LogFC, fill = level)) +
  geom_bar(stat= "identity", position = position_dodge(0.9)) +
  geom_errorbar(limits_1, position = position_dodge(0.9), width = 0.25) +
  labs(y = "Log2 Fold Change in a Non-O157 Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = "italic")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  scale_x_discrete(limits = titles_1) +
  #scale_colour_manual(values = c("turquoise3", "orange", "purple"))
  scale_fill_grey(start = 0.2, end = 0.9)
 
             
## what are the family/genus/species full names??
## METAGENOMESEQ RESULTS
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

## DESEQ2 RESULTS
prausnitzii <- subset(cowonlytaxa, Rank7 %in% "s__prausnitzii")#all faecalibacterium like above!
atopobium <- subset(cowonlytaxa, Rank6 %in% "g__Atopobium") # nothing to the species level
corynebacterium <- subset(cowonlytaxa, Rank6 %in% "g__Corynebacterium") #many diff spp
acetobacteraceae <- subset(cowonlytaxa, Rank5 %in% "f__Acetobacteraceae") #same as above for metagenomeseq
BS11 <- subset(cowonlytaxa, Rank5 %in% "f__BS11") #new cleanup reference means found in the final OTU 
# picking step when using open referencepicking
pirellulaceae <- subset(cowonlytaxa, Rank5 %in% "f__Pirellulaceae") # same as above, different genuses in this




## Want to try to turn a couple of tables into dot plots for easier understanding
table_6 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 13)
all_plot <- ggplot(data=table_6, aes(x = Levels, y = OR, ymin = Lower_CL, ymax = Upper_CL, 
                         color= Alpha_Diversity)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange", "purple")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  scale_x_discrete(labels=c("O157 Present", "O157 Present", 
                            "O157 Present", "Ever shed",
                            "Ever shed", "Ever shed",
                            "Intermittent","Intermittent",
                            "Intermittent", "Multi-day",
                            "Multi-day", "Multi-day"))+
  coord_flip()


table_7 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 10)
pathotype_plot <- ggplot(data=table_7, aes(x = Levels, y = OR, ymin = Lower_CL, ymax = Upper_CL, 
                         color= Alpha_Diversity)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange", "purple")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  scale_x_discrete(labels=c("O157 present", "O157 present", 
                            "O157 present"), name="Pathotype (Sample Level)")+
  
  coord_flip()

table_8 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 11)
ever_plot <- ggplot(data=table_8, aes(x = Levels, y = OR, ymin = Lower_CL, ymax = Upper_CL, 
                         color= Alpha_Diversity)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange", "purple")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  scale_x_discrete(labels=c("Ever shed", "Ever shed", 
                            "Ever shed"), name="Ever vs Never (Cow Level)")+
  
  coord_flip()


table_9 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 12)
pattern_plot <- ggplot(data=table_9, aes(x = Levels, y = OR, ymin = Lower_CL, ymax = Upper_CL, 
                         color= Alpha_Diversity)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange", "purple")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  scale_x_discrete(labels=c("Intermittent", "Intermittent", 
                            "Intermittent", "Multi-Day",
                            "Multi-Day", "Multi-Day"), 
                              name="Pattern (Cow Level)")+
  
  coord_flip()

table_new <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 16)
table_new$Levels <- as.factor(table_new$Levels)
ggplot(data=table_new, aes(x = name, y = OR, ymin = Lower_CL, ymax = Upper_CL, 
                           color= Levels)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  scale_y_continuous(limits = c(0.05, 1.25))+
  scale_x_reverse()+
  #scale_x_discrete(name="Pattern (Cow Level)", 
  #                labels= c("Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
  #              "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment",
  #             "Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
  #            "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment")) + 
  coord_flip()







###################### TRYING STUFF #################
table_10 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 14)
pattern_2nd_plot <- ggplot(data=table_10, aes(x = Adjustment, y = OR, ymin = Lower_CL, ymax = Upper_CL,
                                              color = Levels)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange", "purple")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  scale_x_discrete(labels=c("Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
                            "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment",
                            "Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
                            "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment"), 
                            name="Pattern (Cow Level)")+
  scale_y_continuous(limits = c(0.25, 1.75)) +
  
  coord_flip()

table_15 <- read.xlsx("excel sheets/Test output for alpha div models.xlsx", 15)
try <- ggplot(data=table_15, aes(x = Adjustment, y = OR, ymin = Lower_CL, ymax = Upper_CL,
                                              color = Levels)) +
  geom_point(size=3)+
  geom_errorbar(width = 0.2)+
  scale_colour_manual(values = c("turquoise3", "orange")) +
  geom_hline(aes(yintercept=1),color="red") + 
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_classic()+
  #scale_x_discrete(labels=c("Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
                            #"Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment",
                           # "Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
                           # "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment"), 
                  # name="Pattern (Cow Level)")+
  scale_y_continuous(limits = c(0.25, 1.75)) +
  
  coord_flip()

  

pattern_plot + scale_x_discrete(labels=c("Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
                            "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment",
                             "Crude", "Parity", "Treatment", "Farm", "Treatment + Parity",
                             "Treatment + Farm", "Parity + Farm", "Parity + Farm + Treatment"), 
                             name="Pattern (Cow Level)") 


                  
  

