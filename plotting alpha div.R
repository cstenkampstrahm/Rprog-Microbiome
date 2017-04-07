### variance = avg of squared differences from mean, std deviation is square
### root of the variance
library("xlsx")
alpha_div <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
library(ggplot2)
library(tidyverse)

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
