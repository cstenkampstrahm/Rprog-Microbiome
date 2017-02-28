### variance = avg of squared differences from mean, std deviation is square
### root of the variance
library("xlsx")
alpha_div <- read.xlsx("Cow_map_wrichnshansnevennormed.xlsx", 1)
library(ggplot2)
library(tidyverse)

#Richness
path_richness <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathrichavg = mean(Obs_richness), pathrichsd = sd(Obs_richness), 
            pathrichmin = min(Obs_richness),pathrichmax = max(Obs_richness), 
            pathnval = n()) %>% 
  ungroup()

evnev_richness <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevrichavg = mean(Obs_richness), evnevrichsd = sd(Obs_richness), 
            evnevrichmin = min(Obs_richness),evnevrichmax = max(Obs_richness), 
            evnevnval = n()) %>% 
  ungroup()
  
patt_richness <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattrichavg = mean(Obs_richness), pattrichsd = sd(Obs_richness), 
            pattrichmin = min(Obs_richness),pattrichmax = max(Obs_richness), 
            pattnval = n())

#Shannons
path_shann <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathshanavg = mean(Shannon), pathshannsd = sd(Shannon), 
            pathshannmin = min(Shannon),pathshannmax = max(Shannon), 
            pathnval = n()) %>% 
  ungroup()

evnev_shann <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevshannavg = mean(Shannon), evnevshannsd = sd(Shannon), 
            evnevshannmin = min(Shannon),evnevshannmax = max(Shannon), 
            evnevnval = n()) %>% 
  ungroup()

patt_shann <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattshannavg = mean(Shannon), pattshannsd = sd(Shannon), 
            pattshannmin = min(Shannon),pattshannmax = max(Shannon), 
            pattnval = n())

# Evenness
path_even <- alpha_div %>% group_by(Pathotype_1) %>%
  summarise(pathevenavg = mean(Evenness), pathevensd = sd(Evenness), 
            pathevennmin = min(Evenness),pathevenmax = max(Evenness), 
            pathnval = n()) %>% 
  ungroup()

evnev_even <- alpha_div %>% group_by(EvNev_1) %>%
  summarise(evnevevenavg = mean(Evenness), evnevevensd = sd(Evenness), 
            evnevevenmin = min(Evenness),evnevevenmax = max(Evenness), 
            evnevnval = n()) %>% 
  ungroup()

patt_even <- alpha_div %>% group_by(Pattern_1) %>%
  summarise(pattevenavg = mean(Evenness), pattevensd = sd(Evenness), 
            pattevenmin = min(Evenness),pattevenmax = max(Evenness), 
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
