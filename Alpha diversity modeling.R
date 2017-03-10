# Want to run some models using richness, shannons, evenness as the outcome of 
# interest, and then I have the type of bacteria they are shedding as the input 
# value (independent variable). But need to control for farm/anima/day. how best
# to do this? Generalized linear mixed model, random effects model, linear mixed 
# model etc (diff names for the same thing....). will use a random effects model
# with e.coli status as the fixed effect, the outcome alpha diversity metric 
# as the dependent variable, and the random effects as farm, animal and day

library("xlsx")
library("lme4")

# removed the one sample with NA values (7A)

alpha_div <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormed.xlsx", 1)
# need to check for normality of Normrich Normshann Normeven. Must convert to 
# class numeric

library("tidyverse")
alpha_div <- mutate(alpha_div, "Normrich" = as.numeric(Normrich), 
       "Normshann" = as.numeric(Normshann),
       "Normeven" = as.numeric(Normeven))
class(alpha_div$Normrich)
class(alpha_div$Normshann)
class(alpha_div$Normeven)
shapiro.test(alpha_div$Normrich)
shapiro.test(alpha_div$Normshann)
shapiro.test(alpha_div$Normeven)

# Got output:
# data:  alpha_div$Normrich
# W = 0.94881, p-value = 1.808e-06

# data:  alpha_div$Normshann
# W = 0.96648, p-value = 0.0001263

# data:  alpha_div$Normeven
# W = 0.98961, p-value = 0.1666

# Going to check the qqplot, looks like all values but the evenness are NOT 
# normally distributed, esp richness. shannons really hard to tell....
qqnorm(alpha_div$Normeven, ylab= "Evenness")
qqnorm(alpha_div$Normrich, ylab= "Richness")
qqnorm(alpha_div$Normshann, ylab= "Shannons")

# will try to log transform??
alpha_div$Normrich1 <-log10(alpha_div$Normrich)
shapiro.test(alpha_div$Normrich1)

# data:  alpha_div$Normrich1
# W = 0.79537, p-value = 2.699e-15

alpha_div$Normshann1 <- log10(alpha_div$Normshann)
shapiro.test(alpha_div$Normshann1)

# data:  alpha_div$Normshann1
# W = 0.94952, p-value = 2.11e-06

# Discussed the issues of normality with Ryan Gan on 2/8/27. Showed qq plots
# etc. Because plots look VERY close to normal... and the log transforming does
# not change the shapiro wilk outcome, it shows that the distributions are robust
# and we can conclude that the normality is good enough for our model assumptions
# will continue to model using the Normshann, Normeven, Normrich as fixed effects
# Further discussed how to put the random effects in for the pattern and ev nev 
# metrics of interest. Probably only stratify by farm and cow, since every day
# we have classified as the same with an individual, yet we have still sampled them
# five repeated times on different farms.


#EVENNESS 

###Pathotype
m1.lme4.even <- lmer(Normeven ~ Pathotype_1 + (1|Individual_animal) +
                     (1|Day) + (1|Farm), data = alpha_div)
summary(m1.lme4.even)
anova(m1.lme4.even)

# calculating the null model so can compare max liklihood estimates and pvals
# added REML=FALSE to above too
m1.lme4.even <- lmer(Normeven ~ Pathotype_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m1.lme4.evennull <- lmer(Normeven ~ (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)

anova(m1.lme4.evennull,m1.lme4.even)

#Got the following output:
#Data: alpha_div
#Models:
 # m1.lme4.evennull: Normeven ~ (1 | Individual_animal) + (1 | Day) + (1 | Farm)
#m1.lme4.even: Normeven ~ Pathotype_1 + (1 | Individual_animal) + (1 | Day) + 
#  m1.lme4.even:     (1 | Farm)
#Df     AIC     BIC logLik deviance Chisq Chi Df Pr(>Chisq)
#m1.lme4.evennull  5 -880.83 -864.44 445.41  -890.83                        
#m1.lme4.even      6 -878.83 -859.16 445.41  -890.83 2e-04      1     0.9877

#Pattern
m2.lme4.even <- lmer(Normeven ~ Pattern_1 + (1|Individual_animal) + (1|Farm), 
                     data = alpha_div)
summary(m2.lme4.even)
anova(m2.lme4.even)

# calc the null
m2.lme4.even <- lmer(Normeven ~ Pattern_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m2.lme4.evennull <- lmer(Normeven ~ (1|Individual_animal) +
                           (1|Farm), data = alpha_div, REML=FALSE)
anova(m2.lme4.evennull,m2.lme4.even)

#EvNev
m3.lme4.even <- lmer(Normeven ~ EvNev_1 + (1|Individual_animal) + (1|Farm), 
                     data = alpha_div)
summary(m3.lme4.even)
anova(m3.lme4.even)
#calc the null
m3.lme4.even <- lmer(Normeven ~ EvNev_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m3.lme4.evennull <- lmer(Normeven ~ (1|Individual_animal) +
                           (1|Farm), data = alpha_div, REML=FALSE)
anova(m3.lme4.evennull,m3.lme4.even)


#RICHNESS 
#Pathotype
m1.lme4.rich <- lmer(Normrich ~ Pathotype_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div)
summary(m1.lme4.rich)
anova(m1.lme4.rich)
#calc the null
m1.lme4.rich <- lmer(Normrich ~ Pathotype_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m1.lme4.richnull <- lmer(Normrich ~ (1|Individual_animal) +
                           (1|Farm), data = alpha_div, REML=FALSE)

anova(m1.lme4.richnull,m1.lme4.rich)
#Pattern
m2.lme4.rich <- lmer(Normrich ~ Pattern_1 + (1|Individual_animal) + (1|Farm), 
                     data = alpha_div)
summary(m2.lme4.rich)
anova(m2.lme4.rich)
#calc the null
m2.lme4.rich <- lmer(Normrich ~ Pattern_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m2.lme4.richnull <- lmer(Normrich ~ (1|Individual_animal) +
                           (1|Farm), data = alpha_div, REML=FALSE)

anova(m2.lme4.richnull,m2.lme4.rich)
#EvNev
m3.lme4.rich <- lmer(Normrich ~ EvNev_1 + (1|Individual_animal) + (1|Farm), 
                     data = alpha_div)
summary(m3.lme4.rich)
anova(m3.lme4.rich)
#calc the null
m3.lme4.rich <- lmer(Normrich ~ EvNev_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m3.lme4.richnull <- lmer(Normrich ~ (1|Individual_animal) +
                           (1|Farm), data = alpha_div, REML=FALSE)

anova(m3.lme4.richnull,m3.lme4.rich)


#SHANNONS 
#Pathotype
m1.lme4.shann <- lmer(Normshann ~ Pathotype_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div)
summary(m1.lme4.shann)
anova(m1.lme4.shann)
#calc the null
m1.lme4.shann <- lmer(Normshann ~ Pathotype_1 + (1|Individual_animal) +
                       (1|Farm), data = alpha_div, REML=FALSE)
m1.lme4.shannnull <- lmer(Normshann ~ (1|Individual_animal) +
                           (1|Farm), data = alpha_div, REML=FALSE)

anova(m1.lme4.shannnull,m1.lme4.shann)

#Pattern
m2.lme4.shann <- lmer(Normshann ~ Pattern_1 + (1|Individual_animal) + (1|Farm), 
                     data = alpha_div)
summary(m2.lme4.shann)
anova(m2.lme4.shann)
#calc the null
m2.lme4.shann <- lmer(Normshann ~ Pattern_1 + (1|Individual_animal) +
                        (1|Farm), data = alpha_div, REML=FALSE)
m2.lme4.shannnull <- lmer(Normshann ~ (1|Individual_animal) +
                            (1|Farm), data = alpha_div, REML=FALSE)

anova(m2.lme4.shannnull,m2.lme4.shann)

#EvNev
m3.lme4.shann <- lmer(Normshann ~ EvNev_1 + (1|Individual_animal) + (1|Farm), 
                     data = alpha_div)
summary(m3.lme4.shann)
anova(m3.lme4.shann)
#calc the null
m3.lme4.shann <- lmer(Normshann ~ EvNev_1 + (1|Individual_animal) +
                        (1|Farm), data = alpha_div, REML=FALSE)
m3.lme4.shannnull <- lmer(Normshann ~ (1|Individual_animal) +
                            (1|Farm), data = alpha_div, REML=FALSE)

anova(m3.lme4.shannnull,m3.lme4.shann)


#Output (how do I know if it's significant? look at an f distribution table I imagine)

# going to look into evaluating p-values

library(pbkrtest)
PBmodcomp(m1.lme4.even)

# Running model with different package to see if output similar/same:
library(nlme)
m1.nlme.even <- lme(Normeven ~ Pathotype_1, random = ~ 1|c(Individual_animal, 
              Day, Farm), data = alpha_div)
# code above is not working.....

## Discussed these models in class during epi meeting on 2/16/17- Sheryl suggests just doing 
# a simple linear regression:One problem is the hierarchy of the random effects model 
#when my outcome (for pattern and ev nev) is clustered at the cow level. Sheryl is 
#saying to try a simple linear regression, making dummy variables for the 0, 1, 2 category 
#(so they are all either 0 or 1) and then run a regression using ‘cow’ as an additional 
#variable in the model to control for it. But first should just do the simple linear with 
#the dummy variables and see if it is significant. Then do the same thing with the ever never
#outcome.  Maggie Clark in the meeting thinks that the number of cows is going to be small 
#enough that adding the variable ‘cow’ might not make a difference, including the clustering 
#variable for ‘farm’ as well. For the Sample level outcome can do the same linear regression 
#for pathotype…..
library(ggplot2)
# made the dummy variables in excel in the cowmapwrichnshansnevennormed sheet. Path_1.0, Path_1.1,
# for the values for 0,1, for Pathotype_1. And Pattern_1.0, Pattern_1.1, Pattern_1.2 
# for Pattern_1 (0,1,2). Also made Evnev_1.0, Evnev_1.1 for EvNev_1
# How does the call 'lm' look at my variables? Do I need to make sure they are classified as factors?
# all my variables of interest (individual animal, farm, pattern, pathotype, evernever) should be factors,
# but I want to make sure the models interpret them correctly, so what should they be?
# They should all be factors. And the factors with two variables (0 and 1) just run as individual 
# factors not with the design variables. The pattern (0, 1, 2) can run as a factor too and R will
# use the lowest as the comparison/contrast/baseline variable for the estimates. 

# Ran all the first code from the top again. 

# changing variables to factors:

alpha_div <- alpha_div %>% mutate(Path_1.0 = as.factor(Path_1.0),
                                  Path_1.1 = as.factor(Path_1.1),
                                  Pattern_1.0 = as.factor(Pattern_1.0),
                                  Pattern_1.1 = as.factor(Pattern_1.1),
                                  Pattern_1.2 = as.factor(Pattern_1.2),
                                  Evnev_1.0 = as.factor(Evnev_1.0),
                                  Evnev_1.1 = as.factor(Evnev_1.1),
                                  Individual_animal = as.factor(Individual_animal),
                                  Farm = as.factor(Farm),
                                  Pattern_1 = as.factor(Pattern_1),
                                  EvNev_1 = as.factor(EvNev_1),
                                  Pathotype_1 = as.factor(Pathotype_1))
# Models want to build:

# Normeven, Pattern_1.0 + Pattern_1.1 + Pattern_1.2

Even_pattern1 <- lm(Normeven ~ Pattern_1, alpha_div)
summary(Even_pattern1)
anova(Even_pattern1)
plot_patteven <- alpha_div %>% select(Pattern_1, Normeven) %>%
    ggplot(aes(x=Pattern_1, y=Normeven)) + geom_point()
# now with cow as a factor
Even_pattern2 <- lm(Normeven ~ Pattern_1 +
                      Individual_animal, alpha_div)
summary(Even_pattern2)
anova(Even_pattern2)

# now with farm as a factor 
Even_pattern3 <- lm(Normeven ~ Pattern_1 + 
                      Farm, alpha_div)
summary(Even_pattern3)
anova(Even_pattern3)

# now with both as factors

Even_pattern4 <- lm(Normeven ~ Pattern_1 + 
                      Farm + Individual_animal, alpha_div)
summary(Even_pattern4)
anova(Even_pattern4)

# Normshann, Pattern_1.0 + Pattern_1.1 + Pattern_1.2
Shann_pattern1 <- lm(Normshann ~ Pattern_1, alpha_div)
summary(Shann_pattern1)
anova(Shann_pattern1)
plot_pattshann <- alpha_div %>% select(Pattern_1, Normshann) %>%
  ggplot(aes(x=Pattern_1, y=Normshann)) + geom_point()
# now with cow as a factor
Shann_pattern2 <- lm(Normshann ~ Pattern_1 + 
                      Individual_animal, alpha_div)
summary(Shann_pattern2)
anova(Shann_pattern2)

# now with farm as a factor 
Shann_pattern3 <- lm(Normshann ~ Pattern_1 + 
                      Farm, alpha_div)
summary(Shann_pattern3)
anova(Shann_pattern3)

# now with both as factors
Shann_pattern4 <- lm(Normshann ~ Pattern_1 +
                      Farm + Individual_animal, alpha_div)
summary(Shann_pattern4)
anova(Shann_pattern4)

# Normrich, Pattern_1.0 + Pattern_1.1 + Pattern_1.2
Rich_pattern1 <- lm(Normrich ~ Pattern_1, alpha_div)
summary(Rich_pattern1)
anova(Rich_pattern1)
plot_pattrich <- alpha_div %>% select(Pattern_1, Normrich) %>%
  ggplot(aes(x=Pattern_1, y=Normrich)) + geom_point()

# now with cow as a factor
Rich_pattern2 <- lm(Normrich ~ Pattern_1 +
                       Individual_animal, alpha_div)
summary(Rich_pattern2)
anova(Rich_pattern2)

# now with farm as a factor 
Rich_pattern3 <- lm(Normrich ~ Pattern_1 +
                       Farm, alpha_div)
summary(Rich_pattern3)
anova(Rich_pattern3)

# now with both as factors
Rich_pattern4 <- lm(Normrich ~ Pattern_1 +
                       Farm + Individual_animal, alpha_div)
summary(Rich_pattern4)
anova(Rich_pattern4)

# Normeven, Path_1.1 + Path_1.0
even_path1 <- lm(Normeven ~ Pathotype_1, alpha_div)
summary(even_path1)
anova(even_path1)
plot_patheven <- alpha_div %>% select(Pathotype_1, Normeven) %>%
  ggplot(aes(x=Pathotype_1, y=Normeven)) + geom_point()
plot_patheven2 <- alpha_div %>% select(Pathotype_2, Normeven) %>%
  ggplot(aes(x=Pathotype_2, y=Normeven)) + geom_point()
# now with cow as a factor
even_path2 <- lm(Normeven ~ Pathotype_1 +
                    Individual_animal, alpha_div)
summary(even_path2)
anova(even_path2)

# now with farm as a factor 
even_path3 <- lm(Normeven ~ Pathotype_1 +
                    Farm, alpha_div)
summary(even_path3)
anova(even_path3)

# now with both as factors
even_path4 <- lm(Normeven ~ Pathotype_1 +
                    Farm + Individual_animal, alpha_div)
summary(even_path4)
anova(even_path4)

# Normshann, Path_1.1 + Path_1.0
shann_path1 <- lm(Normshann ~ Pathotype_1, alpha_div)
summary(shann_path1)
anova(shann_path1)
plot_pathshann <- alpha_div %>% select(Pathotype_1, Normshann) %>%
  ggplot(aes(x=Pathotype_1, y=Normshann)) + geom_point()
plot_pathshann2 <- alpha_div %>% select(Pathotype_2, Normshann) %>%
  ggplot(aes(x=Pathotype_2, y=Normshann)) + geom_point()
# now with cow as a factor
shann_path2 <- lm(Normshann ~ Pathotype_1 + 
                   Individual_animal, alpha_div)
summary(shann_path2)
anova(shann_path2)

# now with farm as a factor 
shann_path3 <- lm(Normshann ~ Pathotype_1 +
                   Farm, alpha_div)
summary(shann_path3)
anova(shann_path3)

# now with both as factors
shann_path4 <- lm(Normshann ~ Pathotype_1 +
                   Farm + Individual_animal, alpha_div)
summary(shann_path4)
anova(shann_path4)

# Normrich, Path_1.1 + Path_1.0
rich_path1 <- lm(Normrich ~ Pathotype_1, alpha_div)
summary(rich_path1)
anova(rich_path1)
plot_pathrich <- alpha_div %>% select(Pathotype_1, Normrich) %>%
  ggplot(aes(x=Pathotype_1, y=Normrich)) + geom_point()
plot_pathrich2 <- alpha_div %>% select(Pathotype_2, Normrich) %>%
  ggplot(aes(x=Pathotype_2, y=Normrich)) + geom_point()
# now with cow as a factor
rich_path2 <- lm(Normrich ~ Pathotype_1 + 
                    Individual_animal, alpha_div)
summary(rich_path2)
anova(rich_path2)

# now with farm as a factor 
rich_path3 <- lm(Normrich ~ Pathotype_1 + 
                    Farm, alpha_div)
summary(rich_path3)
anova(rich_path3)

# now with both as factors
rich_path4 <- lm(Normrich ~ Pathotype_1 + 
                    Farm + Individual_animal, alpha_div)
summary(rich_path4)
anova(rich_path4)


# Normeven, Evnev_1.0 + Evnev_1.1
even_evnev1 <- lm(Normeven ~ EvNev_1, alpha_div)
summary(even_evnev1)
anova(even_evnev1)
plot_evneveven <- alpha_div %>% select(EvNev_1, Normeven) %>%
  ggplot(aes(x=EvNev_1, y=Normeven)) + geom_point()

# now with cow as a factor
even_evnev2 <- lm(Normeven ~ EvNev_1 + 
                     Individual_animal, alpha_div)
summary(even_evnev2)
anova(even_evnev2)

# now with farm as a factor 
even_evnev3 <- lm(Normeven ~ EvNev_1 + 
                     Farm, alpha_div)
summary(even_evnev3)
anova(even_evnev3)

# now with both as factors
even_evnev4 <- lm(Normeven ~ EvNev_1 + 
                     Farm + Individual_animal, alpha_div)
summary(even_evnev4)
anova(even_evnev4)


# Normshann, Evnev_1.0 + Evnev_1.1
Shann_evnev1 <- lm(Normshann ~ EvNev_1, alpha_div)
summary(Shann_evnev1)
anova(Shann_evnev1)
plot_evnevshann <- alpha_div %>% select(EvNev_1, Normshann) %>%
  ggplot(aes(x=EvNev_1, y=Normshann)) + geom_point()

# now with cow as a factor
Shann_evnev2 <- lm(Normshann ~ EvNev_1 + 
                    Individual_animal, alpha_div)
summary(Shann_evnev2)
anova(Shann_evnev2)

# now with farm as a factor 
Shann_evnev3 <- lm(Normshann ~ EvNev_1 + 
                    Farm, alpha_div)
summary(Shann_evnev3)
anova(Shann_evnev3)

# now with both as factors
Shann_evnev4 <- lm(Normshann ~ EvNev_1 + 
                    Farm + Individual_animal, alpha_div)
summary(Shann_evnev4)
anova(Shann_evnev4)

# Normrich, Evnev_1.0 + Evnev_1.1
Rich_evnev1 <- lm(Normrich ~ EvNev_1, alpha_div)
summary(Rich_evnev1)
anova(Rich_evnev1)
plot_evnevrich <- alpha_div %>% select(EvNev_1, Normrich) %>%
  ggplot(aes(x=EvNev_1, y=Normrich)) + geom_point()

# now with cow as a factor
Rich_evnev2 <- lm(Normrich ~ EvNev_1 + 
                      Individual_animal, alpha_div)
summary(Rich_evnev2)
anova(Rich_evnev2)

# now with farm as a factor 
Rich_evnev3 <- lm(Normrich ~ EvNev_1 + 
                      Farm, alpha_div)
summary(Rich_evnev3)
anova(Rich_evnev3)

# now with both as factors
Rich_evnev4 <- lm(Normrich ~ EvNev_1 + 
                      Farm + Individual_animal, alpha_div)
summary(Rich_evnev4)
anova(Rich_evnev4)

## seeing if the 'controlling' variables are significantly associated with any 
## alpha div measures
#Day
Rich_day1 <- lm(Normrich ~ Day, alpha_div)
summary(Rich_day1)
anova(Rich_day1)
Shann_day1 <- lm(Normshann ~ Day, alpha_div)
summary(Shann_day1)
anova(Shann_day1)
Even_day1 <- lm(Normeven ~ Day, alpha_div)
summary(Even_day1)
anova(Even_day1)

#Farm
Rich_farm1 <- lm(Normrich ~ Farm, alpha_div)
summary(Rich_farm1)
anova(Rich_farm1)
Shann_farm1 <- lm(Normshann ~ Farm, alpha_div)
summary(Shann_farm1)
anova(Shann_farm1)
Even_farm1 <- lm(Normeven ~ Farm, alpha_div)
summary(Even_farm1)
anova(Even_farm1)
 #Individual_animal
Rich_animal1 <- lm(Normrich ~ Individual_animal, alpha_div)
summary(Rich_animal1)
anova(Rich_animal1)
Shann_animal1 <- lm(Normshann ~ Individual_animal, alpha_div)
summary(Shann_animal1)
anova(Shann_animal1)
Even_animal1 <- lm(Normeven ~ Individual_animal, alpha_div)
summary(Even_animal1)
anova(Even_animal1)

######### try other things
Even_pattern1 <- lm(Normeven ~ Pattern_1, alpha_div)
summary(Even_pattern1)
anova(Even_pattern1)

even_evnev1 <- lm(Normeven ~ EvNev_1, alpha_div)
summary(even_evnev1)
anova(even_evnev1)

even_path1 <- lm(Normeven ~ Pathotype_1, alpha_div)
summary(even_path1)
anova(even_path1)

even_path2 <-lm(Normeven ~ Pathotype_1 + Individual_animal, alpha_div)
summary(even_path2)
anova(even_path2)
