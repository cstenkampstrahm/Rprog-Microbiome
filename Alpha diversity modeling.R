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

alpha_div <- read.xlsx("Cow_map_wrichnshansnevennormed.xlsx", 1)
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