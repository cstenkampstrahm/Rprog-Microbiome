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
# normally distributed
qqnorm(alpha_div$Normeven, ylab= "Evenness")
qqnorm(alpha_div$Normrich, ylab= "Richness")
qqnorm(alpha_div$Normshann, ylab= "Shannons")

# will try to log transform??
alpha_div$Normrich1 <-log10(10*alpha_div$Normrich)
shapiro.test(alpha_div$Normrich1)

# data:  alpha_div$Normrich1
# W = 0.79537, p-value = 2.699e-15

alpha_div$Normshann1 <- log10(10*alpha_div$Normshann)
shapiro.test(alpha_div$Normshann1)

# data:  alpha_div$Normshann1
# W = 0.94952, p-value = 2.11e-06

m1.lme4.even <- lmer(Normeven ~ Pathotype_1 + (1|Individual_animal) +
                     (1|Day) + (1|Farm), data = alpha_div)
summary(m1.lme4.even)
anova(m1.lme4.even)
#Output (how do I know if it's significant? look at an f distribution table I imagine)
#Analysis of Variance Table
#Df     Sum Sq    Mean Sq F value
#Pathotype_1  1 2.5998e-07 2.5998e-07   6e-04
# going to look into evaluating p-values
library(pbkrtest)
PBmodcomp(m1.lme4.even)

# Running model with different package to see if output similar/same:
library(nlme)
m1.nlme.even <- lme(Normeven ~ Pathotype_1, random = ~ 1|c(Individual_animal, 
              Day, Farm), data = alpha_div)
# code above is not working.....