### Final models and code used for epidemiologic modeling of cow 16s data

##### FROM parsed down alpha div. R
library("xlsx")
library("lme4")
library("nnet") 
library("gee")
library("tidyverse")
alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
avg_alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)

alpha_div_scaled <- mutate(alpha_div_scaled, "Normrich" = as.numeric(Normrich), 
                           "Normshann" = as.numeric(Normshann),
                           "Normeven" = as.numeric(Normeven),
                           "DIM" = as.numeric(DIM),
                           "Parity_1" = as.factor(Parity_1),
                           "Disease" = as.factor(Disease),
                           "Farm" = as.factor(Farm),
                           "Treatment" = as.factor(Treatment))
avg_alpha_div_scaled <- mutate(avg_alpha_div_scaled, "Avgrich" = as.numeric(Avgrich), 
                               "Avgshann" = as.numeric(Avgshann),
                               "Avgeven" = as.numeric(Avgeven),
                               "DIM" = as.numeric(DIM),
                               "Parity_1" = as.factor(Parity_1),
                               "Disease" = as.factor(Disease),
                               "Farm" = as.factor(Farm),
                               "Treatment" = as.factor(Treatment))
alpha_div_scaled <- alpha_div_scaled %>% mutate(Individual_animal = as.factor(Individual_animal),
                                                Farm = as.factor(Farm),
                                                Pattern_1 = as.factor(Pattern_1),
                                                EvNev_1 = as.factor(EvNev_1),
                                                Pathotype_1 = as.factor(Pathotype_1))
avg_alpha_div_scaled <- avg_alpha_div_scaled %>% mutate(Individual_animal = as.factor(Individual_animal),
                                                        Farm = as.factor(Farm),
                                                        Pattern_1 = as.factor(Pattern_1),
                                                        EvNev_1 = as.factor(EvNev_1),
                                                        Pathotype_1 = as.factor(Pathotype_1))
#####PATHOTYPE
# richness
m1.1 <- glmer(Pathotype_1 ~ Normrich + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.1)

m1.1Parity <-  update(m1.1, .~. + Parity_1)
summary(m1.1Parity)
m1.1treat <- update(m1.1, .~. + Treatment)
summary(m1.1treat)
m1.1farm <- update(m1.1, .~. + Farm)
summary(m1.1farm)
m1.1treatparity <- update(m1.1, .~. + Treatment + Parity_1)
summary(m1.1treatparity)
m1.1treatfarm <- update(m1.1, .~. + Treatment + Farm)
summary(m1.1treatfarm)
m1.1both <- update(m1.1, .~. + Farm + Parity_1)
summary(m1.1both)
m1.1all <- update(m1.1, .~. + Farm + Parity_1 + Treatment)
summary(m1.1all)

# shannons
m1.2 <- glmer(Pathotype_1 ~ Normshann + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))

summary(m1.2)

m1.2Parity <-  update(m1.2, .~. + Parity_1)
summary(m1.2Parity)
m1.2treat <- update(m1.2, .~. + Treatment)
summary(m1.2treat)
m1.2farm <- update(m1.2, .~. + Farm)
summary(m1.2farm)
m1.2treatparity <- update(m1.2, .~. + Treatment + Parity_1)
summary(m1.2treatparity)
m1.2treatfarm <- update(m1.2, .~. + Treatment + Farm)
summary(m1.2treatfarm)
m1.2both <- update(m1.2, .~. + Farm + Parity_1)
summary(m1.2both)
m1.2all <- update(m1.2, .~. + Farm + Parity_1 + Treatment)
summary(m1.2all)

# evenness
m1.3 <- glmer(Pathotype_1 ~ Normeven + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.3)

m1.3Parity <-  update(m1.3, .~. + Parity_1)
summary(m1.3Parity)
m1.3treat <- update(m1.3, .~. + Treatment)
summary(m1.3treat)
m1.3farm <- update(m1.3, .~. + Farm)
summary(m1.3farm)
m1.3treatparity <- update(m1.3, .~. + Treatment + Parity_1)
summary(m1.3treatparity)
m1.3treatfarm <- update(m1.3, .~. + Treatment + Farm)
summary(m1.3treatfarm)
m1.3both <- update(m1.3, .~. + Farm + Parity_1)
summary(m1.3both)
m1.3all <- update(m1.3, .~. + Farm + Parity_1 + Treatment)
summary(m1.3all)

#### Ever_Never
# richness
m2.1 <- glm(EvNev_1 ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1)

m2.1Parity <-  update(m2.1, .~. + Parity_1)
summary(m2.1Parity)
m2.1treat <- update(m2.1, .~. + Treatment)
summary(m2.1treat)
m2.1farm <- update(m2.1, .~. + Farm)
summary(m2.1farm)
m2.1treatparity <- update(m2.1, .~. + Treatment + Parity_1)
summary(m2.1treatparity)
m2.1treatfarm <- update(m2.1, .~. + Treatment + Farm)
summary(m2.1treatfarm)
m2.1both <- update(m2.1, .~. + Farm + Parity_1)
summary(m2.1both)
m2.1all <- update(m2.1, .~. + Farm + Parity_1 + Treatment)
summary(m2.1all)

# shannons
m2.2 <- glm(EvNev_1 ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(m2.2)

m2.2Parity <-  update(m2.2, .~. + Parity_1)
summary(m2.2Parity)
m2.2treat <- update(m2.2, .~. + Treatment)
summary(m2.2treat)
m2.2farm <- update(m2.2, .~. + Farm)
summary(m2.2farm)
m2.2treatparity <- update(m2.2, .~. + Treatment + Parity_1)
summary(m2.2treatparity)
m2.2treatfarm <- update(m2.2, .~. + Treatment + Farm)
summary(m2.2treatfarm)
m2.2both <- update(m2.2, .~. + Farm + Parity_1)
summary(m2.2both)
m2.2all <- update(m2.2, .~. + Farm + Parity_1 + Treatment)
summary(m2.2all)

# evenness
m2.3 <- glm(EvNev_1 ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(m2.3)

m2.3Parity <-  update(m2.3, .~. + Parity_1)
summary(m2.3Parity)
m2.3treat <- update(m2.3, .~. + Treatment)
summary(m2.3treat)
m2.3farm <- update(m2.3, .~. + Farm)
summary(m2.3farm)
m2.3treatparity <- update(m2.3, .~. + Treatment + Parity_1)
summary(m2.3treatparity)
m2.3treatfarm <- update(m2.3, .~. + Treatment + Farm)
summary(m2.3treatfarm)
m2.3both <- update(m2.3, .~. + Farm + Parity_1)
summary(m2.3both)
m2.3all <- update(m2.3, .~. + Farm + Parity_1 + Treatment)
summary(m2.3all)

#### Pattern

#richness
m3.1 <- multinom(Pattern_1 ~ Avgrich, data = avg_alpha_div_scaled)
summary(m3.1)
z <- summary(m3.1)$coefficients/summary(m3.1)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p

m3.1Parity <-  update(m3.1, .~. + Parity_1)
summary(m3.1Parity)
m3.1treat <- update(m3.1, .~. + Treatment)
summary(m3.1treat)
m3.1farm <- update(m3.1, .~. + Farm)
summary(m3.1farm)
m3.1treatparity <- update(m3.1, .~. + Treatment + Parity_1)
summary(m3.1treatparity)
m3.1treatfarm <- update(m3.1, .~. + Treatment + Farm)
summary(m3.1treatfarm)
m3.1both <- update(m3.1, .~. + Farm + Parity_1)
summary(m3.1both)
m3.1all <- update(m3.1, .~. + Farm + Parity_1 + Treatment)
summary(m3.1all)
 
# shannons
m3.2 <- multinom(Pattern_1 ~ Avgshann, data = avg_alpha_div_scaled)
summary(m3.2)
z <- summary(m3.2)$coefficients/summary(m3.2)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.2Parity <-  update(m3.2, .~. + Parity_1)
summary(m3.2Parity)
m3.2treat <- update(m3.2, .~. + Treatment)
summary(m3.2treat)
m3.2farm <- update(m3.2, .~. + Farm)
summary(m3.2farm)
m3.2treatparity <- update(m3.2, .~. + Treatment + Parity_1)
summary(m3.2treatparity)
m3.2treatfarm <- update(m3.2, .~. + Treatment + Farm)
summary(m3.2treatfarm)
m3.2both <- update(m3.2, .~. + Farm + Parity_1)
summary(m3.2both)
m3.2all <- update(m3.2, .~. + Farm + Parity_1 + Treatment)
summary(m3.2all)

# evenness
m3.3 <- multinom(Pattern_1 ~ Avgeven, data = avg_alpha_div_scaled)
summary(m3.3)
z <- summary(m3.3)$coefficients/summary(m3.3)$standard.errors
z
p <- (1-pnorm(abs(z),0,1)) * 2
p

m3.3Parity <-  update(m3.3, .~. + Parity_1)
summary(m3.3Parity)
m3.3treat <- update(m3.3, .~. + Treatment)
summary(m3.3treat)
m3.3farm <- update(m3.3, .~. + Farm)
summary(m3.3farm)
m3.3treatparity <- update(m3.3, .~. + Treatment + Parity_1)
summary(m3.3treatparity)
m3.3treatfarm <- update(m3.3, .~. + Treatment + Farm)
summary(m3.3treatfarm)
m3.3both <- update(m3.3, .~. + Farm + Parity_1)
summary(m3.3both)
m3.3all <- update(m3.3, .~. + Farm + Parity_1 + Treatment)
summary(m3.3all)

###### making table of confounders by O157 metrics
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
fisher.test(table_1$Pathotype_1, table_1$Parity_1)  
fisher.test(table_1$Parity_1, table_1$Pathotype_1)

# Treatment
pattern_treament <- table_2 %>% group_by(Pattern_1, Treatment) %>% summarise(treatnval = n())
fisher.test(table_2$Pattern_1, table_2$Treatment)

evnev_treatment <- table_2 %>% group_by(EvNev_1, Treatment) %>% summarise(treatnval = n())
fisher.test(table_2$EvNev_1, table_2$Treatment)

pathotype_treatment <- table_1 %>% group_by(Pathotype_1, Treatment) %>% summarise(treatnval = n())
fisher.test(table_1$Pathotype_1, table_1$Treatment)  

### FROM Alpha_div_confounders.R

# DIM 
# (mixed linear with cow as random for DIM, and mixed logistic with cow as random for DIM_1)

confmod3 <- lmer(DIM ~ Normrich + (1|Individual_animal), data = alpha_div_scaled)
anova(confmod3, nullmodDIM)
summary(confmod3)

confmod6 <- lmer(DIM ~ Normshann + (1|Individual_animal), data = alpha_div_scaled)
anova(confmod6, nullmodDIM)
summary(confmod6)

confmod8 <- lmer(DIM ~ Normeven + (1|Individual_animal), data = alpha_div_scaled)
anova(confmod8, nullmodDIM)
summary(confmod8)

# Parity_1

confmod10 <- multinom(Parity_1 ~ Avgrich, data = avg_alpha_div_scaled)
summary(confmod10)
z <- summary(confmod10)$coefficients/summary(confmod10)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 

confmod12 <- multinom(Parity_1 ~ Avgshann, data = avg_alpha_div_scaled)
summary(confmod12)
z <- summary(confmod12)$coefficients/summary(confmod12)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 

confmod13 <- multinom(Parity_1 ~ Avgeven, data = avg_alpha_div_scaled)
summary(confmod13)
z <- summary(confmod13)$coefficients/summary(confmod13)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 

# Disease
confmod19 <- glm(Disease ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod19)
confmod20 <- glm(Disease ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(confmod20)
confmod21 <- glm(Disease ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(confmod21)

# Farm
confmod23 <- glm(Farm ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod23)
confmod24 <- glm(Farm ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(confmod24)
confmod25 <- glm(Farm ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(confmod25)

# Treatment
confmod26 <- multinom(Treatment ~ Avgrich, data = avg_alpha_div_scaled)
summary(confmod26)
z <- summary(confmod26)$coefficients/summary(confmod26)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 

confmod28 <- multinom(Treatment ~ Avgshann, data = avg_alpha_div_scaled)
summary(confmod28)
z <- summary(confmod28)$coefficients/summary(confmod28)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 

confmod29 <- multinom(Treatment ~ Avgeven, data = avg_alpha_div_scaled)
summary(confmod29)
z <- summary(confmod29)$coefficients/summary(confmod29)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 

## to get more meaningful values, can calculate the IQR of diversity values 
# (between 1st and 3rd quartile) and multiply the B coefficients by that:
UQR_normrich <- IQR(alpha_div_scaled$Normrich)
IQR_normrich  
# 1255
IQR_avgrich <- IQR(avg_alpha_div_scaled$Avgrich)
IQR_avgrich 
# 797.65
IQR_normshann <- IQR(alpha_div_scaled$Normshann)
IQR_normshann
# 0.4804946
IQR_avgshann <- IQR(avg_alpha_div_scaled$Avgshann)
IQR_avgshann
# 0.3716254
IQR_normeven  <- IQR(alpha_div_scaled$Normeven)
IQR_normeven
# 0.0366165
IQR_avgeven <- IQR(avg_alpha_div_scaled$Avgeven)
IQR_avgeven
# 0.03348146

# the only model that is significant is one that looks for association between 
# multiday shedder and richness. Look at spread of avgrich values to make sure 
# no outliers:
avgrichspread <- avg_alpha_div_scaled %>% group_by(Pattern_1) %>% summarize(mean(Avgrich), 
                                                                            min(Avgrich), 
                                                                            max(Avgrich))
avgrichspread
# # A tibble: 3 × 4
#       Pattern_1 `mean(Avgrich)` `min(Avgrich)` `max(Avgrich)`
#       <fctr>        <dbl>          <dbl>          <dbl>
#1         0        4186.431         3359.8         5180.0
#2         1        4252.625         2877.4         7644.4
#3         2        3842.000         2440.0         5000.4

# others
avgshannspread <- avg_alpha_div_scaled %>% group_by(Pattern_1) %>% summarize(mean(Avgshann), 
                                                                             min(Avgshann), 
                                                                             max(Avgshann))
avgshannspread

avgevenspread <- avg_alpha_div_scaled %>% group_by(Pattern_1) %>% summarize(mean(Avgeven), 
                                                                            min(Avgeven), 
                                                                            max(Avgeven))
avgevenspread

