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
                           "Scaledrich" = as.numeric(Scaledrich),
                           "DIM" = as.numeric(DIM),
                           "Parity_1" = as.factor(Parity_1),
                           "Disease" = as.factor(Disease),
                           "Farm" = as.factor(Farm),
                           "Treatment" = as.factor(Treatment))
avg_alpha_div_scaled <- mutate(avg_alpha_div_scaled, "Avgrich" = as.numeric(Avgrich), 
                               "Avgshann" = as.numeric(Avgshann),
                               "Avgeven" = as.numeric(Avgeven),
                               "Avgscaledrich" = as.numeric(Avgscaledrich),
                               "Newscaleavgrich" = as.numeric(Newscaleavgrich),
                               "DIM" = as.numeric(DIM),
                               "Parity_1" = as.factor(Parity_1),
                               "Disease" = as.factor(Disease),
                               "Farm" = as.factor(Farm),
                               "Treatment" = as.factor(Treatment))

alpha_div_scaled <- alpha_div_scaled %>% mutate(Path_1.0 = as.factor(Path_1.0),
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

avg_alpha_div_scaled <- avg_alpha_div_scaled %>% mutate(Path_1.0 = as.factor(Path_1.0),
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


#PATHOTYPE
m1.1.0 <- glmer(Pathotype_1 ~ Scaledrich + (1|Individual_animal), data = alpha_div_scaled,
                family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.1.0)
m1.1.0Parity <-  update(m1.1.0, .~. + Parity_1)
summary(m1.1.0Parity)
m1.1.0farm <- update(m1.1.0, .~. + Farm)
summary(m1.1.0farm)
m1.1.0both <- update(m1.1.0, .~. + Farm + Parity_1)
summary(m1.1.0both)

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

#EverNever
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


m2.1.1 <- glm(EvNev_1 ~ Newscaleavgrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1.1)

m2.1.1Parity <-  update(m2.1.1, .~. + Parity_1)
summary(m2.1.1Parity)
m2.1.1farm <- update(m2.1.1, .~. + Farm)
summary(m2.1.1farm)
m2.1.1both <- update(m2.1.1, .~. + Farm + Parity_1)
summary(m2.1.1both)

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

#PATTERN
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


m3.1.0 <- multinom(Pattern_1 ~ Newscaleavgrich, data = avg_alpha_div_scaled)
summary(m3.1.0)
z <- summary(m3.1.0)$coefficients/summary(m3.1.0)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.1.0Parity <-  update(m3.1.0, .~. + Parity_1)
summary(m3.1.0Parity)
m3.1.0farm <- update(m3.1.0, .~. + Farm)
summary(m3.1.0farm)
m3.1.0both <- update(m3.1.0, .~. + Farm + Parity_1)
summary(m3.1.0both)


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

# looking at 3x6 table of farm, pattern and parity. for outcome descriptive purposes

threebysix <- alpha_div_scaled %>% group_by(Pattern_1, Farm, Parity_1) %>% summarize(n = n())



threebysix1 <- avg_alpha_div_scaled %>% group_by(Pattern_1, Farm, Parity_1) %>% summarize(n = n())

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
fisher.test(table_1$Pathotype_1, table_1$Parity_1)  
fisher.test(table_1$Parity_1, table_1$Pathotype_1)

# Treatment
pattern_treament <- table_2 %>% group_by(Pattern_1, Treatment) %>% summarise(treatnval = n())
fisher.test(table_2$Pattern_1, table_2$Treatment)

evnev_treatment <- table_2 %>% group_by(EvNev_1, Treatment) %>% summarise(treatnval = n())
fisher.test(table_2$EvNev_1, table_2$Treatment)

pathotype_treatment <- table_1 %>% group_by(Pathotype_1, Treatment) %>% summarise(treatnval = n())
fisher.test(table_1$Pathotype_1, table_1$Treatment)  

