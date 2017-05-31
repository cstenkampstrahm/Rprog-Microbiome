# Going to omit cow #30 from dataset to see if alpha diversity models have diff
# outcomes

library("xlsx")
library("lme4")
library("nnet") 
library("gee") 
library("tidyverse")

avg_alpha_div <- read.xlsx("excel sheets/Cow_map_w_no30.xlsx", 2)
alpha_div <- read.xlsx("excel sheets/Cow_map_w_no30.xlsx", 1)
alpha_div <- mutate(alpha_div, "Normrich" = as.numeric(Normrich), 
                    "Normshann" = as.numeric(Normshann),
                    "Normeven" = as.numeric(Normeven))
avg_alpha_div <- mutate(avg_alpha_div, "Avgrich" = as.numeric(Avgrich), 
                        "Avgshann" = as.numeric(Avgshann),
                        "Avgeven" = as.numeric(Avgeven))
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
avg_alpha_div <- avg_alpha_div %>% mutate(Path_1.0 = as.factor(Path_1.0),
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


#PATTERN
m3.1 <- multinom(Pattern_1 ~ Avgrich, data = avg_alpha_div)
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

#PATHOTYPE
m1.1 <- glmer(Pathotype_1 ~ Normrich + (1|Individual_animal), data = alpha_div,
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


#EverNever
m2.1 <- glm(EvNev_1 ~ Avgrich, data = avg_alpha_div, family = binomial)
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
