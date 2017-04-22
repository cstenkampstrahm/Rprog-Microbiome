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
                           "Farm" = as.factor(Farm))
avg_alpha_div_scaled <- mutate(avg_alpha_div_scaled, "Avgrich" = as.numeric(Avgrich), 
                               "Avgshann" = as.numeric(Avgshann),
                               "Avgeven" = as.numeric(Avgeven),
                               "Avgscaledrich" = as.numeric(Avgscaledrich),
                               "Newscaleavgrich" = as.numeric(Newscaleavgrich),
                               "DIM" = as.numeric(DIM),
                               "Parity_1" = as.factor(Parity_1),
                               "Disease" = as.factor(Disease),
                               "Farm" = as.factor(Farm))

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
m1.1farm <- update(m1.1, .~. + Farm)
summary(m1.1farm)
m1.1both <- update(m1.1, .~. + Farm + Parity_1)
summary(m1.1both)


m1.2 <- glmer(Pathotype_1 ~ Normshann + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))

summary(m1.2)

m1.2Parity <-  update(m1.2, .~. + Parity_1)
summary(m1.2Parity)
m1.2farm <- update(m1.2, .~. + Farm)
summary(m1.2farm)
m1.2both <- update(m1.2, .~. + Farm + Parity_1)
summary(m1.2both)

m1.3 <- glmer(Pathotype_1 ~ Normeven + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.3)

m1.3Parity <-  update(m1.3, .~. + Parity_1)
summary(m1.3Parity)
m1.3farm <- update(m1.3, .~. + Farm)
summary(m1.3farm)
m1.3both <- update(m1.3, .~. + Farm + Parity_1)
summary(m1.3both)

#EverNever
m2.1 <- glm(EvNev_1 ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1)

m2.1Parity <-  update(m2.1, .~. + Parity_1)
summary(m2.1Parity)
m2.1farm <- update(m2.1, .~. + Farm)
summary(m2.1farm)
m2.1both <- update(m2.1, .~. + Farm + Parity_1)
summary(m2.1both)


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
m2.2farm <- update(m2.2, .~. + Farm)
summary(m2.2farm)
m2.2both <- update(m2.2, .~. + Farm + Parity_1)
summary(m2.2both)

m2.3 <- glm(EvNev_1 ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(m2.3)

m2.3Parity <-  update(m2.3, .~. + Parity_1)
summary(m2.3Parity)
m2.3farm <- update(m2.3, .~. + Farm)
summary(m2.3farm)
m2.3both <- update(m2.3, .~. + Farm + Parity_1)
summary(m2.3both)

#PATTERN
m3.1 <- multinom(Pattern_1 ~ Avgrich, data = avg_alpha_div_scaled)
summary(m3.1)
z <- summary(m3.1)$coefficients/summary(m3.1)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p

m3.1Parity <-  update(m3.1, .~. + Parity_1)
summary(m3.1Parity)
m3.1farm <- update(m3.1, .~. + Farm)
summary(m3.1farm)
m3.1both <- update(m3.1, .~. + Farm + Parity_1)
summary(m3.1both)


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


m3.1.00 <- multinom(Pattern_1, data = avg_alpha_div_scaled)

m3.2 <- multinom(Pattern_1 ~ Avgshann, data = avg_alpha_div_scaled)
summary(m3.2)
z <- summary(m3.2)$coefficients/summary(m3.2)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.2Parity <-  update(m3.2, .~. + Parity_1)
summary(m3.2Parity)
m3.2farm <- update(m3.2, .~. + Farm)
summary(m3.2farm)
m3.2both <- update(m3.2, .~. + Farm + Parity_1)
summary(m3.2both)

m3.3 <- multinom(Pattern_1 ~Avgeven, data = avg_alpha_div_scaled)
summary(m3.3)
z <- summary(m3.3)$coefficients/summary(m3.3)$standard.errors
z
p <- (1-pnorm(abs(z),0,1)) * 2
p

m3.3Parity <-  update(m3.3, .~. + Parity_1)
summary(m3.3Parity)
m3.3farm <- update(m3.3, .~. + Farm)
summary(m3.3farm)
m3.3both <- update(m3.3, .~. + Farm + Parity_1)
summary(m3.3both)

# looking at 3x6 table of farm, pattern and parity. for outcome descriptive purposes

threebysix <- alpha_div_scaled %>% group_by(Pattern_1, Farm, Parity_1) %>% summarize(n = n())



threebysix1 <- avg_alpha_div_scaled %>% group_by(Pattern_1, Farm, Parity_1) %>% summarize(n = n())


