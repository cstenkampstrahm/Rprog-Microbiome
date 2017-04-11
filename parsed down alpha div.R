alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
avg_alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)

alpha_div_scaled <- mutate(alpha_div_scaled, "Normrich" = as.numeric(Normrich), 
                           "Normshann" = as.numeric(Normshann),
                           "Normeven" = as.numeric(Normeven),
                           "Scaledrich" = as.numeric(Scaledrich))
avg_alpha_div_scaled <- mutate(avg_alpha_div_scaled, "Avgrich" = as.numeric(Avgrich), 
                               "Avgshann" = as.numeric(Avgshann),
                               "Avgeven" = as.numeric(Avgeven),
                               "Avgscaledrich" = as.numeric(Avgscaledrich),
                               "Newscaleavgrich" = as.numeric(Newscaleavgrich))

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

m1.2 <- glmer(Pathotype_1 ~ Normshann + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))

summary(m1.2)

m1.3 <- glmer(Pathotype_1 ~ Normeven + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.3)

#EverNever
m2.1 <- glm(EvNev_1 ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1)

m2.1.1 <- glm(EvNev_1 ~ Newscaleavgrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1.1)

m2.2 <- glm(EvNev_1 ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(m2.2)

m2.3 <- glm(EvNev_1 ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(m2.3)

#PATTERN
m3.1 <- multinom(Pattern_1 ~ Avgrich, data = avg_alpha_div_scaled)
summary(m3.1)
z <- summary(m3.1)$coefficients/summary(m3.1)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p


m3.1.0 <- multinom(Pattern_1 ~ Newscaleavgrich, data = avg_alpha_div_scaled)
summary(m3.1.0)
z <- summary(m3.1.0)$coefficients/summary(m3.1.0)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.1.00 <- multinom(Pattern_1, data = avg_alpha_div_scaled)

m3.2 <- multinom(Pattern_1 ~ Avgshann, data = avg_alpha_div_scaled)
summary(m3.2)
z <- summary(m3.2)$coefficients/summary(m3.2)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.3 <- multinom(Pattern_1 ~Avgeven, data = avg_alpha_div_scaled)
summary(m3.3)
z <- summary(m3.3)$coefficients/summary(m3.3)$standard.errors
z
p <- (1-pnorm(abs(z),0,1)) * 2
p









