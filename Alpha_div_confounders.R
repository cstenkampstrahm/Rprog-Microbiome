##### Alpha div looking at potential confounding variables: ####

alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
avg_alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)

alpha_div_scaled <- mutate(alpha_div_scaled, "Normrich" = as.numeric(Normrich), 
                           "Normshann" = as.numeric(Normshann),
                           "Normeven" = as.numeric(Normeven),
                           "Scaledrich" = as.numeric(Scaledrich),
                           "DIM" = as.numeric(DIM),
                           "DIM_1" = as.factor(DIM_1),
                           "Parity" = as.factor(Parity),
                          "Disease" = as.factor(Disease),
                           "Farm" = as.factor(Farm))
avg_alpha_div_scaled <- mutate(avg_alpha_div_scaled, "Avgrich" = as.numeric(Avgrich), 
                               "Avgshann" = as.numeric(Avgshann),
                               "Avgeven" = as.numeric(Avgeven),
                               "Avgscaledrich" = as.numeric(Avgscaledrich),
                               "Newscaleavgrich" = as.numeric(Newscaleavgrich),
                               "DIM" = as.numeric(DIM),
                               "DIM_1" = as.factor(DIM_1),
                               "Parity" = as.factor(Parity),
                               "Parity_1" = as.factor(Parity_1),
                               "Disease" = as.factor(Disease),
                               "Farm" = as.factor(Farm))

# DIM 
# (mixed linear with cow as random for DIM, and mixed logistic with cow as random for DIM_1)
confmod1 <- lmer(DIM ~ Scaledrich + (1|Individual_animal), data = alpha_div_scaled)
nullmodDIM <- lmer(DIM ~ (1|Individual_animal), data = alpha_div_scaled)
anova(confmod1, nullmodDIM)
summary(confmod1)

confmod3 <- lmer(DIM ~ Normrich + (1|Individual_animal), data = alpha_div_scaled)
anova(confmod3, nullmodDIM)
summary(confmod3)
# got error when not using scaled that some variables are on very different scales again

confmod4 <- glmer(DIM_1 ~ Scaledrich + (1|Individual_animal), data = alpha_div_scaled, 
                  family = binomial)
summary(confmod4)

confmod5 <- glmer(DIM_1 ~ Normrich + (1|Individual_animal), data = alpha_div_scaled, 
                  family = binomial)
summary(confmod5)
# got error when not using scaled that some variables are on very different scales
confmod6 <- lmer(DIM ~ Normshann + (1|Individual_animal), data = alpha_div_scaled)
anova(confmod6, nullmodDIM)
summary(confmod6)

confmod7 <- glmer(DIM_1 ~ Normshann + (1|Individual_animal), data = alpha_div_scaled, 
                  family = binomial)
summary(confmod7)

confmod8 <- lmer(DIM ~ Normeven + (1|Individual_animal), data = alpha_div_scaled)
anova(confmod8, nullmodDIM)
summary(confmod8)

confmod9 <- glmer(DIM_1 ~ Normeven + (1|Individual_animal), data = alpha_div_scaled, 
                  family = binomial)
summary(confmod9)


# Parity
summary(avg_alpha_div_scaled$Parity)
#  1  2  3  4  5  6 
# 19  7  8  3  1  2 

# Parity_1

confmod10 <- multinom(Parity_1 ~ Avgrich, data = avg_alpha_div_scaled)
summary(confmod10)
z <- summary(confmod10)$coefficients/summary(confmod10)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p 


confmod11 <- multinom(Parity_1 ~ Newscaleavgrich, data = avg_alpha_div_scaled)
summary(confmod11)
z <- summary(confmod11)$coefficients/summary(confmod11)$standard.errors
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

# Parity_2

confmod14 <- glm(Parity_2 ~ Avgscaledrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod14)
confmod15 <- glm(Parity_2 ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod15)
confmod16 <- glm(Parity_2 ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(confmod16)
confmod17 <- glm(Parity_2 ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(confmod17)

# Disease

confmod18 <- glm(Disease ~ Newscaleavgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod18)
confmod19 <- glm(Disease ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod19)
confmod20 <- glm(Disease ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(confmod20)
confmod21 <- glm(Disease ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(confmod21)

# Farm

confmod22 <- glm(Farm ~ Newscaleavgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod22)
confmod23 <- glm(Farm ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(confmod23)
confmod24 <- glm(Farm ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(confmod24)
confmod25 <- glm(Farm ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(confmod25)



