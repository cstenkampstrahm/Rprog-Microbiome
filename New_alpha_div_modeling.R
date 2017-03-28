## Met with Sheryl on 3/21- decided tht we would use the alpha diversity measures
# as the independent variable, and the O157 metrics as the dependent ones. Will use
# multinomial regression with pattern, logistic reression with the others. For contr
# for animal, going to use the average value for the cow level variables. Going to 
# control for the effect of animal with the pathotype metric using generalized estimating 
# estimating equations (or can continue to use random effects....)
# made a new excel sheet (2nd tab in existing cow_map_wrichnshnsnevennormed.xlsl) 
# with averages by cow for evenness, pathotype and evnev

library("xlsx")
library("lme4")
library("nnet") # multinomial modeling
library("gee") # if using generalized estimating equations
library("tidyverse")

avg_alpha_div <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormed.xlsx", 2)
alpha_div <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormed.xlsx", 1)
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

# Pathotype_1, Pattern_1, EvNev_1
# Decided to code model names as follows: Pathotype outcome will be 1, EvNev
# will be 2 and Pattern will be 3. For independent vars, richness will be
# 1, shannons 2 and evenness 3. So 1.1 will be pathotype with richness, for
# example

## Pathotype

m1.1 <- glmer(Pathotype_1 ~ Normrich + (1|Individual_animal), data = alpha_div,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.1)

####### am getting lots of errors with this richness one, see below
######## errors with richness and pathotype:
#Warning messages:
 # 1: Some predictor variables are on very different scales: consider rescaling 
#2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
 #                 Model failed to converge with max|grad| = 1.1526 (tol = 0.001, component 1)
  #              3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
   #                               Model is nearly unidentifiable: very large eigenvalue
    #                            - Rescale variables?;Model is nearly unidentifiable: large eigenvalue ratio
     #                           - Rescale variables?
# per internet sleuthing of methods and ways to deal with issues in R:
# Options: rescale and center continuous parameters (warning 3)
# Options: check singularity (warning 2)
# Options: increase interations (warning 2)
# Options: change the optimizer (From bobyqa to other)

# trying rescale (scale command first centers (subtr the mean) then scales (divides by SD))
Normrich1 <- scale(alpha_div$Normrich, center = TRUE, scale = TRUE)
Normrich1 <- apply(Normrich1,1,as.numeric)
Normrich1 <- as.data.frame(Normrich1)
write.xlsx(Normrich1, "scaledrich.xlsx")
# added the column of scaled values to the xlsx because could not turn matrix back
# into numeric class and get it merged! Variable is now called Scaledrich in new data frame
# for both avg and non-avg sheets:
alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 1)
avg_alpha_div_scaled <- read.xlsx("excel sheets/Cow_map_wrichnshansnevennormednscaled.xlsx", 2)

alpha_div_scaled <- mutate(alpha_div_scaled, "Normrich" = as.numeric(Normrich), 
                    "Normshann" = as.numeric(Normshann),
                    "Normeven" = as.numeric(Normeven),
                    "Scaledrich" = as.numeric(Scaledrich))

avg_alpha_div_scaled <- mutate(avg_alpha_div_scaled, "Avgrich" = as.numeric(Avgrich), 
                        "Avgshann" = as.numeric(Avgshann),
                        "Avgeven" = as.numeric(Avgeven),
                        "Avgscaledrich" = as.numeric(Avgscaledrich))

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

m1.1 <- glmer(Pathotype_1 ~ Scaledrich + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.1)

# not getting any errors now


m1.2 <- glmer(Pathotype_1 ~ Normshann + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
m1.2.0 <- glmer(Pathotype_1 ~ (1|Individual_animal), data = alpha_div_scaled, family =
                  binomial, control = glmerControl(optimizer = "bobyqa"))
## the above is the general null model for pathotype. do we even need this though? dont think so
summary(m1.2)
anova(m1.2,m1.2.0)


m1.3 <- glmer(Pathotype_1 ~ Normeven + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.3)

## EverNever

m2.1 <- glm(EvNev_1 ~ Avgrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1)
#trying scaled average richness
m2.1.1 <- glm(EvNev_1 ~ Avgscaledrich, data = avg_alpha_div_scaled, family = binomial)
summary(m2.1.1)

m2.2 <- glm(EvNev_1 ~ Avgshann, data = avg_alpha_div_scaled, family = binomial)
summary(m2.2)

m2.3 <- glm(EvNev_1 ~ Avgeven, data = avg_alpha_div_scaled, family = binomial)
summary(m2.3)


## Pattern

m3.1 <- multinom(Pattern_1 ~ Avgrich, data = avg_alpha_div_scaled)
summary(m3.1)
# calculate p values using wald tests (z-tests)
z <- summary(m3.1)$coefficients/summary(m3.1)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2  
p  

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

  