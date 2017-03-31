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

m1.1.0 <- glmer(Pathotype_1 ~ Scaledrich + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.1.0)

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

# is it appropriate to calc the p vals with a normal distribution when
# the avgrich isn't normally distributed? Check normality of avgrich and
# scaled avgrich, and look at model with scaled avg rich

shapiro.test(avg_alpha_div_scaled$Avgrich)
shapiro.test(avg_alpha_div_scaled$Avgscaledrich)
qqnorm(avg_alpha_div_scaled$Avgrich, ylab= "Average Richness")
qqnorm(avg_alpha_div_scaled$Avgscaledrich, ylab= "Average Richness Scaled")

# both have sig p vals, not normally distributed. quantile plots look 
# pretty normal though. 

m3.1.0 <- multinom(Pattern_1 ~ Avgscaledrich, data = avg_alpha_div_scaled)
summary(m3.1.0)
z <- summary(m3.1.0)$coefficients/summary(m3.1.0)$standard.errors
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


## to get more meaningful values, can calculate the IQR of diversity values 
# (between 1st and 3rd quartile) and multiply the B coefficients by that:

IQR_scaledrich <- IQR(alpha_div_scaled$Scaledrich)
IQR_scaledrich
# 1.036271
IQR_normrich <- IQR(alpha_div_scaled$Normrich)
IQR_normrich  
# 1255
IQR_avgrich <- IQR(avg_alpha_div_scaled$Avgrich)
IQR_avgrich 
# 797.65
IQR_avgrichscaled <- IQR(avg_alpha_div_scaled$Avgscaledrich)
IQR_avgrichscaled 
# 0.6586304
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

richspread <- alpha_div_scaled %>% group_by(Pattern_1) %>% summarize(mean(Normrich), 
                                                                 min(Normrich), 
                                                                 max(Normrich))
richspread

# Look at adding different variables to models and see how coefficients change
# Variables of interest are: Disease, Parity, DIM, Farm
# m1.1 richness and pathotype (cow randome)
# m1.1.0 scaled richness and pathotype (cow random)
# m1.2 shannons and pathotype (cow random)
# m 1.3 evenness and pathotype (cow random)

alpha_div_scaled <- alpha_div_scaled %>% mutate("Parity" = as.factor(Parity),
                                                "DIM" = as.numeric(DIM),
                                                "Disease" = as.factor(Disease),
                                                "Farm" = as.factor(Farm))
avg_alpha_div_scaled <- avg_alpha_div_scaled %>% mutate("Parity" = as.factor(Parity),
                                                "DIM" = as.numeric(DIM),
                                                "Disease" = as.factor(Disease),
                                                "Farm" = as.factor(Farm))
# run model again with the right scaled data set
m1.1 <- glmer(Pathotype_1 ~ Normrich + (1|Individual_animal), data = alpha_div_scaled,
              family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(m1.1)

##PATHOTYPE
m1.1farm <- update(m1.1, .~. + Farm)
summary(m1.1farm)
m1.1Parity <-  update(m1.1, .~. + Parity)
summary(m1.1Parity)
m1.1Disease <- update(m1.1, .~. + Disease)
summary(m1.1Disease)
m1.1DIM <-  update(m1.1, .~. + DIM)
summary(m1.1DIM)
  
m1.1.0farm <- update(m1.1.0, .~. + Farm)
summary(m1.1.0farm)
m1.1.0Parity <-  update(m1.1.0, .~. + Parity)
summary(m1.1.0Parity)
m1.1.0Disease <- update(m1.1.0, .~. + Disease)
summary(m1.1.0Disease)
m1.1.0DIM <-  update(m1.1.0, .~. + DIM)
summary(m1.1.0DIM)

m1.2farm <- update(m1.2, .~. + Farm)
summary(m1.2farm)
m1.2Parity <-  update(m1.2, .~. + Parity)
summary(m1.2Parity)
m1.2Disease <- update(m1.2, .~. + Disease)
summary(m1.2Disease)
m1.2DIM <-  update(m1.2, .~. + DIM)
summary(m1.2DIM)

m1.3farm <- update(m1.3, .~. + Farm)
summary(m1.3farm)
m1.3Parity <-  update(m1.3, .~. + Parity)
summary(m1.3Parity)
m1.3Disease <- update(m1.3, .~. + Disease)
summary(m1.3Disease)
m1.3DIM <-  update(m1.3, .~. + DIM)
summary(m1.3DIM)

##EVERNEVER

m2.1farm <- update(m2.1, .~. + Farm)
summary(m2.1farm)
m2.1Parity <-  update(m2.1, .~. + Parity)
summary(m2.1Parity)
m2.1Disease <- update(m2.1, .~. + Disease)
summary(m2.1Disease)
m2.1DIM <-  update(m2.1, .~. + DIM)
summary(m2.1DIM)

m2.2farm <- update(m2.2, .~. + Farm)
summary(m2.2farm)
m2.2Parity <-  update(m2.2, .~. + Parity)
summary(m2.2Parity)
m2.2Disease <- update(m2.2, .~. + Disease)
summary(m2.2Disease)
m2.2DIM <-  update(m2.2, .~. + DIM)
summary(m2.2DIM)

m2.3farm <- update(m2.3, .~. + Farm)
summary(m2.3farm)
m2.3Parity <-  update(m2.3, .~. + Parity)
summary(m2.3Parity)
m2.3Disease <- update(m2.3, .~. + Disease)
summary(m2.3Disease)
m2.3DIM <-  update(m2.3, .~. + DIM)
summary(m2.3DIM)

##PATTERN
m3.1farm <- update(m3.1, .~. + Farm)
summary(m3.1farm)
z <- summary(m3.1farm)$coefficients/summary(m3.1farm)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.1Parity <-  update(m3.1, .~. + Parity)
summary(m3.1Parity)
z <- summary(m3.1Parity)$coefficients/summary(m3.1Parity)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.1Disease <- update(m3.1, .~. + Disease)
summary(m3.1Disease)
z <- summary(m3.1Disease)$coefficients/summary(m3.1Disease)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.1DIM <-  update(m3.1, .~. + DIM)
summary(m3.1DIM)
z <- summary(m3.1DIM)$coefficients/summary(m3.1DIM)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.1.0farm <- update(m3.1.0, .~. + Farm)
summary(m3.1.0farm)
z <- summary(m3.1.0farm)$coefficients/summary(m3.1.0farm)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.1.0Parity <-  update(m3.1.0, .~. + Parity)
summary(m3.1.0Parity)
z <- summary(m3.1.0Parity)$coefficients/summary(m3.1.0Parity)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.1.0Disease <- update(m3.1.0, .~. + Disease)
summary(m3.1.0Disease)
z <- summary(m3.1.0Disease)$coefficients/summary(m3.1.0Disease)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.1.0DIM <-  update(m3.1.0, .~. + DIM)
summary(m3.1.0DIM)
z <- summary(m3.1.0DIM)$coefficients/summary(m3.1.0DIM)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.2farm <- update(m3.2, .~. + Farm)
summary(m3.2farm)
z <- summary(m3.2farm)$coefficients/summary(m3.2farm)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.2Parity <-  update(m3.2, .~. + Parity)
summary(m3.2Parity)
z <- summary(m3.2Parity)$coefficients/summary(m3.2Parity)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.2Disease <- update(m3.2, .~. + Disease)
summary(m3.2Disease)
z <- summary(m3.2Disease)$coefficients/summary(m3.2Disease)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.2DIM <-  update(m3.2, .~. + DIM)
summary(m3.2DIM)
z <- summary(m3.2DIM)$coefficients/summary(m3.2DIM)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p

m3.3farm <- update(m3.3, .~. + Farm)
summary(m3.3farm)
z <- summary(m3.3farm)$coefficients/summary(m3.3farm)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.3Parity <-  update(m3.3, .~. + Parity)
summary(m3.3Parity)
z <- summary(m3.3Parity)$coefficients/summary(m3.3Parity)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.3Disease <- update(m3.3, .~. + Disease)
summary(m3.3Disease)
z <- summary(m3.3Disease)$coefficients/summary(m3.3Disease)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
m3.3DIM <-  update(m3.3, .~. + DIM)
summary(m3.3DIM)
z <- summary(m3.3DIM)$coefficients/summary(m3.3DIM)$standard.errors
z
p <- (1-pnorm(abs(z), 0, 1)) * 2
p
