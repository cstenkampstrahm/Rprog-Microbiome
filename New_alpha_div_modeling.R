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
library("nnet")
library("gee")
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

## Pathotype

m1 <- glmer(Pathotype_1 ~ Normrich + (1|Individual_animal), data = alpha_div,
            family = binomial, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)



