#Module 3 Project Report - Natalie Martin & Sofia Rugova 

#Setup: loading libraries & data
library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")

#1 - Establish the anole.log data tibble
anole.log.tibble <- anole %>% 
  left_join(anole.eco, by = "Species") %>%
  filter(!Ecomorph %in% c("U", "CH")) %>% 
  na.omit() %>% 
  mutate_at(c("SVL", "HTotal", "PH", "ArbPD"), log)

#2 - Construct two linear models to assess effect of perch diameter and height 
#Effect of perch diameter: 
pd.lm <- lm(HTotal~SVL + ArbPD, data = anole.log.tibble)
#Effect of perch height 
ph.lm <- lm(HTotal~SVL + PH, data = anole.log.tibble)

#3 - Plotting the residuals of the linear models 
#Mutating the tibble to include the residuals from both models 
anole.log.tibble <- anole.log.tibble %>% 
  mutate(residuals.pd = residuals(pd.lm), residuals.ph = residuals(ph.lm))

#Plotting perch diameter residuals 
ggplot(anole.log.tibble, aes(x=Ecomorph2, y=residuals.pd)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="point", size=3) + 
  labs(y= "Perch Diameter Residuals")

#Plotting perch height residuals 
ggplot(anole.log.tibble, aes(x=Ecomorph2, y=residuals.ph)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="point", size=3) + 
  labs(y= "Perch Height Residuals")

#4 - Constructing phylogenetic least squares models
anole.tree <- read.tree("anole.tre")

#Hindlimb-SVL relationship + perch height 
pgls.ph <- gls(HTotal~SVL + PH, 
               correlation = corBrownian(1, phy = anole.tree, form = ~Species),
               data = anole.log.tibble, 
               method = "ML")
#Hindlimb-SVL relationship + perch diameter 
pgls.pd <- gls(HTotal~SVL + ArbPD, 
               correlation = corBrownian(1, phy = anole.tree, form = ~Species),
               data = anole.log.tibble, 
               method = "ML")
#Hindlimb-SVL relationship + perch height + perch diameter 
pgls.ph.pd <- gls(HTotal~SVL + PH + ArbPD, 
               correlation = corBrownian(1, phy = anole.tree, form = ~Species),
               data = anole.log.tibble, 
               method = "ML")

#5 - Assess the fit of the models using AICc and AICw
PGLSmodassess <- MuMIn::AICc(pgls.ph, pgls.pd, pgls.ph.pd)
aicw(PGLSmodassess$AICc)
#When evaluting AIC scores, the lowest score is suggestive of the best model. 
#The lowest scoring model was the model that accounted for both perch height and perch diameter. 
#This suggests that both of the covariates are significant predictors of hindlimb length in a phylogenetic context. 

#6 - Produce a plot that concisely visualizes the best fitting PGLS model
anole.log.tibble <- anole.log.tibble %>% 
  mutate(residuals.hindlimb = residuals(pgls.ph.pd))
anole.log.tibble %>% 
  dplyr::select(Ecomorph2, residuals.pd, residuals.ph, residuals.hindlimb) %>% 
  pivot_longer(cols=c("residuals.pd", "residuals.ph", "residuals.hindlimb")) %>%
  ggplot(aes(x=Ecomorph2, y=value)) +
  geom_boxplot()+
  stat_summary(fun=mean, geom = "point", size=3)+
  facet_grid(name~.,scales = "free_y")+
  ylab("residual")
