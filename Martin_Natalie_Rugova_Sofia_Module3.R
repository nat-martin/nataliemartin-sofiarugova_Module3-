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
pd.lm <- lm(HTotal~SVL + ArbPD, anole.log.tibble)
#Effect of perch height 
ph.lm <- lm(HTotal~SVL + PH, anole.log.tibble)

#3 - Plotting the residuals of the linear models 
#Mutating the tibble to include the residuals from both models 
anole.log.tibble <- anole.log.tibble %>% 
  mutate(residuals.pd = residuals(pd.lm))

anole.log.tibble <- anole.log.tibble %>% 
  mutate(residuals.ph = residuals(ph.lm))

#Plotting perch diameter residuals 
ggplot(anole.log.tibble, aes(x=Ecomorph2, y=residuals.pd)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="point", size=3)

#Plotting perch height residuals 
ggplot(anole.log.tibble, aes(x=Ecomorph2, y=residuals.ph)) + 
  geom_boxplot() + 
  stat_summary(fun=mean, geom="point", size=3)

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
anole.phylo.aic <- MuMIn::AICc(pgls.ph, pgls.pd, pgls.ph.pd)
aicw(anole.phylo.aic$AICc)
#When evaluating AIC scores, the lowest fit score (AICc) or highest w score (AICw) is suggestive of the best model. 
#The lowest fit and highest AICw scoring model was the model that accounted for both perch height and perch diameter. 
#This suggests that both of the covariates are significant predictors of hindlimb length in a phylogenetic context. 

#6 - Produce a plot that concisely visualizes the best fitting PGLS model
anole.log.tibble <- anole.log.tibble %>% 
  mutate(residuals.hindlimb = residuals(pgls.ph.pd))

anole.pglsHDgraph <- ggplot(anole.log.tibble) +
  geom_point(aes(x = ArbPD + PH, y = residuals.hindlimb, color = Ecomorph2), size = 4) 
anole.pglsHDgraph + scale_color_manual(values = c("#9f4bf5","#efa2df","#ff3568","#f4b52a","#2bd6ff","#3a3aff"))