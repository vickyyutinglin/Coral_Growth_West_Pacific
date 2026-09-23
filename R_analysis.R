#' ---
#' title: "Taxon- and size-specific latitudinal variation in coral growth across the western Pacific"
#' author: "XXX"
#' date: "2025-12-24"
#' ---
#' 
#' Data and script to replicate analyses in XXX et al. (XXX). If used in full or in part, please cite the original publication: 
#' 
#' XXX (Under Revision in Ecography) Taxon- and size-specific latitudinal variation in coral growth across the western Pacific. DOI: XXX
#' 
#' Raw data also available on [Dryad](XXX)
#' 

#####  Packages
rm(list = ls()) 
library(lme4)#
library(lmerTest)#
library(dplyr)#
library(ggplot2)#
library(ggrepel)#
library(stats)#
library(tidyr)#
library(utils)#
library(vegan)#
library(stats) #
library(glmmTMB)
library(DHARMa)#
library(sjPlot)#
library(rcompanion)#


##### Environmental data analysis
env_14 <- read.csv('Data//env.csv', header = T)
Site <- env_14[,2]
Lat_GP <- env_14[,5]
env_env <- env_14[-c(1:5)]
env_env$Lat_GP <- Lat_GP

# nmds
set.seed(123)
nmds.fct <- metaMDS(env_env[,1:16], distance = "bray")  
nmds.fct$stress # check the stress value of nMDS
stressplot(nmds.fct) # Appendix G

en  <-  envfit(nmds.fct, env_env, permutations = 999, na.rm = TRUE)
site.scores.fct <- as.data.frame(scores(nmds.fct)$sites) # extract site scores
site.scores.fct <- data.frame(env_14$X, site.scores.fct) 
colnames(site.scores.fct) <- c('no','NMDS1','NMDS2')
en_coord_cont <- as.data.frame(scores(en, "vectors"))[c(1:6,8:15),] * 0.18 

ggplot(data = site.scores.fct) +
  geom_point(aes(x = NMDS1, y = NMDS2, color = Lat_GP),size = 8, alpha = 0.5) + 
  stat_ellipse(geom = "polygon",level = 0.95, aes(x = NMDS1, y = NMDS2, group = Lat_GP, fill =Lat_GP), alpha = 0.3)+
  geom_text(mapping = aes(x = NMDS1, y = NMDS2, label = no), size = 4) +
  scale_color_manual(values = c('#008FF8FF','#E23B0EFF'))+
  scale_fill_manual(values = c('#008FF8FF','#E23B0EFF'))+
  theme_minimal() +
  coord_fixed() +
  theme(legend.position = "right",text = element_text(size = 12))+
  xlim(c(-0.25,0.25))+
  ylim(c(-0.20,0.20)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, alpha = 0.5, colour = "grey30") +
  geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2, label = row.names(en_coord_cont)), size = 5) 

##### Data preparation
coral_growth <- read.csv('Data\\coral_growth.csv', header = T, sep = ',', stringsAsFactors = F)
coral_growth <- coral_growth %>% 
  filter((Area_1 != 0 & Area_2 != 0& Area_3 != 0)|
           (Area_1 != 0 & Area_2 != 0& Area_3 == 0)|
           (Area_1 != 0 & Area_2 != 0& (is.na(Area_3)))|
           (Area_1 == 0 & Area_2 != 0& Area_3 != 0)|
           ((is.na(Area_1)) & Area_2 != 0& Area_3 != 0))%>% 
  mutate(Area_stan_2 = Area_1 + (((Area_2 - Area_1) / Interval_days_12) *365),
         Area_stan_3 = Area_2 + (((Area_3 - Area_2) / Interval_days_23) *365),
         Net_horizental_growth_12 = Area_stan_2 - Area_1,
         Net_horizental_growth_23 = Area_stan_3 - Area_2,
         log_Area_1 = log10(Area_1 ),
         log_Area_stan_2 = log10(Area_stan_2 ),
         log_Area_2 = log10(Area_2 ),
         log_Area_stan_3 = log10(Area_stan_3 ))

Acropora.r.table <- coral_growth%>%
  filter(Genus == 'Acropora', 
         Morpho == 'table')
Acropora.r <- lm.table.r(Acropora.r.table) 
Acropora.p <- lm.table.p(Acropora.r.table) 

Pocillopora.r.table <- coral_growth%>%
  filter(Genus == 'Pocillopora', 
         Morpho == 'bushy')
Pocillopora.r <- lm.table.r(Pocillopora.r.table)
Pocillopora.p <- lm.table.p(Pocillopora.r.table)

Dipsastraea.r.table <- coral_growth%>%
  filter(Genus == 'Dipsastraea', 
         Morpho == 'encrusting', )
Dipsastraea.r <- lm.table.r(Dipsastraea.r.table)
Dipsastraea.p <- lm.table.p(Dipsastraea.r.table)

Favites.r.table <- coral_growth%>%
  filter(Genus == 'Favites', 
         Morpho == 'encrusting')
Favites.r <- lm.table.r(Favites.r.table)
Favites.p <- lm.table.p(Favites.r.table)

Porites.r.table <- coral_growth%>%
  filter(Genus == 'Porites', 
         Morpho == 'encrusting')
Porites.r <- lm.table.r(Porites.r.table)
Porites.p <- lm.table.p(Porites.r.table)

###### Potential growth #####
### Acropora ###
lmm.Acropora <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Acropora.p)
summary(lmm.Acropora)
lmm.Acropora.simres <- simulateResiduals(fittedModel = lmm.Acropora, plot = T)
testCategorical(lmm.Acropora.simres, Acropora.p$Lat_reg)

# adding weight to fixed effect in LMM and choose the best model with AIC
variance.com(Acropora.p) # model 3 is the best
lmm.Acropora.3 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Location/Site) + (1 | CoralID) + (1 | Sampling_year),
  dispformula = ~ Lat_reg + log_Area_year1,   # models residual variance as function of Lat_reg
  family = gaussian(),
  data = Acropora.p
)
lmm.Acropora.3.simres <- simulateResiduals(fittedModel = lmm.Acropora.3, plot = T)
testCategorical(lmm.Acropora.3,Acropora.p$Lat_reg,plot=T)
summary(lmm.Acropora.3)

Predicted <- predict(lmm.Acropora.3, type="response")
Residuals <- residuals(lmm.Acropora.3)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Acropora.3, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF'))+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")

### Pocillopora ###
lmm.Pocillopora <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Pocillopora.p)
summary(lmm.Pocillopora)
lmm.Pocillopora.simres <- simulateResiduals(fittedModel = lmm.Pocillopora, plot = T)
testCategorical(lmm.Pocillopora.simres, Pocillopora.p$Lat_reg)

# adding weight to fixed effect in LMM and choose the best model with AIC
variance.com(Pocillopora.p) # model 4 is the best
lmm.Pocillopora.4 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Location/Site) + (1 | CoralID) + (1 | Sampling_year),
  dispformula = ~  log_Area_year1,  
  family = gaussian(),
  data = Pocillopora.p
)
lmm.Pocillopora.4.simres <- simulateResiduals(fittedModel = lmm.Pocillopora.4, plot = T)
testCategorical(lmm.Pocillopora.4,Pocillopora.p$Lat_reg,plot=T)
summary(lmm.Pocillopora.4)

Predicted <- predict(lmm.Pocillopora.4, type="response")
Residuals <- residuals(lmm.Pocillopora.4)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Pocillopora.4, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF')) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")


### Porites ###
lmm.Porites <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Porites.p)
summary(lmm.Porites)
lmm.Porites.simres <- simulateResiduals(fittedModel = lmm.Porites, plot = T)
testCategorical(lmm.Porites.simres, Porites.p$Lat_reg)

# adding weight to fixed effect in LMM and choose the best model with AIC
variance.com(Porites.p) # model 3 is the best
lmm.Porites.3 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Location/Site) + (1 | CoralID) + (1 | Sampling_year),
  dispformula = ~ Lat_reg + log_Area_year1,   # models residual variance as function of Lat_reg
  family = gaussian(),
  data = Porites.p
)
lmm.Porites.3.simres <- simulateResiduals(fittedModel = lmm.Porites.3, plot = T)
testCategorical(lmm.Porites.3,Porites.p$Lat_reg,plot=T)
summary(lmm.Porites.3)

Predicted <- predict(lmm.Porites.3, type="response")
Residuals <- residuals(lmm.Porites.3)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Porites.3, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF')) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")

### Dipsastraea ###
lmm.Dipsastraea <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Dipsastraea.p)
summary(lmm.Dipsastraea)
lmm.Dipsastraea.simres <- simulateResiduals(fittedModel = lmm.Dipsastraea, plot = T)
testCategorical(lmm.Dipsastraea.simres, Dipsastraea.p$Lat_reg)

# check the normality in LMM and choose the best model with AIC
normality.com(Dipsastraea.p, 'no') # model 5 is the best
lmm.Dipsastraea.5 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | CoralID) + (1 | Sampling_year),
  family = gaussian(),
  data = Dipsastraea.p
)
lmm.Dipsastraea.5.simres <- simulateResiduals(fittedModel = lmm.Dipsastraea.5, plot = T)
testCategorical(lmm.Dipsastraea.5.simres, Dipsastraea.p$Lat_reg)
summary(lmm.Dipsastraea.5)


variance.com.2(Dipsastraea.p,log_Area_year2 ~ log_Area_year1 * Lat_reg  + (1 | CoralID) + (1 | Sampling_year))
lmm.Dipsastraea.4 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | CoralID) + (1 | Sampling_year),
  family = gaussian(),
  dispformula = ~  log_Area_year1,   # models residual variance as function of Lat_reg
  data = Dipsastraea.p
)
lmm.Dipsastraea.4.simres <- simulateResiduals(fittedModel = lmm.Dipsastraea.4, plot = T)
testCategorical(lmm.Dipsastraea.4.simres, Dipsastraea.p$Lat_reg)
summary(lmm.Dipsastraea.4)

Predicted <- predict(lmm.Dipsastraea.4, type="response")
Residuals <- residuals(lmm.Dipsastraea.4)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Dipsastraea.4, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF')) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")

### Favites ###
lmm.Favites <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Favites.p)
summary(lmm.Favites)
lmm.Favites.simres <- simulateResiduals(fittedModel = lmm.Favites, plot = T)
testCategorical(lmm.Favites.simres, Favites.p$Lat_reg)
Predicted <- predict(lmm.Favites, type="response")
Residuals <- residuals(lmm.Favites)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")
plot_model(lmm.Favites, type = "int", show.data = T,
           colors = c('#008FF8FF','#E23B0EFF')) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")


###### Realized growth #####
### Acropora ###
lmm.Acropora.r <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Acropora.r)
summary(lmm.Acropora.r)
lmm.Acropora.r.simres <- simulateResiduals(fittedModel = lmm.Acropora.r, plot = T)
testCategorical(lmm.Acropora.r.simres, Acropora.r$Lat_reg)

# adding weight to fixed effect in LMM and choose the best model with AIC
variance.com(Acropora.r) # model 3 is the best
lmm.Acropora.r.3 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | CoralID) + (1 | Sampling_year) + (1|Location/Site),
  dispformula = ~ Lat_reg + log_Area_year1,  
  family = gaussian(),
  data = Acropora.r
)
lmm.Acropora.r.3.simres <- simulateResiduals(fittedModel = lmm.Acropora.r.3, plot =T)
testCategorical(lmm.Acropora.r.3.simres, Acropora.r$Lat_reg)
summary(lmm.Acropora.r.3)

Predicted <- predict(lmm.Acropora.r.3, type="response")
Residuals <- residuals(lmm.Acropora.r.3)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Acropora.r.3, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF'))+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")

### Pocillopora ###
lmm.Pocillopora.r <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Pocillopora.r)
summary(lmm.Pocillopora.r)
lmm.Pocillopora.r.simres <- simulateResiduals(fittedModel = lmm.Pocillopora.r, plot = T)
testCategorical(lmm.Pocillopora.r.simres, Pocillopora.r$Lat_reg)

# check the normality in LMM and choose the best model with AIC
normality.com(Pocillopora.r, 'no') # model 5 is the best

lmm.Pocillopora.r.14 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Location),
  family = gaussian(),
  data = Pocillopora.r
)
lmm.Pocillopora.14.simres <- simulateResiduals(fittedModel = lmm.Pocillopora.r.14, plot = T)
testCategorical(lmm.Pocillopora.14.simres, Pocillopora.r$Lat_reg)
summary(lmm.Pocillopora.r.14)

Predicted <- predict(lmm.Pocillopora.r.14, type="response")
Residuals <- residuals(lmm.Pocillopora.r.14)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Pocillopora.r.14, type = "int", show.data = T,
           colors = c('#008FF8FF','#E23B0EFF'))+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")

### Porites ### 
lmm.Porites.r <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Porites.r)
summary(lmm.Porites.r)
lmm.Porites.r.simres <- simulateResiduals(fittedModel = lmm.Porites.r, plot = T)
testCategorical(lmm.Porites.r.simres, Porites.r$Lat_reg)

# check the normality in LMM and choose the best model with AIC
normality.com(Porites.r, 'no')
# none of models meet the normality assumption, so I choose the one with lowest AIC - model 5
lmm.Porites.r.18 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Site)  ,
  family = gaussian(),
  data = Porites.r
)
summary(lmm.Porites.r.18)
lmm.Porites.r.18.simres <- simulateResiduals(fittedModel = lmm.Porites.r.18, plot = T)
testCategorical(lmm.Porites.r.18.simres, Porites.r$Lat_reg)
Predicted <- predict(lmm.Porites.r.18, type="response")
Residuals <- residuals(lmm.Porites.r.18)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Porites.r.18, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF'))+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")


### Dipsastraea ###
lmm.Dipsastraea.r <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Dipsastraea.r)
summary(lmm.Dipsastraea.r)
lmm.Dipsastraea.r.simres <- simulateResiduals(fittedModel = lmm.Dipsastraea.r, plot = T)
testCategorical(lmm.Dipsastraea.r.simres, Dipsastraea.r$Lat_reg)

# adding weight to fixed effect in LMM and choose the best model with AIC
variance.com(Dipsastraea.r) # model 3 is the best
lmm.Dipsastraea.r.2 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Location/Site) + (1 | CoralID) + (1 | Sampling_year),
  dispformula = ~  Lat_reg ,  
  family = gaussian(),
  data = Dipsastraea.r
)
lmm.Dipsastraea.r.2.simres <- simulateResiduals(fittedModel = lmm.Dipsastraea.r.2, plot = T)
testCategorical(lmm.Dipsastraea.r.2.simres, Dipsastraea.r$Lat_reg)
summary(lmm.Dipsastraea.r.2)

# check the normality in LMM and choose the best model with AIC
normality.com(Dipsastraea.r, ~ Lat_reg)
# none of models meet the normality assumption, so I choose the one with lowest AIC - model 13
lmm.Dipsastraea.r.13 <- glmmTMB(
  log_Area_year2 ~ log_Area_year1 * Lat_reg 
  + (1 | Location) + (1 | Sampling_year) ,
  dispformula = ~ Lat_reg + log_Area_year1,   
  family = gaussian(),
  data = Dipsastraea.r
)
lmm.Dipsastraea.r.13.simres <- simulateResiduals(fittedModel = lmm.Dipsastraea.r.13, n = 1000, plot = T)
testCategorical(lmm.Dipsastraea.r.13.simres, Dipsastraea.r$Lat_reg)
summary(lmm.Dipsastraea.r.13)
redicted <- predict(lmm.Dipsastraea.r.13, type="response")
Residuals <- residuals(lmm.Dipsastraea.r.13)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Dipsastraea.r.13, type = "int", show.data = T,
           colors = c('#E23B0EFF','#008FF8FF'))+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")


### Favites ###
lmm.Favites.r <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = Favites.r)
summary(lmm.Favites.r)
lmm.Favites.r.simres <- simulateResiduals(fittedModel = lmm.Favites.r, plot = T)
testCategorical(lmm.Favites.r.simres, Favites.r$Lat_reg)

Predicted <- predict(lmm.Favites.r, type="response")
Residuals <- residuals(lmm.Favites.r)
efronRSquared(residual = Residuals, 
              predicted = Predicted, 
              statistic = "EfronRSquared")

plot_model(lmm.Favites.r, type = "int", show.data = T,
           colors = c('#008FF8FF','#E23B0EFF'))+ 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = log10(5), linetype = "dashed", color = "red")

