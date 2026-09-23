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

### run this script before running the formal analysis
### functions to compile the growth data frame
lm.table.r <- function(coral.table){
  a <- coral.table[ coral.table$Area_1 != 0,]
  info <- a[,1:10]
  log_Area_year1 <- a$log_Area_1
  log_Area_year2 <- a$log_Area_stan_2
  Area_year1 <- a$Area_1
  Area_year2 <- a$Area_stan_2
  Sampling_year <- rep('y12', length(log_Area_year1))
  lm.table1 <- data.frame(info, Sampling_year, Area_year1, Area_year2, log_Area_year1, log_Area_year2)
  
  b <- coral.table[ coral.table$Area_3 != 0,]
  info <- b[,1:10]
  log_Area_year1 <- b$log_Area_2
  log_Area_year2 <- b$log_Area_stan_3
  Area_year1 <- b$Area_2
  Area_year2 <- b$Area_stan_3
  Sampling_year <- rep('y23', length(log_Area_year1))
  lm.table2 <- data.frame(info, Sampling_year, Area_year1, Area_year2,log_Area_year1, log_Area_year2)
  
  lm.table <- rbind(lm.table1,lm.table2)
  lm.table <-na.omit(lm.table)
  lm.table
}

lm.table.p <- function(coral.table){
  a <- coral.table[coral.table$Net_horizental_growth_12>0 & coral.table$Area_1 != 0,]
  info <- a[,1:10]
  log_Area_year1 <- a$log_Area_1
  log_Area_year2 <- a$log_Area_stan_2
  Area_year1 <- a$Area_1
  Area_year2 <- a$Area_stan_2
  Sampling_year <- rep('y12', length(log_Area_year1))
  lm.table1 <- data.frame(info, Sampling_year, Area_year1, Area_year2, log_Area_year1, log_Area_year2)
  
  b <- coral.table[coral.table$Net_horizental_growth_23>0 & coral.table$Area_3 != 0,]
  info <- b[,1:10]
  log_Area_year1 <- b$log_Area_2
  log_Area_year2 <- b$log_Area_stan_3
  Area_year1 <- b$Area_2
  Area_year2 <- b$Area_stan_3
  Sampling_year <- rep('y23', length(log_Area_year1))
  lm.table2 <- data.frame(info, Sampling_year, Area_year1, Area_year2,log_Area_year1, log_Area_year2)
  
  lm.table <- rbind(lm.table1,lm.table2)
  lm.table <-na.omit(lm.table)
  lm.table
}

### functions to check the assumptions of LMM - heterogeneity and normality
variance.com <- function(coral.growth){
  lmm.1 <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = coral.growth)
  lmm.simres <- simulateResiduals(fittedModel = lmm.1)
  lmm.simres.1 <- testCategorical(lmm.simres, coral.growth$Lat_reg)
  title(main = "Model 1")
  
  lmm.2 <- glmmTMB(
    log_Area_year2 ~ log_Area_year1 * Lat_reg 
    + (1 | CoralID) + (1 | Sampling_year)+ (1 | Location/Site) ,
    dispformula = ~  Lat_reg,   # models residual variance as function of Lat_reg
    family =  gaussian(),
    data = coral.growth
  )
  lmm.2.simres <- simulateResiduals(fittedModel = lmm.2)
  lmm.simres.2 <- testCategorical(lmm.2.simres, coral.growth$Lat_reg)
  title(main = "Model 2")
  
  lmm.3 <- glmmTMB(
    log_Area_year2 ~ log_Area_year1 * Lat_reg 
    + (1 | Location/Site) + (1 | CoralID) + (1 | Sampling_year),
    dispformula = ~ Lat_reg + log_Area_year1,   # models residual variance as function of Lat_reg
    family = gaussian(),
    data = coral.growth
  )
  lmm.3.simres <- simulateResiduals(fittedModel = lmm.3)
  lmm.simres.3 <- testCategorical(lmm.3.simres, coral.growth$Lat_reg)
  title(main = "Model 3")
  
  lmm.4 <- glmmTMB(
    log_Area_year2 ~ log_Area_year1 * Lat_reg 
    + (1 | Location/Site) + (1 | CoralID) + (1 | Sampling_year),
    dispformula = ~  log_Area_year1,   # models residual variance as function of Lat_reg
    family = gaussian(),
    data = coral.growth
  )
  lmm.4.simres <- simulateResiduals(fittedModel = lmm.4)
  lmm.simres.4 <- testCategorical(lmm.4.simres, coral.growth$Lat_reg)
  title(main = "Model 4")
  
  print(AIC(lmm.1,lmm.2,lmm.3,lmm.4))
  print(c(lmm.simres.1,lmm.simres.2,lmm.simres.3,lmm.simres.4))
}




variance.com.2 <- function(coral.growth, effects){
  lmm.1 <- lmer(log_Area_year2 ~ log_Area_year1*Lat_reg + (1|Location/Site) + (1|CoralID)+ (1|Sampling_year), REML= F, data = coral.growth)
  lmm.simres <- simulateResiduals(fittedModel = lmm.1)
  lmm.simres.1 <- testCategorical(lmm.simres, coral.growth$Lat_reg)
  title(main = "Model 1")
  
  lmm.2 <- glmmTMB(
    effects ,
    dispformula = ~  Lat_reg,   # models residual variance as function of Lat_reg
    family =  gaussian(),
    data = coral.growth
  )
  lmm.2.simres <- simulateResiduals(fittedModel = lmm.2)
  lmm.simres.2 <- testCategorical(lmm.2.simres, coral.growth$Lat_reg)
  title(main = "Model 2")
  
  lmm.3 <- glmmTMB(
    effects,
    dispformula = ~ Lat_reg + log_Area_year1,   # models residual variance as function of Lat_reg
    family = gaussian(),
    data = coral.growth
  )
  lmm.3.simres <- simulateResiduals(fittedModel = lmm.3)
  lmm.simres.3 <- testCategorical(lmm.3.simres, coral.growth$Lat_reg)
  title(main = "Model 3")
  
  lmm.4 <- glmmTMB(
    effects,
    dispformula = ~  log_Area_year1,   # models residual variance as function of Lat_reg
    family = gaussian(),
    data = coral.growth
  )
  lmm.4.simres <- simulateResiduals(fittedModel = lmm.4)
  lmm.simres.4 <- testCategorical(lmm.4.simres, coral.growth$Lat_reg)
  title(main = "Model 4")
  
  print(AIC(lmm.1,lmm.2,lmm.3,lmm.4))
  print(c(lmm.simres.1,lmm.simres.2,lmm.simres.3,lmm.simres.4))
}

normality.com <- function(coral.growth,weight){
  if(weight == 'no'){
    lmm.original <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | CoralID) + (1 | Sampling_year) + (1|Location/Site),
      family = gaussian(),
      data = coral.growth
    )
    lmm.original.simres <- simulateResiduals(fittedModel = lmm.original, plot = T)
    lmm.original.simres <- testCategorical(lmm.original.simres, coral.growth$Lat_reg)
    title(main = "Model original")
    
    lmm.5 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | CoralID) + (1 | Sampling_year),
      family = gaussian(),
      data = coral.growth
    )
    lmm.5.simres <- simulateResiduals(fittedModel = lmm.5, plot = T)
    title(main = "Model 5")
    
    lmm.6 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location/Site)  + (1 | Sampling_year),
      family = gaussian(),
      data = coral.growth
    )
    lmm.6.simres <- simulateResiduals(fittedModel = lmm.6, plot = T)
    title(main = "Model 6")
    
    lmm.7 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location/Site) + (1 | CoralID) ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.7.simres <- simulateResiduals(fittedModel = lmm.7, plot = T)
    title(main = "Model 7")
    
    lmm.8 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location/Site),
      family = gaussian(),
      data = coral.growth
    )
    lmm.8.simres <- simulateResiduals(fittedModel = lmm.8, plot = T)
    title(main = "Model 8")
    
    lmm.9 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | CoralID)  ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.9.simres <- simulateResiduals(fittedModel = lmm.9, plot = T)
    title(main = "Model 9")
    
    lmm.10 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Sampling_year)  ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.10.simres <- simulateResiduals(fittedModel = lmm.10, plot = T)
    title(main = "Model 10")
    
    lmm.11 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location) + (1 | CoralID) + (1 | Sampling_year),
      family = gaussian(),
      data = coral.growth
    )
    lmm.11.simres <- simulateResiduals(fittedModel = lmm.11, plot = T)
    title(main = "Model 11")
    
    lmm.12 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location) + (1 | CoralID) ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.12.simres <- simulateResiduals(fittedModel = lmm.12, plot = T)
    title(main = "Model 12")
    
    lmm.13 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location) + (1 | Sampling_year) ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.13.simres <- simulateResiduals(fittedModel = lmm.13, plot = T)
    title(main = "Model 13")
    
    lmm.14 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location),
      family = gaussian(),
      data = coral.growth
    )
    lmm.14.simres <- simulateResiduals(fittedModel = lmm.14, plot = T)
    title(main = "Model 14")
    
    lmm.15 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | CoralID) + (1 | Sampling_year),
      family = gaussian(),
      data = coral.growth
    )
    lmm.15.simres <- simulateResiduals(fittedModel = lmm.15, plot = T)
    title(main = "Model 15")
    
    lmm.16 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | CoralID) ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.16.simres <- simulateResiduals(fittedModel = lmm.16, plot = T)
    title(main = "Model 16")
    
    lmm.17 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | Sampling_year) ,
      data = coral.growth
    )
    lmm.17.simres <- simulateResiduals(fittedModel = lmm.17, plot = T)
    title(main = "Model 17")
    
    lmm.18 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.18.simres <- simulateResiduals(fittedModel = lmm.18, plot = T)
    title(main = "Model 18")
    
    lmm.19 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | Location) ,
      family = gaussian(),
      data = coral.growth
    )
    lmm.19.simres <- simulateResiduals(fittedModel = lmm.19, plot = T)
    title(main = "Model 19")
    
    AIC(lmm.original,lmm.5,lmm.6,lmm.7,lmm.8,lmm.9,
        lmm.10,lmm.11,lmm.12,lmm.13,lmm.14,lmm.15,
        lmm.16,lmm.17,lmm.18,lmm.19)
    
  }
  else{
    lmm.original <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | CoralID) + (1 | Sampling_year) + (1|Location/Site),
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.original.simres <- simulateResiduals(fittedModel = lmm.original, plot = T)
    title(main = "Model original")
    
    lmm.5 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | CoralID) + (1 | Sampling_year),
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.5.simres <- simulateResiduals(fittedModel = lmm.5, plot = T)
    title(main = "Model 5")
    
    lmm.6 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location/Site)  + (1 | Sampling_year),
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.6.simres <- simulateResiduals(fittedModel = lmm.6, plot = T)
    title(main = "Model 6")
    
    lmm.7 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location/Site) + (1 | CoralID) ,
      dispformula = weight,  
      family = gaussian(),
      data = coral.growth
    )
    lmm.7.simres <- simulateResiduals(fittedModel = lmm.7, plot = T)
    title(main = "Model 7")
    
    lmm.8 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location/Site)  ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.8.simres <- simulateResiduals(fittedModel = lmm.8, plot = T)
    title(main = "Model 8")
    
    lmm.9 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | CoralID)  ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.9.simres <- simulateResiduals(fittedModel = lmm.9, plot = T)
    title(main = "Model 9")
    
    lmm.10 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Sampling_year)  ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.10.simres <- simulateResiduals(fittedModel = lmm.10, plot = T)
    title(main = "Model 10")
    
    lmm.11 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location) + (1 | CoralID) + (1 | Sampling_year),
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.11.simres <- simulateResiduals(fittedModel = lmm.11, plot = T)
    title(main = "Model 11")
    
    lmm.12 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location) + (1 | CoralID) ,
      dispformula = weight,  
      family = gaussian(),
      data = coral.growth
    )
    lmm.12.simres <- simulateResiduals(fittedModel = lmm.12, plot = T)
    title(main = "Model 12")
    
    lmm.13 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location) + (1 | Sampling_year) ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.13.simres <- simulateResiduals(fittedModel = lmm.13, plot = T)
    title(main = "Model 13")
    
    lmm.14 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Location)  ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.14.simres <- simulateResiduals(fittedModel = lmm.14, plot = T)
    title(main = "Model 14")
    
    lmm.15 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | CoralID) + (1 | Sampling_year),
      dispformula = weight,  
      family = gaussian(),
      data = coral.growth
    )
    lmm.15.simres <- simulateResiduals(fittedModel = lmm.15, plot = T)
    title(main = "Model 15")
    
    lmm.16 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | CoralID) ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.16.simres <- simulateResiduals(fittedModel = lmm.16, plot = T)
    title(main = "Model 16")
    
    lmm.17 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | Sampling_year) ,
      dispformula = weight,   
      data = coral.growth
    )
    lmm.17.simres <- simulateResiduals(fittedModel = lmm.17, plot = T)
    title(main = "Model 17")
    
    lmm.18 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site)  ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.18.simres <- simulateResiduals(fittedModel = lmm.18, plot = T)
    title(main = "Model 18")
    
    lmm.19 <- glmmTMB(
      log_Area_year2 ~ log_Area_year1 * Lat_reg 
      + (1 | Site) + (1 | Location) ,
      dispformula = weight,   
      family = gaussian(),
      data = coral.growth
    )
    lmm.19.simres <- simulateResiduals(fittedModel = lmm.19, plot = T)
    title(main = "Model 19")
    
    
    AIC(lmm.original,lmm.5,lmm.6,lmm.7,lmm.8,lmm.9,
        lmm.10,lmm.11,lmm.12,lmm.13,lmm.14,lmm.15,
        lmm.16,lmm.17,lmm.18,lmm.19)
  }
}


