# GPP SEMs with "corrected" data
# LMC Nov 2022

## LOAD LIBRARIES
library(tidyverse)
library(ggthemes)
library(gridExtra)
library(GGally)
library(grid)
library(scales)
library(lavaan)
library(semPlot)
library(broom)

## THEME ##################################################
theme_set(theme_few(base_size = 16, base_family = "Helvetica"))
windowsFonts(Helvetica=windowsFont("TT Helvetica"))

# LOAD DATA
dat2017 <- read.csv("01_Data/06_IceChan_2017_dataforSEM.csv") %>% 
  select(-X)

#### CONCEPTUAL FRAMEWORK ####

## Starting with full model-- based on hypothesized framework in manuscript
# Model modification criteria: Run model, check summary statistics and modification indices (MIs)
# Remove any non-significant paths (p > 0.05) and add MIs (MI > 3) that are ecologically relevant (and not redundant with existing paths in model)
# Re-run model, re-check summary stats and MIs, re-adjust model, repeat until model fit indices indicate well-fitting model and all paths sig (p < 0.05)

# N and N:P are numerical, but "technically" categorical b/c only two levels. Using the numeric values-- the standardized coefficients are the same as when 
# N is treated as a factor in the model. Using the numeric value lets me interpret a one unit increase in N in the "estimate" column

#### FULL MODEL #######
## Ran model:
## NS paths: Nfixers ~ Nfix, GPP ~ Nfix + Nfixers
## highest MIs (that make ecological sense): Nfix ~~ GPP, Nfixers ~~ GPP

path_full = "
  #regression model
Nfixers ~ temp + N + Nfix
nonNfixers ~ N + Nup
Nfix ~ temp + N
Nup ~ N + temp
GPP ~ Nup + Nfix + Nfixers + nonNfixers + temp
  #residual correlation
"
# Fit model
fit_full = lavaan(path_full, data = dat2017,
                        auto.var=TRUE, auto.fix.first=TRUE,
                        auto.cov.lv.x=TRUE)
summary(fit_full, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_full)
semPaths(fit_full, "model1", "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
mi <- modindices(fit_full)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Removed NS paths from inital full model: Nfixers ~ Nfix, GPP ~ Nfix + Nfixers
## Added MIs: Nfix ~~ GPP, Nfixers ~~ GPP

## NS paths: Nfix ~~ GPP (p = 0.126) and Nfixers ~~ GPP (p = 0.222)
## highest MIs (that make sense): GPP ~ Nfix + Nfixers

## Model fit indices good-- all paths p < 0.05 (except for covariances)
## GPP ~~ Nfixers + Nfix NS and MI suggests GPP ~ Nfixers + Nfix, so going to try that
path_fulla = "
  #regression model
Nfixers ~ temp + N
nonNfixers ~ N + Nup
Nfix ~ temp + N
Nup ~ N + temp
GPP ~ Nup + nonNfixers + temp
  #residual correlation
GPP ~~ Nfix + Nfixers
"

fit_fulla = lavaan(path_fulla, data = dat2017,
                  auto.var=TRUE, auto.fix.first=TRUE,
                  auto.cov.lv.x=TRUE)
summary(fit_fulla, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fulla)
semPaths(fit_fulla, "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
mi <- modindices(fit_fulla)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Removed GPP ~~ Nfixers + Nfix
## Added GPP ~ Nfixers + Nfix (based on MI)

## GPP ~ Nfixers + Nfix NS, so going back to GPP ~~ Nfixers + Nfix

path_fullb = "
  #regression model
Nfixers ~ temp + N
nonNfixers ~ N + Nup
Nfix ~ temp + N
Nup ~ N + temp
GPP ~ Nup + Nfix + Nfixers + nonNfixers + temp
  #residual correlation
"

fit_fullb = lavaan(path_fullb, data = dat2017,
                   auto.var=TRUE, auto.fix.first=TRUE,
                   auto.cov.lv.x=TRUE)
summary(fit_fullb, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullb)
semPaths(fit_fullb, "model1", "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
inspect(fit_fullb, 'r2') # obtain parameter R2
fitMeasures(fit_fullb)
mi <- modindices(fit_fullb)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Removed GPP ~ Nfixers + Nfix, added back GPP ~~ Nfixers + Nfix
# Tried adding Nfixers ~~ Nfix, but NS (and not suggested by MIs)

## Model fit indices good, all paths sig
path_fullc = "
  #regression model
Nfixers ~ temp + N
nonNfixers ~ N + Nup
Nfix ~ temp + N
Nup ~ temp + N
GPP ~ temp + Nup + nonNfixers
  #residual correlation
GPP ~~ Nfix + Nfixers
"

fit_fullc = lavaan(path_fullc, data = dat2017,
                   auto.var=TRUE, auto.fix.first=TRUE,
                   auto.cov.lv.x=TRUE)
summary(fit_fullc, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullc)
semPaths(fit_fullc, "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
inspect(fit_fullc, 'r2') # obtain parameter R2
fitMeasures(fit_fullc)
mi <- modindices(fit_fullc)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)



##### TESTING N:P USING FULL/ORIGINAL MODEL
### Using same criteria for model modification as N-only model

## GPP ~ Nfix + Nfixers and Nfixers ~ Nfix NS
## None of the MIs make much sense to add right now
path_fullratio = "
  #regression model
Nfixers ~ temp + NPratio + Nfix
nonNfixers ~ NPratio + Nup
Nfix ~ temp + NPratio
Nup ~ NPratio + temp
GPP ~ Nup + Nfix + Nfixers + nonNfixers + temp
  #residual correlation
"

fit_fullratio = lavaan(path_fullratio, data = dat2017,
                  auto.var=TRUE, auto.fix.first=TRUE,
                  auto.cov.lv.x=TRUE)
summary(fit_fullratio, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullratio)
mi <- modindices(fit_fullratio)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Removed NS paths: GPP ~ Nfix + Nfixers and Nfixers ~ Nfix
## All paths sig, but model fit not good. 
## Only one MI seems semi-reasonable: Nfixers ~~ nonNfixers
## Could make sense, as you'd expect a general neg correlation between the two
path_fullratioa = "
  #regression model
Nfixers ~ temp + NPratio
nonNfixers ~ NPratio + Nup
Nfix ~ temp + NPratio
Nup ~ temp + NPratio
GPP ~ temp + Nup + nonNfixers
  #residual correlation
"

fit_fullratioa = lavaan(path_fullratioa, data = dat2017,
                       auto.var=TRUE, auto.fix.first=TRUE,
                       auto.cov.lv.x=TRUE)
summary(fit_fullratioa, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullratioa)
inspect(fit_fullratioa, 'r2')
mi <- modindices(fit_fullratioa)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Added Nfixers ~~ nonNfixers (MI)
## Could argue this isn't a reasonable model w/ the Nfixers ~~ nonNfixers path
path_fullratiob = "
# regression model
Nfixers ~ temp + NPratio
nonNfixers ~ NPratio + Nup
Nfix ~ temp + NPratio
Nup ~ temp + NPratio
GPP ~ temp + Nup + nonNfixers
# residual correlation
Nfixers ~~ nonNfixers
"

fit_fullratiob = lavaan(path_fullratiob, data = dat2017,
                        auto.var=TRUE, auto.fix.first=TRUE,
                        auto.cov.lv.x=TRUE)
summary(fit_fullratiob, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullratiob)
inspect(fit_fullratiob, 'r2')
mi <- modindices(fit_fullratiob)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

### This is the same "best" model as N-only
### Running it w/ N:P doesn't result in a good model fit.
path_fullratioc = "
  #regression model
Nfixers ~ temp + NPratio
nonNfixers ~ NPratio + Nup
Nfix ~ temp + NPratio
Nup ~ NPratio + temp
GPP ~ Nup + nonNfixers + temp
  #residual correlation
GPP ~~ Nfix + Nfixers
"

fit_fullratioc = lavaan(path_fullratioc, data = dat2017,
                        auto.var=TRUE, auto.fix.first=TRUE,
                        auto.cov.lv.x=TRUE)
summary(fit_fullratioc, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullratioc)
inspect(fit_fullratioc, 'r2')
mi <- modindices(fit_fullratioc)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

##### COMPARE N and N:P MODELS
cbind(m1=inspect(fit_fullc, 'fit.measures'), m2=inspect(fit_fullratioc, 'fit.measures'))
anova(fit_fullc, fit_fullratioc)

### Testing P using full model

## P terms all NS
## Looking at other NS terms, this model should converge to same final model as N-only
path_fullP = "
  #regression model
Nfixers ~ temp + N + P + Nfix
nonNfixers ~ N + P + Nup
Nfix ~ temp + N + P
Nup ~ temp + N + P
GPP ~ temp + Nup + Nfix + Nfixers + nonNfixers
  #residual correlation
"

fit_fullP = lavaan(path_fullP, data = dat2017,
                  auto.var=TRUE, auto.fix.first=TRUE,
                  auto.cov.lv.x=TRUE)
summary(fit_fullP, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullP)
inspect(fit_fullP, 'fit.measures')
semPaths(fit_fullP, "model1", "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
mi <- modindices(fit_fullP)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Removing NS paths: Nfixers ~ Nfix + P, nonNfixers ~ P, Nfix ~ P, GPP ~ Nfix + Nfixers
## Model now contains no P terms
## Model fit good, Nup ~ P p = 0.064
path_fullPa = "
  #regression model
Nfixers ~ temp + N
nonNfixers ~ N + Nup
Nfix ~ temp + N
Nup ~ N + temp
GPP ~ Nup + nonNfixers + temp
  #residual correlation
"

fit_fullPa = lavaan(path_fullPa, data = dat2017,
                   auto.var=TRUE, auto.fix.first=TRUE,
                   auto.cov.lv.x=TRUE)
summary(fit_fullPa, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullPa)
inspect(fit_fullPa, 'fit.measures')
semPaths(fit_fullPa, "model1", "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
mi <- modindices(fit_fullPa)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)

## Fitting final "N-only" model w/ P term
## Model fit indices good, but all P terms NS
path_fullPc = "
  #regression model
Nfixers ~ temp + N + P + Nfix
nonNfixers ~ N + P + Nup
Nfix ~ temp + N + P
Nup ~ N + P + temp
GPP ~ Nup + nonNfixers + temp
  #residual correlation
# nonNfixers ~~ Nup
GPP ~~ Nfixers + Nfix
"

fit_fullPc = lavaan(path_fullPc, data = dat2017,
                   auto.var=TRUE, auto.fix.first=TRUE,
                   auto.cov.lv.x=TRUE)
summary(fit_fullPc, standardized = TRUE, fit.measures = TRUE, modindices = T)
varTable(fit_fullPc)
inspect(fit_fullPc, 'fit.measures')
semPaths(fit_fullPc, "model1", "std", bifactor = "g", layout = "tree3", residuals = FALSE, exoCov = FALSE)
mi <- modindices(fit_fullPc)  
head(mi[order(mi$mi, decreasing=TRUE), ], 10)
