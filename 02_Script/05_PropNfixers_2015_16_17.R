# Proportion N fixers for 2015, 2016, 2017
# JMH, Nov 2022, Jul 2024

# Libraries ----
# general
library(tidyverse)
# baysian zero-inflated beta regression
library(brms)

# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::install_cmdstan()

library(tidybayes)
library(broom)
library(broom.mixed)
library(emmeans)
library(grid)
library(egg)

# Data ----
chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv", row.names = 1) %>% 
  mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
         MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate)) %>% 
  # broom.mixed hates underscores
  select(Year, channel, MetDate2, MeanPre2wksTemp, NuM = N_uM, PuM = P_uM, propNfixer)

chan15 <- chan %>% 
  filter(Year == "2015") %>% 
  # only data for second date
  filter(MetDate2 == "2015-07-27") %>% 
  # convert 1's to a number just smaller than one
  mutate(propNfixerAdj = ifelse(propNfixer == 1, 0.99, propNfixer))

chan16 <- chan %>% 
  filter(Year == "2016") %>% 
  filter(MetDate2 == "2016-08-02") %>% 
  # too many 1's so modeling the proportion that are NOT N fixers
  mutate(propNOT_Nfixer = 1- propNfixer)

chan17 <- chan %>% 
  filter(Year == "2017") %>% 
  filter(MetDate2 == "2017-07-22") %>% 
  # convert 1's to a number just smaller than one
  mutate(propNfixerAdj = ifelse(propNfixer == 1, 0.99, propNfixer),
         NP = as.factor(round(NuM/PuM,1)),
         NuM_f = as.factor(NuM))


# 2015 ----
# about a third of the values are zero (i.e., no N fixers)
chan15 %>% 
  count(propNfixer == 0) %>% 
  mutate(prop = n/sum(n))

# There are 4 ones
# I converted these to 0.99 in propNfixer_adj
chan15 %>% 
  count(propNfixer == 1) %>% 
  mutate(prop = n/sum(n))



## Models ----
# https://www.andrewheiss.com/blog/2021/11/08/beta-regression-guide/#super-fancy-detailed-model-with-lots-of-moving-parts-just-for-fun

# Look at the default priors
ZIBR_2015_Formula <- bf(
  # models the mean - coefs need to be plogis()
  propNfixerAdj ~ NuM * MeanPre2wksTemp,
  # models the precision - coefs need to be exponentiated
  phi ~ NuM * MeanPre2wksTemp,
  # models zero/not-zero process
  zi ~ NuM * MeanPre2wksTemp)

# set priors
# this greatly helps convergence
ZIBR_2015_priors <- get_prior(ZIBR_2015_Formula, 
                              data = chan15,
                              family = zero_inflated_beta())

priorBLAH <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
               set_prior("normal(0, 1)", class = "b"),
               set_prior("normal(0, 1)", class = "b", dpar = "phi"),
               set_prior("normal(0, 1)", class = "b", dpar = "zi"))

### ZIBR_2015_1muX_phiX_ziX ----
# takes a ~5 mins
ZIBR_2015_1muX_phiX_ziX <- brm(
  bf(
    # models the mean - coefs need to be plogis()
    propNfixerAdj ~ NuM * MeanPre2wksTemp,
    # models the precision - coefs need to be exponentiated
    phi ~ NuM * MeanPre2wksTemp,
    # models zero/not-zero process
    zi ~ NuM * MeanPre2wksTemp),
  data = chan15,
  prior = priorBLAH,
  init = "0",
  control = list(adapt_delta = 0.999, max_treedepth = 30),
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4,
  backend = "cmdstanr"
)

ZIBR_2015_1muX_phiX_ziX.loo <- LOO(ZIBR_2015_1muX_phiX_ziX)
tidy(ZIBR_2015_1muX_phiX_ziX, effects = "fixed")
# Mean interaction is not significant so remove to achieve simpler model

### ZIBR_2015_1muA_phiX_ziX ----
# this is considerably faster
ZIBR_2015_1muA_phiX_ziX <- brm(
  bf(
    # models the mean - coefs need to be plogis()
    # removing the interactions
    propNfixerAdj ~ NuM + MeanPre2wksTemp,
    # models the precision - coefs need to be exponentiated
    phi ~ NuM * MeanPre2wksTemp,
    # models zero/not-zero process
    zi ~ NuM * MeanPre2wksTemp),
  data = chan15,
  prior = priorBLAH,
  init = "0",
  control = list(adapt_delta = 0.9999, max_treedepth = 30),
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4,
  backend = "cmdstanr"
)

ZIBR_2015_1muA_phiX_ziX.loo <- LOO(ZIBR_2015_1muA_phiX_ziX)
tidy(ZIBR_2015_1muA_phiX_ziX, effects = "fixed")

function_test <- function(arg_)

library(sjPlot)
tab_model(ZIBR_2015_1muA_phiX_ziX)
## Plot ----
### Coefs ----
  ZIBR_2015_1muA_phiX_ziX_pos_beta <- ZIBR_2015_1muA_phiX_ziX %>% 
    gather_draws(`b_.*`, regex = TRUE) %>% 
    mutate(component = ifelse(str_detect(.variable, "phi_"), "Precision", 
                              ifelse(str_detect(.variable, "zi_"), "ZeroNoneZero", "Mean")),
           intercept = str_detect(.variable, "Intercept"))
  
  ggplot(ZIBR_2015_1muA_phiX_ziX_pos_beta, aes(x = .value, y = fct_rev(.variable), fill = component)) +
    geom_vline(xintercept = 0) +
    stat_halfeye(aes(slab_alpha = intercept), 
                 .width = c(0.8, 0.95), point_interval = "median_hdi") +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    scale_slab_alpha_discrete(range = c(1, 0.4)) +
    guides(fill = "none", slab_alpha = "none") +
    labs(x = "Coefficient", y = "Variable",
         caption = "80% and 95% credible intervals shown in black") +
    facet_wrap(vars(component), ncol = 1, scales = "free") 

### Response to Temp ----
  # takes a bit
  ZIBR_2015_1muA_phiX_ziX_pos_beta_Temp <- ZIBR_2015_1muA_phiX_ziX %>% 
    epred_draws(newdata = expand_grid(NuM = c(0.11, 1.9, 3.68, 7.25, 10.82, 14.4),
                                      MeanPre2wksTemp = seq(9,24, by = 0.5))) %>% 
    mutate(NcompInt = ifelse(NuM == 0.11 | NuM == 3.68,"Comp","NotComp"),
           NcompInt = as.factor(NcompInt),
           NuM = as.factor(NuM),
           NuM = fct_recode(NuM, "0 µM-N" = "0.11",
                          "1.8 µM-N" = "1.9",
                          "3.6 µM-N" = "3.68",
                          "7.1 µM-N" = "7.25",
                          "10.7 µM-N" = "10.82",
                          "14.3 µM-N" = "14.4"),
           NuM = fct_relevel(NuM, "0 µM-N",
                          "1.8 µM-N",
                          "3.6 µM-N",
                          "7.1 µM-N",
                          "10.7 µM-N",
                          "14.3 µM-N"))
  
  ColorVals2015 <- c("black", "#FED976", "#339900", "#FD8D3C", "#F03B20", "#BD0026")
  
 p.PerNfixers2015 <-  ggplot(ZIBR_2015_1muA_phiX_ziX_pos_beta_Temp, 
           aes(y = .epred, x = MeanPre2wksTemp, color = NuM, fill = NuM, size = NcompInt)) +
           geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
           geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
          stat_lineribbon(.width = c(0.95), alpha = 0.1) + #
          stat_lineribbon(.width = c(0), alpha = 1) + #
          scale_color_manual(values = ColorVals2015, name = "Treatment")+
          scale_fill_manual(values = ColorVals2015, name = "Treatment") +
          scale_size_manual(values = c(1.2,0.5)) +
          ylab("") +
          xlab("Temperature (°C)") +
          scale_y_continuous(limits = c(0,1.25), breaks = c(0, 0.25, 0.5, 0.75, 1.00, 1.25)) +
          scale_x_continuous(limits = c(5,25), breaks = c(5, 10, 15, 20, 25)) +
          theme(panel.background = element_rect(fill = "white", color = "white"),
                panel.border = element_rect(color = "black", fill = "NA", size = 1),
                axis.title.x = element_text(size = 22),
                axis.title.y = element_text(size = 60),
                axis.text = element_text(size = 18),
                axis.line = element_line(color = "black", size = 1),
                plot.background = element_rect(fill = "white", color =  "white"),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 24, face = "bold"),
                strip.text.y = element_blank(),
                legend.title = element_text(size = 24, face = "bold"),
                legend.text = element_text(size = 20),
                legend.key.height = unit(1,"cm"),
                legend.position = c(0.15, 0.63),
                legend.background = element_rect(fill = "transparent")) +
          guides(size = "none") 

### Response to N ----
  ZIBR_2015_1muA_phiX_ziX_pos_beta_N <- ZIBR_2015_1muA_phiX_ziX %>% 
    epred_draws(newdata = expand_grid(NuM = seq(0.1,14.5, by = 0.1),
                                      MeanPre2wksTemp = c(9.2, 12.5, 15.8, 19.7, 24.5)))
  
ggplot(ZIBR_2015_1muA_phiX_ziX_pos_beta_N %>% 
         mutate(MeanPre2wksTempF= as.factor(MeanPre2wksTemp)), aes(y = .epred, x = NuM,
                                                                   color = MeanPre2wksTempF,
                                                                   fill = MeanPre2wksTempF)) +
  stat_lineribbon(.width = c(0.95), alpha = 0.5)


### Percent change with N and temp ----

ZIBR_2015_1muA_phiX_ziX_pos_beta_Temp_PerChange <- ZIBR_2015_1muA_phiX_ziX %>% 
  epred_draws(newdata = expand_grid(NuM = c(0.11, 3.68),
                                    MeanPre2wksTemp = c(11,21))) %>% 
  group_by(NuM, MeanPre2wksTemp) %>% 
  summarise(PerNfixers = mean(.epred, na.rm = T)) %>% 
  mutate(Ntemp = paste0(NuM, "_", MeanPre2wksTemp)) %>% 
  pivot_wider(names_from = Ntemp, values_from = PerNfixers)

# Temp, low N
(0.586-0.225)/0.225*100

# Temp, high N
(0.264 - 0.0890)/0.0890*100
# N, cold
(0.0890 - 0.225)/0.225*100
# N, warm
(0.264 - 0.586)/0.586*100

# 2016 ----
# 13 zeros - these are 100% N fixers
chan16 %>% 
  count(propNOT_Nfixer == 0) %>% 
  mutate(prop = n/sum(n))

# No 100% non-n fixers
chan16 %>% 
  count(propNOT_Nfixer == 1) %>% 
  mutate(prop = n/sum(n))

## ZIBR_2016_1muX_phiX_ziX ----
  # set priors
  # this greatly helps convergence
  # Look at the default priors
  ZIBR_2016_Formula <- bf(
    # models the mean - coefs need to be plogis()
    propNOT_Nfixer ~ PuM * MeanPre2wksTemp,
    # models the precision - coefs need to be exponentiated
    phi ~ PuM * MeanPre2wksTemp,
    # models zero/not-zero process
    zi ~ PuM * MeanPre2wksTemp)
  

  ZIBR_2016_priors <- get_prior(ZIBR_2016_Formula, 
                                data = chan16,
                                family = zero_inflated_beta())
  
  priorBLAH <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
                 set_prior("normal(0, 1)", class = "b"),
                 set_prior("normal(0, 1)", class = "b", dpar = "phi"),
                 set_prior("normal(0, 1)", class = "b", dpar = "zi"))



ZIBR_2016_muX_phiX_ziX <- brm(
  bf(
    # models the mean - coefs need to be plogis()
    propNOT_Nfixer ~ PuM * MeanPre2wksTemp,
    # models the precision - coefs need to be exponentiated
    phi ~ PuM * MeanPre2wksTemp,
    # models zero/not-zero process
    zi ~ PuM * MeanPre2wksTemp),
  data = chan16,
  prior = priorBLAH,
  init = "0",
  control = list(adapt_delta = 0.999, max_treedepth = 30),
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4,
  backend = "cmdstanr"
)

ZIBR_2016_muX_phiX_ziX.loo <- LOO(ZIBR_2016_muX_phiX_ziX)
tidy(ZIBR_2016_muX_phiX_ziX, effects = "fixed")
# no interactions are significant, so turn off

## ZIBR_2016_muA_phiA_ziA ----
ZIBR_2016_muA_phiA_ziA <- brm(
  bf(
    # models the mean - coefs need to be plogis()
    propNOT_Nfixer ~ PuM + MeanPre2wksTemp,
    # models the precision - coefs need to be exponentiated
    phi ~ PuM + MeanPre2wksTemp,
    # models zero/not-zero process
    zi ~ PuM + MeanPre2wksTemp),
  data = chan16,
  prior = priorBLAH,
  init = "0",
  control = list(adapt_delta = 0.999, max_treedepth = 30),
  family = zero_inflated_beta(),
  chains = 4, iter = 4000, warmup = 1000,
  cores = 4,
  backend = "cmdstanr"
)

ZIBR_2016_muA_phiA_ziA.loo <- LOO(ZIBR_2016_muA_phiA_ziA)
tidy(ZIBR_2016_muA_phiA_ziA, effects = "fixed")
# SOME P EFFECT ON ZI
# temperature influences the zi parameter
tab_model(ZIBR_2016_muA_phiA_ziA)

ZIBR_2016_muA_phiA_ziA_2 <- add_criterion(ZIBR_2016_muA_phiA_ziA, "waic")
ZIBR_2016_muX_phiX_ziX_2 <- add_criterion(ZIBR_2016_muX_phiX_ziX, "waic")
loo_compare(ZIBR_2016_muA_phiA_ziA_2, ZIBR_2016_muX_phiX_ziX_2, criterion = "waic")
## Plot ----
  ### Coefs ----
ZIBR_2016_muA_phiA_ziA_pos_beta <- ZIBR_2016_muA_phiA_ziA %>% 
    gather_draws(`b_.*`, regex = TRUE) %>% 
    mutate(component = ifelse(str_detect(.variable, "phi_"), "Precision", 
                              ifelse(str_detect(.variable, "zi_"), "ZeroNoneZero", "Mean")),
           intercept = str_detect(.variable, "Intercept"))
  
  ggplot(ZIBR_2016_muA_phiA_ziA_pos_beta, aes(x = .value, y = fct_rev(.variable), fill = component)) +
    geom_vline(xintercept = 0) +
    stat_halfeye(aes(slab_alpha = intercept), 
                 .width = c(0.8, 0.95), point_interval = "median_hdi") +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    scale_slab_alpha_discrete(range = c(1, 0.4)) +
    guides(fill = "none", slab_alpha = "none") +
    labs(x = "Coefficient", y = "Variable",
         caption = "80% and 95% credible intervals shown in black") +
    facet_wrap(vars(component), ncol = 1, scales = "free") 
  
  ### Response to Temp ----
  # takes a bit
  ZIBR_2016_1muA_phiX_ziX_pos_beta_Temp <- ZIBR_2016_muA_phiA_ziA %>% 
    epred_draws(newdata = expand_grid(PuM = c(1.97, 6.81, 3.59, 1.17, 0.36, 5.20),
                                      MeanPre2wksTemp = seq(9,24, by = 0.5))) %>% 
              mutate(.epred2 = 1 - .epred,
                     PuM = as.factor(PuM),
                     PuM = fct_recode(PuM, "0 µM-P" = "0.36",
                                            "0.8 µM-P" = "1.17",
                                      "1.6 µM-P" = "1.97",
                                      "3.6 µM-P" = "3.59",
                                      "5.8 µM-P" = "5.2",
                                      "6.5 µM-P" = "6.81"),
                     PuM = fct_relevel(PuM, "0 µM-P",
                                      "0.8 µM-P",
                                      "1.6 µM-P",
                                      "3.6 µM-P",
                                      "5.8 µM-P",
                                      "6.5 µM-P")) %>% 
    mutate(NcompInt = ifelse(PuM == "0 µM-P" | PuM == "3.6 µM-P","Comp","NotComp"),
           NcompInt = as.factor(NcompInt))
  
  ColorVals2016 <- c("black", "#C7E9B4", "#7FCDBB", "#339900", "#2C7FB8", "#253494")
  
  p.PerNfixer2016 <- ggplot(ZIBR_2016_1muA_phiX_ziX_pos_beta_Temp %>% 
                              mutate(PuMF = as.factor(PuM)), 
                            aes(y = .epred2, x = MeanPre2wksTemp, color = PuMF, fill = PuMF, size = NcompInt)) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
    stat_lineribbon(.width = c(0.95), alpha = 0.1) +
    stat_lineribbon(.width = c(0), alpha = 1) +
    scale_color_manual(values = ColorVals2016, name = "Treatment")+
    scale_fill_manual(values = ColorVals2016, name = "Treatment") +
    ylab("") +
    xlab("Temperature (°C)") +
    scale_y_continuous(limits = c(0,1.25), breaks = c(0, 0.25, 0.5, 0.75, 1.00, 1.25)) +
    scale_x_continuous(limits = c(5,25), breaks = c(5, 10, 15, 20, 25)) +
    scale_size_manual(values = c(1.2,0.5)) +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 1),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 60),
          axis.text = element_text(size = 18),
          axis.line = element_line(color = "black", size = 1),
          plot.background = element_rect(fill = "white", color =  "white"),
          strip.background = element_blank(),
          # strip.background = element_rect(fill = "grey", color =  "white"),
          strip.text.x = element_text(size = 24, face = "bold"),
          # strip.text.y = element_text(size = 28, face = "bold"),
          strip.text.y = element_blank(),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 20),
          legend.key.height = unit(1,"cm"),
          legend.position = c(0.15, 0.4),
          legend.background = element_rect(fill = "transparent")) +
    guides(size = "none") 
  
  ### Response to N ----
  ZIBR_2016_muA_phiA_ziA_pos_beta_N <- ZIBR_2016_muA_phiA_ziA %>% 
    epred_draws(newdata = expand_grid(PuM = seq(1.97,5.2, by = 0.1),
                                      MeanPre2wksTemp = c(11.0, 13.8, 17.3, 20.7, 25.5)))%>% 
    mutate(.epred2 = 1 - .epred)
  
  ggplot(ZIBR_2016_muA_phiA_ziA_pos_beta_N %>% 
           mutate(MeanPre2wksTempF= as.factor(MeanPre2wksTemp)), aes(y = .epred2, x = PuM,
                                                                     color = MeanPre2wksTempF,
                                                                     fill = MeanPre2wksTempF)) +
    stat_lineribbon(.width = c(0.95), alpha = 0.5)
  
  
  ZIBR_2016_1muA_phiX_ziX_pos_beta_Temp_PerChange <- ZIBR_2016_muA_phiA_ziA  %>% 
    epred_draws(newdata = expand_grid(PuM = c(0.11, 3.68),
                                      MeanPre2wksTemp = c(11,21))) %>% 
    mutate(.epred2 = 1 - .epred) %>% 
    group_by(PuM, MeanPre2wksTemp) %>% 
    summarise(PerNfixers = mean(.epred2, na.rm = T)) %>% 
    mutate(Ptemp = paste0(PuM, "_", MeanPre2wksTemp)) 
  
  # Temp, low N
  (0.921-0.955)/0.955*100
  # Temp, high N
  (0.921 - 0.982)/0.982*100
  # N, cold
  (0.982 - 0.955)/0.955*100
  # N, warm
  (0.921 - 0.921)/0.921*100

  # 2017 ----
  # about a third of the values are zero (i.e., no N fixers)
  chan17 %>% 
    count(propNfixer == 0) %>% 
    mutate(prop = n/sum(n))
  
  # There are 4 ones
  # I converted these to 0.99 in propNfixer_adj
  chan17 %>% 
    count(propNfixer == 1) %>% 
    mutate(prop = n/sum(n))
  
  
  # Look at the default priors
  ZIBR_2017_Formula <- bf(
    # models the mean - coefs need to be plogis()
    propNfixerAdj ~ NP * MeanPre2wksTemp,
    # models the precision - coefs need to be exponentiated
    phi ~ NP * MeanPre2wksTemp,
    # models zero/not-zero process
    zi ~ NP * MeanPre2wksTemp)
  
  # set priors
  # this greatly helps convergence
  ZIBR_2017_priors <- get_prior(ZIBR_2017_Formula, 
                                data = chan17,
                                family = zero_inflated_beta())
  
  priorBLAH <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
                 set_prior("normal(0, 1)", class = "b"),
                 set_prior("normal(0, 1)", class = "b", dpar = "phi"),
                 set_prior("normal(0, 1)", class = "b", dpar = "zi"))
  
  # approach: Fit N and N:P model and compare with LOO
  # Need to take a different approach with these bayesian models
  
  ### ZIBR_2017_1muNP_phiNP_ziNP ----
  ZIBR_2017_1muNP_phiNP_ziNP <- brm(
    bf(
      # models the mean - coefs need to be plogis()
      propNfixerAdj ~ NP * MeanPre2wksTemp,
      # models the precision - coefs need to be exponentiated
      # switched to additive b/c CIs offered no support
      phi ~ NP + MeanPre2wksTemp,
      # models zero/not-zero process
      # switched to additive b/c CIs offered no support
      zi ~ NP + MeanPre2wksTemp),
    data = chan17,
    prior = priorBLAH,
    init = "0",
    control = list(adapt_delta = 0.999, max_treedepth = 30),
    family = zero_inflated_beta(),
    chains = 4, iter = 4000, warmup = 1000,
    cores = 4,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE)
  )
  
  ZIBR_2017_1muNP_phiNP_ziNP.loo <- LOO(ZIBR_2017_1muNP_phiNP_ziNP, moment_match = TRUE)
  tidy(ZIBR_2017_1muNP_phiNP_ziNP, effects = "fixed")
  
  ### ZIBR_2017_1muN_phiN_ziN ----
  ZIBR_2017_1muN_phiN_ziN <- brm(
    bf(
      # models the mean - coefs need to be plogis()
      propNfixerAdj ~ NuM_f * MeanPre2wksTemp,
      # models the precision - coefs need to be exponentiated
      # switched to additive b/c CIs offered no support
      phi ~ NuM_f + MeanPre2wksTemp,
      # models zero/not-zero process
      # switched to additive b/c CIs offered no support
      zi ~ NuM_f + MeanPre2wksTemp),
    data = chan17,
    prior = priorBLAH,
    init = "0",
    control = list(adapt_delta = 0.999, max_treedepth = 30),
    family = zero_inflated_beta(),
    chains = 4, iter = 4000, warmup = 1000,
    cores = 4,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE)
  )
  
  ZIBR_2017_1muN_phiN_ziN.loo <- LOO(ZIBR_2017_1muN_phiN_ziN, moment_match = TRUE)
  tidy(ZIBR_2017_1muN_phiN_ziN, effects = "fixed")
  
  # elpd_diff is larger than se difference indicating that the N model has better predictive performance
  # https://mc-stan.org/loo/articles/loo2-example.html#comparing-the-models-on-expected-log-predictive-density-1
  loo::loo_compare(ZIBR_2017_1muNP_phiNP_ziNP.loo, ZIBR_2017_1muN_phiN_ziN.loo)
  tab_model(ZIBR_2017_1muN_phiN_ziN)
  
## Plot ----
  ### Coefs ----
  ZIBR_2017_1muN_phiN_ziN_pos_beta <- ZIBR_2017_1muN_phiN_ziN %>% 
    gather_draws(`b_.*`, regex = TRUE) %>% 
    mutate(component = ifelse(str_detect(.variable, "phi_"), "Precision", 
                              ifelse(str_detect(.variable, "zi_"), "ZeroNoneZero", "Mean")),
           intercept = str_detect(.variable, "Intercept"))
  
  ggplot(ZIBR_2017_1muN_phiN_ziN_pos_beta, aes(x = .value, y = fct_rev(.variable), fill = component)) +
    geom_vline(xintercept = 0) +
    stat_halfeye(aes(slab_alpha = intercept), 
                 .width = c(0.8, 0.95), point_interval = "median_hdi") +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    scale_slab_alpha_discrete(range = c(1, 0.4)) +
    guides(fill = "none", slab_alpha = "none") +
    labs(x = "Coefficient", y = "Variable",
         caption = "80% and 95% credible intervals shown in black") +
    facet_wrap(vars(component), ncol = 1, scales = "free") 
  
  ### Response to Temp ----
  # takes a bit
  ZIBR_2017_1muN_phiN_ziN_pos_beta_Temp <- ZIBR_2017_1muN_phiN_ziN %>% 
    epred_draws(newdata = expand_grid(NuM_f = c("0.11", "3.68"),
                                      MeanPre2wksTemp = seq(7.9, 21, by = 0.5))) %>% 
    mutate(NcompInt = ifelse(NuM_f == 0 | NuM_f == 3.6,"Comp","NotComp"),
           NcompInt = as.factor(NcompInt),
           NuM_f = as.factor(NuM_f),
           NuM_f = fct_recode(NuM_f, "0 µM-N" = "0.11",
                              "3.6 µM-N" = "3.68"),
           NuM_f = fct_relevel(NuM_f, "0 µM-N",
                              "3.6 µM-N"))
  
  ColorVals2017 <- c("black", "#339900")
  
  p.PerNfixer2017 <- ggplot(ZIBR_2017_1muN_phiN_ziN_pos_beta_Temp, 
                            aes(y = .epred, x = MeanPre2wksTemp, color = NuM_f, fill = NuM_f,size = NcompInt)) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
    stat_lineribbon(.width = c(0.95), alpha = 0.1) +
    stat_lineribbon(.width = c(0), alpha = 1) +
    scale_color_manual(values = ColorVals2017, name = "Treatment")+
    scale_fill_manual(values = ColorVals2017, name = "Treatment") +
    ylab(expression(paste("Proportion ", N[2],"-fixers"))) +
    xlab("Temperature (°C)") +
    scale_y_continuous(limits = c(0,1.25), breaks = c(0, 0.25, 0.5, 0.75, 1.00, 1.25)) +
    scale_x_continuous(limits = c(5,25), breaks = c(5, 10, 15, 20, 25)) +
    scale_size_manual(values = c(1.2,0.5)) +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 1),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          axis.text = element_text(size = 18),
          axis.line = element_line(color = "black", size = 1),
          plot.background = element_rect(fill = "white", color =  "white"),
          strip.background = element_blank(),
          # strip.background = element_rect(fill = "grey", color =  "white"),
          strip.text.x = element_text(size = 24, face = "bold"),
          # strip.text.y = element_text(size = 28, face = "bold"),
          strip.text.y = element_blank(),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 20),
          legend.key.height = unit(1,"cm"),
          legend.position = c(0.15, 0.75),
          legend.background = element_rect(fill = "transparent")) +
    guides(size = "none") 
  
  ### Response to N ----
  ZIBR_2017_1muN_phiN_ziN_pos_beta_N <- ZIBR_2017_1muN_phiN_ziN %>% 
    epred_draws(newdata = expand_grid(NuM_f = c("0.11", "3.68"),
                                      MeanPre2wksTemp = seq(7.9, 19.8, by = 0.5)))
  
  ggplot(ZIBR_2017_1muN_phiN_ziN_pos_beta_N %>% 
           mutate(MeanPre2wksTempF= as.factor(MeanPre2wksTemp)), aes(y = .epred, x = NuM_f,
                                                                     color = MeanPre2wksTempF,
                                                                     fill = MeanPre2wksTempF)) +
          stat_lineribbon(.width = c(0.95), alpha = 0.5) 
  
  # Print ----
  p1.gtf <- gtable_frame(ggplotGrob(p.PerNfixers2015), width = unit(1, "null"), height = unit(1, "null"))
  p2.gtf <- gtable_frame(ggplotGrob(p.PerNfixer2016), width = unit(1, "null"), height = unit(1, "null"))
  # p3.gtf <- gtable_frame(ggplotGrob(p.PerNfixer2017), width = unit(1, "null"), height = unit(1, "null"))
  
  p123.gtf <- gtable_frame(gtable_rbind(p1.gtf, p2.gtf), width = unit(3,"null"), height = unit(1,"null"))
  
  png("05_Figures4MS/05_Fig4_2015to16_PerNfixer2.png", units = "in", height = 12, width = 8, res = 300)
  grid.newpage()
  grid.draw(p123.gtf)
  grid.text(expression(paste("Proportion ", N[2],"-fixers")), x = unit(0.04,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 36, face = "bold"), rot = 90)
  grid.text("a) N-only", x = unit(0.31,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 28, fontface = "bold"))
  grid.text("b) P-only", x = unit(0.31,"npc"), y = unit(0.47,"npc"), gp=gpar(fontsize = 28, fontface = "bold"))
  # grid.text("c) N + P", x = unit(0.31,"npc"), y = unit(0.315,"npc"), gp=gpar(fontsize = 28, fontface = "bold"))
  dev.off()
  
  
  p3.gtf <- gtable_frame(ggplotGrob(p.PerNfixer2017), width = unit(1, "null"), height = unit(1, "null"))
  

  
  png("05_Figures4MS/05_FigS10_2017_PerNfixer2.png", units = "in", height = 6, width = 8, res = 300)
  p.PerNfixer2017
  dev.off()

  
# Table ----
  Tab_2017 <- tidy(ZIBR_2017_1muN_phiN_ziN, effects = "fixed") %>% 
    mutate(Year = "2017",
           component = ifelse(grepl("phi", term), "phi",component))
  
  # NOTE: PROPORTION NON-N2 FIXERS UNLIKE OTHER MODELS
  Tab_2016 <- tidy(ZIBR_2016_muA_phiA_ziA, effects = "fixed")%>% 
    mutate(Year = "2016",
           component = ifelse(grepl("phi", term), "phi",component))
  
  Tab_2015 <- tidy(ZIBR_2015_1muA_phiX_ziX, effects = "fixed")%>% 
    mutate(Year = "2015",
           component = ifelse(grepl("phi", term), "phi",component))
  
  Tab_PorNfixer <- rbind(Tab_2015, Tab_2016, Tab_2017) %>% 
    mutate(estimate_t = case_when(component %in% c("cond", "zi") ~ inv.logit(estimate),
                                  component == "phi" ~ exp(estimate)),
           conf.low_t = case_when(component %in% c("cond", "zi") ~ inv.logit(conf.low),
                                  component == "phi" ~ exp(conf.low)),
           conf.high_t = case_when(component %in% c("cond", "zi") ~ inv.logit(conf.high),
                                   component == "phi" ~ exp(conf.high))) %>% 
    select(Year, component, term, estimate_t, conf.low_t, conf.high_t) %>% 
    mutate(term2 = case_when(term == "(Intercept)" ~ "Intercept",
                             term == "phi_(Intercept)" ~ "Intercept",
                             term == "NuM" ~ "Nitrogen",
                             term == "phi_NuM" ~ "Nitrogen",
                             term == "phi_NuM:MeanPre2wksTemp" ~ "NitrogenXTemperature",
                             term == "NuM:MeanPre2wksTemp" ~ "NitrogenXTemperature",
                             term == "PuM" ~ "Phosphorus",
                             term == "MeanPre2wksTemp" ~ "Temperature",
                             term == "phi_PuM" ~ "Phosphorus",
                             term == "phi_MeanPre2wksTemp" ~ "Temperature",
                             term == "NuM_f3.68" ~ "Nitrogen",
                             term == "NuM_f3.68:MeanPre2wksTemp" ~ "NitrogenXTemperature",
                             term == "phi_NuM_f3.68" ~ "Nitrogen"),
           value = paste0(round(estimate_t,2)," (", round(conf.low_t,2),"-", round(conf.high_t,2),")")) %>% 
    select(Year, component, term2, value) %>% 
    arrange(Year, component, term2) %>% 
    mutate(term3 = case_when(term2 == "Nitrogen" ~ "Nutrient",
                             term2 == "Phosphorus" ~ "Nutrient",
                             term2 == "NitrogenXTemperature" ~ "NutrientXTemperature",
                             !(term2 %in% c("Nitrogen", "Phosphorus", "NitrogenXTemperature")) ~ term2)) %>% 
    pivot_wider(id_cols = c(component, term3), names_from = "Year", values_from = "value") %>% 
    arrange(component, term3)
  
  write.csv(Tab_PorNfixer, "04_Tables4MS/05_PropNfixers_Table.csv")

# save image ----
  #NOTE: RDAT is not in Repo due to its size. Available upon request.
  OutOfRepoPath <- "/Users/hood.211/Library/CloudStorage/OneDrive-TheOhioStateUniversity/IcelandChannels_MS/20211016_EcoMongrphs/R_code_20221008_REDO/TooLarge4SumbissionUpload/"
  # save.image(paste0(OutOfRepoPath,"05_PropNfix_2015_2016_Rdat"))
  # load(paste0(OutOfRepoPath,"05_PropNfix_2015_2016_Rdat"))


