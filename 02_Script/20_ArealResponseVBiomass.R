# relationships between areal responses and biomass
# JMH, Nov 2022

# Libraries ----
# general
library(tidyverse)

# tables
library(modelsummary)

# Data ----
chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv", row.names = 1) %>% 
  mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
         # Adding 2015 minimum to Nuptake to allow modeling
         NUp_uM_N_m2_hr = NUp_uM_N_m2_hr + 202,
         TotAssim_uM_N_m2_hr = NUp_uM_N_m2_hr + Nfix_uM_N_m2h) %>% 
  select(MetDate:StartDate,
         NEP_uM_C_m2h, R_uM_C_m2h, GPP_uM_C_m2h, Met_gAFDM_m2,
         NUp_uM_N_m2_hr, Pup_gAFDM_m2,
         Nfix_uM_N_m2h, BMS_N2fixers_g_m2, 
         TotAssim_uM_N_m2_hr) %>% 
  rowwise() %>% 
  mutate(AvgBMnoNFix = mean(c(Pup_gAFDM_m2, Met_gAFDM_m2), na.rm = T)) %>% 
  ungroup()

met <- chan %>% 
  select(Year, sample_event, channel, tempF, N_uM, P_uM, NPratio,
         NEP_uM_C_m2h, R_uM_C_m2h, GPP_uM_C_m2h, Met_gAFDM_m2) %>% 
  pivot_longer(cols = NEP_uM_C_m2h:GPP_uM_C_m2h) %>% 
  mutate(logValue = log(value + 0.01),
         logBM = log(Met_gAFDM_m2 + 0.01)) %>% 
  select(Year:NPratio, name, logBM, logValue)

Nstuff <- chan %>% 
  select(Year, sample_event, channel, tempF, N_uM, P_uM, NPratio,
         NUp_uM_N_m2_hr, Pup_gAFDM_m2,
         Nfix_uM_N_m2h, BMS_N2fixers_g_m2, 
         TotAssim_uM_N_m2_hr,
         AvgBMnoNFix) %>% 
  mutate(BMS_N2fixers_g_m2_New = ifelse(BMS_N2fixers_g_m2 == 0, AvgBMnoNFix, BMS_N2fixers_g_m2)) %>% 
  pivot_longer(cols = c(NUp_uM_N_m2_hr, Nfix_uM_N_m2h, TotAssim_uM_N_m2_hr)) %>% 
  rowwise() %>% 
  mutate(AvgBM = mean(c(Pup_gAFDM_m2, BMS_N2fixers_g_m2_New), na.rm = T)) %>% 
  ungroup() %>% 
  mutate(logBM = log(ifelse(name == "NUp_uM_N_m2_hr", Pup_gAFDM_m2,
                            ifelse(name == "Nfix_uM_N_m2h", BMS_N2fixers_g_m2_New,
                                   ifelse(name == "TotAssim_uM_N_m2_hr", AvgBM, AvgBM))) + 0.01),
         logValue = log(value + 0.01)) %>% 
  select(Year:NPratio, name,  logBM, logValue)


chan2 <- rbind(met, Nstuff) %>% 
  mutate(NameYr = paste0(name, "_B_", Year)) %>% 
  filter(complete.cases(.))

ggplot(chan2, aes(y = logValue, x = logBM)) +
  geom_point()+
  stat_smooth(method = "lm") +
  facet_grid(name ~ Year)

chan2.lm <- chan2 %>% 
  group_by(name, Year) %>% 
  mutate(N_uMf = as.character(N_uM)) %>% 
  nest() %>% 
  mutate(ValueBM_fit = map(data, ~lm(logValue ~ logBM, data = .x))) %>% 
  # remove N fix
  filter(name != "Nfix_uM_N_m2h")%>% 
  mutate(ID = paste0(name, " | ", Year))


chan2.lm_PuMN_2015 <- chan2 %>% 
  # filter(name == "Nfix_uM_N_m2h") %>%
  filter(Year == "2015") %>% 
  group_by(name, Year) %>% 
  mutate(N_uMf = as.character(N_uM)) %>% 
  nest() %>% 
  mutate(ValueBM_fit = map(data, ~lm(logValue ~ logBM * N_uMf, data = .x))) %>% 
  mutate(ID = paste0(name, " | ", Year, "_interaction"))

chan2.lm_PuMN_2016 <- chan2 %>% 
  filter(name == "Nfix_uM_N_m2h") %>%
  filter(Year == "2016") %>% 
  group_by(name, Year) %>% 
  nest() %>% 
  mutate(ValueBM_fit = map(data, ~lm(logValue ~ logBM, data = .x)))  %>% 
  mutate(ID = paste0(name, " | ", Year, "_interaction"))

chan2.lm_PuMN_2017 <- chan2 %>% 
  # filter(name == "Nfix_uM_N_m2h") %>%
  filter(Year == "2017") %>% 
  group_by(name, Year) %>% 
  mutate(N_uMf = as.character(N_uM)) %>% 
  nest() %>% 
  mutate(ValueBM_fit = map(data, ~lm(logValue ~ logBM * N_uMf, data = .x)))  %>% 
  mutate(ID = paste0(name, " | ", Year, "_interaction"))

chan2.lm2 <- rbind(chan2.lm, chan2.lm_PuMN_2015, chan2.lm_PuMN_2016, chan2.lm_PuMN_2017)


# Table ----
Mods <- chan2.lm2$ValueBM_fit
names(Mods) <- unique(chan2.lm2$ID)

# make CSV then convert to table
modelsummary(Mods, 
             stars = TRUE,
             statistic = NULL,
             # coef_map = cm,
             output = "04_Tables4MS/20_ArealResponseByBiomass.csv")

# Figure ----
# https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html

FluxLabels <- c("NEP", "ER", "GPP", "N uptake", "N fixation", "Total N Assim.")
names(FluxLabels) <-  c("NEP_uM_C_m2h", "R_uM_C_m2h", "GPP_uM_C_m2h", "NUp_uM_N_m2_hr", "Nfix_uM_N_m2h", "TotAssim_uM_N_m2_hr")

ExpLabels <- c("N-only experiment", "P-only experiment", "N+P experiment")
names(ExpLabels) <- c("2015", "2016", "2017")

png("05_Figures4MS/20_FluxVersusBiomass_FigS8.jpg", units = "in", height = 17, width = 15, res = 300)
ggplot(chan2 %>% 
         mutate(name = as.factor(name),
                name = fct_relevel(name,
                                   "GPP_uM_C_m2h", "R_uM_C_m2h", "NEP_uM_C_m2h", "NUp_uM_N_m2_hr", "Nfix_uM_N_m2h", "TotAssim_uM_N_m2_hr")), 
       aes(y = logValue, x = logBM)) +
  geom_point()+
  stat_smooth(method = "lm") +
  facet_grid(name ~ Year, 
             labeller = labeller(name = FluxLabels,
                                 Year = ExpLabels),
             scales = "free_y") +
  ylab(expression(paste(log[e]," transformed flux"))) +
  xlab(expression(paste(log[e], " Biofilm biomass (g AFDM ",m^-2,")")))+
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(color = "black", fill = "NA", linewidth = 1),
        panel.spacing.x = unit(15, "pt"),
        panel.spacing.y = unit(15, "pt"),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        axis.text = element_text(size = 22),
        axis.line = element_line(color = "black", linewidth = 1),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color =  "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24, face = "bold"),
        legend.key.height = unit(1,"cm")) 
dev.off()
