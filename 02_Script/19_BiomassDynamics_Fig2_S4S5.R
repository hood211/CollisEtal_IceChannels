# JMH 2023
# Generates Fig. 2, S4, and S5
# update 2 Jul 24

library(tidyverse)

# Data ---- 

## smp dates ----
  D1_2015_met <- as.Date("2015-07-15", format = "%Y-%m-%d")
  D2_2015_met <- as.Date("2015-07-27", format = "%Y-%m-%d")
  D1_2015_up <- as.Date("2015-07-20", format = "%Y-%m-%d")
  D2_2015_up <- as.Date("2015-07-28", format = "%Y-%m-%d")
  # these are split over 2 days, took second
  D1_2015_nfix <- as.Date("2015-07-17", format = "%Y-%m-%d")
  D2_2015_nfix <- as.Date("2015-07-30", format = "%Y-%m-%d")
  
  D1_2016_met <- as.Date("2016-07-27", format = "%Y-%m-%d")
  D2_2016_met <- as.Date("2016-08-02", format = "%Y-%m-%d")
  D1_2016_up <- as.Date("2016-07-28", format = "%Y-%m-%d") 
  D2_2016_up <- as.Date("2016-08-03", format = "%Y-%m-%d")
  # these are split over 2 days, took second
  D1_2016_nfix <- as.Date("2016-07-26", format = "%Y-%m-%d")
  D2_2016_nfix <- as.Date("2016-08-01", format = "%Y-%m-%d")
  
  D1_2017_met <- as.Date("2017-07-11", format = "%Y-%m-%d")
  D2_2017_met <- as.Date("2017-07-22", format = "%Y-%m-%d")
  D1_2017_up <- as.Date("2017-07-12", format = "%Y-%m-%d")
  D2_2017_up <- as.Date("2017-07-23", format = "%Y-%m-%d")
  # these are split over 2 days, took second
  D1_2017_nfix <- as.Date("2017-07-15", format = "%Y-%m-%d")
  D2_2017_nfix <- as.Date("2017-07-21", format = "%Y-%m-%d")

BMsmp2015 <- read.csv("01_Data/2015ChannelAFDM_cleaned.csv", row.names = 1) %>% 
              pivot_longer(cols = c(Met_gAFDM_m2:Nfix_gAFDM_m2), names_to = "Response", values_to = "gAFDM_m2") %>% 
              mutate(MetDate2 = case_when(MetDate == "D1" & Response == "Met_gAFDM_m2" ~ D1_2015_met,
                                          MetDate == "D2" & Response == "Met_gAFDM_m2" ~ D2_2015_met,
                                          MetDate == "D1" & (Response == "Nup_gAFDM_m2" | Response == "Pup_gAFDM_m2") ~ D1_2015_up,
                                          MetDate == "D2" & (Response == "Nup_gAFDM_m2" | Response == "Pup_gAFDM_m2") ~ D2_2015_up,
                                          MetDate == "D1" & Response == "Nfix_gAFDM_m2" ~ D1_2015_nfix,
                                          MetDate == "D2" & Response == "Nfix_gAFDM_m2" ~ D2_2015_nfix))

BMsmp2016 <- read.csv("01_Data/2016ChannelAFDM_cleaned.csv", row.names = 1) %>% 
  pivot_longer(cols = c(Met_gAFDM_m2:Nfix_gAFDM_m2), names_to = "Response", values_to = "gAFDM_m2") %>% 
  mutate(MetDate2 = case_when(ExpDay2 == "D1" & Response == "Met_gAFDM_m2" ~ D1_2016_met,
                              ExpDay2 == "D2" & Response == "Met_gAFDM_m2" ~ D2_2016_met,
                              ExpDay2 == "D1" & (Response == "Nup_gAFDM_m2" | Response == "Pup_gAFDM_m2") ~ D1_2016_up,
                              ExpDay2 == "D2" & (Response == "Nup_gAFDM_m2" | Response == "Pup_gAFDM_m2") ~ D2_2016_up,
                              ExpDay2 == "D1" & Response == "Nfix_gAFDM_m2" ~ D1_2016_nfix,
                              ExpDay2 == "D2" & Response == "Nfix_gAFDM_m2" ~ D2_2016_nfix))

# no n fix afdm on this date
BMsmp2017 <- read.csv("01_Data/2017ChannelAFDM_cleaned.csv", row.names = 1)%>% 
  pivot_longer(cols = c(Met_gAFDM_m2:Pup_gAFDM_m2), names_to = "Response", values_to = "gAFDM_m2") %>% 
  mutate(MetDate2 = case_when(ExpDay2 == "D1" & Response == "Met_gAFDM_m2" ~ D1_2017_met,
                              ExpDay2 == "D2" & Response == "Met_gAFDM_m2" ~ D2_2017_met,
                              ExpDay2 == "D1" & (Response == "Nup_gAFDM_m2" | Response == "Pup_gAFDM_m2") ~ D1_2017_up,
                              ExpDay2 == "D2" & (Response == "Nup_gAFDM_m2" | Response == "Pup_gAFDM_m2") ~ D2_2017_up))

## time series
BioDynExpCodes <- read.csv("01_Data/00_DefinBMts.csv", row.names = 1) %>% 
                      select(Y, channel, TempF, nitrogen_ugL, phosphorus_ugL) %>% 
                      distinct()

ExpDays2015_date <- seq(as.Date("2015-06-04"), as.Date("2015-08-01"), by = "day")
ExpDays2015_days <- seq(seq(1,59, by = 1))
ExpDays2015 <- as.data.frame(cbind(ExpDays2015_date, ExpDays2015_days)) %>% 
                  mutate(ExpDays2015_date = as.Date(ExpDays2015_date))
names(ExpDays2015) <- c("date", "ExpDays")


ExpDays2016_date <- seq(as.Date("2016-06-03"), as.Date("2016-08-03"), by = "day")
ExpDays2016_days <- seq(seq(1,62, by = 1))
ExpDays2016 <- as.data.frame(cbind(ExpDays2016_date, ExpDays2016_days)) %>% 
  mutate(ExpDays2016_date = as.Date(ExpDays2016_date))
names(ExpDays2016) <- c("date", "ExpDays")

ExpDays2017_date <- seq(as.Date("2017-05-23"), as.Date("2017-07-31"), by = "day")
ExpDays2017_days <- seq(seq(1,70, by = 1))
ExpDays2017 <- as.data.frame(cbind(ExpDays2017_date, ExpDays2017_days)) %>% 
  mutate(ExpDays2017_date = as.Date(ExpDays2017_date))
names(ExpDays2017) <- c("date", "ExpDays")

ExpDays <- rbind(ExpDays2015, ExpDays2016, ExpDays2017)

BioDyn <- read.csv("01_Data/00_DefinBMts.csv", row.names = 1) %>% 
  select(-c(ID, Y_chan, date_chan, afdm_g_m2.la, TempF:phosphorus_ugL, TempK, Y, ExpDay)) %>% 
  filter(!is.na(afdm_g_m2)) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>% 
  full_join(BMsmp2015 %>% 
              select(-MetDate) %>% 
              rename(gAFDM_m2_2015 = gAFDM_m2), by = c("date" = "MetDate2", "channel" = "channel")) %>% 
  mutate(Response = ifelse(is.na(Response), "TS", Response)) %>% 
  full_join(BMsmp2016 %>% 
              select(-ExpDay2, -Year) %>% 
              rename(gAFDM_m2_2016 = gAFDM_m2), by = c("date" = "MetDate2", "channel" = "channel", "Response")) %>% 
  full_join(BMsmp2017 %>% 
              select(-ExpDay2, -Year) %>% 
              rename(gAFDM_m2_2017 = gAFDM_m2), by = c("date" = "MetDate2", "channel" = "channel", "Response")) %>% 
  rowwise() %>% 
  # avg any times when there are AFDM from ts and sampling
  mutate(gAFDM_m2_F = mean(c(afdm_g_m2, gAFDM_m2_2015, gAFDM_m2_2016, gAFDM_m2_2017), na.rm = T),
         Y = as.numeric(strftime(date, format = "%Y"))) %>% 
  select(Y, date, channel, Response, gAFDM_m2_F) %>% 
  left_join(ExpDays, by = "date") %>% 
  left_join(BioDynExpCodes, by = c("Y", "channel")) %>% 
  mutate(Nf = as.factor(as.character(nitrogen_ugL)),
         Pf = as.factor(as.character(phosphorus_ugL)),
         date = as.Date(date)) 

# 2015 -----
bio15 <- BioDyn %>% 
  filter(Y == "2015") %>% 
  mutate(TempFval = ifelse(TempF == "A", "8°C", 
                           ifelse(TempF == "B", "12°C",
                                  ifelse(TempF == "C", "16°C",
                                         ifelse(TempF == "D", "20°C",
                                                ifelse(TempF == "E", "24°C","BLAH" )))))) %>% 
  mutate(TempFval = as.factor(TempFval),
    TempFval = fct_relevel(TempFval, c("8°C","12°C","16°C","20°C","24°C"))) %>% 
  mutate(Nf = fct_relevel(Nf, c("0", "25", "50", "100", "150", "200")),
         NuM = case_when(nitrogen_ugL == 0 ~ 0,
                         nitrogen_ugL == 25 ~ 1.8,
                         nitrogen_ugL == 50 ~ 3.6,
                         nitrogen_ugL == 100 ~ 7.1,
                         nitrogen_ugL == 150 ~ 10.7,
                         nitrogen_ugL == 200 ~ 14.3),
         NuM = as.factor(NuM),
         NuM = fct_relevel(NuM, "0", "1.8", "3.6", "7.1", "10.7", "14.3"),
         date_col = case_when(date == "2015-07-13" | date == "2015-07-14" | date == "2015-07-15" | date == "2015-07-16" | date == "2015-07-17" | date == "2015-07-18" ~ "MD1", 
                                date == "2015-07-27" | date == "2015-07-28" | date == "2015-07-29" | date == "2015-07-30" ~ "MD2",
                              TRUE ~ "TS"),
         date_col = as.factor(date_col)) %>% 
  droplevels()


Nlabels_2015 <- c("0 µM-N", "1.8 µM-N", "3.6 µM-N", "7.1 µM-N", "10.7 µM-N", "14.3 µM-N")
names(Nlabels_2015) <-  c("0", "1.8", "3.6", "7.1", "10.7", "14.3")

png("05_Figures4MS/19_Fig1_BMtimeSeries2015.png", units = "in", height = 13, width = 15, res = 300)
ggplot() +
  geom_line(data = bio15, aes(y = gAFDM_m2_F, x = ExpDays), color = "grey20") +
  geom_point(data = bio15, aes(y = gAFDM_m2_F, x = ExpDays, fill = date_col),
             shape = 21, size = 4) +
  facet_grid(NuM ~ TempFval, labeller = labeller(NuM = Nlabels_2015)) +
  xlab("Day of experiment") +
  ylab(expression(paste("Biofilm biomass (g AFDM ",m^-2,")"))) +
  scale_fill_manual(values = c("grey40", "violet", "violetred3"), 
                    name = "Measurement type") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(color = "black", fill = "NA", linewidth = 1),
        panel.spacing.x = unit(15, "pt"),
        panel.spacing.y = unit(15, "pt"),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        axis.text = element_text(size = 22),
        axis.line = element_line(color = "black", linewidth = 1),
        plot.background = element_rect(fill = "white", color =  "white"),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24, face = "bold"),
        legend.key.height = unit(1,"cm")) 
dev.off()

# 2016 ----
bio16 <- BioDyn %>% 
  filter(Y == "2016") %>% 
  mutate(TempFval = ifelse(TempF == "A", "10°C", 
                           ifelse(TempF == "B", "13°C",
                                  ifelse(TempF == "C", "16°C",
                                         ifelse(TempF == "D", "20°C",
                                                ifelse(TempF == "E", "24°C","BLAH" )))))) %>% 
  mutate(TempFval = as.factor(TempFval),
         TempFval = fct_relevel(TempFval, c("10°C","13°C","16°C","20°C","24°C"))) %>% 
  mutate(Pf = fct_relevel(Pf, c("0", "25", "50", "100", "150", "200")),
         PuM = case_when(phosphorus_ugL == 0 ~ 0,
                         phosphorus_ugL == 25 ~ 0.8,
                         phosphorus_ugL == 50 ~ 1.6,
                         phosphorus_ugL == 100 ~ 3.6,
                         phosphorus_ugL == 150 ~ 4.8,
                         phosphorus_ugL == 200 ~ 6.5),
         PuM = as.factor(PuM),
         PuM = fct_relevel(PuM, "0", "0.8", "1.6", "3.6", "4.8", "6.5"),
         date_col = case_when(date %in% c(D1_2016_met, D1_2016_up, D1_2016_nfix) ~ "MD1",
                              date %in% c(D2_2016_met, D2_2016_up, D2_2016_nfix) ~ "MD2",
                              TRUE ~ "TS"),
         date_col = as.factor(date_col)) %>% 
  droplevels()

Plabels_2016 <- c("0 µM-N", "0.8 µM-P", "1.6 µM-P", "3.6 µM-P", "4.8 µM-P", "6.5 µM-P")
names(Plabels_2016) <-  c("0", "0.8", "1.6", "3.6", "4.8", "6.5")

png("05_Figures4MS/19_FigS5_BMtimeSeries2016.png", units = "in", height = 12, width = 15, res = 300)
ggplot() +
  geom_line(data = bio16, aes(y = gAFDM_m2_F, x = ExpDays), color = "grey20") +
  geom_point(data = bio16, aes(y = gAFDM_m2_F, x = ExpDays, fill = date_col),
             shape = 21, size = 4) +
  facet_grid(PuM ~ TempFval,
             labeller = labeller(PuM = Plabels_2016)) +
  xlab("Day of experiment") +
  ylab(expression(paste("Biofilm biomass (g AFDM ",m^-2,")"))) +
  scale_fill_manual(values = c("grey40", "violet", "violetred3"), 
                    name = "Measurement type") +
  # scale_y_continuous(sec.axis = sec_axis(trans = ~.))+
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(color = "black", fill = "NA", size = 1),
        panel.spacing.x = unit(15, "pt"),
        panel.spacing.y = unit(15, "pt"),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        axis.text = element_text(size = 22),
        axis.line = element_line(color = "black", size = 1),
        plot.background = element_rect(fill = "white", color =  "white"),
        strip.background = element_blank(),
        # strip.background = element_rect(fill = "grey", color =  "white"),
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        # strip.placement = "outside",
        legend.position = "top",
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24, face = "bold"),
        legend.key.height = unit(1,"cm")) 
dev.off()


# 2017 ----
bio17 <- BioDyn %>% 
  filter(Y == "2017")%>% 
  mutate(TempFval = ifelse(TempF == "A", "8°C", 
                           ifelse(TempF == "B", "11°C",
                                  ifelse(TempF == "C", "13°C",
                                         ifelse(TempF == "D", "17°C",
                                                ifelse(TempF == "E", "19°C","BLAH" )))))) %>% 
  mutate(TempFval = as.factor(TempFval),
         TempFval = fct_relevel(TempFval, c("8°C","11°C","13°C","17°C","19°C"))) %>% 
  mutate(NandP = paste0(Nf," ", Pf)) %>% 
  mutate(NandP2 = ifelse(NandP == "0 0", "0 N, 0 P",
                         ifelse(NandP == "50 0", "50 N, 0 P",
                                ifelse(NandP == "50 111", "50 N, 111 P", "BLAH")))) %>% 
  mutate(NandP_uM = case_when(NandP2 == "0 N, 0 P" ~ "0 N, 0 P",
                              NandP2 == "50 N, 0 P" ~ "3.6 N, 0 P",
                              NandP2 == "50 N, 111 P" ~ "3.6 N, 3.6 P"),
         NandP_uM = as.factor(NandP_uM),
         NandP_uM = fct_relevel(NandP_uM, "0 N, 0 P","3.6 N, 0 P","3.6 N, 3.6 P"),
         date_col = case_when(date %in% c(D1_2017_met, D1_2017_up, D1_2017_nfix) ~ "MD1",
                              date %in% c(D2_2017_met, D2_2017_up, D2_2017_nfix) ~ "MD2",
                              TRUE ~ "TS"),
         date_col = as.factor(date_col),
         Rep = as.factor(ifelse(channel <= 15, "R1", "R2"))) %>%
  # remove 2 big outliers
  filter(!(date == "2017-07-12" & channel == "12" & Response == "Nup_gAFDM_m2")) %>% 
  filter(!(date == "2017-07-10" & channel == "12" )) %>% 
  droplevels() 

NPlabels_2017 <- c("0 µM-N, 0 µM-P","3.6 µM-N, 0 µM-P","3.6 µM-N, 3.6 µM-P")
names(NPlabels_2017) <-  c("0 N, 0 P","3.6 N, 0 P","3.6 N, 3.6 P")

png("05_Figures4MS/19_FigS6_BMtimeSeries2017.png", units = "in", height = 12, width = 15, res = 300)
ggplot() +
  geom_line(data = bio17, aes(y = gAFDM_m2_F, x = ExpDays, linetype = Rep), color = "grey20") +
  geom_point(data = bio17, aes(y = gAFDM_m2_F, x = ExpDays, fill = date_col, shape = Rep),
              size = 4) +
  facet_grid(NandP_uM ~ TempFval,
             labeller = labeller(NandP_uM = NPlabels_2017)) +
  xlab("Day of experiment") +
  ylab(expression(paste("Biofilm biomass (g AFDM ",m^-2,")"))) +
  scale_fill_manual(values = c("grey40", "violet", "violetred3"), 
                    name = "Measurement type") +
  scale_shape_manual(values = c(21,22)) +
  guides(fill = guide_legend(override.aes = list(shape=21, size = 6))) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(color = "black", fill = "NA", size = 1),
        panel.spacing.x = unit(15, "pt"),
        panel.spacing.y = unit(15, "pt"),
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        axis.text = element_text(size = 22),
        axis.line = element_line(color = "black", size = 1),
        plot.background = element_rect(fill = "white", color =  "white"),
        strip.background = element_blank(),
        # strip.background = element_rect(fill = "grey", color =  "white"),
        strip.text.x = element_text(size = 25, face = "bold"),
        strip.text.y = element_text(size = 25, face = "bold"),
        legend.position = "top",
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 24, face = "bold"),
        legend.key.height = unit(1,"cm")) 
dev.off()

# save/load
# save.image("02b_Script_SavedImages/19_Fig2_S2S3_BiomassDynamics_Rdat")
# load("02b_Script_SavedImages/19_Fig2_S2S3_BiomassDynamics_Rdat")

