library(tidyverse)

biovol <- read.csv("01_Data/PaulaCommComp/NfixerBiovol_2015.csv") %>% 
  mutate(Biovol_Unit_AFDM = Biovolume.unitCm2/(Nfix_gAFDM_m2*10000))

# NO RELATIONSHIP BETWEEN PAULAS DATA AND "OURS"
ggplot(biovol, aes(y = AFDMmgcm2*10000*(1/1000), x = Nfix_gAFDM_m2)) +
  geom_point() +
  ylab("Paula AFDM (g AFDM/m2)") +
  xlab("'Lyndsie' AFDM (g AFDM/m2)")

# PRETTY STRONG RELATIONSHIP
summary(lm(Biovolume.unitCm2 ~ Nfix_gAFDM_m2, data = biovol))
ggplot(biovol, aes(y = Biovolume.unitCm2, x = Nfix_gAFDM_m2)) +
  geom_point(shape = 21, fill = "lightgreen", size = 5) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()+
  xlab(expression(paste("AFDM (g ",m^-1,")"))) +
  ylab(expression(paste("Biovolume (grid intersections ", cm^-1,")")))+
  annotate("text", x = 15, y = 1e3, label = "R2 = 0.55, P < 0.001", hjust = 0, size = 5) +
  stat_smooth(method = "lm")+
  theme(axis.title = element_text(size = 18))

summary(lm(log10(Biovolume.unitCm2) ~ log10(Nfix_gAFDM_m2), data = biovol))

ggplot(biovol, aes(y = TempC_Avg, x = Nut.Conc, color = Biovol_Unit_AFDM)) +
  geom_point()

ggplot(biovol, aes(y = Biovol_Unit_AFDM, x = Nut.Conc, color = as.factor(Temp.Code))) +
  geom_point() +
  geom_line() +
  scale_y_log10()

ggplot(biovol, aes(y = Biovol_Unit_AFDM, x = TempC_Avg, color = as.factor(Nut.Conc))) +
  geom_point() +
  geom_line() +
  scale_y_log10()

ggplot(biovol, aes(y = Biovol_Unit_AFDM/Nfix_gAFDM_m2, x = TempC_Avg, color = as.factor(Nut.Conc))) +
  geom_point() 

ggplot(biovol, aes(y = Biovol_Unit_AFDM/Nfix_gAFDM_m2, x = Nut.Conc)) +
  geom_point() 
  
