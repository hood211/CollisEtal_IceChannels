# Ordination of channel biofilm composition
# JMH, 5 June 23, Jul 2024

library(tidyverse)
library(vegan)
library(BiodiversityR)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

# data ----
  # biofilm composition at finest resolution
  # D0 mislabeled as C10 - changed in csv
  # 2 C29's in 2016, one was C5 - changed in csv

  BF_f <-  read.csv("01_Data/ChannelAlgae_finestRes.csv") %>% 
    mutate(year = as.character(year),
           ChanNum = as.character(ChanNum))
  
  # biofilm composition at coarsest resolution
  BF_c <-  read.csv("01_Data/ChannelAlgae_CourseRes.csv")

  # taxon table
  BG_taxa <- read.csv("01_Data/ChannelAlgae_TaxaCode_Matrix.csv")
  
  # General Channel data/info
  SecondSamplingDates <- c("7/27/2015", "8/2/2016", "7/22/2017")
  
  ChanInfo <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv", row.names = 1) %>% 
    # all biofilm comm data is for second measurement day
    filter(MetDate %in% SecondSamplingDates) %>% 
    select(Year, chan, tempF, MeanPre2wksTemp:P_cat) %>% 
    mutate(Year = as.character(Year),
           chan = as.character(chan))

# prepare data ----
  BF_f2 <- ChanInfo %>% 
            full_join(BF_f, by = c("Year" = "year", "chan" = "ChanNum"))
  
  BF_f2_2015  <- BF_f2 %>% 
    filter(Year == "2015") %>% 
    # two bottles Paula said delete this one
    filter(chan != "3b")
  
  BF_f2_2016  <- BF_f2 %>% 
    filter(Year == "2016") 
  
  BF_f2_2017  <- BF_f2 %>% 
    filter(Year == "2017") 

# RDA  ----
  
  ## 2015 ----
  # https://rpubs.com/Roeland-KINDT/706490
  # https://chrischizinski.github.io/SNR_R_Group/2016-08-10-Data-Transformations
  # https://ordnews.colostate.narkive.com/lMWF502c/1593-log-sqrt-and-other-transformation-with-bray-curtis-dissimilarity
  # https://stats.stackexchange.com/questions/296884/should-functional-groups-abundances-be-transformed-before-using-rda
  BF_f2_2015com <- BF_f2_2015 %>% 
    select(EPIADN:CROCOC) %>% 
    # square root data, which are already relative abundance = Hellinger
    mutate(across(EPIADN:CROCOC, sqrt)) %>% 
    # remove columns with all NA's
    purrr::discard(~all(is.na(.))) %>% 
    # remove colums with all zeros
    select(where(~ any(. != 0)))
  
  BF_f2_2015env <- BF_f2_2015 %>% 
    select(Year:Iceland)%>% 
    # coding in a interaction between T and N so that I can plot
    mutate(TempNitInt = MeanPre2wksTemp * N_tr_uM)

  RDA_fine0 <- rda(BF_f2_2015com ~ MeanPre2wksTemp + N_tr_uM + MeanPre2wksTemp:N_tr_uM, data = BF_f2_2015env)  
  ordistep(RDA_fine0, direction = "both") # interaction is best
  RDA_fine <- rda(BF_f2_2015com ~ MeanPre2wksTemp + N_tr_uM + MeanPre2wksTemp:N_tr_uM, data = BF_f2_2015env)  
  
  # there are three rda axes here
  anova(RDA_fine, permutations = 9999)
  anova(RDA_fine, by = "terms", permutations = 9999)
  anova(RDA_fine, by = "axis", permutations = 9999)
  RsquareAdj(RDA_fine)
  # BiodiversityR::nested.npmanova(BF_f2_2015com ~ MeanPre2wksTemp + N_tr_uM, data = BF_f2_2015env, method="euclidean", permutations = 999)
  # https://www.worldagroforestry.org/output/tree-diversity-analysis
  
  
  # https://rpubs.com/Roeland-KINDT/694016
  
  
  ### RDA axes 1 and 2 ----
  RDA_fine_op <-  ordiplot(RDA_fine, choices = c(1,2))
  sites.long2 <- sites.long(RDA_fine_op, env.data=BF_f2_2015env)
  species.long2 <- species.long(RDA_fine_op)
  axis.long2 <- axis.long(RDA_fine, choices=c(1, 2))
  
  spec.envfit <- envfit(RDA_fine_op, env=BF_f2_2015com)
  spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
  species.long2 <- species.long(RDA_fine_op, spec.data=spec.data.envfit)
  # species that explain > 60% of variation
  species.long3 <- species.long2[species.long2$r >= 0.5, ]
  # species.long3 <- species.long2[species.long2$p < 0.05, ]
  vectors.envfit <- envfit(RDA_fine_op, env=BF_f2_2015env[,c("MeanPre2wksTemp", "N_tr_uM", "TempNitInt")])
  vectors.long3 <- vectorfit.long(vectors.envfit)
  
  species.long4 <- species.long3 %>% 
    mutate(labels = case_when(labels == "DIAMES" ~ "ODOMES",
                              labels == "HANNEA" ~ "HANARC",
                              labels == "MELOSI" ~ "MELVAR",
                              labels == 'MERIDO' ~ "MERCIR",
                              labels == "NITPAL" ~ "NITSP5",
                              labels == "NOSSPO" ~ "NOSSPO"))
  
  vectors.long4 <- vectors.long3 %>% 
                    mutate(vector = ifelse(vector == "MeanPre2wksTemp", "Temperature",
                                           ifelse(vector == "N_tr_uM", "µM-N",
                                              ifelse(vector == "TempNitInt", "Temp x µM-N", "BLAH"))))
  
  #### 2015 plot ----
  # https://community.rstudio.com/t/feature-discussion-layers-aware-of-each-other-in-ggplot2/108156
  geom_text_repel2 <- function(...) {
    layer <- ggrepel::geom_text_repel(...)
    layer$ggrepel <- TRUE
    class(layer) <- c("ggrepel", class(layer))
    return(layer)
  }
  
  ggplot_add.ggrepel <- function(object, plot, object_name) {
    if (any(do.call(c, lapply(plot$layer, function(x) x$ggrepel)))) {
      warning(
        "There is more than one ggrepel layers. ",
        "This may cause overlap of labels"
      )
    }
    # Optionally, one may modify `object` here.
    NextMethod("ggplot_add")
  }
  
  p_rda12a <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey50", linetype = 1, linewidth = 1) +
    geom_hline(yintercept = c(0), color = "grey50", linetype = 1, linewidth = 1) +  
    xlab("RDA1 (22%; P < 0.001)") +
    ylab("RDA2 (10%; P < 0.001)") +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) + 
    geom_point(data=sites.long2 %>% 
                 mutate(Temperature = as.factor(round(MeanPre2wksTemp,0))),
               aes(x=axis1, y=axis2, fill=Temperature, size=N_uM),
               shape = 21) +
    geom_segment(data=vectors.long4,
                 aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
                 colour="black", linewidth=1, arrow=arrow(length = unit(2,"mm")), alpha = 1) +
    geom_text_repel2(data=vectors.long4, 
                    aes(x=axis1*1.25, y=axis2*1.25-0.25, label=vector),
                    colour="black", fontface = "bold", size = 5,
                    point.size = NA,
                    box.padding = 0.5) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    scale_size_binned(breaks = c(0.11, 0.19, 3.68, 7.25, 10.82, 14.4)) +
    guides(size = guide_legend("µM-N"),
           fill = guide_legend(override.aes=list(shape=21, size = 6), "Temperature (°C)")) +
    annotate(geom = "text", x = -2.5, y = 4,
             label = "a) N-only", size = 7, fontface = "bold", hjust = 0) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position = "none")
  
  p_rda12b <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey50", linetype = 1, linewidth = 1) +
    geom_hline(yintercept = c(0), color = "grey50", linetype = 1, linewidth = 1) +  
    xlab("RDA1 (22%; P < 0.001)") +
    ylab("RDA2 (10%; P < 0.001)") +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) + 
    geom_point(data=sites.long2 %>% 
                 mutate(Temperature = as.factor(round(MeanPre2wksTemp,0))),
               aes(x=axis1, y=axis2, fill=Temperature, size=N_uM),
               shape = 21) +
    geom_segment(data=species.long4,
                 aes(x=0, y=0, xend=axis1*1, yend=axis2*1),
                 colour="red", linewidth=1, arrow=arrow(length = unit(2,"mm")), alpha = 0.7) +
    geom_text_repel2(data=species.long4,
                     aes(x=axis1*1.25, y=axis2*1.25, label=labels),
                     colour="red", fontface = "bold", size = 5,
                     point.size = NA,
                     box.padding = 0.5) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    scale_size_binned(breaks = c(0.11, 0.19, 3.68, 7.25, 10.82, 14.4)) +
    guides(size = guide_legend("µM-N"),
           fill = guide_legend(override.aes=list(shape=21, size = 6), "Temperature (°C)")) +
    annotate(geom = "text", x = -2.5, y = 4,
             label = "b) N-only", size = 7, fontface = "bold", hjust = 0) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 12, face = "bold"))
  
   ### RDA axes 1 and 3 ----
  RDA_fine_op_rda13 <-  ordiplot(RDA_fine, choices = c(1,3))
  sites.long2_rda13 <- sites.long(RDA_fine_op_rda13, env.data=BF_f2_2015env)
  species.long2_rda13 <- species.long(RDA_fine_op_rda13)
  axis.long2_rda13 <- axis.long(RDA_fine, choices=c(1,3))
  
  spec.envfit_rda13 <- envfit(RDA_fine_op_rda13, env=BF_f2_2015com)
  spec.data.envfit_rda13 <- data.frame(r=spec.envfit_rda13$vectors$r, p=spec.envfit_rda13$vectors$pvals)
  species.long2_rda13 <- species.long(RDA_fine_op_rda13, spec.data=spec.data.envfit_rda13)
  # species with p values < 0.05
  species.long3_rda13  <- species.long2_rda13[species.long2_rda13$p < 0.05, ]
  vectors.envfit_rda13 <- envfit(RDA_fine_op_rda13, env=BF_f2_2015env[,c("MeanPre2wksTemp", "N_tr_uM")])
  vectors.long3_rda13 <- vectorfit.long(vectors.envfit_rda13)
  
  p_rda13 <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axis.long2_rda13[1, "label"]) +
    ylab(axis.long2_rda13[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_point(data=sites.long2_rda13,
               aes(x=axis1, y=axis2, colour=tempF, shape=as.factor(N_uM)),
               size=5) +
    geom_point(data=species.long2_rda13,
               aes(x=axis1, y=axis2)) +
    geom_segment(data=species.long3_rda13,
                 aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
                 colour="red", size=0.7, arrow=arrow()) +
    geom_text_repel(data=species.long3_rda13,
                    aes(x=axis1*2, y=axis2*2, label=labels),
                    colour="red") +
    geom_segment(data=vectors.long3_rda13,
                 aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
                 colour="blue", size=0.7, arrow=arrow()) +
    geom_text_repel(data=vectors.long3_rda13, 
                    aes(x=axis1*1, y=axis2*1, label=vector),
                    colour="black") +
    # BioR.theme +
    coord_fixed(ratio=1) +
    scale_color_brewer(palette = "RdBu", direction = -1)
  
  
  ### RDA-PCA axes 1 and 2 ----
  RDA_fine_op_pca <-  ordiplot(RDA_fine, choices = c(4,5))
  sites.long2_pca <- sites.long(RDA_fine_op_pca, env.data=BF_f2_2015env)
  species.long2_pca <- species.long(RDA_fine_op_pca)
  axis.long2_pca <- axis.long(RDA_fine, choices=c(4,5))
  
  spec.envfit_pca <- envfit(RDA_fine_op_pca, env=BF_f2_2015com)
  spec.data.envfit_pca <- data.frame(r=spec.envfit_pca$vectors$r, p=spec.envfit_pca$vectors$pvals)
  species.long2_pca <- species.long(RDA_fine_op_pca, spec.data=spec.data.envfit_pca)
  # species with p values < 0.05
  species.long3_pca  <- species.long2_pca[species.long2_pca$p < 0.05, ]
  vectors.envfit_pca <- envfit(RDA_fine_op_pca, env=BF_f2_2015env[,c("MeanPre2wksTemp", "N_tr_uM")])
  vectors.long3_pca <- vectorfit.long(vectors.envfit_pca)
  
  p_rda_pca12 <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axis.long2_pca[1, "label"]) +
    ylab(axis.long2_pca[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_point(data=sites.long2_pca,
               aes(x=axis1, y=axis2, colour=tempF, shape=as.factor(N_uM)),
               size=5) +
    geom_point(data=species.long2_pca,
               aes(x=axis1, y=axis2)) +
    geom_segment(data=species.long3_pca,
                 aes(x=0, y=0, xend=axis1*2, yend=axis2*2),
                 colour="red", size=0.7, arrow=arrow()) +
    geom_text_repel(data=species.long3_pca,
                    aes(x=axis1*2, y=axis2*2, label=labels),
                    colour="red") +
    geom_segment(data=vectors.long3_pca,
                 aes(x=0, y=0, xend=axis1*1, yend=axis2*1), 
                 colour="blue", size=0.7, arrow=arrow()) +
    geom_text_repel(data=vectors.long3_pca, 
                    aes(x=axis1*1, y=axis2*1, label=vector),
                    colour="black") +
    # BioR.theme +
    coord_fixed(ratio=1) +
    scale_color_brewer(palette = "RdBu", direction = -1)

  grid.arrange(p_rda12, p_rda13, p_rda_pca12, nrow = 3, ncol = 1)  
  
  ## 2016 ----
  BF_f2_2016com <- BF_f2_2016 %>% 
    select(EPIADN:CROCOC) %>% 
    # square root data, which are already relative abundance = Hellinger
    mutate(across(EPIADN:CROCOC, sqrt)) %>% 
    # remove columns with all NA's
    purrr::discard(~all(is.na(.))) %>% 
    # remove colums with all zeros
    select(where(~ any(. != 0)))
  
  BF_f2_2016env <- BF_f2_2016 %>% 
    select(Year:Iceland)
  
  RDA_fine_2016 <- rda(BF_f2_2016com ~ MeanPre2wksTemp + P_tr_uM + MeanPre2wksTemp:P_tr_uM, data = BF_f2_2016env)  
  ordistep(RDA_fine_2016, direction = "both")
  RDA_fine_2016f <- rda(BF_f2_2016com ~ MeanPre2wksTemp, data = BF_f2_2016env)  
  # there are three rda axes here
  anova(RDA_fine_2016f, permutations = 9999)
  anova(RDA_fine_2016f, by = "terms", permutations = 9999)
  anova(RDA_fine_2016f, by = "axis", permutations = 9999)
  RsquareAdj(RDA_fine_2016f)
  # BiodiversityR::nested.npmanova(BF_f2_2015com ~ MeanPre2wksTemp + N_tr_uM, data = BF_f2_2015env, method="euclidean", permutations = 999)
  # https://www.worldagroforestry.org/output/tree-diversity-analysis
  
  
  # https://rpubs.com/Roeland-KINDT/694016
  
  
  ## RDA axes 1 and 2 ----
  RDA_fine_op_2016 <-  ordiplot(RDA_fine_2016f, choices = c(1,2))
  sites.long2_2016 <- sites.long(RDA_fine_op_2016, env.data=BF_f2_2016env)
  species.long2_2016 <- species.long(RDA_fine_op_2016)
  axis.long2_2016 <- axis.long(RDA_fine_2016f, choices=c(1, 2))
  
  spec.envfit_2016 <- envfit(RDA_fine_op_2016, env=BF_f2_2016com)
  spec.data.envfit_2016 <- data.frame(r=spec.envfit_2016$vectors$r, p=spec.envfit_2016$vectors$pvals)
  species.long2_2016 <- species.long(RDA_fine_op_2016, spec.data = spec.data.envfit_2016)
  # species that explain > 50% of variation
  species.long3_2016 <- species.long2_2016[species.long2_2016$r >= 0.5, ]
  vectors.envfit_2016 <- envfit(RDA_fine_op_2016, env=BF_f2_2016env[,c("MeanPre2wksTemp")])
  vectors.long3_2016 <- vectorfit.long(vectors.envfit_2016)
  
  species.long4_2016 <- species.long3_2016 %>% 
    mutate(labels = case_when(labels == "RHOGIB" ~ "EPIGIB",
                              labels == "NOSSPO" ~ "NOSSPO",
                              labels == "ANSTR" ~ "ANASP1"))
  
  vectors.long4_2016 <- vectors.long3_2016 %>% 
    mutate(vector = ifelse(vector == "1", "Temperature", "BLAH"))
  
  #### 2016 plot ----
  p_rda12_2016a <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey50", linetype = 1, size = 1) +
    geom_hline(yintercept = c(0), color = "grey50", linetype = 1, size = 1) +  
    xlab("RDA1 (9%, P = 0.015)") +
    ylab(axis.long2_2016[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-2,3)) + 
    geom_point(data=sites.long2_2016 %>% 
                 mutate(Temperature = as.factor(round(MeanPre2wksTemp,0))),
               aes(x=axis1, y=axis2, fill=Temperature),
               shape = 21, size=5) +
    geom_segment(data=vectors.long4_2016,
                 aes(x=0, y=0, xend=axis1*1.2, yend=axis2*1), 
                 colour="black", size=1, arrow=arrow(length = unit(2,"mm")), alpha = 1) +
    geom_text_repel2(data=vectors.long4_2016, 
                    aes(x=axis1*1, y=axis2*1.2+0.25, label=vector),
                    colour="black", fontface = "bold", size = 5,
                    point.size = NA,
                    box.padding = 0.5) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    scale_size_binned(breaks = c(0.11, 0.19, 3.68, 7.25, 10.82, 14.4)) +
    annotate(geom = "text", x = -3, y = 3,
             label = "c) P-only", size = 7, fontface = "bold", hjust = 0) +
    guides(size = guide_legend("µM-N"),
           fill = guide_legend("Temperature (°C)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position = "none")
  
  p_rda12_2016b <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey50", linetype = 1, size = 1) +
    geom_hline(yintercept = c(0), color = "grey50", linetype = 1, size = 1) +  
    xlab("RDA1 (9%, P = 0.015)") +
    ylab(axis.long2_2016[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL), limits = c(-2,3)) + 
    geom_point(data=sites.long2_2016 %>% 
                 mutate(Temperature = as.factor(round(MeanPre2wksTemp,0))),
               aes(x=axis1, y=axis2, fill=Temperature),
               shape = 21, size=5) +
    geom_segment(data=species.long4_2016,
                 aes(x=0, y=0, xend=axis1*1, yend=axis2*1),
                 colour="red", size=1, arrow=arrow(length = unit(2,"mm")), alpha = 0.7) +
    geom_text_repel2(data=species.long4_2016,
                     aes(x=axis1*1.25, y=axis2*1.25, label=labels),
                     colour="red", fontface = "bold", size = 5,
                     point.size = NA,
                     box.padding = 0.5) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    scale_size_binned(breaks = c(0.11, 0.19, 3.68, 7.25, 10.82, 14.4)) +
    annotate(geom = "text", x = -3, y = 3,
             label = "d) P-only", size = 7, fontface = "bold", hjust = 0) +
    guides(size = guide_legend("µM-N"),
           fill = guide_legend("Temperature (°C)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 12, face = "bold"))

  ## 2017 ----
  BF_f2_2017com <- BF_f2_2017 %>% 
    select(EPIADN:CROCOC) %>% 
    # square root data, which are already relative abundance = Hellinger
    mutate(across(EPIADN:CROCOC, sqrt)) %>% 
    # remove columns with all NA's
    purrr::discard(~all(is.na(.))) %>% 
    # remove colums with all zeros
    select(where(~ any(. != 0)))
  
  BF_f2_2017env <- BF_f2_2017 %>% 
    select(Year:Iceland) %>% 
    # coding in a interaction between T and N so that I can plot
    mutate(TempNitInt = MeanPre2wksTemp * N_tr_uM)
  
  RDA_fine_2017 <- rda(BF_f2_2017com ~ MeanPre2wksTemp * NPratio + MeanPre2wksTemp * N_tr_uM , data = BF_f2_2017env)  
  ordistep(RDA_fine_2017, direction = "both") # best interaction bw Temp and N
  
  RDA_fine_2017f <- rda(BF_f2_2017com ~ MeanPre2wksTemp + N_tr_uM + MeanPre2wksTemp:N_tr_uM, data = BF_f2_2017env)  
  # there are three rda axes here
  anova(RDA_fine_2017f, permutations = 9999)
  anova(RDA_fine_2017f, by = "terms", permutations = 9999)
  anova(RDA_fine_2017f, by = "axis", permutations = 9999)
  RsquareAdj(RDA_fine_2017f)
  # BiodiversityR::nested.npmanova(BF_f2_2015com ~ MeanPre2wksTemp + N_tr_uM, data = BF_f2_2015env, method="euclidean", permutations = 999)
  # https://www.worldagroforestry.org/output/tree-diversity-analysis
  
  
  # https://rpubs.com/Roeland-KINDT/694016
  
  
  ## RDA axes 1 and 2 ----
  RDA_fine_op_2017 <-  ordiplot(RDA_fine_2017f, choices = c(1,2))
  sites.long2_2017 <- sites.long(RDA_fine_op_2017, env.data=BF_f2_2017env)
  species.long2_2017 <- species.long(RDA_fine_op_2017)
  axis.long2_2017 <- axis.long(RDA_fine_2017f, choices=c(1, 2))
  
  spec.envfit_2017 <- envfit(RDA_fine_op_2017, env=BF_f2_2017com)
  spec.data.envfit_2017 <- data.frame(r=spec.envfit_2017$vectors$r, p=spec.envfit_2017$vectors$pvals)
  species.long2_2017 <- species.long(RDA_fine_op_2017, spec.data = spec.data.envfit_2017)
  # species that explain > 50% of variation
  species.long3_2017 <- species.long2_2017[species.long2_2017$r >= 0.5, ]
  vectors.envfit_2017 <- envfit(RDA_fine_op_2017, env=BF_f2_2017env[,c("MeanPre2wksTemp", "N_tr_uM", "TempNitInt")])
  vectors.long3_2017 <- vectorfit.long(vectors.envfit_2017)
  
  species.long4_2017 <- species.long3_2017 %>% 
    mutate(labels = case_when(labels == "DIAMES" ~ "ODOMES",
                              labels == "FRADMD" ~ "FRASP3",
                              labels == "FRALNG" ~ "FRASP4",
                              labels == "FRAPCH" ~ "FRASP5",
                              labels == "MELOSI" ~ "MELVAR",
                              labels == "NITPAL" ~ "NITSP5",
                              labels == "NIT086" ~ "NITSP6",
                              labels == "ULNULN" ~ "ULNULN",
                              labels == "NOSSPO" ~ "NOSSPO"))
  
  vectors.long4_2017 <- vectors.long3_2017 %>% 
    mutate(vector = ifelse(vector == "MeanPre2wksTemp", "Temperature",
                           ifelse(vector == "N_tr_uM", "µM-N", 
                                  ifelse(vector == "TempNitInt", "Temp x µM-N", "BLAH"))))
  
  #### 2017 plot ----
  p_rda12_2017a <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey50", linetype = 1, size = 1) +
    geom_hline(yintercept = c(0), color = "grey50", linetype = 1, size = 1) +  
    xlab("RDA1 (29%, P < 0.001)") +
    ylab("RDA2 (7%, P = 0.012)") +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) + 
    geom_point(data=sites.long2_2017 %>% 
                 mutate(Temperature = as.factor(round(MeanPre2wksTemp,0))),
               aes(x=axis1, y=axis2, fill=Temperature, size = as.factor(N_uM)),
               shape = 21) +
    geom_segment(data=vectors.long4_2017,
                 aes(x=0, y=0, xend=axis1*3, yend=axis2*3), 
                 colour="black", size=1, arrow=arrow(length = unit(2,"mm")), alpha = 1) +
    geom_text_repel2(data=vectors.long4_2017, 
                    aes(x=axis1*3.2, y=axis2*3.2 - 0.25, label=vector),
                    colour="black", fontface = "bold", size = 5,
                    point.size = NA,
                    box.padding = 0.5) +
    scale_size_manual(values = c(2,5), guide = guide_legend("µM-N")) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    annotate(geom = "text", x = -2.7, y = 3.8,
             label = "a)", size = 7, fontface = "bold", hjust = 0) +
    guides(fill=guide_legend(override.aes=list(shape=21, size = 6), "Temperature (°C)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position = "none")
  
  p_rda12_2017b <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey50", linetype = 1, size = 1) +
    geom_hline(yintercept = c(0), color = "grey50", linetype = 1, size = 1) +  
    xlab("RDA1 (29%, P < 0.001)") +
    ylab("RDA2 (7%, P = 0.012)") +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) + 
    geom_point(data=sites.long2_2017 %>% 
                 mutate(Temperature = as.factor(round(MeanPre2wksTemp,0))),
               aes(x=axis1, y=axis2, fill=Temperature, size = as.factor(N_uM)),
               shape = 21) +
    geom_segment(data=species.long4_2017,
                 aes(x=0, y=0, xend=axis1*1, yend=axis2*1),
                 colour="red", size=1, arrow=arrow(length = unit(2,"mm")), alpha = 0.7) +
    geom_text_repel2(data=species.long4_2017,
                     aes(x=axis1*1.1, y=axis2*1.1, label=labels),
                     colour="red", fontface = "bold", size = 5,
                     point.size = NA,
                     box.padding = 0.5) +
    scale_size_manual(values = c(2,5), guide = guide_legend("µM-N")) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    annotate(geom = "text", x = -2.7, y = 3.8,
             label = "b)", size = 7, fontface = "bold", hjust = 0) +
    guides(fill=guide_legend(override.aes=list(shape=21, size = 6), "Temperature (°C)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 12, face = "bold"))

  png("05_Figures4MS/18_Fig3_Channels_Nexperiment_BiofilmComp_2015_16.png", 
      units = "in", height = 10, width = 11, res = 300)
  plot_grid(p_rda12a, p_rda12b, 
            p_rda12_2016a, p_rda12_2016b, 
            ncol = 2,
            nrow = 2,
            rel_widths = c(0.75,1))
dev.off()  

png("05_Figures4MS/18_FigS12_Channels_Nexperiment_BiofilmComp_2017.png", 
    units = "in", height = 5, width = 11, res = 300)
plot_grid(p_rda12_2017a, p_rda12_2017b,
          ncol = 2,
          nrow = 1,
          rel_widths = c(0.75,1))
dev.off()  


# Rel biomass plots ----
# New taxa table
# taxon table
BG_taxa <- read.csv("01_Data/ChannelAlgae_TaxaCode_Matrix_JMH.csv") %>% 
  select(MajorGroups3, FixerNonFixer, MinorGroupJMH, AlgalCodeName)

# clean up data
BF_f3 <- BF_f2 %>% 
  # two bottles Paula said delete this one
  filter(!(chan == "3b" & Year == "2015")) %>% 
  pivot_longer(cols = EPIADN:CROCOC, names_to = "Taxa", values_to = "RelBM") %>% 
  left_join(BG_taxa, by = c("Taxa" = "AlgalCodeName")) %>% 
  group_by(Year, chan, tempF, N_uM, P_uM, NPratio, MajorGroups3, FixerNonFixer, MinorGroupJMH) %>% 
  summarize(RelBM = sum(RelBM, na.rm = T)) 

# identify rare taxa - 
# rare = taxa <10% rel BM in a yehar
Group2Other <- BF_f3 %>% 
  group_by(Year, MajorGroups3, MinorGroupJMH) %>% 
  summarize(max_RelBM = max(RelBM)) %>% 
  filter(max_RelBM <= 10) %>% 
  arrange(Year, MajorGroups3, MinorGroupJMH)

# List of rare taxa, plus unID, grouped by fixer/nonfixer
Other_FX_2015 <- Group2Other[Group2Other$MajorGroups3 %in% c("CYFX", "DIAFX") &
                               Group2Other$Year == "2015",]$MinorGroupJMH
Other_FX_2016 <- Group2Other[Group2Other$MajorGroups3 %in% c("CYFX", "DIAFX") &
                               Group2Other$Year == "2016",]$MinorGroupJMH
Other_FX_2017 <- Group2Other[Group2Other$MajorGroups3 %in% c("CYFX", "DIAFX") &
                               Group2Other$Year == "2017",]$MinorGroupJMH
Other_NFX_2015 <- Group2Other[Group2Other$MajorGroups3 %in% c("CYNFX", "DIANFX",
                                                              "CRPT", "DINFLG", "GREEN", "UNID", "YLGR") &
                                Group2Other$Year == "2015",]$MinorGroupJMH

Other_NFX_2016 <- Group2Other[Group2Other$MajorGroups3 %in% c("CYNFX", "DIANFX",
                                                              "CRPT", "DINFLG", "GREEN", "UNID", "YLGR") &
                                Group2Other$Year == "2016",]$MinorGroupJMH

Other_NFX_2017 <- Group2Other[Group2Other$MajorGroups3 %in% c("CYNFX", "DIANFX",
                                                              "CRPT", "DINFLG", "GREEN", "UNID", "YLGR") &
                                Group2Other$Year == "2017",]$MinorGroupJMH

# add other tazxa
Other_NFX_2015 <- c(Other_NFX_2015, "Green_unID")

BF_f4 <- BF_f3 %>% 
  mutate(MinorGroupJMH2 = case_when(MinorGroupJMH %in% Other_FX_2015 & Year == "2015" ~ "Other_F",
                                    MinorGroupJMH %in% Other_FX_2016 & Year == "2016" ~ "Other_F",
                                    MinorGroupJMH %in% Other_FX_2017 & Year == "2017" ~ "Other_F",
                                    MinorGroupJMH %in% Other_NFX_2015 & Year == "2015" ~ "Other_NF",
                                    MinorGroupJMH %in% Other_NFX_2016 & Year == "2016" ~ "Other_NF",
                                    MinorGroupJMH %in% Other_NFX_2017 & Year == "2017" ~ "Other_NF",
                                    .default = as.character(MinorGroupJMH)),
         MinorGroupJMH2 = case_when(MinorGroupJMH2 == "Hannea" ~ "Hannaea",
                                    MinorGroupJMH2 == "Diatoma" ~ "Odontidium",
                                    MinorGroupJMH2 == "Synedra" ~ "Ulnaria",
                                    MinorGroupJMH2 == "Nitschia" ~ "Nitzschia",
                                    .default = as.character(MinorGroupJMH2))) 


# 2015 ----
BF_f4_2015 <- BF_f4 %>% 
  filter(Year == "2015") %>% 
  mutate(MinorGroupJMH2 = as.factor(MinorGroupJMH2),
         MinorGroupJMH2 = fct_relevel(MinorGroupJMH2, 
                                      TaxaOrder2015)) %>% 
  mutate(TempFval = ifelse(tempF == "A", "8°C", 
                           ifelse(tempF == "B", "12°C",
                                  ifelse(tempF == "C", "16°C",
                                         ifelse(tempF == "D", "20°C",
                                                ifelse(tempF == "E", "24°C","BLAH" )))))) %>% 
  mutate(TempFval = as.factor(TempFval),
         TempFval = fct_relevel(TempFval, c("8°C","12°C","16°C","20°C","24°C"))) %>% 
  mutate(N_uM = case_when(N_uM == 0.11 ~ 0,
                          N_uM == 1.90 ~ 1.8,
                          N_uM == 3.68 ~ 3.6,
                          N_uM == 7.25 ~ 7.1,
                          N_uM == 10.82 ~ 10.7,
                          N_uM == 14.40 ~ 14.3),
         N_uM = as.factor(N_uM),
         N_uM = fct_relevel(N_uM, "0", "1.8", "3.6", "7.1", "10.7", "14.3")) 

# labels  
Nlabels_2015 <- c("0 µM-N", "1.8 µM-N", "3.6 µM-N", "7.1 µM-N", "10.7 µM-N", "14.3 µM-N")
names(Nlabels_2015) <-  c("0", "1.8", "3.6", "7.1", "10.7", "14.3")

# organize sp
as.vector(BF_f4 %>% 
            filter(Year == "2015") %>% 
            group_by(MinorGroupJMH2) %>% 
            summarize(RelBM = mean(RelBM)) %>% 
            arrange(desc(RelBM)) %>% 
            select(MinorGroupJMH2))


TaxaOrder2015 <- c("Melosira", "Fragilaria", "Odontidium", "Ulnaria", "Meridion", "Hannaea", "Nitzschia", "Other_NF", "Nostoc", "Anabaena", "Epithemia", "Other_F")

ColorsNfix_2015 <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
ColorsNonNfix_2015 <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec")
Colors2015 <- c(ColorsNonNfix_2015, ColorsNfix_2015)

# plot
png("05_Figures4MS/18_AutotrophRelBMPlots2015_FigS7.jpg", units = "in", height = 13, width = 15, res = 300)
ggplot(BF_f4_2015, 
       aes(x = TempFval, y = RelBM, fill = MinorGroupJMH2)) +
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.75) +
  scale_fill_manual(values = Colors2015,
                    labels = c(expression(italic("Melosira")), 
                               expression(italic("Fragilaria")), 
                               expression(italic("Odontidium")), 
                               expression(italic("Ulnaria")), 
                               expression(italic("Meridion")), 
                               expression(italic("Hannea")), 
                               expression(italic("Nitzschia")), 
                               "Other_NF", 
                               expression(italic("Nostoc")), 
                               expression(italic("Anabaena")), 
                               expression(italic("Epithemia")), 
                               "Other_F")) +
  facet_grid(N_uM ~., labeller = labeller(N_uM = Nlabels_2015)) +
  ylab("Relative biomass (%)") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "Taxa group")) +
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

# 2016 ----
TaxaOrder2016 <- c("Nostoc", "Epithemia", "Anabaena", "Tolypothrix", "Calothrix", "Cylindrospermum", "Other_F", "Mougeotia", "Meridion", "Ulnaria", "Other_NF")

BF_f4_2016 <- BF_f4 %>% 
  filter(Year == "2016") %>% 
  mutate(MinorGroupJMH2 = as.factor(MinorGroupJMH2),
         MinorGroupJMH2 = fct_relevel(MinorGroupJMH2, 
                                      TaxaOrder2016)) %>% 
  mutate(TempFval = ifelse(tempF == "A", "10°C", 
                           ifelse(tempF == "B", "13°C",
                                  ifelse(tempF == "C", "16°C",
                                         ifelse(tempF == "D", "20°C",
                                                ifelse(tempF == "E", "24°C","BLAH" )))))) %>% 
  mutate(TempFval = as.factor(TempFval),
         TempFval = fct_relevel(TempFval, c("10°C","13°C","16°C","20°C","24°C"))) %>% 
  mutate(P_uM = case_when(P_uM == 0.36 ~ 0,
                          P_uM == 1.17 ~ 0.8,
                          P_uM == 1.97 ~ 1.6,
                          P_uM == 3.59 ~ 3.6,
                          P_uM == 5.20 ~ 4.8,
                          P_uM == 6.81 ~ 6.5),
         P_uM = as.factor(P_uM),
         P_uM = fct_relevel(P_uM, "0", "0.8", "1.6", "3.6", "4.8", "6.5")) 

# organize sp
as.vector(BF_f4 %>% 
            filter(Year == "2016") %>% 
            group_by(MinorGroupJMH2) %>% 
            summarize(RelBM = mean(RelBM)) %>% 
            arrange(desc(RelBM)) %>% 
            select(MinorGroupJMH2))





ColorsNfix_2016 <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
ColorsNonNfix_2016 <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4")
Colors2016 <- c(ColorsNfix_2016, ColorsNonNfix_2016)

# labels  
Plabels_2016 <- c("0 µM-N", "0.8 µM-P", "1.6 µM-P", "3.6 µM-P", "4.8 µM-P", "6.5 µM-P")
names(Plabels_2016) <-  c("0", "0.8", "1.6", "3.6", "4.8", "6.5")

# plot
png("05_Figures4MS/18_AutotrophRelBMPlots2016_FigS9.jpg", units = "in", height = 13, width = 15, res = 300)
ggplot(BF_f4_2016, 
       aes(x = TempFval, y = RelBM, fill = MinorGroupJMH2)) +
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.75) +
  scale_color_manual(values = c("black", "grey90")) +
  scale_fill_manual(values = Colors2016,
                    labels = c(expression(italic("Nostoc")), 
                               expression(italic("Epithemia")), 
                               expression(italic("Anabaena")), 
                               expression(italic("Tolypothrix")), 
                               expression(italic("Calothrix")), 
                               expression(italic("Cylindrospermum")), 
                               "Other_F", 
                               expression(italic("Mougeotia")), 
                               expression(italic("Meridion")), 
                               expression(italic("Ulnaria")), 
                               "Other_NF")) +
  facet_grid(P_uM ~., labeller = labeller(P_uM = Plabels_2016)) +
  ylab("Relative biomass (%)") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "Taxa group")) +
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


# 2017 ----
TaxaOrder2017 <- c("Nostoc", "Anabaena", "Epithemia", "CyanoNonF_unID", "Other_F", #5
                   "Melosira", "Odontidium", "Ulnaria", "Fragilaria", "Nitzschia", "Meridion", "DiatomNonF_unID",
                   "Green_unID", "CyanoNonF_unID", "Achnanthes", "Hannaea", "Other_NF") #11

BF_f4_2017 <- BF_f4 %>% 
  filter(Year == "2017") %>% 
  mutate(MinorGroupJMH2 = as.factor(MinorGroupJMH2),
         MinorGroupJMH2 = fct_relevel(MinorGroupJMH2, 
                                      TaxaOrder2017)) %>% 
  mutate(TempFval = ifelse(tempF == "A", "10°C", 
                           ifelse(tempF == "B", "13°C",
                                  ifelse(tempF == "C", "16°C",
                                         ifelse(tempF == "D", "20°C",
                                                ifelse(tempF == "E", "24°C","BLAH" )))))) %>% 
  mutate(TempFval = as.factor(TempFval),
         TempFval = fct_relevel(TempFval, c("10°C","13°C","16°C","20°C","24°C"))) %>% 
  mutate(trtLabel = case_when(NPratio == "0.31" ~ "0 N, 0 P",
                              NPratio == "0.93" ~ "3.6 N, 3.6 P",
                              NPratio == "10.22" ~ "3.6 N, 0 P"),
         trtLabel = as.factor(trtLabel),
         trtLabel = fct_relevel(trtLabel, "0 N, 0 P","3.6 N, 0 P","3.6 N, 3.6 P"))


# organize sp
as.vector(BF_f4 %>% 
            filter(Year == "2017") %>% 
            group_by(MinorGroupJMH2) %>% 
            summarize(RelBM = mean(RelBM)) %>% 
            arrange(desc(RelBM)) %>% 
            select(MinorGroupJMH2))







ColorsNfix_2017 <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
ColorsNonNfix_2017 <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec",
                        "#d9d9d9",
                        "#bc80bd",
                        "#ccebc5")
Colors2017 <- c(ColorsNfix_2017, ColorsNonNfix_2017)



# labels  
NPlabels_2017 <- c("0 µM-N, 0 µM-P","3.6 µM-N, 0 µM-P","3.6 µM-N, 3.6 µM-P")
names(NPlabels_2017) <-  c("0 N, 0 P","3.6 N, 0 P","3.6 N, 3.6 P")

# plot
png("05_Figures4MS/18_AutotrophRelBMPlots2017_FigS14.jpg", units = "in", height = 13, width = 15, res = 300)
ggplot(BF_f4_2017, 
       aes(x = TempFval, y = RelBM, fill = MinorGroupJMH2)) +
  geom_bar(position = "stack", stat = "identity", color = "black", linewidth = 0.75) +
  scale_color_manual(values = c("black", "grey90")) +
  scale_fill_manual(values = Colors2017,
                    labels = c(c(expression(italic("Nostoc")), 
                                 expression(italic("Anabaena")), 
                                 expression(italic("Epithemia")), 
                                 "Other_F", #5
                                 expression(italic("Melosira")), 
                                 expression(italic("Odontidium")), 
                                 expression(italic("Ulnaria")), 
                                 expression(italic("Fragilaria")), 
                                 expression(italic("Nitzschia")), 
                                 expression(italic("Meridion")), 
                                 "Diatom-unidentified",
                                 "Green-unidentified",
                                 "Cyano-unidentified",
                                 expression(italic("Achnanthes")), 
                                 expression(italic("Hannaea")), 
                                 "Other_NF"))) +
  facet_grid(trtLabel ~., labeller = labeller(trtLabel = NPlabels_2017)) +
  ylab("Relative biomass (%)") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "Taxa group")) +
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


# save/load ----
# save.image("02b_Script_SavedImages/18_BiofilmCompositionOrdination_Rdat")
# load("02b_Script_SavedImages/18_BiofilmCompositionOrdination_Rdat")
