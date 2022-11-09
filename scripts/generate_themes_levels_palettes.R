
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Save cancer levels
cancer_levels_12 <- 
  c('AML','CLL','LYMPH','MYEL','LUNGC','CRC','GLIOM','PRC','BRC', 'CVX','ENDC','OVC')

cancers_12_mapping <- 
  data.frame(GROUP = cancer_levels_12,
             Cancer = c("AML",
                        "CLL",
                        "DLBCL",
                        "Myeloma",
                        "Lung",
                        "Colorectal",
                        "Glioma",
                        "Prostate",
                        "Breast",
                        "Cervical",
                        "Endometrial",
                        "Ovarian")) 


levels <- list("cancers_12" = cancer_levels_12,
               "cancers_12_mapping" = cancers_12_mapping)

saveRDS(levels, "data/processed/others/levels.rds")


# Save palettes
group_pal <- 
  c("AML" = "#A6CEE3",
    "CLL" = "#2271B5",
    "LYMPH" = "#08585A",
    "MYEL" = "#66C2A5",
    "CRC" = "#B89B74",
    "LUNGC" = "#ADC74F",
    "GLIOM" = "#FFD321",
    "BRC" = "#E8A29A",  
    "CVX" =  "#9E0142", 
    "ENDC" = "#B195AE",
    "OVC" = "#603479",
    "PRC" = "#E7662B")

group_pal_comp <- 
  group_pal %>% 
  enframe("GROUP","color") %>% 
  left_join(cancers_12_mapping, by = "GROUP") %>% 
  mutate(Cancer = factor(Cancer, cancers_12_mapping$Cancer)) %>% 
  arrange(Cancer) %>% 
  select(Cancer,color) %>% 
  deframe()

sex_pal <- 
  c("Male" = "#04D2AE",
    "Female" = "#F46B97") 

stage_pal <- c("0"= "#8FCC90", 
               "1" = "#FEC00E", 
               "2" = "#F17025", 
               "3" = "#C12026", 
               "4" = "#751012", 
               "Unknown" = "grey66")

heatmap_pal <-
  brewer.pal(9, name = "YlOrRd") %>% 
  colorRampPalette()

auc_pal <-
  c("all proteins" = "#20639B",
    "protein panel" = "#3CAEA3",
    "top 3" = "#F6D55C", 
    "top 1" = "#ED553B")

palettes <- list("alt_group" = group_pal,
                 "alt_group_complete" = group_pal_comp,
                 "sex" = sex_pal,
                 "stage" = stage_pal,
                 "heatmap" = heatmap_pal,
                 "auc" = auc_pal)

saveRDS(palettes, "data/processed/others/palettes.rds")


# Save ggplot themes
theme_main <- theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.spacing = unit(0.2, "lines"), 
                    panel.background=element_rect(fill="white"),
                    panel.border = element_blank(),
                    plot.title = element_text(face = "bold",
                                              size = rel(1), hjust = 0.5),
                    plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.ticks = element_line(),
                    axis.ticks.length = unit(.25, "cm"),
                    axis.line = element_line(size = 0.5),
                    axis.text = element_text(size = rel(1), color = 'black'),
                    legend.key = element_blank(),
                    legend.position = "right",
                    legend.text = element_text(size=rel(0.8)),
                    legend.key.size= unit(0.7, "cm"),
                    legend.title = element_text(size=rel(1)),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="grey90",fill="grey90"),
                    strip.text = element_text(face="bold"))


theme_simple <- 
  theme_main +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        strip.background = element_rect(color="white", fill="white"))

plot_themes <- list("main" = theme_main,
                    "simple" = theme_simple)

saveRDS(plot_themes, "data/processed/others/plot_themes.rds")

