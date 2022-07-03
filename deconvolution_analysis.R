library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)

#stacked barplots
palette <- c("#fce94f", "#eecc00", "#c4a000", "#fcaf3e", "#ce5c00", 
             "#e9b96e", "#c17d11", "#8ae234", "#73d216", "#4e9a06", 
             "#729fcf", "#204a87", "#ad7fa8", "#75507b", "#ef3535", 
             "#cc0000", "#a40000", "#eeeeec", "#babdb6", "#bfc1bb", 
             "#888a85", "#2e3436")

tcga_cibersort <- read_csv("tcga_cibersort.csv", col_names = T)
tcga_barplot <- tcga_cibersort %>%
  select(-c(`P-value`, `Correlation`, `RMSE`)) %>%
  melt(value.name = "Relative Percent",
       variable.name = "Cell Type",
       id = "Mixture")

ggplot(tcga_barplot, aes(fill = `Cell Type`,
                         y = `Relative Percent`,
                         x = `Mixture`)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   margin = margin(-15,1,1,1)),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        panel.background = element_blank()) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(values = palette)

geo_cibersort <- read_csv("geo_cibersort.csv", col_names = T)
geo_barplot <- geo_cibersort %>%
  select(-c(`P-value`, `Correlation`, `RMSE`)) %>%
  melt(value.name = "Relative Percent",
       variable.name = "Cell Type",
       id = "Mixture")

ggplot(geo_barplot, aes(fill = `Cell Type`,
                        y = `Relative Percent`,
                        x = `Mixture`)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   margin = margin(-15,1,1,1)),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        panel.background = element_blank()) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(values = palette)

geo_17_cibersort <- read_csv("geo_17_cibersort.csv", col_names = T)
geo_17_barplot <- geo_17_cibersort %>%
  select(-c(`P-value`, `Correlation`, `RMSE`)) %>%
  melt(value.name = "Relative Percent",
       variable.name = "Cell Type",
       id = "Mixture")

ggplot(geo_17_barplot, aes(fill = `Cell Type`,
                           y = `Relative Percent`,
                           x = `Mixture`)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   margin = margin(-15,1,1,1)),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        panel.background = element_blank()) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(values = palette)

#dataset all samples-datasets
tcga_geo <- dplyr::bind_rows(list(TCGA=tcga_cibersort,
                                  GEO=geo_cibersort),
                             .id = 'Dataset')

rm(list=ls()[! ls() %in% c("geo_cibersort","tcga_cibersort", "tcga_geo")])

dataset_var_val <- tcga_geo %>%
  select(-`P-value`, -RMSE, -Correlation, -Mixture) %>%
  melt(id.var = "Dataset", value.name = "Relative Percent",
       variable.name = "Cell Type")

#facet boxplot
ggplot(data = dataset_var_val, aes(x = `Cell Type`,
                                   y = `Relative Percent`)) + 
  geom_boxplot(aes(fill = Dataset, color = Dataset), outlier.size = 0.5) +
  scale_color_manual(values = c("#cb4740", "#76a4dc")) +
  scale_fill_manual(values = c("#e49d99", "#caddf1")) +
  facet_wrap(~`Cell Type`, scales="free", ncol = 4) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = NULL, y = NULL) +
  stat_compare_means(aes(group = Dataset),
                     method = "t.test",
                     label = "p.format",
                     vjust = -2.0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.position = c(0.88, 0.07),
        strip.text.x = element_text(face = "bold.italic"))
