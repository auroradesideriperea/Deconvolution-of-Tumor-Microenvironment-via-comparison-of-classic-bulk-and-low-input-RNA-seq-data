library(tidyverse)
library(ggplot2)
library(reshape2)
library(readxl)
library(splitstackshape)
library(ggpubr)
library(pheatmap)

signatures <- read_xlsx("immune-checkpoint inhibitors-related inflammation signatures.xlsx",
                        skip = 1) %>%
  select("signature" = "Name", "gene_name" = "Genes") %>%
  head(n = 8) %>%
  cSplit("gene_name", sep = ",", direction = "long") %>%
  mutate("gene_name" = gsub("\\(\\w+\\)", "", gene_name))

#tcga dataset
tcga_mixture <- read_tsv("tcga_mixture.tsv", col_names = T)

tcga_signatures <- tcga_mixture %>%
  inner_join(signatures, by = "gene_name")%>%
  relocate(signature, .before = gene_name)

tcga_signatures_gm <- tcga_signatures %>%
  select(-"gene_name") %>%
  group_by(signature) %>%
  summarise_each(funs(exp(mean(log(.[.>0])))))

write_tsv(tcga_signatures_gm, file = "tcga_signatures_gm.tsv")

#geo dataset
geo_mixture <- read_tsv("geo_mixture.tsv", col_names = T)

geo_signatures <- geo_mixture %>%
  inner_join(signatures, by = "gene_name")%>%
  relocate(signature, .before = gene_name)

geo_signatures_gm <- geo_signatures %>%
  select(-"gene_name") %>%
  group_by(signature) %>%
  summarise_each(funs(exp(mean(log(.[.>0])))))

write_tsv(geo_signatures_gm, file = "geo_signatures_gm.tsv")

#boxplot
rm(list=ls()[! ls() %in% c("geo_signatures_gm","tcga_signatures_gm")])

signatures_tcga_geo <- dplyr::bind_rows(list(TCGA=tcga_signatures_gm,
                                             GEO=geo_signatures_gm),
                                        .id = 'Dataset')

signatures_dataset_var_val <- signatures_tcga_geo %>%
  melt(id.var = c("Dataset", "signature"),
       value.name = "Geometric mean",
       variable.name = "Sample",
       na.rm = T)

ggplot(data = signatures_dataset_var_val, aes(x = `signature`,
                                              y = `Geometric mean`)) +
  geom_boxplot(aes(color = Dataset, fill = Dataset), outlier.size = 0.5) +
  scale_color_manual(values = c("#cb4740", "#76a4dc")) +
  scale_fill_manual(values = c("#e49d99", "#caddf1")) +
  facet_wrap(~`signature`, scales="free", ncol = 4) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  labs(x = NULL, y = NULL) +
  stat_compare_means(aes(group = Dataset),
                     method = "t.test",
                     label = "p.format",
                     vjust = -2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.position = "bottom",
        strip.text.x = element_text(face = "bold.italic"))

#base R
signatures_heatmap <- cbind(tcga_signatures_gm, geo_signatures_gm) %>%
  column_to_rownames("signature") %>%
  select(-signature) %>%
  data.matrix()

pheatmap(signatures_heatmap,
        scale= "column",
        fontsize_row = 13,
        cluster_cols=F)

#ggplot2
# ggplot(data = signatures_dataset_var_val,
#        mapping = aes(x = `signature`,
#                      y = `Sample`,
#                      fill = `Geometric mean`)) +
#   geom_tile(color = "white") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,
#                                    margin(1,1,1,1))) +
#   labs(x = NULL, y = NULL) +
#   scale_fill_gradient(low = "yellow", high = "red")