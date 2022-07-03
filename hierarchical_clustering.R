library(tidyverse)
library(ggplot2)
library(ggdendro)

geo_mixture <- read_tsv("geo_mixture.tsv", col_names = T)

geo_dendro <- geo_mixture %>%
  column_to_rownames("gene_name") %>% #making gene_name column rownames
  t %>% #switching rownames and colnames
  dist %>% #creating distance matrix
  hclust %>% # hierarchical clustering
  as.dendrogram %>%
  dendro_data 

geo_dendro[["labels"]] <- geo_dendro[["labels"]] %>%
  mutate(patient = gsub("_.*", "", label))

# palette17 = c("#00B2FF", "#78FF69", "#DB011C", "#00007A", "#FF28C0",
#               "#FFF104", "#018511", "#DB0068", "#4E1980", "#E68700",
#               "#EBBAFA", "#8A7138", "#00DBB2", "#851901", "#496AE6",
#               "#FA8A90", "#000000")

palette = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
            "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF",
            "#7E6148FF", "#B09C85FF", "#008B45FF", "#631879FF",
            "#5F559BFF", "#A20056FF", "#808180FF", "#1B1919FF",
            "#BB0021FF")

ggplot(segment(geo_dendro)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = label(geo_dendro),
            aes(x, y, label = label, hjust = 0, color = patient),
            size = 5) +
  scale_color_manual(values =  palette, name = "Patient") +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0)) +
  theme_dendro() +
  theme(legend.position="none")
