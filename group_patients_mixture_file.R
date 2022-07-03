library(tidyverse)
library(reshape2)

geo_mixture <- read_tsv("geo_mixture.tsv", col_names = T)

geo_mixture_17 <- geo_mixture %>%
  melt(variable.name = "sample") %>%
  mutate(sample = gsub("_.*", "", sample)) %>%
  group_by(sample, gene_name) %>%
  summarise_each(funs(mean)) %>%
  pivot_wider(names_from = sample)

write_tsv(geo_mixture_17, file = "geo_mixture_17.tsv") #saving mixture file