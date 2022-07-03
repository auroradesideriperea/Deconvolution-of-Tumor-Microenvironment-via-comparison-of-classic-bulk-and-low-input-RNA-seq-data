library(tidyverse)
library(plyr)

search_folder <- "./primary tumor"
expression.files <- list.files(path = search_folder,
                               recursive = T,
                               pattern = "rna_seq.augmented_star_gene_counts.tsv",
                               full.names = T)


basename(expression.files)

sample_sheet <- read_tsv("./primary tumor/gdc_sample_sheet.2022-04-04.tsv",
                         show_col_types = FALSE)
head(sample_sheet)

read_tsv_filename <- function(filename){
  # file_path = basename(filename)
  file_path = sample_sheet$`Sample ID`[sample_sheet$`File Name`==basename(filename)]
  # read count matrix, select only columns of interest
  d <- read_tsv(filename, skip=6, col_names=F, col_select=c(X1,X2,X3,X4), show_col_types = FALSE) %>%
    dplyr::rename("gene_id"=X1,"gene_name"=X2,"gene_type"=X3,!!file_path:=X4)
  d
}

import.list <- llply(expression.files, read_tsv_filename)

df <- reduce(import.list, inner_join,
             by = c("gene_id","gene_name","gene_type"))

final_df <- df %>%
  select(contains(c("gene_id","gene_name","gene_type","01A"))) %>%
  filter(gene_type == "protein_coding") #keeping only protein coding

#setting readcounts dataframe
tcga_readcounts <- final_df %>%
  select(-"gene_id", -"gene_type") %>%
  group_by(gene_name) %>%
  summarise_each(funs(sum)) %>% #summing expression same genes
  remove_rownames %>%
  column_to_rownames(var="gene_name")