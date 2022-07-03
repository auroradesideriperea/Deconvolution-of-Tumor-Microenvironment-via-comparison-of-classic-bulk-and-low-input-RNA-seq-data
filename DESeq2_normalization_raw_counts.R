library(tidyverse)
library(DESeq2)

#TCGA
tcga_sampleinfo <- DataFrame(condition = gsub("_[0-9]+", "", names(tcga_readcounts)),
                             row.names = names(tcga_readcounts) ) #DEseq coldata

tcga_DESeq.ds <- DESeqDataSetFromMatrix(countData = tcga_readcounts,
                                        colData = tcga_sampleinfo,
                                        design = ~ condition) #creating DEseq object
dim(tcga_DESeq.ds) #dimensions: 19938

tcga_keep_genes <- rowSums(counts(tcga_DESeq.ds)) > 0
tcga_DESeq.ds <- tcga_DESeq.ds[ tcga_keep_genes, ] #removing non-expressed genes
dim(tcga_DESeq.ds) #new dimensions: 19164

tcga_DESeq.ds <- estimateSizeFactors(tcga_DESeq.ds)
tcga_normalized_readcounts <- counts(tcga_DESeq.ds, normalized=TRUE)
tcga_normalized_readcounts #normalizing by "median ratio method" 

tcga_mixture <- tibble::rownames_to_column(as.data.frame(tcga_normalized_readcounts),
                                           "gene_name") #rownames as first column

write_tsv(tcga_mixture, file = "tcga_mixture.tsv") #saving mixture file

#GEO
geo_readcounts <- read.table(file = "GSE132107_counts.txt",
                             header = TRUE) #import raw counts
geo_sampleinfo <- DataFrame(condition = gsub("_[0-9]+", "", names(geo_readcounts)),
                            row.names = names(geo_readcounts) ) #DEseq coldata

geo_DESeq.ds <- DESeqDataSetFromMatrix(countData = geo_readcounts,
                                       colData = geo_sampleinfo,
                                       design = ~ condition) #creating DEseq object
dim(geo_DESeq.ds) #dimensions: 23710

keep_genes <- rowSums(counts(geo_DESeq.ds)) > 0
geo_DESeq.ds <- geo_DESeq.ds[ keep_genes, ] #removing non-expressed genes
dim(geo_DESeq.ds) #new dimensions: 22663

geo_DESeq.ds <- estimateSizeFactors(geo_DESeq.ds)
geo_normalized_readcounts <- counts(geo_DESeq.ds, normalized=TRUE)
geo_normalized_readcounts #normalizing by "median ratio method" 

geo_mixture <- tibble::rownames_to_column(as.data.frame(geo_normalized_readcounts),
                                          "gene_name") #rownames as first column

write_tsv(geo_mixture, file = "geo_mixture.tsv") #saving mixture file 39 samples