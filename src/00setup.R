if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DOSE")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(DOSE)
library(enrichplot)
library(clusterProfiler)
# Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers
library(org.Hs.eg.db)
library(ggplot2)
library(ggnewscale)
library(ggupset)

data(geneList)
head(geneList)

de <- names(geneList)[abs(geneList) > 2]
head(de)


 