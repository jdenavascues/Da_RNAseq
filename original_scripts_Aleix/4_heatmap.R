library(dplyr)
library(purrr)
library("pheatmap")
library(cetcolor)
library(genefilter)

setwd('/Users/aleix/Documents/bioinformatics/RNAseq/output/')
basedir <- getwd()

#Read all tables in a list. Datasets are in the 'differential_expression' folder in txt format.
data_list <- list.files("differential_expression", pattern = "txt")
datasets <- pmap(list(file = paste("differential_expression",data_list[-3], sep = "/"), header=TRUE, sep="\t"), read.table)
names(datasets) <- c("DaDa", "Da", "Sc")

# Filter for DE genes
# Need to find how to do absolute using map
DE_genes <- datasets %>% 
  map(filter, padj <= 0.05) %>%
  map(filter, log2FoldChange >= 1 | log2FoldChange <= -1) 
DE_unique_genes <- unique(unlist(lapply(DE_genes, '[[', "ensemblGeneID")))

tpmNormalisedCounts <- read.table(file="tpm_normalised_data.txt", header=TRUE, sep="\t")[-c(7:9)]
colnames(tpmNormalisedCounts) <- c("Da:Da1", "Da:Da2", "Da:Da3", "Da1", "Da2", "Da3", "Sc1", "Sc2", "Sc3", "WT1", "WT2", "WT3")
tpmNormalisedCounts$EnsemblGeneID <- rownames(tpmNormalisedCounts)
DE_tpm <- tpmNormalisedCounts[match(DE_unique_genes, tpmNormalisedCounts$EnsemblGeneID),]
DE_tpm <- DE_tpm[-13]

###
means_DE <- data.frame(matrix(nrow=length(DE_unique_genes), ncol = 0))

for (i in 1:(length(colnames(DE_tpm))/3)) {
  
  means_DE <- cbind(means_DE, (rowMeans(DE_tpm[(3*(i-1)+1):(3*i)])))
}
colnames(means_DE) <- c("Da:Da", "Da", "Sc", "Wt")
###

Zscore_DE <- data.frame(matrix(nrow=length(DE_unique_genes), ncol = 0))

# choose to use all samples or the mean of the replicates for each sample
var1 <- DE_tpm



## Create Heatmap

# create color palet
col_pal <- cet_pal(n = 256, name = "cbd1", alpha = 1)

# define metrics for clustering
drows1 <- "manhattan"
dcols1 <- "manhattan"

hm.parameters <- list(var1, 
                      color = col_pal,
                      #cellwidth = 15, cellheight = 12, 
                      scale = "row", #Does the z-score
                      treeheight_row = 90,
                      cellwidth = 18,
                      kmeans_k = NA,
                      cutree_rows = 6,
                      show_rownames = F, show_colnames = T,
                      main = "Fold change 4 mean ward.D2",
                      clustering_method = "ward.D2",
                      cluster_rows = TRUE, cluster_cols = TRUE,
                      clustering_distance_rows = drows1, 
                      clustering_distance_cols = dcols1)

# prepare path and file.name for output
filename <- "Fold change 4_m_w2_6.pdf"
outfile <- paste(basedir,"heatmap", filename, sep="/")

# To draw the heat map on screen 
do.call("pheatmap", hm.parameters)

# To draw to file 
#do.call("pheatmap", c(hm.parameters, filename=outfile))

##############################

col_pal <- cet_pal(n = 256, name = "cbd1", alpha = 1)

trial <- dutta_R4_UP[c("ensemblGeneID", "TYPE")]
trial <- na.omit(trial)
trial <- trial[trial$TYPE != "NON",]
rownames(trial) <- trial$ensemblGeneID

trial <- data.frame(trial)
trial$TYPE <- as.factor(trial$TYPE)


short_list <- DE_genes %>% map(subset,ensemblGeneID %in% trial$ensemblGeneID) %>% 
  map(arrange, desc(abs(log2FoldChange))) %>% map(top_n, 100, abs(log2FoldChange)) %>% map(column_to_rownames, var="ensemblGeneID")

norm_terms <- map(short_list, select, starts_with("norm"))

carac <- c(color = col_pal, cluster_rows = TRUE, show_rownames = FALSE, annotation_row = select(trial, TYPE), scale = "row")

walk(norm_terms, pheatmap,
         color = col_pal, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,,
         annotation_row = select(trial, TYPE), 
        scale = "row")
