library(rGREAT)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
library(enrichplot)
library(GOSemSim)
library(pheatmap)
library(tidyr)
library(KEGGREST)


file_paths <- list(
  chg = "/shares/grossniklaus.botinst.uzh/eharputluoglu/meth1001_code_DKT/data_figs/chg.txt.gz",
  cpg = "/shares/grossniklaus.botinst.uzh/eharputluoglu/meth1001_code_DKT/data_figs/cpg.txt.gz",
  chh = "/shares/grossniklaus.botinst.uzh/eharputluoglu/meth1001_code_DKT/data_figs/chh.txt.gz"
)
output_base_dir <- "/shares/grossniklaus.botinst.uzh/eharputluoglu/meth1001_code_DKT/analysis/18_GO_heatmaps"

# Create the base output directory if it doesn't exist
if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)


run_GO <- function(gene_list, temperature) {
  if (length(gene_list) == 0) return(NULL)
  
  ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "ALL", 
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
  
  return(as.data.frame(ego))
}

classify_go <- function(ego_df, temperature, output_dir) {
  if (is.null(ego_df) || nrow(ego_df) == 0) return(NULL)
  
  ego_df$Category <- factor(ego_df$ONTOLOGY, levels = c("BP", "MF", "CC"))
  
  p <- ggplot(ego_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("GO Classification for", temperature),
         x = "GO Term", y = "-log10(p.adjust)", fill = "Category") +
    theme_minimal()
  
  filename <- paste0("GO_Classification_", temperature, ".png")
  ggsave(file.path(output_dir, filename), p, width = 15, height = 10)
  
  return(p)
}

analyze_TE_effects <- function(data, temperature) {
  data$TE_associated <- ifelse(data$class == "transposable_element", "TE", "Non-TE")
  
  p <- ggplot(data, aes(x = TE_associated, y = get(paste0("JSD_bit_", temperature)), fill = TE_associated)) +
    geom_boxplot() +
    labs(title = paste("TE Influence on JSD at", temperature), x = "Gene Type", y = "JSD") +
    theme_minimal()
  
  return(p)
}

correlate_methylation_jsd <- function(data, temperature) {
  jsd_col <- paste0("JSD_bit_", temperature)
  meth_col <- paste0("meth_", temperature)
  
  correlation <- cor(data[[meth_col]], data[[jsd_col]], use = "complete.obs")
  return(correlation)
}

identify_kegg <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  
  kegg_res <- enrichKEGG(
    gene = gene_list,
    organism = "ath", 
    keyType = "kegg",
    pvalueCutoff = 0.05
  )
  
  return(as.data.frame(kegg_res))
}

cluster_genes <- function(data) {
  jsd_data <- data[, c("JSD_bit_10C", "JSD_bit_16C", "JSD_bit_22C")]
  jsd_data <- data.frame(lapply(jsd_data, as.numeric))
  jsd_matrix <- as.matrix(jsd_data)
  valid_rows <- complete.cases(jsd_data) & rowSums(is.finite(jsd_matrix)) == ncol(jsd_data) & !apply(jsd_data, 1, function(row) any(is.na(row)))
  jsd_data <- jsd_data[valid_rows, ]
  if (nrow(jsd_data) < 2) {
    warning("Not enough data for clustering after removing problematic values.")
    return(NULL)
  }
  dist_matrix <- dist(jsd_data, method = "euclidean")
  if (any(!is.finite(dist_matrix))) {
    warning("Distance matrix contains NA/NaN/Inf values.")
    return(NULL)
  }
  clustering <- hclust(dist_matrix, method = "ward.D2")
  return(clustering)
}

visualize_heatmap <- function(data, output_dir, context_name) {
  jsd_cols <- c("JSD_bit_10C", "JSD_bit_16C", "JSD_bit_22C")
  for (col in jsd_cols) {
    data[[col]] <- as.numeric(data[[col]])
  }
  valid_rows <- complete.cases(data[, jsd_cols]) & rowSums(is.finite(as.matrix(data[, jsd_cols]))) == length(jsd_cols)
  data <- data[valid_rows, ]
  zero_var_cols <- sapply(jsd_cols, function(col) var(data[[col]], na.rm = TRUE) == 0)
  if (any(zero_var_cols)) {
    warning(paste("Columns with zero variance:", paste(jsd_cols[zero_var_cols], collapse = ", ")))
    warning("Heatmap generation skipped due to zero variance.")
    return()
  }
  if (nrow(data) > 1){
    pheatmap(as.matrix(data[, jsd_cols]),
             cluster_rows = TRUE, cluster_cols = TRUE,
             main = paste("Methylation Heatmap -", context_name))
    ggsave(file.path(output_dir, paste0("Heatmap_", context_name, ".png")))
  }else{
    warning("Not enough data to create heatmap.")
  }
}



for (context in names(file_paths)) {
  data <- read.table(file_paths[[context]], header = TRUE, sep = "\t")
  output_dir <- file.path(output_base_dir, context)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  temperatures <- c("10C", "16C", "22C")
  for (temp in temperatures) {
    jsd_col <- paste0("JSD_bit_", temp)
    genes <- unique(data$nearestTSS.gene_id[!is.na(data[[jsd_col]])])
    ego_df <- run_GO(genes, temp)
    if (!is.null(ego_df) && nrow(ego_df) > 0) {
      classify_go(ego_df, temp, output_dir)
    }
    kegg_df <- identify_kegg(genes)
    if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
      write.csv(kegg_df, file.path(output_dir, paste0("KEGG_Pathways_", temp, ".csv")))
    }
    te_effect <- analyze_TE_effects(data, temp)
    correlation <- correlate_methylation_jsd(data, temp)
  }
  clustering <- cluster_genes(data)
  if (!is.null(clustering)) {
    png(file.path(output_dir, paste0("Clustering_Dendrogram_", context, ".png"))) # Save clustering plot
    plot(clustering, main = paste("Clustering Dendrogram -", context))
    dev.off()
  }
  visualize_heatmap(data, output_dir, context)
  if (context == names(file_paths)[1]) {
    all_data <- data[, c("ID", "JSD_bit_10C", "JSD_bit_16C", "JSD_bit_22C")]
  } else {
    all_data <- rbind(all_data, data[, c("ID", "JSD_bit_10C", "JSD_bit_16C", "JSD_bit_22C")])
  }
}

all_output_dir <- output_base_dir
visualize_heatmap(all_data, all_output_dir, "Combined_Contexts")