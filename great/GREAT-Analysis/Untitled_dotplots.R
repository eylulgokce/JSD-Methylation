library(GenomicRanges)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.At.tair.db)
library(viridis)
library(patchwork)

# Paths
file_paths <- list(
  chg = "/Users/eylul/Desktop/JSD-Methylation/data_figs/chg.txt.gz",
  cpg = "/Users/eylul/Desktop/JSD-Methylation/data_figs/cpg.txt.gz",
  chh = "/Users/eylul/Desktop/JSD-Methylation/data_figs/chh.txt.gz"
)

# GO enrichment function
run_GO <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  
  ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "ALL", 
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.05
  )
  
  return(as.data.frame(ego))
}

# Dotplot generatorß
plot_go_dotplot <- function(ego_df, title_text) {
  ego_result <- suppressWarnings(
    new("enrichResult", result = ego_df, pvalueCutoff = 0.1, 
        pAdjustMethod = "BH", qvalueCutoff = 0.05)
  )
  
  p <- dotplot(ego_result, showCategory = 10, title = title_text) +
    scale_color_viridis() +
    theme(axis.text.y = element_text(size = 5))
  return(p)
}

# Main loop
for (context in names(file_paths)) {
  data <- read.table(file_paths[[context]], header = TRUE, sep = "\t")
  data$class[is.na(data$class)] <- "unknown"
  
  region_classes <- sort(unique(data$class))
  temperatures <- c("10C", "16C", "22C")
  
  blank_plot <- ggplot() + theme_void() + ggtitle("No data")
  plot_grid <- list()
  
  for (region in region_classes) {
    for (temp in temperatures) {
      jsd_col <- paste0("JSD_bit_", temp)
      sub_data <- data[data$class == region & !is.na(data[[jsd_col]]), ]
      
      genes <- unique(sub_data$nearestTSS.gene_id)
      genes <- genes[genes != "" & !is.na(genes)]
      
      ego_df <- run_GO(genes)
      
      if (!is.null(ego_df) && nrow(ego_df) > 0) {
        p <- plot_go_dotplot(ego_df, paste(context, temp, region))
        plot_grid[[paste(region, temp)]] <- p
      } else {
        plot_grid[[paste(region, temp)]] <- blank_plot
      }
    }
  }
  
  # Build final layout: rows = region classes, cols = temperatures
  ordered_plots <- list()
  for (region in region_classes) {
    for (temp in temperatures) {
      key <- paste(region, temp)
      ordered_plots[[length(ordered_plots) + 1]] <- plot_grid[[key]]
    }
  }
  
  # Combine into one figure
  final_plot <- wrap_plots(ordered_plots, ncol = 3)
  
  # Save to PDF
  output_dir <- file.path("GO_Enrichment_Plots", context)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_path <- file.path(output_dir, paste0("GO_combined_", context, ".pdf"))
  
  ggsave(
    output_path,
    final_plot,
    width = 3 * 5,       # 3 columns × 5 inches
    height = length(region_classes) * 3.5,  # Rows × 3.5 inches
    units = "in",
    limitsize = FALSE
  )
}
