library(GenomicRanges)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.At.tair.db)
library(viridis)
library(patchwork)

# --- Paths ---
file_paths <- list(
  chg = "/Users/eylul/Desktop/JSD-Methylation/data_figs/chg.txt.gz",
  cpg = "/Users/eylul/Desktop/JSD-Methylation/data_figs/cpg.txt.gz",
  chh = "/Users/eylul/Desktop/JSD-Methylation/data_figs/chh.txt.gz"
)

# --- GO Function ---
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

# --- Dotplot Wrapper ---
plot_go_dotplot <- function(ego_df, title_text) {
  ego_result <- suppressWarnings(
    new("enrichResult", result = ego_df, pvalueCutoff = 0.1, 
        pAdjustMethod = "BH", qvalueCutoff = 0.05)
  )
  
  p <- dotplot(ego_result, showCategory = 20, title = title_text) +
    scale_color_viridis()
  return(p)
}

# --- Main Loop ---
for (context in names(file_paths)) {
  data <- read.table(file_paths[[context]], header = TRUE, sep = "\t")
  data$class[is.na(data$class)] <- "unknown"
  region_classes <- sort(unique(data$class))
  temperatures <- c("10C", "16C", "22C")
  
  plot_list <- list()
  
  for (region in region_classes) {
    for (temp in temperatures) {
      jsd_col <- paste0("JSD_bit_", temp)
      sub_data <- data[data$class == region & !is.na(data[[jsd_col]]), ]
      
      genes <- unique(sub_data$nearestTSS.gene_id)
      genes <- genes[genes != "" & !is.na(genes)]
      
      ego_df <- run_GO(genes)
      
      if (!is.null(ego_df) && nrow(ego_df) > 0) {
        p <- plot_go_dotplot(ego_df, paste(region, temp))
        plot_list[[paste(region, temp)]] <- p
      } else {
        plot_list[[paste(region, temp)]] <- NULL
      }
    }
  }
  
  # Layout as 3 cols × 7 rows
  plots_matrix <- matrix(plot_list, nrow = length(region_classes), ncol = length(temperatures), byrow = TRUE)
  
  # Convert NULLs to empty plots
  blank_plot <- ggplot() + theme_void()
  plots_matrix[is.na(plots_matrix)] <- list(blank_plot)
  plots_filled <- lapply(plots_matrix, function(x) if (is.null(x)) blank_plot else x)
  
  # Combine plots
  all_plots <- wrap_plots(plots_filled, ncol = 3)
  
  # Save combined figure
  output_dir <- file.path("GO_Enrichment_Plots", context)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_path <- file.path(output_dir, paste0("GO_combined_", context, ".pdf"))
  
  ggsave(
    output_path,
    all_plots,
    width = 3 * 5,       # 3 columns × 5 inches
    height = 7 * 3.5,    # 7 rows × 3.5 inches
    units = "in",
    limitsize = FALSE
  )
}
