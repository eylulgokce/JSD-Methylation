# Load necessary libraries
library(rGREAT)
library(GenomicRanges)
library(ggplot2)
library(enrichplot)
library(dplyr)
library(readr)
library(writexl)

# Function to perform rGREAT analysis
run_great <- function(gr, db, txdb = "TxDb.Athaliana.BioMart.plantsmart51"){
  l <- sum(gr$qval <= 0.1)
  dmr <- NULL
  if (l >= 0) {
    dmr <- gr[gr$qval <= 0.1]
  } else {
    dmr <- gr[gr$pval <= 0.05]
  }
  
  res <- data.frame()
  try(
    res <- great(gr = dmr,
                 background = gr,
                 min_gene_set_size = 5, 
                 verbose = T,
                 gene_sets = db,
                 tss_source = txdb)
  )
  
  if (is(res, "data.frame")) {
    return(res)
  } else {
    return(res@table)
  }
}

# Function to generate dot plot for GO analysis
generate_dot_plot <- function(go_results, context, go_category) {
  top_go_terms <- go_results %>%
    arrange(pvalue) %>%
    head(20)
  
  plot <- ggplot(top_go_terms, aes(x = reorder(Description, -pvalue), y = -log10(pvalue), size = Count, color = pvalue)) +
    geom_point() +
    coord_flip() +
    theme_minimal() +
    labs(title = paste("Top 20 GO Terms for", context, "-", go_category),
         x = "GO Term", y = "-log10(p-value)") +
    scale_size_continuous(range = c(3, 10)) +
    scale_color_gradient(low = "blue", high = "red")
  
  return(plot)
}

# Process data for CpG context
process_context_data <- function(file_paths, context, txdb) {
  context_data <- lapply(file_paths, function(file_path) {
    data <- read.table(file_path, header = TRUE, sep = "\t")
    data <- data[!is.na(data$nearestTSS.gene_id), ]
    
    # For each temperature, we perform rGREAT analysis
    results <- lapply(c("10C", "16C", "22C"), function(temp) {
      jsd_col <- paste0("JSD_bit_", temp)
      genes <- unique(data$nearestTSS.gene_id[!is.na(data[[jsd_col]])])
      
      # Create GRanges object for rGREAT
      gr <- GRanges(
        seqnames = data$chromosome,
        ranges = IRanges(start = data$location, end = data$location),
        gene_id = genes,
        qval = data[[paste0("JSD_bit_", temp)]]
      )
      
      great_res <- run_great(gr = gr, db = "GO:BP", txdb = txdb)
      return(great_res)
    })
    
    return(results)
  })
  
  return(context_data)
}

# Read CpG data
cpg_file_paths <- list.files(path = "input", pattern = "CpG", full.names = T)
names(cpg_file_paths) <- gsub(pattern = "CpG_|\\.rds", replacement = "", x = basename(cpg_file_paths))

# Process CpG data
cpg_results <- process_context_data(cpg_file_paths, "CpG", txdb)

# Plot results for CpG
for (i in seq_along(cpg_results)) {
  for (temp in c("10C", "16C", "22C")) {
    plot <- generate_dot_plot(cpg_results[[i]][[temp]], "CpG", paste0("GO:BP - ", temp))
    print(plot)  # Display the plot
  }
}

# Process data for CHG context
chg_file_paths <- list.files(path = "input", pattern = "CHG", full.names = T)
names(chg_file_paths) <- gsub(pattern = "CHG_|\\.rds", replacement = "", x = basename(chg_file_paths))

# Process CHG data
chg_results <- process_context_data(chg_file_paths, "CHG", txdb)

# Plot results for CHG
for (i in seq_along(chg_results)) {
  for (temp in c("10C", "16C", "22C")) {
    plot <- generate_dot_plot(chg_results[[i]][[temp]], "CHG", paste0("GO:BP - ", temp))
    print(plot)  # Display the plot
  }
}

# Save results for CpG and CHG analysis as Excel files
save_results_to_excel <- function(results, context) {
  for (i in seq_along(results)) {
    result_file <- paste0("./output/", context, "/", context, "_", names(results)[i], ".xlsx")
    write_xlsx(results[[i]], path = result_file, col_names = TRUE, format_headers = TRUE)
  }
}

# Save the CpG and CHG results
save_results_to_excel(cpg_results, "CpG")
save_results_to_excel(chg_results, "CHG")
