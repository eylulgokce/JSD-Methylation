library(rGREAT)
library(dplyr)
library(DT)
library(writexl)
library(KEGGREST)
library(GenomicRanges)
library(data.table)
library(TxDb.Athaliana.BioMart.plantsmart51)
txdb <- TxDb.Athaliana.BioMart.plantsmart51

## Reactome pathways
reactome <- "https://plantreactome.gramene.org/download/current/gene_ids_by_pathway_and_species.tab"
react <- data.frame(data.table::fread(input = reactome, header = F, nThread = 16))
rdb <- react[grep(pattern = "^R-ATH", x = react$V1),]
reactome_pathways <- split(rdb$V4, paste(rdb$V1, rdb$V2, sep = ": "))

## KEGG pathways
kg <- keggList("organism")
pathway2gene = keggLink("pathway", "ath")
pathwayName <- keggList("pathway","ath") 
df1 <- data.frame(gene = gsub("ath:", "", names(pathway2gene)),
                  pathID = gsub("path:", "", pathway2gene))
df2 <- data.frame(pathID = gsub("path:", "", names(pathwayName)),
                  name = pathwayName)
df_kegg <- merge(df2, df1)
kegg_pathways = split(df_kegg$gene, paste(df_kegg$pathID, df_kegg$name, sep = ": "))

run_great <- function(gr, background_gr = NULL, db, txdb = "TxDb.Athaliana.BioMart.plantsmart51"){
  
  dmr <- gr[!is.na(gr$JSD_bit)]
  
  cat(paste("Using", length(dmr), "MSC regions\n"))
  
  if(length(dmr) < 5) {
    cat("Not enough regions for analysis\n")
    return(data.frame())
  }
  
  if(is.null(background_gr)) {
    background_gr <- dmr  
  } else {
    cat(paste("Using background with", length(background_gr), "regions\n"))
  }
  
  res <- data.frame()
  try({
    res <- great(gr = dmr,
                 background = background_gr,
                 min_gene_set_size = 5, 
                 verbose = TRUE,
                 gene_sets = db,
                 tss_source = txdb)
  })
  
  if(is(res, "data.frame")){
    return(res)
  } else{
    return(res@table)
  }
}

convert_to_dmr_format <- function(data, temp_col, context_type, width = 200) {
  data_valid <- data[!is.na(data[[temp_col]]), ]
  
  if(nrow(data_valid) == 0) {
    cat(paste("No data available for", temp_col, "in", context_type, "\n"))
    return(NULL)
  }
  
  gr <- GRanges(
    seqnames = data_valid$chromosome,
    ranges = IRanges(start = pmax(1, data_valid$location - width/2), 
                     end = data_valid$location + width/2),
    JSD_bit = data_valid[[temp_col]],
    gene_id = data_valid$nearestTSS.gene_id,
    class = data_valid$class
  )
  
  gr$name <- paste("region", 1:length(gr), sep = "_")
  
  cat(paste("Created", length(gr), "regions for", temp_col, "in", context_type, "\n"))
  cat(paste("JSD range:", round(min(gr$JSD_bit, na.rm = TRUE), 3), 
            "to", round(max(gr$JSD_bit, na.rm = TRUE), 3), "\n"))
  
  return(gr)
}

create_genome_background <- function(txdb, tile_width = 200) {
  
  seqinfo_txdb <- seqinfo(txdb)
  chr_lengths <- seqlengths(seqinfo_txdb)
  
  background_ranges <- GRangesList()
  
  for(chr in names(chr_lengths)) {
    if(!is.na(chr_lengths[chr])) {
      chr_tiles <- GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = seq(1, chr_lengths[chr], by = tile_width),
          width = tile_width
        )
      )
      background_ranges[[chr]] <- chr_tiles
    }
  }
  
  background_gr <- unlist(background_ranges)
  background_gr$name <- paste("bg_region", 1:length(background_gr), sep = "_")
  
  return(background_gr)
}

setwd("/shares/grossniklaus.botinst.uzh/eharputluoglu/Meth_heatmaps/en_heatmaps")

contexts <- c("cpg", "chg", "chh")
temperatures <- c("10C", "16C", "22C")
data_list <- list()

for (context in contexts) {
  file_path <- paste0(context, ".txt.gz")
  if (file.exists(file_path)) {
    data_list[[context]] <- fread(file_path)
    cat(paste("Loaded", context, "data with", nrow(data_list[[context]]), "rows\n"))
  } else {
    cat(paste("File not found:", file_path, "\n"))
  }
}

dir.create("output_2", showWarnings = FALSE)
for(context in contexts) {
  dir.create(paste0("output_2/", toupper(context)), showWarnings = FALSE, recursive = TRUE)
}

# genome_background <- create_genome_background(txdb)
# cat(paste("Created genome background with", length(genome_background), "regions\n"))

for(context in contexts) {
  if(is.null(data_list[[context]])) next
  
  cat(paste("\n=== Processing", toupper(context), "context ===\n"))
  
  context_data <- list()
  
  for(temp in temperatures) {
    temp_col <- paste0("JSD_bit_", temp)
    
    if(temp_col %in% colnames(data_list[[context]])) {
      gr_temp <- convert_to_dmr_format(data_list[[context]], temp_col, context)
      
      if(!is.null(gr_temp) && length(gr_temp) >= 5) {
        context_data[[temp]] <- list(dmr = gr_temp)
        cat(paste("Prepared", temp, "data with", length(gr_temp), "regions\n"))
      }
    }
  }
  
  if(length(context_data) > 0) {
    res_context <- list(
      BP = lapply(context_data, function(x) run_great(gr = x$dmr, db = "GO:BP")),
      MF = lapply(context_data, function(x) run_great(gr = x$dmr, db = "GO:MF")),
      CC = lapply(context_data, function(x) run_great(gr = x$dmr, db = "GO:CC")),
      KEGG = lapply(context_data, function(x) run_great(gr = x$dmr, db = kegg_pathways)),
      Reactome = lapply(context_data, function(x) run_great(gr = x$dmr, db = reactome_pathways))
    )
    
    for(i in seq_along(res_context)){
      n <- paste0("./output_2/", toupper(context), "/", toupper(context), "_", names(res_context)[i], ".xlsx")
      write_xlsx(x = res_context[[i]], path = n, col_names = TRUE, format_headers = TRUE)
      cat(paste("Saved:", n, "\n"))
    }
  }
}

cat("\n=== Analysis by Genomic Region and Temperature ===\n")

for(context in contexts) {
  if(is.null(data_list[[context]])) next
  
  unique_regions <- unique(data_list[[context]]$class[!is.na(data_list[[context]]$class)])
  unique_regions <- unique_regions[unique_regions != ""]
  
  for(region in unique_regions) {
    cat(paste("Processing", region, "in", context, "\n"))
    
    for(temp in temperatures) {
      temp_col <- paste0("JSD_bit_", temp)
      
      if(temp_col %in% colnames(data_list[[context]])) {
        region_temp_data <- data_list[[context]][
          data_list[[context]]$class == region & 
            !is.na(data_list[[context]]$class) &
            !is.na(data_list[[context]][[temp_col]]), ]
        
        if(nrow(region_temp_data) >= 5) {
          cat(paste("  ", temp, ":", nrow(region_temp_data), "regions\n"))
          
          gr_region_temp <- convert_to_dmr_format(
            region_temp_data, temp_col, paste(context, region, temp, sep = "_")
          )
          
          if(!is.null(gr_region_temp)) {
            res_region_temp <- list(
              BP = run_great(gr = gr_region_temp, db = "GO:BP"),
              MF = run_great(gr = gr_region_temp, db = "GO:MF"),
              CC = run_great(gr = gr_region_temp, db = "GO:CC"),
              KEGG = run_great(gr = gr_region_temp, db = kegg_pathways),
              Reactome = run_great(gr = gr_region_temp, db = reactome_pathways)
            )
            
            for(analysis_type in names(res_region_temp)) {
              if(nrow(res_region_temp[[analysis_type]]) > 0) {
                output_2_file <- paste0("output_2/", toupper(context), "/", 
                                      toupper(context), "_", region, "_", temp, "_", 
                                      analysis_type, ".xlsx")
                write_xlsx(x = list(results = res_region_temp[[analysis_type]]), 
                           path = output_2_file, col_names = TRUE, format_headers = TRUE)
                cat(paste("    Saved:", basename(output_2_file), "\n"))
              }
            }
          }
        } else {
          cat(paste("  ", temp, ": insufficient data (", nrow(region_temp_data), " regions)\n"))
        }
      }
    }
  }
}

cat("\n=== Analysis Summary ===\n")
output_2_files <- list.files("output_2", pattern = "*.xlsx", recursive = TRUE, full.names = TRUE)
cat(paste("Total output_2 files created:", length(output_2_files), "\n"))

for(context in contexts) {
  if(!is.null(data_list[[context]])) {
    cat(paste("\n=== Summary for", toupper(context), "===\n"))
    cat(paste("Total regions:", nrow(data_list[[context]]), "\n"))
    
    for(temp in temperatures) {
      temp_col <- paste0("JSD_bit_", temp)
      if(temp_col %in% colnames(data_list[[context]])) {
        valid_data <- data_list[[context]][[temp_col]][!is.na(data_list[[context]][[temp_col]])]
        if(length(valid_data) > 0) {
          cat(paste(temp, "- Valid regions:", length(valid_data), 
                    "Mean JSD:", round(mean(valid_data, na.rm = TRUE), 3),
                    "Range:", round(min(valid_data), 3), "-", round(max(valid_data), 3), "\n"))
        }
      }
    }
    
    cat("Genomic regions distribution:\n")
    region_counts <- table(data_list[[context]]$class, useNA = "always")
    print(region_counts)
    
    # Temperature x Region 
    cat("\nData availability by temperature and region:\n")
    for(temp in temperatures) {
      temp_col <- paste0("JSD_bit_", temp)
      if(temp_col %in% colnames(data_list[[context]])) {
        temp_data <- data_list[[context]][!is.na(data_list[[context]][[temp_col]]), ]
        if(nrow(temp_data) > 0) {
          cat(paste(temp, "by region:\n"))
          temp_region_table <- table(temp_data$class, useNA = "always")
          print(temp_region_table)
        }
      }
    }
  }
}

cat("\n=== Analysis Complete ===\n")