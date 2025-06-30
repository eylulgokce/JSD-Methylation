#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom data.table fread rbindlist
#' @importFrom utils read.table
NULL

#' Get Regions for Processing
#'
#' Get stepsize and list of regions for tabix-indexed files.
#'
#' @param tabixfiles Character vector of tabix-indexed file paths
#' @param chrom Character string specifying the chromosome
#' @param exp_numsites Expected number of sites for chunking (default: 1000)
#' @return A list containing progress percentages and regions
#' @export
get_regions <- function(tabixfiles, chrom = NULL, exp_numsites = 1000) {

  if (is.null(chrom)) {
    stop("Chromosome must be specified")
  }

  sup_position <- supremum_position(tabixfiles, chrom)

  if (is.null(sup_position)) {
    message("Skipping because chromosome is missing.")
    return(list(progress = NULL, regions = NULL))
  }

  sup_numsites <- supremum_numsites(tabixfiles, chrom)

  if (is.null(sup_numsites) || sup_numsites == 0) {
    message("Skipping because there are no entries.")
    return(list(progress = NULL, regions = NULL))
  }

  step <- ceiling(sup_position / sup_numsites * exp_numsites)

  if (step < sup_position) {
    stepsize <- step
  } else {
    stepsize <- sup_position
  }

  pos_start <- seq(0, sup_position - 1, by = stepsize + 1)
  pos_end <- c(seq(stepsize, sup_position - 1, by = stepsize + 1), sup_position)

  # Ensure pos_end has same length as pos_start
  if (length(pos_end) > length(pos_start)) {
    pos_end <- pos_end[1:length(pos_start)]
  }

  progress <- round(100 * pos_end / sup_position, 1)

  regions <- data.frame(
    chrom = rep(chrom, length(pos_start)),
    start = pos_start,
    end = pos_end,
    stringsAsFactors = FALSE
  )

  return(list(progress = progress, regions = regions))
}

#' Get Data from Tabix Files
#'
#' Combines tabix-indexed genome position files for specified regions.
#'
#' @param files Character vector of file paths
#' @param labels Character vector of labels for the files
#' @param data_columns List of column specifications for each file
#' @param regions Data frame with regions to query
#' @param join Character string specifying join type (default: "outer")
#' @param preset Character string specifying file format preset (default: "bed")
#' @return A list of data frames, one for each region
#' @export
get_data <- function(files, labels = NULL, data_columns = NULL, regions = NULL,
                     join = "outer", preset = "bed") {

  # Check input arguments
  if (is.null(labels)) {
    keys <- paste0("unit_", seq_along(files))
  } else if (length(labels) == length(files)) {
    keys <- labels
  } else {
    stop("Number of files and labels must match!")
  }

  if (is.null(data_columns)) {
    stop("The list of data_columns must have at least one entry!")
  } else if (length(data_columns) == length(files)) {
    # Each file has its own column specification
  } else if (length(data_columns) == 1) {
    # Replicate single column specification for all files
    data_columns <- rep(data_columns, length(files))
  } else {
    stop("Either supply a single entry in data_columns or the number of entries must match the number of files!")
  }

  # Define column indices based on preset
  if (preset == "bed") {
    index_cols <- c(1, 2, 3)  # chrom, start, end
    index_names <- c("chrom", "start", "end")
    index_types <- c("character", "integer", "integer")
  }  else {
    stop("Unknown preset: ", preset)
  }

  # Process each region
  result_list <- list()

  for (i in 1:nrow(regions)) {
    region <- regions[i, ]
    query <- paste0(region$chrom, ":", region$start, "-", region$end)

    # Read data from each file for this region
    file_data_list <- list()

    for (j in seq_along(files)) {
      file_path <- files[j]

      # Use system command to call tabix
      cmd <- paste("tabix", shQuote(file_path), shQuote(query))

      tryCatch({
        # Execute tabix command and read result
        tabix_output <- system(cmd, intern = TRUE)

        if (length(tabix_output) > 0) {
          # Parse the tabix output
          con <- textConnection(tabix_output)

          # Determine which columns to read
          all_cols <- c(index_cols, data_columns[[j]])
          all_names <- c(index_names, paste0("col_", data_columns[[j]]))

          # Read the data
          df <- read.table(con, sep = "\t", header = FALSE,
                           comment.char = "#", stringsAsFactors = FALSE)
          close(con)

          if (ncol(df) >= max(all_cols)) {
            # Select only the required columns
            df <- df[, all_cols, drop = FALSE]
            colnames(df) <- all_names

            # Set proper column types
            df$chrom <- as.character(df$chrom)
            df$start <- as.integer(df$start)
            df$end <- as.integer(df$end)

            # Set row names using genomic coordinates
            rownames(df) <- paste0(df$chrom, ":", df$start, "-", df$end)

            file_data_list[[j]] <- df
          } else {
            file_data_list[[j]] <- NULL
          }
        } else {
          file_data_list[[j]] <- NULL
        }
      }, error = function(e) {
        warning("Error reading file ", file_path, " for region ", query, ": ", e$message)
        file_data_list[[j]] <- NULL
      })
    }

    # Merge data from all files for this region
    valid_data <- file_data_list[!sapply(file_data_list, is.null)]

    if (length(valid_data) > 0) {
      if (length(valid_data) == 1) {
        merged_df <- valid_data[[1]]
        # Add sampling unit information to column names
        data_cols <- setdiff(colnames(merged_df), c("chrom", "start", "end"))
        colnames(merged_df)[colnames(merged_df) %in% data_cols] <-
          paste0(keys[which(!sapply(file_data_list, is.null))][1], ".", data_cols)
      } else {
        # Merge multiple data frames
        merged_df <- valid_data[[1]]
        base_coords <- merged_df[, c("chrom", "start", "end")]

        # Rename data columns to include sampling unit
        data_cols <- setdiff(colnames(merged_df), c("chrom", "start", "end"))
        colnames(merged_df)[colnames(merged_df) %in% data_cols] <-
          paste0(keys[which(!sapply(file_data_list, is.null))][1], ".", data_cols)

        # Merge with remaining files
        for (k in 2:length(valid_data)) {
          df_k <- valid_data[[k]]
          unit_idx <- which(!sapply(file_data_list, is.null))[k]

          # Rename data columns
          data_cols_k <- setdiff(colnames(df_k), c("chrom", "start", "end"))
          colnames(df_k)[colnames(df_k) %in% data_cols_k] <-
            paste0(keys[unit_idx], ".", data_cols_k)

          # Merge based on genomic coordinates
          if (join == "outer") {
            merged_df <- merge(merged_df, df_k,
                               by = c("chrom", "start", "end"),
                               all = TRUE)
          } else {
            merged_df <- merge(merged_df, df_k,
                               by = c("chrom", "start", "end"))
          }
        }
      }

      result_list[[i]] <- merged_df
    } else {
      # Return empty data frame with correct structure
      result_list[[i]] <- data.frame()
    }
  }

  return(result_list)
}

#' Get Maximum Number of Sites
#'
#' Return the least upper bound for the number of covered sites.
#'
#' @param tabixfiles Character vector of tabix-indexed file paths
#' @param chrom Character string specifying the chromosome
#' @return Integer representing the maximum number of sites, or NULL if no data
#' @export
supremum_numsites <- function(tabixfiles, chrom) {

  sites <- c()

  for (f in tabixfiles) {
    cmd <- paste("tabix", shQuote(f), shQuote(chrom), "| wc -l")

    tryCatch({
      site_count <- as.integer(system(cmd, intern = TRUE))
      if (!is.na(site_count)) {
        sites <- c(sites, site_count)
      }
    }, error = function(e) {
      # Continue with next file if this one fails
    })
  }

  if (length(sites) > 0) {
    return(max(sites))
  } else {
    return(NULL)
  }
}

#' Get Maximum Position
#'
#' Return the least upper bound for the chromosome end coordinate.
#'
#' @param tabixfiles Character vector of tabix-indexed file paths
#' @param chrom Character string specifying the chromosome
#' @return Integer representing the maximum position, or NULL if no data
#' @export
supremum_position <- function(tabixfiles, chrom) {

  end_coordinates <- c()

  for (f in tabixfiles) {
    cmd <- paste("tabix", shQuote(f), shQuote(chrom), "| tail -1 | cut -f3")

    tryCatch({
      base_position <- as.integer(system(cmd, intern = TRUE))
      if (!is.na(base_position)) {
        end_coordinates <- c(end_coordinates, base_position)
      }
    }, error = function(e) {
      # Continue with next file if this one fails
    })
  }

  if (length(end_coordinates) > 0) {
    return(max(end_coordinates))
  } else {
    return(NULL)
  }
}

