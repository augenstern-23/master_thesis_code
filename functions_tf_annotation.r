library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(stringr)
library(data.table)

# define file paths for new files
file_paths_enhancer_cleaned <- list()
file_paths_enhancer_atac <- list()
file_paths_fasta <- list()
file_paths_fimo <- list()

for (tissue in tissues) {
  file_paths_enhancer_atac <- c(file_paths_enhancer_atac, paste("data/enhancer_and_atac/", tissue, "_enhancer_atac.txt", sep = ""))
    file_paths_enhancer_cleaned <- c(file_paths_enhancer_cleaned, paste("data/ENCODE/enhancer_cleaned/",
                                                                        tissue, "_enhancer_cleaned.bed", 
                                                                        sep = ""))
    file_paths_fasta <- c(file_paths_fasta, paste("data/fasta/", tissue, ".fa", sep = ""))
    file_paths_fimo <- c(file_paths_fimo, paste("../../../../pauline_association_plots_data/enhancer_data/fimo/", tissue, sep = ""))
}

# read and clean data frames
read_and_clean <- function(file_paths_enhancer, tissues) {
  cleaned_dataframes <- lapply(file_paths_enhancer, function(file_path) {
    df <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    df <- df[!(df$V10 == 'All cCREs with low DNase Z-scores in a particular cell type are bundled into one “inactive” state for that cell type; the remaining “active” cCREs'), ]
    return(df)
  })
  names(cleaned_dataframes) <- paste(tissues, "enhancer_cleaned", sep = "_")
  output_file_paths <- file.path("data/ENCODE/enhancer_cleaned", paste(tissues, "enhancer_cleaned.bed", sep = "_"))
  
  lapply(seq_along(cleaned_dataframes), function(i) {
    write.table(cleaned_dataframes[[i]], file = output_file_paths[i], row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
  })
  
  list2env(cleaned_dataframes, envir = .GlobalEnv)
}


# bwtool command, combine enhancer & atac
extract_bwtool_data <- function(tissue, enhancer_bed_file, bigwig_file, output_file) {
  system(
    paste("bwtool extract bed", enhancer_bed_file, bigwig_file, output_file),
    intern = TRUE
  )
}

# compute mean of atac and update column 8
calculate_mean_and_update <- function(file_paths) {
  result <- lapply(file_paths, function(file_path) {
    df <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    df$V8 <- sapply(strsplit(df$V8, ","), function(x) {
      numeric_values <- as.numeric(x)
      if (all(!is.na(numeric_values))) {
        return(mean(numeric_values))
      } else {
        return(NA)
      }
    })
    return(df)
  })
  return(result)
}

# plot histograms
plot_histograms <- function(files_enhancer_atac_mean, tissues, bins = 30) {
  for (i in seq_along(files_enhancer_atac_mean)) {
    df <- files_enhancer_atac_mean[[i]]
    tissue <- tissues[i]

    hist(df[[8]], main = paste("Histogram for", tissue, " - Row Count:", nrow(df)),
         xlab = colnames(df)[8], breaks = bins)
    cat("Dataframe:", i, "Tissue:", tissue, "\n")
  }
}

# plot logarithmic histograms
plot_histograms_log <- function(files_enhancer_atac_mean, tissues, bins = 30) {
  for (i in seq_along(files_enhancer_atac_mean)) {
    df <- files_enhancer_atac_mean[[i]]
    tissue <- tissues[i]

    hist(log(df[[8]]), main = paste("Histogram for", tissue, " - Row Count:", nrow(df)),
         xlab = colnames(df)[8], breaks = bins)
    cat("Dataframe:", i, "Tissue:", tissue, "\n")
  }
}

# transform atac value in V8 --> log, z-values, + shift_value
folder_path_bed <- "data/bed/"
file_paths_bed <- list()
transform_and_save_dataframes <- function(dataframes, tissues, shift_value, bins = 30) {
  transformed_dataframes <- lapply(seq_along(dataframes), function(i) {
    df <- dataframes[[i]]
    # Add pseudocount +1 before log to avoid infinite values
    transformed <- scale(log(df$V8 + 1)) + shift_value
    df$V9 <- transformed
    hist(transformed, main = paste("Histogram for", tissues[[i]], " - Row Count:", nrow(df)),
         xlab = colnames(df)[8], breaks = bins)
    filename <- paste0(folder_path_bed, tissues[[i]], "_shift_by_", shift_value, ".bed")
    file_paths_bed <<- append(file_paths_bed, filename)
    write.table(df, file = filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    return(df)
  })
}

# pick two chromosomes to reduce data size
reduce_chr_function <- function(dataframes, tissues, chr_value1, chr_value2) {
  transformed_dataframes <- lapply(seq_along(dataframes), function(i) {
    df <- dataframes[[i]]
    reduced_df <- df[df$V1 %in% c(chr_value1, chr_value2), ]
    output_file <- paste0("data/bed/", tissues[[i]], "_", chr_value1, "_", chr_value2, ".bed")
    write.table(reduced_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    return(reduced_df)
  })
  return(transformed_dataframes)
}

# getfasta
fasta_function <- function(tissues, sequence_file, file_paths_bed, file_paths_fasta){
    system(
        paste("bedtools getfasta -fi", sequence_file, "-bed", file_paths_bed, "-s -fo", file_paths_fasta), 
        intern = TRUE
    )
}

# fimo 
fimo_tf_annotation <- function(tissues, file_paths_fasta, file_path_motif_file, file_paths_fimo) {
  system(
    paste("fimo --o ", file_paths_fimo, " ", file_path_motif_file, " ", file_paths_fasta, sep=""),
    intern = TRUE
  )
}

collapse_fimo_function <- function(tissues, file_paths_fimo) {
  for (i in seq_along(tissues)) {
    tissue <- tissues[i]
    df_tmp <- read.table(paste(file_paths_fimo[[i]], "/fimo.tsv", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
       
    cat(tissue, " nrow before collapsing: ", nrow(df_tmp), "\n", sep="")

    df_tmp_collapsed <- df_tmp %>%
      group_by(motif_id, motif_alt_id, sequence_name, strand, score, p.value, q.value, matched_sequence) %>%
      summarize(start = min(start), stop = max(stop))

    cat(tissue, " nrow after collapsing: ", nrow(df_tmp_collapsed), ", reduction by ", round((nrow(df_tmp)) / nrow(df_tmp_collapsed), 3), "%", "\n", sep="")

    df_tmp_collapsed <- df_tmp_collapsed[, c("motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score", "p.value", "q.value", "matched_sequence")]

    output_file <- file.path(file_paths_fimo[i], "/", paste(tissue, "_collapsed.tsv", sep = ""))
    write.table(df_tmp_collapsed, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

combine_and_write_tfs <- function(tissues, file_paths, output_file_path) {
    df_combined_tf <- data.frame()
    
  for (i in seq_along(tissues)) { 
    #df_tmp <- read.table(file_paths[i], sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
      
    df_tmp <- fread(file_paths[[i]], select = c('motif_id'), data.table = TRUE)
df_tmp[, motif_id := str_extract(motif_id, "(?<=_)[^_]+(?=_)")]

    # getting the # of all TFs for one tissue
    df_tmp_tf <- table(df_tmp$motif_id)
    df_tmp_tf$tissue <- tissues[i]
      
    # bind together the tissues with TFs
    df_combined_tf <- bind_rows(df_combined_tf, df_tmp_tf)

  }
    
  # arrange columns alphabetically
    df_combined_tf <- df_combined_tf[, order(colnames(df_combined_tf))]
    df_combined_tf <- df_combined_tf %>% column_to_rownames(var ="tissue")
  
    # set row names to dataframe names
    #df_combined_tf <- df_combined_tf %>%
    #select(tissue, everything())
    #rownames(df_combined_tf) <- df_combined_tf[, 1]
    #df_combined_tf[, 1] <- NULL

  write.table(df_combined_tf, output_file_path, sep = "\t", row.names = TRUE, col.names = TRUE)
  
  return(df_combined_tf)
}

  convert_tf_to_gene <- function(df1) {
    df2 <- df_transfac_to_genes
    tissues_rows <- rownames(df1)
    df1 <- t(df1)
    df1 <- data.frame(tf = rownames(df1), df1, row.names = NULL)
    df2_select <- df2 %>% select(tf, gene)
    merged_df <- left_join(df1, df2_select, by = "tf")
    df_filtered <- merged_df %>% 
        separate_rows(gene, sep = ",\\s*") %>% 
        ungroup()
    new_df <- df_filtered %>%
      select(gene, everything(), -tf)
    new_df <- new_df %>%
      filter(!is.na(gene))
    new_df <- as.data.frame(t(new_df))
    colnames(new_df) <- new_df[1, ]
    new_df <- new_df[-1, ]
    row_names <- rownames(new_df)
    
    new_df <- as.data.frame(lapply(new_df, as.numeric))
    rownames(new_df) <- row_names
    return(new_df)
}                                 