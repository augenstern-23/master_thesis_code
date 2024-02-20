# make matrix with number of cage tags for TFs


# installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("APL")
BiocManager::install("TENxPBMCData")

# Setup
library(reticulate)
library(APL)
library(Seurat)
library(ggplot2)
set.seed(1234)
library(dplyr)
library(tidyr)
library(dplyr)
install.packages("data.table")
library(data.table)


# make tissue list
files = list.files("/project/pauline_association_plots/code/01_promoter-celltype_matrix/experiments_single/tsv_all_data/complete/")
excluded_ext <- "_all_data.tsv"
tissues <- sub(paste0(excluded_ext, "$"), "", files)

# initialise empty data frame
matrix <- data.frame()

for (tissue in tissues) {

    # load dataframe with cage tags (all data)
    df_all <- read.table(paste('experiments_single/tsv_all_data/complete/', tissue,'_all_data.tsv', sep=''),
                         header=TRUE, sep="\t")
    # load fimo file
    df_fimo <- read.table(paste('experiments_single/fimo/pv-thres/0_collapsed/', tissue, 
                          '_collapsed.tsv', sep=''), header=TRUE, sep='\t')


    # make new df with unique combinations for TF and coordinates
    new_fimo <- distinct(df_fimo, motif_alt_id, sequence_name)

    # Split the "sequence_name" column into separate columns
    new_fimo <- new_fimo %>%
      separate(sequence_name, into = c("chr", "start_end"), sep = ":", remove = FALSE) %>%
      separate(start_end, into = c("start", "end"), sep = "-") %>%
      mutate(end = sub("\\(.*", "", end))

    # Merge fimo and all_data based on matching values in "chr", "start", and "end" columns
    merged_df <- merge(new_fimo, df_all[, c("chr", "start", "end", "nr_ctss")], 
                       by = c("chr", "start", "end"), all.x = TRUE)

    # Rename the "nr_ctss" column in the merged dataframe
    #colnames(merged_df)[ncol(merged_df)] <- "nr_ctss"

    # delete duplicates
    merged_df <- merged_df[!duplicated(merged_df), ]

    # Calculate the sum of "nr_ctss" values per TF
    sum_df <- aggregate(nr_ctss ~ motif_alt_id, data = merged_df, FUN = sum)
    
    # new_colname <- tissue
    colnames(sum_df)[2] <- tissue
    
    # assign tissues as rownames
    rownames(sum_df) <- sum_df[,1]
    sum_df[,1] <- NULL

    matrix <- bind_rows(matrix, data.frame(t(sum_df)))
    }


write.table(matrix, "/project/pauline_association_plots/code/01_promoter-celltype_matrix/experiments_single/results/matrix_tf_cage.tsv", 
            sep="\t", row.names=FALSE, col.names=TRUE)
