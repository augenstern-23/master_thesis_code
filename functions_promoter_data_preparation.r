install.packages('dplyr')
library('dplyr')
library('magrittr')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("FANTOM3and4CAGE")

library(FANTOM3and4CAGE)

# create single dataframes per tissue and save them
new_dfs <- list()
save_tissue_dfs <- function(dataset){

    for(i in seq_along(dataset)){
      tissue_name <- names(dataset[i])
      experiment_names <- names(dataset[[i]])[4:length(names(dataset[[i]]))]
        print(names(dataset[[i]])[4:length(names(dataset[[i]]))])
      
      for(j in seq_along(experiment_names)){
        new_df <- dataset[[i]][, c("chr", "pos", "strand", experiment_names[j])]
        colnames(new_df)[4] <- paste0(tissue_name, "_", j)
        new_dfs[[paste0(tissue_name, "_", j)]] <- new_df
        colnames(new_df) <- NULL
          
        write.table(new_df, paste("/project/pauline_association_plots/code/01_promoter-celltype_matrix/experiments_single/ctss/", tissue_name, "_", j, ".ctss", sep=""), sep="\t", row.names=FALSE, quote=FALSE) 
      }
    }
    return(new_dfs)
}

# show sum of column named like tissue --> not existent in every df
count_colsum_named_like_tissue <- function(df_list){

    for (i in seq_along(df_list)) {
      df_name <- names(df_list)[i]
      df <- df_list[[i]]
      col_sum <- sum(df[[df_name]])
      cat("Sum of", df_name, "column:", col_sum, "\n")
    } 
}

# get tissue names for change of promoter regions
get_tissue_names <- function(file_path){
    all_files <- list.files(path = file_path)
    excluded_ext <- ".bed"
    tissues <- sub(paste0(excluded_ext, "$"), "", all_files)
    return(tissues)
}

# change regions (ctss clusters) to promoter regions and save files
define_promoter_size <- function(bp_upstream, bp_downstream, tissues){
    for (bedfile in tissues) {
    tmp_bed <- read.table(paste(bedfile, sep = ""), header = FALSE, skip = 1, sep="\t", stringsAsFactors=FALSE, quote="")
    tmp_bed$V2 <- tmp_bed$V2 - 200
    tmp_bed$V3 <- tmp_bed$V3 + 50
    write.table(tmp_bed, paste("../promoter_m200_p50/", bedfile, sep = ""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
}
