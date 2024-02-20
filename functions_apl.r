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
library(dplyr)
library(tidyr)
library(dplyr)
library(data.table)
set.seed(1234)


genes_tissuespec <- read.table("/project/pauline_association_plots/code/01_promoter-celltype_matrix/genes_tissuespec.sorted.tsv", header = TRUE, sep="\t")

colnames(genes_tissuespec) <- c("cluster", "size", "clones", "tissue", "symbol", "description", "#libs", "ratio", "entropy")
genes_tissuespec$tissue <- trimws(genes_tissuespec$tissue)
df_transfac_to_genes <- read.table("/project/pauline_association_plots/code/01_promoter-celltype_matrix/enhancer_analysis/data/transfac_to_genesymbol.csv", sep="\t", header=TRUE)


ca_function <- function(df){
    df <- as.data.frame(t(data.matrix(df)))
    ca <- cacomp(obj = as.matrix(df), 
                 top = nrow(df), 
                 python = TRUE, 
                 rm_zeros = FALSE)
    return(ca)
    }

ca_3Dplot_function <- function(df){
    ca <- ca_function(df)
    ca_plot <- ca_3Dplot(ca, type = "ggplot")
    return(ca_plot)
    }


ca_tmp_function <- function(df, dimensions, tissue){
    ca <- ca_function(df)
    df <- as.data.frame(t(data.matrix(df)))
    pick_dims(ca,
              method = "scree_plot")
    pd <- pick_dims(ca,
              mat = as.matrix(df),
              method = "elbow_rule",
              reps = 5, 
              python = TRUE)
    print(pd)
    D <- cacomp_slot(ca, "D")
    expl_inertia <- (D^2/sum(D^2))*100
    sum(expl_inertia[seq_len(pd)])
    ca_dims <- subset_dims(ca, dims = dimensions)
    ca_tmp <- apl_coords(ca_dims, group = tissue)
    ca_tmp <- apl_score(ca_tmp,
                          mat = as.matrix(df),
                          reps = 5,
                          python = TRUE) 
    return(ca_tmp)
    }

df_apl_score <- function(df, dimensions, tissue){
    ca_tmp <- ca_tmp_function(df, dimensions, tissue)
    df_apl_score <- cacomp_slot(ca_tmp, "APL_score")
    return(df_apl_score)
    }
    


apl_function <- function(df, dimensions, tissue){
    ca_tmp <- ca_tmp_function(df, dimensions, tissue)
    df_apl_score <- df_apl_score(df, dimensions, tissue)
    gene_names <- df_apl_score[1:50, ]$Rowname
    print(cat("the top tfs for tissue: '", tissue, "' are : ", "\n", sep=""))
    print(gene_names)
    
    df_top50 <- as.list(gene_names)
    apl_plot <- apl(ca_tmp,
                    show_score = TRUE,
                    type = "ggplot") # type = "plotly" for an interactive plot
    return(apl_plot)
    }

apl_function_plotly <- function(df, dimensions, tissue){
    ca_tmp <- ca_tmp_function(df, dimensions, tissue)
    df_apl_score <- df_apl_score(df, dimensions, tissue)
    gene_names <- df_apl_score[1:50, ]$Rowname
    print(cat("the top tfs for tissue: '", tissue, "' are : ", "\n", sep=""))
    #print(gene_names)
    # Print the genes
    cat(paste(as.list(gene_names), collapse = "\n"))
    df_top50 <- as.list(gene_names)
    apl_plot <- apl(ca_tmp,
                    show_score = TRUE,
                    type = "plotly") # type = "plotly" for an interactive plot
    df_top50
    return(apl_plot)
    }

apl_tf_check <- function(df, dimensions, spec_tissue, genes_tissuespec){
    tf_matrix <- t(data.matrix(df))
    df <- as.data.frame(tf_matrix)
    ca <- cacomp(obj = as.matrix(df),
             top = nrow(df),
             python = TRUE,
             rm_zeros = FALSE)
    pick_dims(ca,
              method = "scree_plot")
    pd <- pick_dims(ca,
              mat = as.matrix(df),
              method = "elbow_rule",
              reps = 5, 
              python = TRUE)
    D <- cacomp_slot(ca, "D")
    expl_inertia <- (D^2/sum(D^2))*100
    ca_dims <- subset_dims(ca, dims = dimensions)
    tissue_tmp <- spec_tissue
    ca_tmp <- apl_coords(ca_dims, group = tissue_tmp)
    ca_tmp <- apl_score(ca_tmp,
                          mat = as.matrix(df),
                          reps = 5,
                          python = TRUE)
    gene_names <- cacomp_slot(ca_tmp, "APL_score")[1:200, ]$Rowname
    df_top50 <- as.list(gsub("\\.", "", gene_names))
    result_list <- lapply(df_top50, function(tf) {
        subset_tissue <- subset(genes_tissuespec[genes_tissuespec$tissue == spec_tissue, ])
        matches <- grepl(tf, subset_tissue$description)
        if (sum(matches) > 0) {
            list(TF = tf, count = sum(matches))
        } else {
            NULL
        }
    })
    result_list <- result_list[sapply(result_list, function(result) !is.null(result))]
    found_count <- sum(sapply(result_list, function(result) result$count > 0))
    cat("Number of TFs found:", found_count, "\n")
    for (result in result_list) {
        cat(result$TF, ": count", result$count, "\n")
    }
}


search_TF <- function(tissue, tf_list) {
    result_list <- lapply(tf_list, function(tf) {
        matches <- subset(genes_tissuespec[genes_tissuespec$tissue == tissue, ], grepl(tf, genes_tissuespec$description))
        list(TF = tf, count = nrow(matches), lines = matches)
    })
    
    for (result in result_list) {
        cat("Results for TF:", result$TF, "\n")
        cat("Count of matches:", result$count, "\n")
        if (result$count > 0) {
            print(result$lines, row.names = FALSE)
        } else {
            cat("No results found\n")
        }
        cat("\n")
    }
}


# convert df with TF names to Gene names
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


 # calculate z-scores per gene and tissue for file: "rna_tissue_consensus" from Human Protein Atlas
z_score_normalization <- function(df) {
  result <- list()
  for (gene in unique(df$Gene.name)) {
    gene_subset <- df[df$Gene.name == gene, ]
    gene_subset$nTPM_log <- log(gene_subset$nTPM+ +1)
    mean_value <- mean(gene_subset$nTPM_log)
    std_dev <- sd(gene_subset$nTPM_log)
    gene_subset$Z_score <- (gene_subset$nTPM_log - mean_value) / std_dev
    #gene_subset$Is_specific <- gene_subset$Z_score > threshold
    result[[gene]] <- gene_subset
  }
  combined_df <- do.call(rbind, result)  # Combine dataframes into one
  return(combined_df)
}
                              
# Calculate the z-score threshold for the highest X percent --> define which genes are tissue specific
# ca 3500 tissue specific genes, 3500/20000 = 0.175 --> take 0.175 as upper percentile

tissuespec_genes_in_rna_tissue_consensus_function <- function(df, upper_percentage){
    threshold <- 1 - upper_percentage
    z_score_threshold <- quantile(df$Z_score, probs = threshold, na.rm = TRUE)
    df$Is_specific <- ifelse(df$Z_score > z_score_threshold, TRUE, FALSE)
    return(df)
}

tissuespec_genes_in_rna_tissue_consensus_function_not_normalized <- function(df, upper_percentage){
    threshold <- 1 - upper_percentage
    nTPM_threshold <- quantile(df$nTPM, probs = threshold, na.rm = TRUE)
    df$Is_specific <- ifelse(df$nTPM > nTPM_threshold, TRUE, FALSE)
    return(df)
}                              

tissue_specificity_test <- function(df_apl, df_tissue_specificity, dimensions, tissue_in_apl, tissue_in_df, number_of_top_genes){
    df_apl_score <- df_apl_score(df_apl, dimensions, tissue_in_apl)
    
    gene_names <- df_apl_score[1:number_of_top_genes, ]$Rowname
    #print(cat("the top genes for tissue: '", tissue_in_apl, "' are : ", "\n", sep=""))
    #print(gene_names)
    list_top_genes <- as.list(gene_names)
    specific_rows <- df_tissue_specificity[df_tissue_specificity$Gene.name %in% list_top_genes & 
                                   df_tissue_specificity$Tissue == tissue_in_df & 
                                   df_tissue_specificity$Is_specific, ]
    
    specific_gene_count <- length(unique(specific_rows$Gene.name))
    gene_percentage <- specific_gene_count/number_of_top_genes*100
    cat("Number of specific genes:", specific_gene_count, "\n")
    cat(gene_percentage, "% of the top ", number_of_top_genes, " genes are specific!", "\n", sep="")
    print(specific_rows)
    return(gene_percentage)
    }

tissue_specificity_test_steps_apl <- function(df_apl, df_tissue_specificity, dimensions, 
                                                tissue_in_apl, tissue_in_df, number_of_top_genes, steps) {
    df_apl_scores <- df_apl_score(df_apl, dimensions, tissue_in_apl)
    result_list <- list()

    for (i in seq(1, number_of_top_genes, by = steps)) {
        gene_names <- df_apl_scores[i:(i + steps-1), ]$Rowname
        specific_rows <- df_tissue_specificity[df_tissue_specificity$Gene.name 
                                               %in% gene_names 
                                               & df_tissue_specificity$Tissue == tissue_in_df 
                                               & df_tissue_specificity$Is_specific, ]
        
        specific_gene_count <- nrow(specific_rows)
        gene_percentage <- specific_gene_count/steps
        result_list[[i]] <- data.frame(Gene_Names = paste(gene_names, collapse = ", "), 
                                       Gene_Percentage = gene_percentage)
        
        #cat("Number of specific genes:", specific_gene_count, "\n")
        #cat(gene_percentage, "% of genes are specific!", "\n", sep="")
        
    }
    result_df <- do.call(rbind, result_list)
    return(result_df)
}

tissue_specificity_test_steps_ordered_list <- function(df_ordered, df_tissue_specificity, 
                                                tissue_in_ordered_list, tissue_in_df, number_of_genes_to_check, steps) {
    df_list <- as.data.frame(t(df))
    ordered_df <- df_list[order(df_list[tissue_in_ordered_list], decreasing=TRUE), ][tissue_in_ordered_list]
    ordered_list <- as.list(row.names(ordered_df))
    result_list <- list()

    for (i in seq(1, number_of_genes_to_check, by = steps)) {
        gene_names <- ordered_list[i:(i + steps-1) ]
        specific_rows <- df_tissue_specificity[df_tissue_specificity$Gene.name 
                                               %in% gene_names 
                                               & df_tissue_specificity$Tissue == tissue_in_df 
                                               & df_tissue_specificity$Is_specific, ]
        specific_gene_count <- nrow(specific_rows)
        gene_percentage <- specific_gene_count/steps
        result_list[[i]] <- data.frame(Gene_Names = paste(gene_names, collapse = ", "), 
                                       Gene_Percentage = gene_percentage)
        
        #cat("Number of specific genes:", specific_gene_count, "\n")
        #cat(gene_percentage, "% of genes are specific!", "\n", sep="")
        
    }
    result_df <- do.call(rbind, result_list)
    return(result_df)
}
                              

# for lists ordered by values (without APL)
tissue_specificity_test_no_apl <- function(list_ordered_by_size, df_tissue_specificity, tissue_in_df, number_of_top_genes){
    list_top_genes <- list_ordered_by_size[1:number_of_top_genes]
    specific_rows <- df_tissue_specificity[df_tissue_specificity$Gene.name %in% list_top_genes & 
                                   df_tissue_specificity$Tissue == tissue_in_df & 
                                   df_tissue_specificity$Is_specific, ]
    
    specific_gene_count <- length(unique(specific_rows$Gene.name))
    gene_percentage <- specific_gene_count/number_of_top_genes*100
    cat("Number of specific genes:", specific_gene_count, "\n")
    cat(gene_percentage, "% of the top ", number_of_top_genes, " genes are specific!", "\n", sep="")
    print(specific_rows)
    
    }

# plot APL-scores & nTPM 
scatter_plot_tissue_specificity <- function(df_for_apl, dimensions, tissue_apl, tissue_zscores, df_zscore){
    df_apl_score <- df_apl_score(df_for_apl, dimensions, tissue_apl)
    filtered_df_zscore <- df_zscore[df_zscore$Tissue == 'cerebellum', ]
    merged_df <- merge(df_apl_score, filtered_df_zscore, by.x = "Rowname", by.y = "Gene.name", 
                   all.x= TRUE)
    plot <- plot(merged_df$Rank, merged_df$Z_score, 
         xlab = "Rank", ylab = "Z_score",
         main = "Scatterplot of Rank vs. Z_score")
    
    fit <- lm(Z_score ~ Rank, data = merged_df)
    abline(fit, col = "red")
    
    return(plot)
    }

hist_plot_tissue_specificity <- function(df_for_apl, dimensions, tissue_apl, tissue_zscores){
    df_apl_score <- df_apl_score(df_for_apl, dimensions, tissue_apl)
    df_zscores <- read.table("/project/pauline_association_plots/code/01_promoter-celltype_matrix/enhancer_analysis/data/rna_tissue_consensus_zscores.tsv", sep="\t", header=TRUE)
    filtered_df_genes_zscore <- df_genes_zscore[df_genes_zscore$Tissue == 'cerebellum', ]
    merged_df <- merge(df_apl_score, filtered_df_genes_zscore, by.x = "Rowname", by.y = "Gene.name", 
                   all.x= TRUE)
    plot <- hist(log(merged_df$Z_score))


    return(plot)
    }


z_normalize_and_shift <- function(matrix) {
  row_means <- rowMeans(matrix)
  row_sds <- apply(matrix, 1, sd)
  normalized_matrix <- t(apply(matrix, 1, function(row) (row - mean(row)) / sd(row)))
  normalized_matrix[!is.finite(normalized_matrix)] <- 0

  shifted_matrix <- normalized_matrix + 1
  
  return(shifted_matrix)
}

has_negative_values <- function(matrix) {
  has_negatives <- any(matrix < 0)
  return(has_negatives)
}