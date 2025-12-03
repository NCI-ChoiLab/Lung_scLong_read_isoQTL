# DEI detection

# Author: Elelta
# Modify: Bolun

# February 19th 2025 - 
# Made corrections to PCT.1 and PCT.2 calculations and saved unfiltered results

# Required Libraries
library(Matrix)
library(DESeq2)
library(writexl)
library(future)
library(future.apply)
library(stringr)
library(BiocParallel)

# Increase R's memory limit (500 GB)
options(java.parameters = "-Xmx1000g") 

# Set up parallel workers
register(MulticoreParam(4))

# Paths to matrix lists
unfiltered_matrix_path <- "/data/Choi_lung/Elelta/BolunTask/isoform_matrices_per_celltype_unfilteredbyisoformsnorindividual.rds"
filtered_matrix_path <- "/data/Choi_lung/Elelta/BolunTask/PostIsoformFiltering/isoform_matrices_filtered_corrected3_1Feb3_2025.rds"
if (!file.exists(unfiltered_matrix_path) || !file.exists(filtered_matrix_path)) 
  stop("Matrix files not found")

# Load unfiltered and filtered matrices
unfiltered_matrices <- readRDS(unfiltered_matrix_path)
filtered_matrices <- readRDS(filtered_matrix_path)

cat(paste("Loaded unfiltered matrix list with", length(unfiltered_matrices), "cell types.\n"))
cat(paste("Loaded filtered matrix list with", length(filtered_matrices), "cell types.\n"))

# Output directory setup
output_dir <- "/data/Choi_lung/Elelta/BolunTask/PostIsoformFiltering/WholeGroupCorrected/Outputfeb19/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir2 <- "/data/Choi_lung/Elelta/BolunTask/PostIsoformFiltering/WholeGroupCorrected/Outputfeb19/unfiltered/"
dir.create(output_dir2, showWarnings = FALSE, recursive = TRUE)

# Process each cell type
for (current_celltype in names(filtered_matrices)) {
  cat("Processing cell type:", current_celltype, "\n")

  tryCatch({
    # Filtered and unfiltered matrices for the current cell type
    filtered_matrix <- filtered_matrices[[current_celltype]]
    if (is.null(filtered_matrix)) stop(paste("Filtered matrix for", current_celltype, "is NULL"))

    # Extract isoforms for the current cell type
    current_isoforms <- rownames(filtered_matrix)
    cat("Number of isoforms in current cell type:", length(current_isoforms), "\n")

    # Combine matrices for the rest of the cell types (excluding current cell type)
    rest_matrix <- do.call(cbind, lapply(setdiff(names(unfiltered_matrices), current_celltype), function(ct) {
      mat <- unfiltered_matrices[[ct]]
      if (!is.null(mat)) {
        mat <- mat[current_isoforms, , drop = FALSE]
        mat[is.na(mat)] <- 0
        return(mat)
      } else {
        return(NULL)
      }
    }))

    if (is.null(rest_matrix) || ncol(rest_matrix) == 0) stop("No columns available in the combined rest matrix")

    cat("Total number of individuals in rest matrix:", ncol(rest_matrix), "\n")
   
    # Remove columns with all-zero values **before** cbind()
    zero_cols <- colSums(rest_matrix) == 0
    if (any(zero_cols)) {
      cat("Removing", sum(zero_cols), "columns with all-zero values from rest_matrix.\n")
      rest_matrix <- rest_matrix[, !zero_cols, drop = FALSE]
    }

    # Combine current and rest matrices
    all_counts <- cbind(filtered_matrix, rest_matrix)
    cat("Final dimension of all_counts matrix:", dim(all_counts), "\n")
    
    # Create colData without incorrect filtering
    colData <- data.frame(
      samples = colnames(all_counts),
      celltype = factor(c(rep("current", ncol(filtered_matrix)), rep("rest", ncol(rest_matrix))), levels = c("rest", "current")),
      stringsAsFactors = FALSE
    )
    
    # Ensure colData aligns with all_counts
    colData <- colData[match(colnames(all_counts), colData$samples), , drop = FALSE]
    
    if (ncol(all_counts) != nrow(colData)) {
      cat("Mismatch detected! all_counts has", ncol(all_counts), "columns but colData has", nrow(colData), "rows.\n")
      cat("All_counts col names:", paste(colnames(all_counts), collapse=", "), "\n")
      cat("ColData samples:", paste(colData$samples, collapse=", "), "\n")
      stop("Mismatch between colData and all_counts dimensions.")
    }
    cat("Final dimension of colData:", dim(colData), "\n")

    # Check for NA values in all_counts
    if (any(is.na(all_counts))) stop("Error: NA values found in all_counts matrix!")

    # Run DESeq2 analysis
    dds <- DESeqDataSetFromMatrix(countData = all_counts, colData = colData, design = ~ celltype)
    dds <- dds[rowSums(counts(dds)) > 0, ]
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds, parallel = TRUE)

    # Calculate PCT.1 and PCT.2 for all genes before filtering
    current_samples <- colData$celltype == "current"
    other_samples <- colData$celltype == "rest"
    pct.1 <- rowSums(all_counts[, current_samples, drop = FALSE] > 0) / sum(current_samples) * 100
    pct.2 <- rowSums(all_counts[, other_samples, drop = FALSE] > 0) / sum(other_samples) * 100

    # Save unfiltered results
    res_unfiltered <- results(dds, contrast = c("celltype", "current", "rest"))
    
    # Assign PCT values to unfiltered results
    res_unfiltered$PCT.1 <- pct.1[rownames(res_unfiltered)]
    res_unfiltered$PCT.2 <- pct.2[rownames(res_unfiltered)]

    saveRDS(res_unfiltered, file = file.path(output_dir2, paste0("differential_isoform_unfiltered_", current_celltype, ".rds")))
    write.csv(as.data.frame(res_unfiltered), file = file.path(output_dir2, paste0("differential_isoform_unfiltered_", current_celltype, ".csv")), row.names = TRUE)

    # Apply filtering
    res_filtered <- res_unfiltered[!is.na(res_unfiltered$padj) & res_unfiltered$padj < 0.05 & (res_unfiltered$log2FoldChange) > 1, ]

    # Assign PCT values to filtered results
    res_filtered$PCT.1 <- pct.1[rownames(res_filtered)]
    res_filtered$PCT.2 <- pct.2[rownames(res_filtered)]

    saveRDS(res_filtered, file = file.path(output_dir, paste0("differential_isoform_res_", current_celltype, "_subgroup.rds")))
    write.csv(res_filtered, file = file.path(output_dir, paste0("differential_isoform_res_", current_celltype, "_subgroup.csv")), row.names = TRUE)

  }, error = function(e) {
    cat("Error processing cell type", current_celltype, ":", e$message, "\n")
  })
}

