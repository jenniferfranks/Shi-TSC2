# ============================================================
# GENERALIZED MONOCLE3 PIPELINE WITH EXTERNAL PYTHON SCRUBLET
# ============================================================

# ============================================================
# 0. USER CONFIGURATION
# ============================================================

project_root <- "/data/hps/assoc/private/franks_lab/user/jfra11/projects/Silica/Shi/TSC2"
setwd(project_root)

# Parent directories containing sample subfolders
parent_dirs <- c(
  "data/TSC2_WS3/outs",
  "data/TSC2_WS4/outs"
)

# Relative path to the h5 file inside each sample folder
target_file_name <- "filtered_feature_bc_matrix.h5"

# Sample map
samplemap_path <- "data/metadata/sample_info.csv"
samplemap_match_col <- "sample_id"
sample_id_col <- "sample_id"
metadata_cols_to_add <- c("model")
cds_sample_col <- "sample"
metadata_numeric_cols <- c( "replicate")

# CDS naming and saving
cds_base_name <- "Shi_TSC2"
cds_save_root <- "data/cds"

# Scrublet configuration
run_scrublet <- TRUE
scrublet_dir <- "scrublet"
scrublet_python <- "/data/hps/assoc/private/franks_lab/user/jfra11/miniforge3/envs/scrublet_env/bin/python"
scrublet_script <- file.path("main/01a_Preprocessing-scrublet.py")
scrublet_expected_doublet_rate <- 0.06
scrublet_score.max <- 0.2

# QC / preprocessing parameters
PCs_to_use <- 100
mitoMax <- 10
cluster.k <- 50
min_genes_expressed <- 100
UMImin <- NULL
UMImax <- NULL

# Save checkpoints
save_after_raw_load <- TRUE
save_after_scrublet <- TRUE
save_after_filtering <- TRUE
save_after_processing <- TRUE

# Seeds
global_seed <- 8248959
cluster_seed <- 3752

# ============================================================
# 1. LIBRARIES
# ============================================================

library(BPCells)
library(monocle3)
library(dplyr)
library(ggplot2)
library(viridis)
library(lme4)
library(gprofiler2)
library(RColorBrewer)
library(splines)
library(monocle3helper)
library(hdf5r)
library(Matrix)

set.seed(global_seed)

# ============================================================
# 2. QUICK SAMPLE MAP PREVIEW
# ============================================================

samplemap <- read.csv(samplemap_path, stringsAsFactors = FALSE)
head(samplemap)
colnames(samplemap)
str(samplemap)

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

make_cds_name <- function(base_name, steps_done = character(0)) {
  if (length(steps_done) == 0) return(base_name)
  paste(c(base_name, steps_done), collapse = "_")
}

save_cds_with_steps <- function(cds, base_name, save_root, steps_done, comment = NULL) {
  out_name <- make_cds_name(base_name, steps_done)
  out_dir <- file.path(save_root, out_name)

  if (!dir.exists(save_root)) {
    dir.create(save_root, recursive = TRUE)
  }

  save_monocle_objects(
    cds,
    directory_path = out_dir,
    comment = comment
  )

  message("Saved CDS to: ", out_dir)
  invisible(out_dir)
}

read_one_h5_as_cds <- function(full_file_path) {
  tryCatch({
    read.cds.cellranger.h5.file(full_file_path)
  }, error = function(e) {
    warning("Error reading ", full_file_path, ": ", e$message)
    return(NULL)
  })
}

get_sample_subdirs <- function(parent_dirs) {
  sub_dirs <- unlist(
    lapply(parent_dirs, function(parent_dir) {
      if (!dir.exists(parent_dir)) {
        warning("Directory does not exist: ", parent_dir)
        return(character(0))
      }

      list.dirs(
        path = parent_dir,
        full.names = TRUE,
        recursive = FALSE
      )
    }),
    use.names = FALSE
  )

  unique(sub_dirs)
}

load_all_sample_cds <- function(parent_dirs, target_file_name) {
  data_list <- list()

  cat("Found", length(parent_dirs), "outs directories to process.\n")
  print(parent_dirs)

  for (folder_path in parent_dirs) {
    if (!dir.exists(folder_path)) {
      warning("Directory does not exist: ", folder_path)
      next
    }

    full_file_path <- file.path(folder_path, target_file_name)
    folder_key <- basename(dirname(folder_path))

    if (file.exists(full_file_path)) {
      message("Processing file for: ", folder_key)

      cds_one <- read_one_h5_as_cds(full_file_path)

      if (!is.null(cds_one)) {
        data_list[[folder_key]] <- cds_one
        message("Successfully loaded and stored data from: ", folder_key)
      }
    } else {
      warning("File not found: ", full_file_path)
    }
  }

  message("Total successful datasets loaded: ", length(data_list))
  print(names(data_list))
  data_list
}

add_sample_metadata <- function(cds,
                                samplemap_path,
                                samplemap_match_col = "sample",
                                sample_id_col = "sample_id",
                                metadata_cols_to_add = NULL,
                                metadata_numeric_cols = NULL,
                                cds_sample_col = "sample") {

  samplemap <- read.csv(samplemap_path, stringsAsFactors = FALSE)
  pd <- pData(cds)

  if (!(sample_id_col %in% colnames(samplemap))) {
    stop("Sample map must contain column: ", sample_id_col)
  }

  if (!(samplemap_match_col %in% colnames(samplemap))) {
    stop("Sample map must contain matching column: ", samplemap_match_col)
  }

  if (!(cds_sample_col %in% colnames(pd))) {
    stop("CDS pData must contain sample column: ", cds_sample_col)
  }

  sample_to_id <- setNames(samplemap[[sample_id_col]], samplemap[[samplemap_match_col]])
  pd[[sample_id_col]] <- sample_to_id[as.character(pd[[cds_sample_col]])]

  if (!is.null(metadata_cols_to_add)) {
    present_cols <- metadata_cols_to_add[metadata_cols_to_add %in% colnames(samplemap)]

    for (col in present_cols) {
      vec <- setNames(samplemap[[col]], samplemap[[sample_id_col]])
      pd[[col]] <- vec[as.character(pd[[sample_id_col]])]

      if (!is.null(metadata_numeric_cols) && col %in% metadata_numeric_cols) {
        pd[[col]] <- suppressWarnings(as.numeric(pd[[col]]))
      }
    }
  }

  pData(cds) <- pd
  cds
}

prepare_scrublet_input <- function(cds,
                                   scrublet_dir,
                                   expected_doublet_rate = 0.06) {

  if (!dir.exists(scrublet_dir)) {
    dir.create(scrublet_dir, recursive = TRUE)
  }

  counts_mat <- t(counts(cds))

  Matrix::writeMM(counts_mat, file.path(scrublet_dir, "counts.mtx"))

  write.csv(
    data.frame(barcode = rownames(counts_mat)),
    file.path(scrublet_dir, "barcodes.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(expected_doublet_rate = expected_doublet_rate),
    file.path(scrublet_dir, "scrublet_params.csv"),
    row.names = FALSE
  )

  message("Wrote Scrublet inputs to: ", scrublet_dir)
}

run_external_scrublet <- function(scrublet_python, scrublet_script, scrublet_dir) {
  if (!file.exists(scrublet_python)) {
    stop("Python executable not found: ", scrublet_python)
  }

  if (!file.exists(scrublet_script)) {
    stop("Scrublet script not found: ", scrublet_script)
  }

  message("Running external Scrublet script...")
  message("Python: ", scrublet_python)
  message("Script: ", scrublet_script)

  result <- system2(
    command = scrublet_python,
    args = c(scrublet_script),
    stdout = TRUE,
    stderr = TRUE
  )

  cat(paste(result, collapse = "\n"), "\n")

  out_file <- file.path(scrublet_dir, "scrublet_results.csv")
  if (!file.exists(out_file)) {
    stop("Scrublet did not produce expected output file: ", out_file)
  }

  invisible(result)
}

add_scrublet_results_to_cds <- function(cds, scrublet_dir) {
  scrublet_results_path <- file.path(scrublet_dir, "scrublet_results.csv")

  if (!file.exists(scrublet_results_path)) {
    stop("Scrublet results file not found: ", scrublet_results_path)
  }

  scrublet_results <- read.csv(scrublet_results_path, stringsAsFactors = FALSE)

  if (!all(c("barcode", "scrublet_score", "scrublet_call") %in% colnames(scrublet_results))) {
    stop("scrublet_results.csv must contain barcode, scrublet_score, and scrublet_call columns")
  }

  pd <- pData(cds)

  score_map <- setNames(scrublet_results$scrublet_score, scrublet_results$barcode)
  call_map  <- setNames(scrublet_results$scrublet_call, scrublet_results$barcode)

  pd$scrublet_score <- as.numeric(score_map[rownames(pd)])
  pd$scrublet_call  <- call_map[rownames(pd)]

  pData(cds) <- pd
  cds
}

apply_basic_qc_filters <- function(cds,
                                   scrublet.score.max = NULL,
                                   min_genes_expressed = NULL,
                                   mitoMax = NULL,
                                   UMImin = NULL,
                                   UMImax = NULL) {

  cds <- cds %>%
    estimate_size_factors() %>%
    detect_genes()

  if (!is.null(scrublet.score.max) && "scrublet_score" %in% colnames(pData(cds))) {
    cds <- cds[, which(pData(cds)$scrublet_score < scrublet.score.max)]
    message("Applied scrublet score filter: < ", scrublet.score.max)
  }

  if (!is.null(UMImin) && "n.umi" %in% colnames(pData(cds))) {
    cds <- cds[, which(pData(cds)$n.umi > UMImin)]
    message("Applied minimum UMI filter: > ", UMImin)
  }

  if (!is.null(UMImax) && "n.umi" %in% colnames(pData(cds))) {
    cds <- cds[, which(pData(cds)$n.umi < UMImax)]
    message("Applied maximum UMI filter: < ", UMImax)
  }

  if (!is.null(min_genes_expressed) && "num_genes_expressed" %in% colnames(pData(cds))) {
    print(summary(pData(cds)$num_genes_expressed))
    cds <- cds[, which(pData(cds)$num_genes_expressed > min_genes_expressed)]
    message("Applied min genes expressed filter: > ", min_genes_expressed)
  }

  if (!is.null(mitoMax)) {
    cds <- calculate_mito(cds)

    if ("perc_mitochondrial_umis" %in% colnames(pData(cds))) {
      print(summary(pData(cds)$perc_mitochondrial_umis))
      cds <- cds[, pData(cds)$perc_mitochondrial_umis < mitoMax]
      message("Applied mitochondrial filter: < ", mitoMax)
    } else {
      warning("perc_mitochondrial_umis not found after calculate_mito()")
    }
  }

  cds
}

process_and_cluster_cds <- function(cds,
                                    PCs_to_use = 100,
                                    cluster.k = 50,
                                    cluster_seed = 3752) {

  cds <- cds %>%
    preprocess_cds(num_dim = PCs_to_use) %>%
    detect_genes() %>%
    reduce_dimension() %>%
    cluster_cells(k = cluster.k, random_seed = cluster_seed)

  pData(cds)$cluster <- clusters(cds)
  cds
}

# ============================================================
# 4. INITIALIZE STEP TRACKER
# ============================================================

steps_done <- character(0)

# ============================================================
# 5. LOAD INDIVIDUAL H5 FILES
# ============================================================

data_list <- load_all_sample_cds(
  parent_dirs = parent_dirs,
  target_file_name = target_file_name
)

if (length(data_list) == 0) {
  stop("No datasets were successfully loaded. Check parent_dirs and target_file_name.")
}

# ============================================================
# 6. COMBINE INTO ONE CDS
# ============================================================

cds <- combine_cds(as.list(data_list))
rm(data_list)

steps_done <- c(steps_done, "raw")

if (save_after_raw_load) {
  save_cds_with_steps(
    cds = cds,
    base_name = cds_base_name,
    save_root = cds_save_root,
    steps_done = steps_done,
    comment = "Read from individual h5 files, no preprocessing performed"
  )
}

print(dim(exprs(cds)))
print(head(pData(cds)))
print(head(fData(cds)))

# ============================================================
# 7. EXTERNAL SCRUBLET
# ============================================================

if (run_scrublet) {
  prepare_scrublet_input(
    cds = cds,
    scrublet_dir = scrublet_dir,
    expected_doublet_rate = scrublet_expected_doublet_rate
  )

  run_external_scrublet(
    scrublet_python = scrublet_python,
    scrublet_script = scrublet_script,
    scrublet_dir = scrublet_dir
  )

  cds <- add_scrublet_results_to_cds(
    cds = cds,
    scrublet_dir = scrublet_dir
  )

  steps_done <- c(steps_done, "scrublet")

  if (save_after_scrublet) {
    save_cds_with_steps(
      cds = cds,
      base_name = cds_base_name,
      save_root = cds_save_root,
      steps_done = steps_done,
      comment = "Scrublet completed externally and imported into pData"
    )
  }
} else {
  message("Skipping scrublet step.")
}

# ============================================================
# 8. ADD SAMPLE METADATA
# ============================================================

cds <- add_sample_metadata(
  cds = cds,
  samplemap_path = samplemap_path,
  samplemap_match_col = samplemap_match_col,
  sample_id_col = sample_id_col,
  metadata_cols_to_add = metadata_cols_to_add,
  metadata_numeric_cols = metadata_numeric_cols,
  cds_sample_col = cds_sample_col
)

steps_done <- c(steps_done, "metadata")

print(colnames(pData(cds)))
print(head(pData(cds)[, intersect(c(sample_id_col, metadata_cols_to_add), colnames(pData(cds))), drop = FALSE]))

# ============================================================
# 9. BASIC QC FILTERING
# ============================================================

cds <- apply_basic_qc_filters(
  cds = cds,
  scrublet.score.max = if (run_scrublet) scrublet_score.max else NULL,
  min_genes_expressed = min_genes_expressed,
  mitoMax = mitoMax,
  UMImin = UMImin,
  UMImax = UMImax
)

steps_done <- c(steps_done, "filtered")

if (save_after_filtering) {
  save_cds_with_steps(
    cds = cds,
    base_name = cds_base_name,
    save_root = cds_save_root,
    steps_done = steps_done,
    comment = "Metadata added and QC filters applied"
  )
}

# ============================================================
# 10. PREPROCESS AND CLUSTER
# ============================================================

cds <- process_and_cluster_cds(
  cds = cds,
  PCs_to_use = PCs_to_use,
  cluster.k = cluster.k,
  cluster_seed = cluster_seed
)

steps_done <- c(steps_done, "processed", "clustered")

if (save_after_processing) {
  save_cds_with_steps(
    cds = cds,
    base_name = cds_base_name,
    save_root = cds_save_root,
    steps_done = steps_done,
    comment = "Filtered, processed, and clustered"
  )
}

# ============================================================
# 11. BASIC PLOTS
# ============================================================
steps_done <- c("raw", "scrublet", "metadata", "filtered", "processed", "clustered")
cds <- load_monocle_objects(file.path(cds_save_root, make_cds_name(cds_base_name, steps_done)))

plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE)

if ("exposure" %in% colnames(pData(cds))) {
  plot_cells(cds, color_cells_by = "exposure")
}

if ("treatment" %in% colnames(pData(cds))) {
  plot_cells(cds, color_cells_by = "treatment")
}

plot_cells(cds, genes = c("Cd68", "Mrc1", "Itgam", "Siglecf"), scale_to_range = FALSE)
plot_cells(cds, genes = c("Cd68", "Cd4", "Itgam", "Siglecf"), scale_to_range = FALSE)
plot_cells(cds, genes = c("Ctsk", "Spp1", "Acp5"), scale_to_range = FALSE)

# ============================================================
# 12. FINAL SUMMARY
# ============================================================

message("Final CDS object name: ", make_cds_name(cds_base_name, steps_done))
message("Final number of cells: ", ncol(cds))
message("Final number of genes: ", nrow(cds))

print(table(pData(cds)$cluster))