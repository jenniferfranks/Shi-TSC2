# ------------------------------------------------------------------------------
# Mixed-effect DEG analysis: GFP+ vs GFP- within each cell type using monocle3
# ------------------------------------------------------------------------------

library(monocle3)
library(dplyr)
library(tibble)
library(stringr)
library(Matrix)

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
main_dir <- "/data/hps/assoc/private/franks_lab/user/jfra11/projects/Silica/Shi/TSC2/"
setwd(main_dir)

cds.stromal <- load_monocle_objects(
  directory_path = "data/cds/Shi_TSC2_stromal-cells"
)

celltype_col <- "subset_cluster_names"
group_col    <- "GFP_status"
sample_col   <- "sample"

out_dir <- file.path(main_dir, "03-GFP_stromalDEG_mixed_models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"),   showWarnings = FALSE)

min_cells_per_group <- 20
min_samples_per_group <- 2
expression_family_to_use <- "mixed-negbinomial"

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

safe_name <- function(x) {
  x <- gsub("[/\\:*?\"<>|[:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

log_message <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
}

flatten_list_columns <- function(df) {
  for (i in seq_along(df)) {
    if (is.list(df[[i]])) {
      df[[i]] <- vapply(
        df[[i]],
        function(x) {
          if (length(x) == 0 || is.null(x)) {
            NA_character_
          } else {
            paste(as.character(unlist(x)), collapse = "; ")
          }
        },
        character(1)
      )
    }
  }
  df
}

write_status_table <- function(df, file) {
  utils::write.csv(df, file = file, row.names = FALSE, quote = TRUE)
}

# ------------------------------------------------------------------------------
# checks
# ------------------------------------------------------------------------------

pd <- as.data.frame(pData(cds.stromal))

stopifnot(celltype_col %in% colnames(pd))
stopifnot(group_col %in% colnames(pd))
stopifnot(sample_col %in% colnames(pd))

# keep only cells with all needed metadata
keep_cells <- !is.na(pd[[celltype_col]]) &
              !is.na(pd[[group_col]]) &
              !is.na(pd[[sample_col]])

cds_use <- cds.stromal[, keep_cells]
pd_use <- as.data.frame(pData(cds_use))

# standardize GFP labels if needed
pd_use[[group_col]] <- as.character(pd_use[[group_col]])
pd_use[[group_col]][pd_use[[group_col]] %in% c("positive", "pos", "GFP positive")] <- "GFP+"
pd_use[[group_col]][pd_use[[group_col]] %in% c("negative", "neg", "GFP negative")] <- "GFP-"

keep_gfp <- pd_use[[group_col]] %in% c("GFP+", "GFP-")
cds_use <- cds_use[, keep_gfp]
pd_use <- pd_use[keep_gfp, , drop = FALSE]

# set factor levels so GFP- is reference
pData(cds_use)[[group_col]] <- factor(pData(cds_use)[[group_col]], levels = c("GFP-", "GFP+"))
pData(cds_use)[[sample_col]] <- as.factor(pData(cds_use)[[sample_col]])
pData(cds_use)[[celltype_col]] <- as.character(pData(cds_use)[[celltype_col]])

celltypes <- sort(unique(pData(cds_use)[[celltype_col]]))

summary_list <- list()

# ------------------------------------------------------------------------------
# loop over cell types
# ------------------------------------------------------------------------------

for (ct in celltypes) {

  ct_safe <- safe_name(ct)
  log_file <- file.path(out_dir, "logs", paste0(ct_safe, "_log.txt"))
  zz <- file(log_file, open = "wt")
  sink(zz, type = "output")
  sink(zz, type = "message")

  tryCatch({

    log_message("====================================================")
    log_message("Starting cell type: ", ct)

    cell_idx <- which(pData(cds_use)[[celltype_col]] == ct)
    cds_ct <- cds_use[, cell_idx]
    pd_ct <- as.data.frame(pData(cds_ct))

    grp_tab <- table(pd_ct[[group_col]])

    n_gfp_neg <- if ("GFP-" %in% names(grp_tab)) unname(grp_tab["GFP-"]) else 0
    n_gfp_pos <- if ("GFP+" %in% names(grp_tab)) unname(grp_tab["GFP+"]) else 0

    n_samples_gfp_neg <- length(unique(pd_ct[[sample_col]][pd_ct[[group_col]] == "GFP-"]))
    n_samples_gfp_pos <- length(unique(pd_ct[[sample_col]][pd_ct[[group_col]] == "GFP+"]))

    log_message("Cells in GFP-: ", n_gfp_neg)
    log_message("Cells in GFP+: ", n_gfp_pos)
    log_message("Samples in GFP-: ", n_samples_gfp_neg)
    log_message("Samples in GFP+: ", n_samples_gfp_pos)

    if (n_gfp_neg < min_cells_per_group || n_gfp_pos < min_cells_per_group) {
      log_message("Skipping: too few cells")
      summary_list[[ct]] <- data.frame(
        cell_type = ct,
        n_cells = ncol(cds_ct),
        n_gfp_neg = n_gfp_neg,
        n_gfp_pos = n_gfp_pos,
        n_samples_gfp_neg = n_samples_gfp_neg,
        n_samples_gfp_pos = n_samples_gfp_pos,
        status = "skipped_too_few_cells",
        n_deg_fdr_0.05 = NA_integer_,
        stringsAsFactors = FALSE
      )
      return(NULL)
    }

    if (n_samples_gfp_neg < min_samples_per_group || n_samples_gfp_pos < min_samples_per_group) {
      log_message("Skipping: too few samples per GFP group")
      summary_list[[ct]] <- data.frame(
        cell_type = ct,
        n_cells = ncol(cds_ct),
        n_gfp_neg = n_gfp_neg,
        n_gfp_pos = n_gfp_pos,
        n_samples_gfp_neg = n_samples_gfp_neg,
        n_samples_gfp_pos = n_samples_gfp_pos,
        status = "skipped_too_few_samples",
        n_deg_fdr_0.05 = NA_integer_,
        stringsAsFactors = FALSE
      )
      return(NULL)
    }

    # optional gene detection filter
    gene_detected <- Matrix::rowSums(exprs(cds_ct) > 0) >= 10
    cds_ct <- cds_ct[gene_detected, ]

    log_message("Genes retained after filtering: ", nrow(cds_ct))

    if (nrow(cds_ct) == 0) {
      log_message("Skipping: no genes after filtering")
      summary_list[[ct]] <- data.frame(
        cell_type = ct,
        n_cells = ncol(cds_ct),
        n_gfp_neg = n_gfp_neg,
        n_gfp_pos = n_gfp_pos,
        n_samples_gfp_neg = n_samples_gfp_neg,
        n_samples_gfp_pos = n_samples_gfp_pos,
        status = "skipped_no_genes",
        n_deg_fdr_0.05 = 0,
        stringsAsFactors = FALSE
      )
      return(NULL)
    }

    log_message("Running fit_models()")

    gene_fits <- monocle3::fit_models(
      cds_ct,
      model_formula_str = paste0("~ ", group_col, " + (1|", sample_col, ")"),
      expression_family = expression_family_to_use,
      clean_model = TRUE,
      cores = 1
    )

    gene_coefs <- monocle3::coefficient_table(gene_fits)

    # --------------------------------------------------------------------------
    # SAFE EXPORT OF FULL COEFFICIENT TABLE
    # --------------------------------------------------------------------------

    tryCatch({
      gene_coefs_export <- as.data.frame(gene_coefs)

      # drop heavy / non-rectangular S4-ish columns before export
      drop_cols <- intersect(
        c("model", "model_summary"),
        colnames(gene_coefs_export)
      )

      if (length(drop_cols) > 0) {
        gene_coefs_export <- gene_coefs_export[, setdiff(colnames(gene_coefs_export), drop_cols), drop = FALSE]
      }

      gene_coefs_export <- flatten_list_columns(gene_coefs_export)

      write_status_table(
        gene_coefs_export,
        file.path(out_dir, "tables", paste0(ct_safe, "_coefficient_table_full.csv"))
      )
      log_message("Wrote full coefficient table")
    }, error = function(e) {
      log_message("Skipping full coefficient table export: ", conditionMessage(e))
    })

    # --------------------------------------------------------------------------
    # USE A PLAIN DATA FRAME FOR DOWNSTREAM FILTERING
    # --------------------------------------------------------------------------

    coef_tab <- tryCatch({
      as.data.frame(gene_coefs)
    }, error = function(e) {
      log_message("Could not coerce coefficient table to data.frame: ", conditionMessage(e))
      NULL
    })

    if (is.null(coef_tab) || nrow(coef_tab) == 0) {
      stop("coefficient_table could not be converted to a usable data.frame")
    }

    log_message("Coefficient table columns: ", paste(colnames(coef_tab), collapse = ", "))

    if (!"term" %in% colnames(coef_tab)) {
      stop("coefficient_table does not contain a 'term' column")
    }

    matching_terms <- unique(coef_tab$term[grepl(group_col, coef_tab$term, fixed = TRUE)])
    log_message("Matching terms for ", group_col, ": ", paste(matching_terms, collapse = ", "))

    coef_tab_sub <- coef_tab %>%
      filter(grepl(group_col, term, fixed = TRUE)) %>%
      mutate(cell_type = ct)

    if (nrow(coef_tab_sub) == 0) {
      stop(paste0("No rows found in coefficient_table for term containing ", group_col))
    }

    # determine the gene identifier column monocle returned
    possible_gene_cols <- c("gene", "gene_id", "id", "feature_id", "gene_short_name")
    gene_col_found <- possible_gene_cols[possible_gene_cols %in% colnames(coef_tab_sub)][1]

    if (is.na(gene_col_found)) {
      stop("Could not find a gene identifier column in coefficient_table")
    }

    coef_tab_sub$gene_id <- as.character(coef_tab_sub[[gene_col_found]])

    # attach gene_short_name if present
    gene_annot <- as.data.frame(rowData(cds_ct))
    gene_annot$gene_id <- rownames(gene_annot)

    if ("gene_short_name" %in% colnames(gene_annot)) {
      coef_tab_sub <- coef_tab_sub %>%
        left_join(
          gene_annot[, c("gene_id", "gene_short_name")],
          by = "gene_id"
        )

      # if monocle returned gene_short_name as the gene_id column, preserve it
      missing_short <- is.na(coef_tab_sub$gene_short_name) | coef_tab_sub$gene_short_name == ""
      coef_tab_sub$gene_short_name[missing_short] <- coef_tab_sub$gene_id[missing_short]
    } else {
      coef_tab_sub$gene_short_name <- coef_tab_sub$gene_id
    }

    # rename only columns that exist
    if ("estimate" %in% colnames(coef_tab_sub)) {
      coef_tab_sub <- coef_tab_sub %>% rename(beta = estimate)
    }

    if ("std_err" %in% colnames(coef_tab_sub)) {
      coef_tab_sub <- coef_tab_sub %>% rename(se = std_err)
    }

    if ("test_val" %in% colnames(coef_tab_sub)) {
      coef_tab_sub <- coef_tab_sub %>% rename(test_stat = test_val)
    }

    # arrange safely
    if ("q_value" %in% colnames(coef_tab_sub) && "beta" %in% colnames(coef_tab_sub)) {
      coef_tab_sub <- coef_tab_sub %>%
        arrange(q_value, desc(abs(beta)))
    } else if ("p_value" %in% colnames(coef_tab_sub) && "beta" %in% colnames(coef_tab_sub)) {
      coef_tab_sub <- coef_tab_sub %>%
        arrange(p_value, desc(abs(beta)))
    }

    coef_tab_sub_export <- flatten_list_columns(as.data.frame(coef_tab_sub))

    out_file <- file.path(out_dir, "tables", paste0(ct_safe, "_GFPpos_vs_GFPneg_mixed_model.csv"))
    write_status_table(coef_tab_sub_export, out_file)

    log_message("Wrote results: ", out_file)

    n_sig <- if ("q_value" %in% colnames(coef_tab_sub)) {
      sum(coef_tab_sub$q_value < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    log_message("Significant genes at q < 0.05: ", n_sig)

    summary_list[[ct]] <- data.frame(
      cell_type = ct,
      n_cells = ncol(cds_ct),
      n_gfp_neg = n_gfp_neg,
      n_gfp_pos = n_gfp_pos,
      n_samples_gfp_neg = n_samples_gfp_neg,
      n_samples_gfp_pos = n_samples_gfp_pos,
      status = "ok",
      n_deg_fdr_0.05 = n_sig,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {

    log_message("ERROR: ", conditionMessage(e))

    summary_list[[ct]] <- data.frame(
      cell_type = ct,
      n_cells = sum(pData(cds_use)[[celltype_col]] == ct, na.rm = TRUE),
      n_gfp_neg = NA_integer_,
      n_gfp_pos = NA_integer_,
      n_samples_gfp_neg = NA_integer_,
      n_samples_gfp_pos = NA_integer_,
      status = paste0("error: ", conditionMessage(e)),
      n_deg_fdr_0.05 = NA_integer_,
      stringsAsFactors = FALSE
    )

  }, finally = {
    sink(type = "message")
    sink(type = "output")
    close(zz)
  })
}

# ------------------------------------------------------------------------------
# combine summary
# ------------------------------------------------------------------------------

summary_df <- bind_rows(summary_list) %>%
  arrange(cell_type)

write.csv(
  summary_df,
  file.path(out_dir, "GFP_mixed_model_summary.csv"),
  row.names = FALSE
)

print(summary_df)