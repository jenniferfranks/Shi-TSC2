## ------------------
# GFP pseudobulk DEG analyses
#
# Jennifer Franks style: monocle + hooke pseudobulk framework
#
# Purpose:
# Run paired pseudobulk DEG testing for GFP+ vs GFP- within each stromal cell type
# using sample_id as the pairing variable and GFP_status as the comparison.
## ------------------

suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(hooke)
  library(gprofiler2)
  library(grid)
})

# ------------------------------------------------------------------------------
# run metadata
# ------------------------------------------------------------------------------
run_date <- format(Sys.Date(), "%Y%m%d")
run_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------
main_dir <- "/data/hps/assoc/private/franks_lab/user/jfra11/projects/Silica/Shi/TSC2/"
deg_root <- file.path(main_dir, "03-GFP_stromalDEG_pseudobulk_hooke")
dir.create(deg_root, recursive = TRUE, showWarnings = FALSE)

setwd(main_dir)

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------
safe_name <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_+.-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

append_log <- function(log_file, ..., sep = "") {
  msg <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste(..., collapse = sep))
  cat(msg, "\n", file = log_file, append = TRUE)
  message(msg)
}

write_status_table <- function(df, outfile) {
  utils::write.csv(df, outfile, row.names = FALSE)
}

flatten_list_columns <- function(df) {
  list_cols <- sapply(df, is.list)
  if (any(list_cols)) {
    df[list_cols] <- lapply(df[list_cols], function(col) {
      vapply(col, function(x) {
        if (length(x) == 0 || all(is.na(x))) {
          NA_character_
        } else {
          paste(as.character(x), collapse = "; ")
        }
      }, character(1))
    })
  }
  df
}

plot_volcano_p_with_q <- function(
  coef_table,
  term_keep = "GFP_for_modelGFP+",
  q_cutoff_1 = 0.01,
  q_cutoff_2 = 0.001,
  label_n_each_side = 8,
  title = NULL
) {
  df <- coef_table %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(
      estimate = as.numeric(estimate),
      p_value = as.numeric(p_value),
      q_value = as.numeric(q_value)
    ) %>%
    filter(!is.na(estimate), !is.na(p_value), !is.na(q_value)) %>%
    filter(term == term_keep) %>%
    mutate(
      neglog10_p = -log10(pmax(p_value, 1e-300)),
      direction_class = case_when(
        estimate > 0 & q_value < q_cutoff_2 ~ "Higher GFP+ (q<0.001)",
        estimate > 0 & q_value < q_cutoff_1 ~ "Higher GFP+ (q<0.01)",
        estimate < 0 & q_value < q_cutoff_2 ~ "Higher GFP- (q<0.001)",
        estimate < 0 & q_value < q_cutoff_1 ~ "Higher GFP- (q<0.01)",
        TRUE ~ "NS"
      )
    )

  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste0(title, "\nNo data"))
    )
  }

  df$direction_class <- factor(
    df$direction_class,
    levels = c(
      "Higher GFP+ (q<0.001)",
      "Higher GFP+ (q<0.01)",
      "Higher GFP- (q<0.001)",
      "Higher GFP- (q<0.01)",
      "NS"
    )
  )

  label_pos <- df %>%
    filter(estimate > 0) %>%
    arrange(q_value, desc(abs(estimate))) %>%
    distinct(gene_short_name, .keep_all = TRUE) %>%
    slice_head(n = label_n_each_side)

  label_neg <- df %>%
    filter(estimate < 0) %>%
    arrange(q_value, desc(abs(estimate))) %>%
    distinct(gene_short_name, .keep_all = TRUE) %>%
    slice_head(n = label_n_each_side)

  label_df <- bind_rows(label_pos, label_neg)

  xmax <- max(abs(df$estimate), na.rm = TRUE) * 1.05
  if (!is.finite(xmax) || xmax <= 0) xmax <- 1

  ggplot(df, aes(x = estimate, y = neglog10_p)) +
    geom_point(aes(color = direction_class), alpha = 0.8, size = 2) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.01), linetype = "dotted", linewidth = 0.4) +
    geom_text_repel(
      data = label_df,
      aes(label = gene_short_name),
      size = 3,
      max.overlaps = Inf
    ) +
    scale_x_continuous(
      limits = c(-xmax, xmax),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(
      values = c(
        "Higher GFP+ (q<0.001)" = "#1B7837",
        "Higher GFP+ (q<0.01)" = "#A6DBA0",
        "Higher GFP- (q<0.001)" = "#B2182B",
        "Higher GFP- (q<0.01)" = "#F4A6A6",
        "NS" = "grey75"
      )
    ) +
    labs(
      title = title,
      x = "Model coefficient",
      y = expression(-log[10](p_value)),
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 9),
      legend.key.width = unit(1.2, "lines"),
      plot.margin = margin(5.5, 20, 5.5, 20)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
}

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------
celltype_col <- "subset_cluster_names"
group_col <- "GFP_status"
sample_col <- "sample_id"

min_cells_total_per_group <- 20
min_cells_for_gene_fraction <- 0.05
min_pseudobulk_groups <- 4
sig_q_cutoff <- 0.01

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
cds <- load_monocle_objects(
  directory_path = file.path(main_dir, "data", "cds", "Shi_TSC2_stromal-cells")
)

pd <- as.data.frame(pData(cds))

if (!celltype_col %in% colnames(pd)) stop(paste(celltype_col, "missing"))
if (!group_col %in% colnames(pd)) stop(paste(group_col, "missing"))
if (!sample_col %in% colnames(pd)) stop(paste(sample_col, "missing"))

# ------------------------------------------------------------------------------
# metadata cleanup
# ------------------------------------------------------------------------------
pd[[group_col]] <- as.character(pd[[group_col]])
pd[[group_col]][pd[[group_col]] %in% c("positive", "pos", "GFP positive")] <- "GFP+"
pd[[group_col]][pd[[group_col]] %in% c("negative", "neg", "GFP negative")] <- "GFP-"

keep_cells <- !is.na(pd[[celltype_col]]) &
              !is.na(pd[[group_col]]) &
              !is.na(pd[[sample_col]]) &
              pd[[group_col]] %in% c("GFP+", "GFP-")

cds <- cds[, keep_cells]
pd <- as.data.frame(pData(cds))

pData(cds)[[group_col]] <- factor(pd[[group_col]], levels = c("GFP-", "GFP+"))
pData(cds)[[sample_col]] <- as.character(pd[[sample_col]])
pData(cds)[[celltype_col]] <- as.character(pd[[celltype_col]])

cell.types <- unique(as.character(pData(cds)[[celltype_col]]))
cell.types <- cell.types[!is.na(cell.types)]
cell.types <- sort(cell.types)

manifest <- vector("list", length(cell.types))

# ------------------------------------------------------------------------------
# main loop
# ------------------------------------------------------------------------------
for (i in seq_along(cell.types)) {

  cell.type <- cell.types[i]
  cell.type.safe <- safe_name(cell.type)

  cell_dir <- file.path(deg_root, cell.type.safe)
  plot_dir <- file.path(cell_dir, "plots")
  result_dir <- file.path(cell_dir, "results")
  object_dir <- file.path(cell_dir, "objects")

  dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(cell_dir, paste0("log_", cell.type.safe, ".txt"))
  cat("", file = log_file)

  append_log(log_file, "====================================================")
  append_log(log_file, "Starting GFP pseudobulk analysis for cell type: ", cell.type)
  append_log(log_file, "Run ID: ", run_id)
  append_log(log_file, "Run started at: ", run_timestamp)

  status_row <- data.frame(
    run_id = run_id,
    run_date = run_date,
    run_timestamp = run_timestamp,
    cell_type = cell.type,
    safe_name = cell.type.safe,
    n_cells_initial = NA_integer_,
    n_cells_gfp_neg = NA_integer_,
    n_cells_gfp_pos = NA_integer_,
    n_samples_total = NA_integer_,
    n_samples_gfp_neg = NA_integer_,
    n_samples_gfp_pos = NA_integer_,
    n_pseudobulk_columns = NA_integer_,
    n_genes_tested = NA_integer_,
    n_sig_q01 = NA_integer_,
    status = "started",
    notes = NA_character_,
    stringsAsFactors = FALSE
  )

  tryCatch({

    keep_cells_ct <- which(pData(cds)[[celltype_col]] == cell.type)
    status_row$n_cells_initial <- length(keep_cells_ct)

    if (length(keep_cells_ct) == 0) stop("No cells found for cell type")

    cds.subset <- cds[, keep_cells_ct]
    pd.subset <- as.data.frame(pData(cds.subset))

    grp_tab <- table(pd.subset[[group_col]])
    samp_tab <- table(pd.subset[[sample_col]], pd.subset[[group_col]])

    status_row$n_cells_gfp_neg <- if ("GFP-" %in% names(grp_tab)) unname(grp_tab["GFP-"]) else 0
    status_row$n_cells_gfp_pos <- if ("GFP+" %in% names(grp_tab)) unname(grp_tab["GFP+"]) else 0
    status_row$n_samples_total <- length(unique(pd.subset[[sample_col]]))
    status_row$n_samples_gfp_neg <- if ("GFP-" %in% colnames(samp_tab)) sum(samp_tab[, "GFP-"] > 0) else 0
    status_row$n_samples_gfp_pos <- if ("GFP+" %in% colnames(samp_tab)) sum(samp_tab[, "GFP+"] > 0) else 0

    append_log(log_file, "Cells in GFP-: ", status_row$n_cells_gfp_neg)
    append_log(log_file, "Cells in GFP+: ", status_row$n_cells_gfp_pos)
    append_log(log_file, "Sample by GFP table:")
    capture.output(print(samp_tab), file = log_file, append = TRUE)

    if (status_row$n_cells_gfp_neg < min_cells_total_per_group ||
        status_row$n_cells_gfp_pos < min_cells_total_per_group) {
      stop("Too few cells in one or both GFP groups")
    }

    cds.subset <- cds.subset %>%
      preprocess_cds() %>%
      estimate_size_factors() %>%
      detect_genes()

    gene_keep <- which(
      fData(cds.subset)$num_cells_expressed >
        min_cells_for_gene_fraction * nrow(pData(cds.subset))
    )
    cds.subset <- cds.subset[gene_keep, ]

    status_row$n_genes_tested <- nrow(cds.subset)
    append_log(log_file, "Genes retained after filtering: ", status_row$n_genes_tested)

    if (nrow(cds.subset) == 0) stop("No genes passed filtering")

    pData(cds.subset)$sample_id.factor <- factor(as.character(pData(cds.subset)[[sample_col]]))
    pData(cds.subset)$GFP_status <- factor(as.character(pData(cds.subset)[[group_col]]), levels = c("GFP-", "GFP+"))

    ccs <- new_cell_count_set(
      cds.subset,
      sample_group = sample_col,
      cell_group = group_col
    )

    pseudobulk_cds <- hooke:::pseudobulk_ccs_for_states(ccs)
    pb_cd <- as.data.frame(colData(pseudobulk_cds))

    append_log(log_file, "Pseudobulk colData columns:")
    append_log(log_file, paste(colnames(pb_cd), collapse = ", "))

    pb_sample_col <- c("sample_id", "sample", "sample_group", "ID")[
      c("sample_id", "sample", "sample_group", "ID") %in% colnames(pb_cd)
    ][1]

    pb_group_col <- c("GFP_status", "cell_group", group_col)[
      c("GFP_status", "cell_group", group_col) %in% colnames(pb_cd)
    ][1]

    if (is.na(pb_sample_col)) {
      stop("Could not find a sample column in pseudobulk_cds colData")
    }
    if (is.na(pb_group_col)) {
      stop("Could not find a GFP/group column in pseudobulk_cds colData")
    }

    append_log(log_file, "Using pseudobulk sample column: ", pb_sample_col)
    append_log(log_file, "Using pseudobulk group column: ", pb_group_col)

    colData(pseudobulk_cds)$sample_for_model <- factor(as.character(pb_cd[[pb_sample_col]]))
    colData(pseudobulk_cds)$GFP_for_model <- factor(as.character(pb_cd[[pb_group_col]]), levels = c("GFP-", "GFP+"))

    append_log(log_file, "Pseudobulk sample by GFP table:")
    capture.output(
      print(table(colData(pseudobulk_cds)$sample_for_model, colData(pseudobulk_cds)$GFP_for_model)),
      file = log_file,
      append = TRUE
    )

    status_row$n_pseudobulk_columns <- ncol(pseudobulk_cds)
    if (status_row$n_pseudobulk_columns < min_pseudobulk_groups) {
      stop("Too few pseudobulk columns")
    }

    pseudo_fit <- fit_models(
      pseudobulk_cds,
      model_formula_str = "~ sample_for_model + GFP_for_model",
      expression_family = "negbinomial",
      weights = colData(pseudobulk_cds)$num_cells_in_group,
      clean_model = TRUE
    )

    pseudo_fit_red <- fit_models(
      pseudobulk_cds,
      model_formula_str = "~ sample_for_model",
      expression_family = "negbinomial",
      weights = colData(pseudobulk_cds)$num_cells_in_group,
      clean_model = TRUE
    )

    compared <- compare_models(pseudo_fit, pseudo_fit_red) %>%
      dplyr::select(gene_short_name, q_value) %>%
      distinct(gene_short_name, .keep_all = TRUE)

    status_row$n_sig_q01 <- sum(compared$q_value < sig_q_cutoff, na.rm = TRUE)

    append_log(log_file, "Genes q<0.01: ", status_row$n_sig_q01)

    write_status_table(
      compared,
      file.path(result_dir, "compared_models.csv")
    )

    # --------------------------------------------------------------------------
    # coefficient table
    # --------------------------------------------------------------------------
    gene_coefs <- coefficient_table(pseudo_fit)

    coef_tab <- as.data.frame(gene_coefs)
    coef_tab <- flatten_list_columns(coef_tab)
    coef_tab$cell_type <- cell.type

    if (!"gene_short_name" %in% colnames(coef_tab)) {
      coef_tab$gene_short_name <- coef_tab$id
    }

    write_status_table(
      coef_tab,
      file.path(result_dir, "coefficient_table.csv")
    )

    coef_tab_gfp <- coef_tab %>%
      filter(term == "GFP_for_modelGFP+")

    coef_tab_gfp <- coef_tab_gfp %>%
      left_join(
        compared %>% rename(model_compare_q_value = q_value),
        by = "gene_short_name"
      )

    write_status_table(
      coef_tab_gfp,
      file.path(result_dir, "GFP_for_modelGFPpos_coefficient_table.csv")
    )

    # --------------------------------------------------------------------------
    # significant genes + g:Profiler
    # --------------------------------------------------------------------------
    sig_genes <- compared %>%
      filter(!is.na(q_value)) %>%
      filter(q_value < sig_q_cutoff) %>%
      distinct(gene_short_name, .keep_all = TRUE) %>%
      arrange(q_value)

    write_status_table(
      sig_genes,
      file.path(result_dir, "significant_genes_q_lt_0.01.csv")
    )

    append_log(log_file, "Significant genes saved: ", nrow(sig_genes))

    sig_gfp_coef <- coef_tab_gfp %>%
      filter(gene_short_name %in% sig_genes$gene_short_name) %>%
      arrange(model_compare_q_value, p_value, desc(abs(as.numeric(estimate))))

    write_status_table(
      sig_gfp_coef,
      file.path(result_dir, "significant_genes_GFP_q_lt_0.01.csv")
    )

    append_log(log_file, "Significant GFP coefficient rows saved: ", nrow(sig_gfp_coef))

    if (nrow(sig_genes) > 0) {

      gp_query <- sig_genes$gene_short_name
      gp_query <- gp_query[!is.na(gp_query) & gp_query != ""]
      gp_query <- unique(gp_query)

      append_log(log_file, "Running g:Profiler on ", length(gp_query), " significant genes")

      gost_res <- tryCatch({
        gost(
          query = gp_query,
          organism = "mmusculus"
        )
      }, error = function(e) {
        append_log(log_file, "g:Profiler failed: ", conditionMessage(e))
        NULL
      })

      if (!is.null(gost_res) && !is.null(gost_res$result) && nrow(gost_res$result) > 0) {
        gp_df <- as.data.frame(gost_res$result, stringsAsFactors = FALSE)
        gp_df <- flatten_list_columns(gp_df)

        write_status_table(
          gp_df,
          file.path(result_dir, "gprofiler_significant_genes_q_lt_0.01.csv")
        )

        append_log(log_file, "g:Profiler results saved: ", nrow(gp_df), " terms")
      } else {
        append_log(log_file, "g:Profiler returned no enriched terms")
      }

    } else {
      append_log(log_file, "No significant genes at q < 0.01, skipping g:Profiler")
    }

    p_volcano <- plot_volcano_p_with_q(
      coef_table = coef_tab,
      term_keep = "GFP_for_modelGFP+",
      q_cutoff_1 = 0.01,
      q_cutoff_2 = 0.001,
      label_n_each_side = 8,
      title = paste0(cell.type, ": pseudobulk GFP+ vs GFP-")
    )

    ggsave(
      file.path(plot_dir, "volcano_GFPpos_vs_GFPneg.png"),
      p_volcano,
      width = 7,
      height = 5,
      dpi = 300
    )

    ggsave(
      file.path(plot_dir, "volcano_GFPpos_vs_GFPneg.pdf"),
      p_volcano,
      width = 7,
      height = 5
    )

    saveRDS(
      list(
        run_id = run_id,
        cell_type = cell.type,
        pseudo_fit = pseudo_fit,
        pseudo_fit_red = pseudo_fit_red,
        compared = compared,
        gene_coefs = gene_coefs,
        coef_tab = coef_tab,
        coef_tab_gfp = coef_tab_gfp,
        pseudobulk_cds = pseudobulk_cds
      ),
      file = file.path(object_dir, "gfp_pseudobulk_fit_objects.rds")
    )

    status_row$status <- "success"

  }, error = function(e) {
    append_log(log_file, "ERROR: ", conditionMessage(e))
    status_row$status <- "failed"
    status_row$notes <- conditionMessage(e)
  })

  manifest[[i]] <- status_row
}

manifest_df <- dplyr::bind_rows(manifest)

write.csv(
  manifest_df,
  file.path(deg_root, paste0("GFP_pseudobulk_manifest_", run_id, ".csv")),
  row.names = FALSE
)

write.csv(
  manifest_df,
  file.path(deg_root, "GFP_pseudobulk_manifest_latest.csv"),
  row.names = FALSE
)

print(manifest_df)

# ------------------------------------------------------------------------------
# sanity check plot for one gene within one cell type
# ------------------------------------------------------------------------------
plot_gene_paired_celltype <- function(
  cds,
  gene,
  celltype,
  celltype_col = "subset_cluster_names",
  sample_col = "sample_id",
  group_col = "GFP_status"
) {

  pd <- as.data.frame(pData(cds))
  fd <- as.data.frame(fData(cds))

  keep_cells <- which(pd[[celltype_col]] == celltype)

  if (length(keep_cells) == 0) {
    stop(paste0("No cells found for cell type: ", celltype))
  }

  cds_sub <- cds[, keep_cells]
  pd_sub <- as.data.frame(pData(cds_sub))

  gene_idx <- which(rownames(fd) == gene | fd$gene_short_name == gene)
  gene_idx <- unique(gene_idx)

  if (length(gene_idx) == 0) {
    stop(paste0("Could not find gene: ", gene))
  }
  if (length(gene_idx) > 1) {
    gene_idx <- gene_idx[1]
    message("Multiple matches found, using first.")
  }

  gene_name <- if ("gene_short_name" %in% colnames(fd)) {
    fd$gene_short_name[gene_idx]
  } else {
    rownames(fd)[gene_idx]
  }

  expr_vec <- as.numeric(exprs(cds_sub)[gene_idx, ])

  df <- pd_sub %>%
    mutate(
      expr = expr_vec,
      sample_plot = .data[[sample_col]],
      group_plot = .data[[group_col]]
    ) %>%
    filter(!is.na(sample_plot), !is.na(group_plot)) %>%
    filter(group_plot %in% c("GFP-", "GFP+")) %>%
    group_by(sample_plot, group_plot) %>%
    summarise(
      mean_expr = mean(expr, na.rm = TRUE),
      median_expr = median(expr, na.rm = TRUE),
      pct_expr = mean(expr > 0, na.rm = TRUE) * 100,
      n_cells = dplyr::n(),
      .groups = "drop"
    )

  print(df)

  p <- ggplot(df, aes(x = group_plot, y = mean_expr, group = sample_plot, color = sample_plot)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    labs(
      title = paste0(
        gene_name, " (", celltype, ")\n",
        "Paired GFP comparison by sample"
      ),
      x = NULL,
      y = "Mean expression"
    ) +
    theme_classic(base_size = 12)

  return(p)
}

# example
plot_gene_paired_celltype(
  cds = cds,
  gene = "Igfbp5",
  celltype = "Adventitial Fibroblasts"
)

