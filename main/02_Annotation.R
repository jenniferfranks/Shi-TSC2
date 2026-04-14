# Mouse whole lung annotation workflow
# Stromal / mesenchymal prioritized, GFP-aware
# Jennifer Franks style: broad annotation first, then iterative subclustering,
# then transfer subset annotations back into the parent cds.
#
# Key biological focus:
# - majority of mesenchymal cells are GFP+
# - identify fibroblast subclusters (AF1, AF2, peribronchial, etc.)
# - quantify GFP+ / GFP- composition within each stromal subtype
#
# Notes:
# - This starts from a processed / clustered monocle cds object.
# - The broad cluster recoding and subtype recoding are placeholders and should be
#   updated after inspecting your own marker plots.
# - GFP status is defined directly from expression of GFP_transgene.
# - Threshold for GFP positivity is currently set to exprs(cds)["GFP_transgene", ] >= 0.1

suppressPackageStartupMessages({
  library(dplyr)
  library(monocle3)
  library(monocle3helper)
  library(Matrix)
  library(tidyr)
  library(SingleCellExperiment)
  library(ggplot2)
})

# ==============================================================================
# 0. SETUP
# ==============================================================================

main_dir <- "/data/hps/assoc/private/franks_lab/user/jfra11/projects/Silica/Shi/TSC2/"
setwd(main_dir)

cds <- load_monocle_objects(
  directory_path = "data/cds/Shi_TSC2_raw_scrublet_metadata_filtered_processed_clustered/"
)

if (!"cell" %in% colnames(pData(cds))) {
  pData(cds)$cell <- rownames(pData(cds))
}

# save original monocle clusters
pData(cds)$global_cluster <- clusters(cds)

# ------------------------------------------------------------------------------
# GFP status from GFP_transgene expression
# ------------------------------------------------------------------------------

if (!"gene_short_name" %in% colnames(rowData(cds))) {
  stop("rowData(cds) must contain gene_short_name")
}

if (!("GFP_transgene" %in% rowData(cds)$gene_short_name) &&
    !("GFP_transgene" %in% rownames(cds))) {
  stop("GFP_transgene not found in rowData(cds)$gene_short_name or rownames(cds)")
}

if ("GFP_transgene" %in% rownames(cds)) {
  gfp_gene_id <- "GFP_transgene"
} else {
  gfp_gene_id <- rownames(cds)[rowData(cds)$gene_short_name == "GFP_transgene"][1]
}

gfp_expr <- exprs(cds)[gfp_gene_id, ]
pData(cds)$GFP_transgene_expr <- as.numeric(gfp_expr)

# binary GFP call
gfp_threshold <- 0.1
pData(cds)$GFP_status <- ifelse(
  pData(cds)$GFP_transgene_expr >= gfp_threshold,
  "GFP+",
  "GFP-"
)
pData(cds)$GFP_status <- factor(pData(cds)$GFP_status, levels = c("GFP-", "GFP+"))

# quick checks
summary(pData(cds)$GFP_transgene_expr)
quantile(
  pData(cds)$GFP_transgene_expr,
  probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1),
  na.rm = TRUE
)
table(pData(cds)$GFP_status, useNA = "always")

# ==============================================================================
# 1. GLOBAL FIRST PASS
# ==============================================================================

pData(cds)
colnames(pData(cds))

if ("Sample.ID" %in% colnames(pData(cds))) table(pData(cds)$Sample.ID)
if ("model" %in% colnames(pData(cds))) table(pData(cds)$model)
table(pData(cds)$GFP_status, useNA = "always")

plot_cells(cds, color_cells_by = "cluster", group_label_size = 8)
if ("sample" %in% colnames(pData(cds))) {
  plot_cells(cds, color_cells_by = "sample", label_cell_groups = FALSE)
}
if ("model" %in% colnames(pData(cds))) {
  plot_cells(cds, color_cells_by = "model", label_cell_groups = FALSE)
}
plot_cells(cds, color_cells_by = "GFP_status", label_cell_groups = FALSE)
plot_cells(cds, color_cells_by = "GFP_transgene_expr", label_cell_groups = FALSE)
plot_cells(cds, color_cells_by = "scrublet_score")
plot_cells(cds, color_cells_by = "num_genes_expressed")

# ------------------------------------------------------------------------------
# Top markers
# ------------------------------------------------------------------------------

marker_test_res <- top_markers(cds, group_cells_by = "global_cluster")

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers$gene_id)

plot_genes_by_group(
  cds,
  top_specific_marker_ids,
  group_cells_by = "global_cluster",
  ordering_type = "maximal_on_diag",
  max.size = 3
)

# ------------------------------------------------------------------------------
# Whole lung marker panel
# ------------------------------------------------------------------------------

plot_genes_by_group(
  cds,
  markers = c(
    # epithelial
    "Epcam", "Ager", "Pdpn", "Hopx", "Sftpc", "Sftpa1", "Scgb1a1", "Cyp2f2", "Foxj1", "Muc5b",
    # endothelial
    "Pecam1", "Cdh5", "Kdr", "Car4", "Gpihbp1", "Ephb4", "Vwf", "Lyve1", "Mmrn1",
    # stromal / mesenchymal
    "Col1a1", "Col1a2", "Dcn", "Pdgfra", "Pi16", "Dpt", "Npnt", "Tcf21",
    "Col15a1", "Inmt", "Gsn", "Acta2", "Tagln", "Myh11", "Rgs5", "Cspg4", "Pdgfrb", "Upk3b", "Msln",
    "Cthrc1", "Postn", "Thy1",
    # immune
    "Ptprc", "Adgre1", "Cd68", "Siglecf", "Marco", "C1qa", "Ccr2", "Ly6g",
    "Flt3", "Clec9a", "Siglech", "Cd3d", "Cd4", "Cd8a", "Ncr1", "Cd79a", "Jchain",
    # gfp
    "GFP_transgene",
    # proliferation
    "Mki67", "Top2a"
  ),
  group_cells_by = "global_cluster",
  ordering_type = "none",
  max.size = 3
)

# ------------------------------------------------------------------------------
# QC by cluster
# ------------------------------------------------------------------------------

pdat <- as.data.frame(pData(cds))

pdat %>%
  group_by(global_cluster) %>%
  summarise(
    mean_scrublet = mean(scrublet_score, na.rm = TRUE),
    median_scrublet = median(scrublet_score, na.rm = TRUE),
    mean_genes = mean(num_genes_expressed, na.rm = TRUE),
    mean_gfp_expr = mean(GFP_transgene_expr, na.rm = TRUE),
    frac_gfp_pos = mean(GFP_status == "GFP+", na.rm = TRUE),
    n = dplyr::n()
  ) %>%
  print(n = 100)

# ==============================================================================
# 2. BROAD ANNOTATION
# ==============================================================================
# Replace this mapping after inspecting your dataset.
# Keep broad labels broad. The point is just to carve out the major compartments.

pData(cds)$broad_annotation <- as.character(pData(cds)$global_cluster)

pData(cds)$broad_annotation <- dplyr::recode(
  pData(cds)$broad_annotation,
  "1" = "Immune",
  "2" = "T cells",
  "3" = "Endothelial",
  "4" = "Alveolar Macrophage",
  "5" = "Epithelial",
  "6" = "NK cells",
  "7" = "Macrophage",
  "8" = "Macrophage",
  "9" = "Immune",
  "10" = "Fibroblast",
  "11" = "DCs",
  "12" = "Fibroblast",
  "13" = "Fibroblast",
  "14" = "Macrophage",
  "15" = "Epithelial",
  "16" = "Doublets?",
  "17" = "Fibroblast",
  "18" = "Macrophage",
  "19" = "Endothelial", 
  "20" = "pDCs",
  "21" = "Epithelial",
  "22" = "Fibroblast",
  "23" = "Fibroblast"
)

plot_cells(cds, color_cells_by = "broad_annotation", label_cell_groups = FALSE)
plot_cells(
  cds,
  color_cells_by = "broad_annotation",
  label_cell_groups = F,
  group_label_size = 6,
  label_groups_by_cluster = FALSE
)

pData(cds)$fine_annotation <- as.character(pData(cds)$broad_annotation)

# ==============================================================================
# 3. OPTIONAL LIGHT CHECKS OF NON-STROMAL COMPARTMENTS
# ==============================================================================
# Keep these light. The main work is stromal.

# ------------------------------------------------------------------------------
# epithelial
# ------------------------------------------------------------------------------
cds.epithelial <- cds[, which(grepl("Epithelial", pData(cds)$broad_annotation))]
if (ncol(cds.epithelial) > 20) {
  cds.epithelial <- cds.epithelial %>%
    preprocess_cds(num_dim = 30) %>%
    reduce_dimension() %>%
    cluster_cells(k = 15, random_seed = 3752)

  plot_cells(cds.epithelial, color_cells_by = "cluster", label_cell_groups = FALSE)

  plot_genes_by_group(
    cds.epithelial,
    markers = c(
      "Epcam", "Ager", "Pdpn", "Hopx", "Sftpc", "Sftpa1",
      "Scgb1a1", "Cyp2f2", "Foxj1", "Muc5b", "Krt5", "Krt14",
      "Ascl1", "Calca", "Mki67", "Top2a"
    ),
    group_cells_by = "cluster",
    ordering_type = "none",
    max.size = 3
  )
}

# ------------------------------------------------------------------------------
# endothelial
# ------------------------------------------------------------------------------
cds.endo <- cds[, which(grepl("Endothelial", pData(cds)$broad_annotation))]
if (ncol(cds.endo) > 20) {
  cds.endo <- cds.endo %>%
    preprocess_cds(num_dim = 30) %>%
    reduce_dimension() %>%
    cluster_cells(k = 15, random_seed = 3752)

  plot_cells(cds.endo, color_cells_by = "cluster", label_cell_groups = FALSE)

  plot_genes_by_group(
    cds.endo,
    markers = c(
      "Pecam1", "Cdh5", "Kdr",
      "Car4", "Gpihbp1",
      "Flt1", "Ramp2",
      "Ephb4", "Vwf",
      "Lyve1", "Mmrn1", "Pdpn",
      "Mki67", "Top2a"
    ),
    group_cells_by = "cluster",
    ordering_type = "none",
    max.size = 3
  )
}

# ------------------------------------------------------------------------------
# immune
# ------------------------------------------------------------------------------
cds.immune <- cds[, which(
  grepl("Myeloid", pData(cds)$broad_annotation) |
    grepl("Lymphoid", pData(cds)$broad_annotation)
)]

if (ncol(cds.immune) > 20) {
  cds.immune <- cds.immune %>%
    preprocess_cds(num_dim = 30) %>%
    reduce_dimension() %>%
    cluster_cells(k = 20, random_seed = 3752)

  plot_cells(cds.immune, color_cells_by = "cluster", label_cell_groups = FALSE)

  plot_genes_by_group(
    cds.immune,
    markers = c(
      "Ptprc", "Adgre1", "Cd68", "Siglecf", "Marco", "C1qa", "Ccr2", "Ly6g",
      "Flt3", "Clec9a", "Siglech", "Cd3d", "Cd4", "Cd8a", "Ncr1", "Cd79a", "Jchain",
      "Mki67", "Top2a"
    ),
    group_cells_by = "cluster",
    ordering_type = "none",
    max.size = 3
  )
}

# ==============================================================================
# 4. STROMAL / MESENCHYMAL COMPARTMENT
# ==============================================================================
# This is the main analysis branch.

cds.stromal <- cds[, which(grepl("Fibroblast", pData(cds)$broad_annotation))]

cds.stromal <- cds.stromal %>%
  preprocess_cds(num_dim = 30) %>%
  reduce_dimension() %>%
  cluster_cells(k = 20, random_seed = 3752)

plot_cells(cds.stromal, color_cells_by = "cluster", label_cell_groups = FALSE)
plot_cells(cds.stromal, color_cells_by = "GFP_status", label_cell_groups = FALSE)
plot_cells(cds.stromal, color_cells_by = "GFP_transgene_expr", label_cell_groups = FALSE)
if ("sample" %in% colnames(pData(cds.stromal))) {
  plot_cells(cds.stromal, color_cells_by = "sample", label_cell_groups = FALSE)
}

marker_test_res <- top_markers(cds.stromal, group_cells_by = "cluster")
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.15) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers$gene_id)

plot_genes_by_group(
  cds.stromal,
  top_specific_marker_ids,
  group_cells_by = "cluster",
  ordering_type = "maximal_on_diag",
  max.size = 3
)

# ------------------------------------------------------------------------------
# Stromal marker panel
# ------------------------------------------------------------------------------

stromal_markers <- c(
  # pan fibroblast / matrix
  "Col1a1", "Col1a2", "Dcn", "Pdgfra", "Col3a1", "Lum",
  # alveolar fibroblast-like
  "Npnt", "Tcf21", "Inmt", "Gsn",
  # adventitial / peribronchial-like
  "Pi16", "Dpt", "Col15a1", "Cxcl14", "Mfap5",
  # activated / myofibroblast-like
  "Cthrc1", "Postn", "Acta2", "Tagln", "Tnc", "Thy1",
  # pericyte / smc
  "Pdgfrb", "Rgs5", "Cspg4",  "Myh11", "Mcam", "Des",
  # mesothelial
  "Upk3b", "Msln", "Krt19", "Krt8",
  # gfp
  "GFP_transgene",
  # proliferation
  "Mki67", "Top2a"
)

plot_genes_by_group(
  cds.stromal,
  stromal_markers,
  group_cells_by = "cluster",
  ordering_type = "none",
  max.size = 3
)

# ------------------------------------------------------------------------------
# First-pass stromal labels
# ------------------------------------------------------------------------------
# Replace after inspection.

pData(cds.stromal)$subset_cluster_names <- as.character(clusters(cds.stromal))
pData(cds.stromal)$subset_cluster_names <- dplyr::recode(
  pData(cds.stromal)$subset_cluster_names,
  "1" = "Alveolar Fibroblasts",
  "2" = "Pericytes/SMCs",
  "3" = "Runx1+ cells, maybe doublets",
  "4" = "Pericytes/SMCs",
  "5" = "Adventitial Fibroblasts",
  "6" = "Myofibroblasts",
  "7" = "Proliferating Fibroblasts",
  "8" = "Alveolar Fibroblasts"
)

plot_cells(cds.stromal, color_cells_by = "subset_cluster_names",
 label_cell_groups = FALSE, cell_size=1.25)

plot_cells(cds.stromal, color_cells_by = "GFP_status",
 label_cell_groups = FALSE, cell_size=1.25)
# ------------------------------------------------------------------------------
# GFP summary within first-pass stromal states
# ------------------------------------------------------------------------------

stromal_gfp_summary <- as.data.frame(pData(cds.stromal)) %>%
  filter(!is.na(subset_cluster_names)) %>%
  group_by(subset_cluster_names, GFP_status) %>%
  summarise(
    n = dplyr::n(),
    mean_GFP_expr = mean(GFP_transgene_expr, na.rm = TRUE),
    median_GFP_expr = median(GFP_transgene_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(subset_cluster_names) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    percent = 100 * prop
  ) %>%
  arrange(subset_cluster_names, GFP_status)

print(stromal_gfp_summary, n = 100)

# ==============================================================================
# 5. FIBROBLAST-FOCUSED RECLUSTERING
# ==============================================================================
# This is where AF1 / AF2 / peribronchial / myofibroblast-style states get resolved.

cds.fib <- cds.stromal[, which(
  pData(cds.stromal)$subset_cluster_names %in% c(
    "Fibroblasts", "Activated Fibroblasts", "Proliferating Fibroblasts"
  )
)]

cds.fib <- cds.fib %>%
  preprocess_cds(num_dim = 30) %>%
  reduce_dimension() %>%
  cluster_cells(k = 20, random_seed = 3752)

plot_cells(cds.fib, color_cells_by = "cluster", label_cell_groups = FALSE)
plot_cells(cds.fib, color_cells_by = "GFP_status", label_cell_groups = FALSE)
plot_cells(cds.fib, color_cells_by = "GFP_transgene_expr", label_cell_groups = FALSE)

marker_test_res <- top_markers(cds.fib, group_cells_by = "cluster")
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.15) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)
top_specific_marker_ids <- unique(top_specific_markers$gene_id)

plot_genes_by_group(
  cds.fib,
  top_specific_marker_ids,
  group_cells_by = "cluster",
  ordering_type = "maximal_on_diag",
  max.size = 3
)

fibro_markers <- c(
  # alveolar fibroblasts
  "Npnt", "Tcf21", "Inmt", "Gsn",
  # adventitial / peribronchial / PI16-like
  "Pi16", "Dpt", "Col15a1", "Cxcl14", "Mfap5",
  # activated / matrix-high / myofibroblast-like
  "Cthrc1", "Postn", "Tnc", "Thy1", "Acta2", "Tagln",
  # general fibroblast
  "Col1a1", "Col1a2", "Dcn", "Pdgfra", "Lum",
  # gfp
  "GFP_transgene",
  # proliferation
  "Mki67", "Top2a"
)

plot_genes_by_group(
  cds.fib,
  fibro_markers,
  group_cells_by = "cluster",
  ordering_type = "none",
  max.size = 3
)

# ------------------------------------------------------------------------------
# Fibroblast subtype labels
# ------------------------------------------------------------------------------
# Replace this after examining the actual marker patterns.

pData(cds.fib)$subset_cluster_names <- as.character(clusters(cds.fib))
pData(cds.fib)$subset_cluster_names <- dplyr::recode(
  pData(cds.fib)$subset_cluster_names,
  "1" = "AF1",
  "2" = "AF2",
  "3" = "Peribronchial Fibroblasts",
  "4" = "Myofibroblasts",
  "5" = "Activated Matrix Fibroblasts",
  "6" = "Proliferating Fibroblasts",
  "7" = "AF1",
  "8" = "AF2"
)

plot_cells(cds.fib, color_cells_by = "subset_cluster_names", label_cell_groups = FALSE)
plot_cells(
  cds.fib,
  color_cells_by = "subset_cluster_names",
  label_cell_groups = TRUE,
  label_groups_by_cluster = FALSE,
  group_label_size = 5
)

# ------------------------------------------------------------------------------
# GFP summary within fibroblast states
# ------------------------------------------------------------------------------

fib_gfp_summary <- as.data.frame(pData(cds.fib)) %>%
  filter(!is.na(subset_cluster_names)) %>%
  group_by(subset_cluster_names, GFP_status) %>%
  summarise(
    n = dplyr::n(),
    mean_GFP_expr = mean(GFP_transgene_expr, na.rm = TRUE),
    median_GFP_expr = median(GFP_transgene_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(subset_cluster_names) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    percent = 100 * prop
  ) %>%
  arrange(subset_cluster_names, GFP_status)

print(fib_gfp_summary, n = 100)

# call each subtype mostly GFP+, mostly GFP-, or mixed
fib_gfp_call_summary <- fib_gfp_summary %>%
  select(subset_cluster_names, GFP_status, prop) %>%
  tidyr::pivot_wider(
    names_from = GFP_status,
    values_from = prop,
    values_fill = 0
  )

if (!"GFP+" %in% colnames(fib_gfp_call_summary)) fib_gfp_call_summary$`GFP+` <- 0
if (!"GFP-" %in% colnames(fib_gfp_call_summary)) fib_gfp_call_summary$`GFP-` <- 0

fib_gfp_call_summary <- fib_gfp_call_summary %>%
  mutate(
    GFP_class = case_when(
      `GFP+` >= 0.80 ~ "Mostly GFP+",
      `GFP-` >= 0.80 ~ "Mostly GFP-",
      TRUE ~ "Mixed"
    )
  )

print(fib_gfp_call_summary, n = 100)

# ------------------------------------------------------------------------------
# visualize GFP composition within fibroblast subtypes
# ------------------------------------------------------------------------------

ggplot(fib_gfp_summary, aes(x = subset_cluster_names, y = percent, fill = GFP_status)) +
  geom_col(position = "fill") +
  coord_flip() +
  ylab("Fraction within subtype") +
  xlab("") +
  theme_classic(base_size = 14)

ggplot(fib_gfp_summary, aes(x = subset_cluster_names, y = n, fill = GFP_status)) +
  geom_col(position = "stack") +
  coord_flip() +
  xlab("") +
  ylab("Cell count") +
  theme_classic(base_size = 14)

# ==============================================================================
# 6. TRANSFER REFINED ANNOTATIONS BACK
# ==============================================================================

cds <- transfer_cell_annotations(
  cds_ref = cds.stromal,
  cds_target = cds,
  annotation_col = "subset_cluster_names",
  new_annotation_col = "fine_annotation",
  cell_id_col = "cell"
)

cds <- transfer_cell_annotations(
  cds_ref = cds.fib,
  cds_target = cds,
  annotation_col = "subset_cluster_names",
  new_annotation_col = "fine_annotation",
  cell_id_col = "cell"
)

plot_cells(cds, color_cells_by = "fine_annotation", label_cell_groups = FALSE)

# ==============================================================================
# 7. FINALIZE LINEAGES
# ==============================================================================

pData(cds)$lineage <- as.character(pData(cds)$fine_annotation)

pData(cds)$lineage <- ifelse(
  pData(cds)$fine_annotation %in% c(
    "AF1", "AF2", "Peribronchial Fibroblasts", "Myofibroblasts",
    "Activated Matrix Fibroblasts", "Proliferating Fibroblasts",
    "Fibroblasts", "Activated Fibroblasts", "Pericytes", "SMCs", "Mesothelial"
  ),
  "Stromal",
  pData(cds)$lineage
)

pData(cds)$lineage <- ifelse(
  pData(cds)$broad_annotation %in% c("Myeloid"),
  "Myeloid",
  pData(cds)$lineage
)

pData(cds)$lineage <- ifelse(
  pData(cds)$broad_annotation %in% c("Lymphoid"),
  "Lymphoid",
  pData(cds)$lineage
)

pData(cds)$lineage <- ifelse(
  pData(cds)$broad_annotation %in% c("Epithelial"),
  "Epithelial",
  pData(cds)$lineage
)

pData(cds)$lineage <- ifelse(
  pData(cds)$broad_annotation %in% c("Endothelial"),
  "Endothelial",
  pData(cds)$lineage
)

pData(cds)$lineage <- factor(
  pData(cds)$lineage,
  levels = c("Myeloid", "Lymphoid", "Epithelial", "Endothelial", "Stromal")
)

plot_cells(cds, color_cells_by = "lineage", label_cell_groups = FALSE)

# ==============================================================================
# 8. FINAL GFP SUMMARIES IN STROMAL / FIBROBLAST STATES
# ==============================================================================

gfp_summary <- as.data.frame(pData(cds)) %>%
  filter(!is.na(fine_annotation)) %>%
  filter(fine_annotation %in% c(
    "AF1", "AF2", "Peribronchial Fibroblasts", "Myofibroblasts",
    "Activated Matrix Fibroblasts", "Proliferating Fibroblasts",
    "Fibroblasts", "Activated Fibroblasts", "Pericytes", "SMCs", "Mesothelial"
  )) %>%
  filter(!is.na(GFP_status)) %>%
  group_by(fine_annotation, GFP_status) %>%
  summarise(
    n = dplyr::n(),
    mean_GFP_expr = mean(GFP_transgene_expr, na.rm = TRUE),
    median_GFP_expr = median(GFP_transgene_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(fine_annotation) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    percent = 100 * prop
  ) %>%
  arrange(fine_annotation, GFP_status)

print(gfp_summary, n = 100)

gfp_summary_wide <- gfp_summary %>%
  select(fine_annotation, GFP_status, n, percent, mean_GFP_expr, median_GFP_expr) %>%
  tidyr::pivot_wider(
    names_from = GFP_status,
    values_from = c(n, percent, mean_GFP_expr, median_GFP_expr),
    values_fill = 0
  )

print(gfp_summary_wide, n = 100)

gfp_call_summary <- gfp_summary %>%
  select(fine_annotation, GFP_status, prop) %>%
  tidyr::pivot_wider(
    names_from = GFP_status,
    values_from = prop,
    values_fill = 0
  )

if (!"GFP+" %in% colnames(gfp_call_summary)) gfp_call_summary$`GFP+` <- 0
if (!"GFP-" %in% colnames(gfp_call_summary)) gfp_call_summary$`GFP-` <- 0

gfp_call_summary <- gfp_call_summary %>%
  mutate(
    GFP_class = case_when(
      `GFP+` >= 0.80 ~ "Mostly GFP+",
      `GFP-` >= 0.80 ~ "Mostly GFP-",
      TRUE ~ "Mixed"
    )
  )

print(gfp_call_summary, n = 100)

# ------------------------------------------------------------------------------
# plot GFP composition by final fibroblast / stromal subtype
# ------------------------------------------------------------------------------

ggplot(gfp_summary, aes(x = fine_annotation, y = percent, fill = GFP_status)) +
  geom_col(position = "fill") +
  coord_flip() +
  ylab("Fraction within subtype") +
  xlab("") +
  theme_classic(base_size = 14)

ggplot(gfp_summary, aes(x = fine_annotation, y = n, fill = GFP_status)) +
  geom_col(position = "stack") +
  coord_flip() +
  xlab("") +
  ylab("Cell count") +
  theme_classic(base_size = 14)

# ==============================================================================
# 9. QUICK SANITY CHECKS
# ==============================================================================

table(pData(cds)$fine_annotation, useNA = "always")
table(pData(cds)$lineage, useNA = "always")
table(pData(cds)$fine_annotation, pData(cds)$GFP_status, useNA = "always")

plot_cells(
  cds,
  genes = c(
    "GFP_transgene",
    "Col1a1", "Dcn", "Pdgfra", "Npnt", "Tcf21", "Pi16", "Dpt",
    "Cthrc1", "Postn", "Acta2", "Rgs5", "Cspg4", "Myh11", "Upk3b"
  ),
  scale_to_range = FALSE
)

plot_genes_by_group(
  cds,
  c(
    "GFP_transgene",
    "Col1a1", "Dcn", "Pdgfra", "Npnt", "Tcf21", "Pi16", "Dpt",
    "Cthrc1", "Postn", "Acta2", "Rgs5", "Cspg4", "Myh11", "Upk3b"
  ),
  group_cells_by = "fine_annotation",
  ordering_type = "none",
  max.size = 3
)

plot_genes_by_group(
  cds,
  "GFP_transgene",
  group_cells_by = "global_cluster",
  ordering_type = "none",
  max.size = 3
)

plot_genes_by_group(
  cds.stromal,
  c(
    "GFP_transgene", "Col1a1", "Dcn", "Pdgfra", "Npnt", "Pi16", "Dpt",
    "Cthrc1", "Postn", "Acta2", "Rgs5", "Cspg4", "Myh11", "Upk3b"
  ),
  group_cells_by = "cluster",
  ordering_type = "none",
  max.size = 3
)

plot_genes_by_group(
  cds.fib,
  c(
    "GFP_transgene", "Npnt", "Tcf21", "Inmt", "Gsn",
    "Pi16", "Dpt", "Col15a1", "Cxcl14",
    "Cthrc1", "Postn", "Tnc", "Acta2", "Tagln",
    "Col1a1", "Col1a2", "Dcn", "Pdgfra"
  ),
  group_cells_by = "subset_cluster_names",
  ordering_type = "none",
  max.size = 3
)

# ==============================================================================
# 10. SAVE
# ==============================================================================

save_monocle_objects(cds, directory_path = "data/cds/Shi_TSC2_annotated_stromal_priority")
save.image("main/TSC2_stromal_annotation_workflow_inprogress.RData")

write.csv(stromal_gfp_summary, "data/cds/Shi_TSC2_stromal_gfp_summary.csv", row.names = FALSE)
write.csv(fib_gfp_summary, "data/cds/Shi_TSC2_fibroblast_gfp_summary.csv", row.names = FALSE)
write.csv(fib_gfp_call_summary, "data/cds/Shi_TSC2_fibroblast_gfp_call_summary.csv", row.names = FALSE)
write.csv(gfp_summary, "data/cds/Shi_TSC2_gfp_summary_by_subtype.csv", row.names = FALSE)
write.csv(gfp_summary_wide, "data/cds/Shi_TSC2_gfp_summary_by_subtype_wide.csv", row.names = FALSE)
write.csv(gfp_call_summary, "data/cds/Shi_TSC2_gfp_call_summary.csv", row.names = FALSE)