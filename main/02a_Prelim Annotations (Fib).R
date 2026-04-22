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
  "1" = "B cells (Cd79a+)",
  "2" = "T cells (Cd3d+)",
  "3" = "Endothelial",
  "4" = "Alveolar Macrophage",
  "5" = "AT1",
  "6" = "NK cells",
  "7" = "Interstitial Macrophage",
  "8" = "Interstitial Macrophage",
  "9" = "Immune",
  "10" = "Fibroblast (Col1a1+",
  "11" = "DCs",
  "12" = "Fibroblast",
  "13" = "Fibroblast (Acta2+)",
  "14" = "Macrophage",
  "15" = "Ciliated Epithelial",
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
  "Pdgfrb", "Rgs5", "Cspg4",  "Myh11", "Mcam", "Des", "Cox4i2",
  # mesothelial
  "Upk3b", "Msln", "Krt19", "Krt8",
  # gfp
  "GFP_transgene",
  # proliferation
  "Mki67", "Top2a",
  #other checks
  "Cd68", "Ptprc", "Epcam", "Pecam1", "Cdh5"
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
  "3" = "Hematopoietic doublets",
  "4" = "Pericytes/SMCs",
  "5" = "Adventitial Fibroblasts",
  "6" = "Myofibroblasts",
  "7" = "Endothelial doublets",
  "8" = "Alveolar Fibroblasts"
)

cds.stromal <- cds.stromal[, which(pData(cds.stromal)$subset_cluster_names != "Hematopoietic doublets" &
  pData(cds.stromal)$subset_cluster_names != "Endothelial doublets")]

cds.stromal <- cds.stromal %>%
  preprocess_cds(num_dim = 30) %>%
  reduce_dimension() %>%
  cluster_cells(k = 20, random_seed = 3752)

plot_cells(cds.stromal, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1.25) + 
  theme_void()
plot_cells(cds.stromal, color_cells_by = "GFP_status", label_cell_groups = FALSE)
plot_cells(cds.stromal, color_cells_by = "GFP_transgene_expr", label_cell_groups = FALSE)
if ("sample" %in% colnames(pData(cds.stromal))) {
  plot_cells(cds.stromal, color_cells_by = "sample", label_cell_groups = FALSE, cell_size = 1.25)
}

pData(cds.stromal)$subset_cluster_names <- as.character(clusters(cds.stromal))
pData(cds.stromal)$subset_cluster_names <- dplyr::recode(
  pData(cds.stromal)$subset_cluster_names,
  "1" = "Alveolar Fibroblasts",
  "2" = "SMCs",
  "3" = "Epithelial doublets",
  "4" = "Pericytes",
  "5" = "Adventitial Fibroblasts",
  "6" = "Myofibroblasts"
)

cds.stromal <- cds.stromal[, which(pData(cds.stromal)$subset_cluster_names != "Epithelial doublets")]


pData(cds.stromal)$subset_cluster_names <- as.character(clusters(cds.stromal))
pData(cds.stromal)$subset_cluster_names <- dplyr::recode(
  pData(cds.stromal)$subset_cluster_names,
  "1" = "Alveolar Fibroblasts",
  "2" = "SMCs",
  "3" = "Pericytes",
  "4" = "Adventitial Fibroblasts",
  "5" = "Myofibroblasts"
)

plot_cells(cds.stromal, color_cells_by = "subset_cluster_names",
 label_cell_groups = F, cell_size=1.25, group_label_size=8) + 
   theme(text=element_text(size=18))

save_monocle_objects(cds.stromal, directory_path = "data/cds/Shi_TSC2_stromal-cells/")
# --------------------------------------------------------------------------------------

plot_cells(cds.stromal, color_cells_by = "GFP_status",
 label_cell_groups = FALSE, cell_size=1.25)+ 
   theme(text=element_text(size=18))

plot_cells(cds.stromal, color_cells_by = "GFP_transgene_expr",
 label_cell_groups = FALSE, cell_size=1.25)+ 
   theme(text=element_text(size=18))
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

plot_genes_by_group(
  cds.stromal,
  group_cells_by = "subset_cluster_names",
  markers = c(
    "GFP_transgene",
    "Pdgfra", "Npnt",       # alveolar fibroblasts
    "Pi16", "Dcn",          # adventitial fibroblasts
    "Pdgfrb", "Cspg4", "Rgs5",  # pericytes
    "Acta2", "Tagln", "Myh11",  # SMCs
    "Col1a1", "Cthrc1", "Fn1"   # myofibroblasts
  ),
  ordering_type = "none"
) + theme(text = element_text(size = 18))




library(dplyr)
stromal_df <- as.data.frame(pData(cds.stromal)) %>%
  filter(!is.na(subset_cluster_names)) %>%
  mutate(
    GFP_status = ifelse(GFP_transgene_expr > 0.1, "GFP+", "GFP-")
  )
stromal_comp <- stromal_df %>%
  group_by(sample_id, subset_cluster_names) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(
    total_cells = sum(n),
    prop = n / total_cells,
    percent = 100 * prop
  )

ggplot(stromal_comp, aes(x = subset_cluster_names, y = percent, fill = sample_id)) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(
    x = "Stromal subtype",
    y = "Percent of stromal cells",
    title = "Relative stromal composition by sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
  )

stromal_gfp_frac <- stromal_df %>%
  group_by(sample_id, subset_cluster_names) %>%
  summarise(
    n_cells = dplyr::n(),
    n_GFP_pos = sum(GFP_transgene_expr > 0.1, na.rm = TRUE),
    frac_GFP_pos = mean(GFP_transgene_expr > 0.1, na.rm = TRUE),
    percent_GFP_pos = 100 * frac_GFP_pos,
    .groups = "drop"
  )

ggplot(stromal_gfp_frac, aes(x = subset_cluster_names, y = percent_GFP_pos, fill = sample_id)) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(
    x = "Stromal subtype",
    y = "Percent GFP+ cells",
    title = "Fraction of GFP+ cells within each stromal subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
  )


stromal_gfp_expr <- stromal_df %>%
  group_by(sample_id, subset_cluster_names) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_GFP = mean(GFP_transgene_expr, na.rm = TRUE),
    median_GFP = median(GFP_transgene_expr, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(stromal_gfp_expr, aes(x = subset_cluster_names, y = mean_GFP, fill = sample_id)) +
  geom_col(position = "dodge") +
  theme_classic() +
  labs(
    x = "Stromal subtype",
    y = "Mean GFP expression",
    title = "Mean GFP expression within each stromal subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
  )

stromal_gfp_stack <- stromal_df %>%
  group_by(sample_id, subset_cluster_names, GFP_status) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(sample_id, subset_cluster_names) %>%
  mutate(
    total = sum(n),
    percent = 100 * n / total
  )

ggplot(stromal_gfp_stack, aes(x = subset_cluster_names, y = percent, fill = GFP_status)) +
  geom_col(position = "fill") +
  facet_wrap(~sample_id) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    x = "Stromal subtype",
    y = "Relative fraction",
    title = "GFP+ and GFP- composition within each stromal subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 16)
  )
table(pData(cds.stromal)$sample)
pData(cds.stromal)


# DEGs between GFP+ and GFP- within each stromal subtype
deg_results <- list()
for (subtype in unique(pData(cds.stromal)$subset_cluster_names)) {
  cat("Running DE analysis for subtype:", subtype, "\n")
  
}
