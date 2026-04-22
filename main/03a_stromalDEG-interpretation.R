 ------------------------------------------------------------------------------
# Mixed-effect DEG analysis: GFP+ vs GFP- within each cell type using monocle3
# ------------------------------------------------------------------------------

library(monocle3)
library(dplyr)
library(tibble)
library(stringr)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(readr)

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

coef_table <- read.csv(
  file.path(out_dir, "tables/Adventitial_Fibroblasts_coefficient_table_full.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)
degs <- coef_table[which(coef_table$term != "(Intercept)" & coef_table$q_value <0.1),] # none are significant
degs <- coef_table[which(coef_table$term != "(Intercept)"& coef_table$p_value <0.01),]


library(gprofiler2)
gostres <- gost(query = degs$gene_short_name, organism = "mmusculus")



df <- as.data.frame(gostres$result, stringsAsFactors = FALSE)
df <- flatten_list_columns(df)
  

write_status_table(df,file.path(out_dir, "AdventitialFibroblast_gprofiler_results.csv"))
write_status_table(degs,file.path(out_dir, "AdventitialFibroblast_deg_results.csv"))

# ------ Alveolar Fibroblasts -------------------------------------------------
coef_table <- read.csv(
  file.path(out_dir, "tables/Alveolar_Fibroblasts_coefficient_table_full.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE)
degs <- coef_table[which(coef_table$term != "(Intercept)" & coef_table$q_value <0.1),] # none are significant
degs <- coef_table[which(coef_table$term != "(Intercept)"& coef_table$p_value <0.01),]


library(gprofiler2)
gostres <- gost(query = degs$gene_short_name, organism = "mmusculus")



df <- as.data.frame(gostres$result, stringsAsFactors = FALSE)
df <- flatten_list_columns(df)
  

write_status_table(df,file.path(out_dir, "AlveolarFibroblast_gprofiler_results.csv"))
write_status_table(degs,file.path(out_dir, "AlveolarFibroblast_deg_results.csv"))


# ---- Function for volcano plot of monocle3 DE results --------------------------------------
plot_volcano_directional <- function(
  coef_table,
  term_keep = "GFP_statusGFP+",
  label_n_each_side = 8,
  p_cutoff_1 = 0.05,
  p_cutoff_2 = 0.01,
  title = NULL
) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  df <- coef_table %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    filter(status == "OK") %>%
    filter(term == term_keep) %>%
    mutate(
      estimate = as.numeric(estimate),
      p_value = as.numeric(p_value)
    ) %>%
    filter(!is.na(p_value), !is.na(estimate)) %>%
    mutate(
      neglog10_p = -log10(p_value),
      direction_class = case_when(
        estimate > 0 & p_value < p_cutoff_2 ~ "Higher in GFP+ (p < 0.01)",
        estimate > 0 & p_value < p_cutoff_1 ~ "Higher in GFP+ (p < 0.05)",
        estimate < 0 & p_value < p_cutoff_2 ~ "Higher in GFP- (p < 0.01)",
        estimate < 0 & p_value < p_cutoff_1 ~ "Higher in GFP- (p < 0.05)",
        TRUE ~ "NS"
      )
    )

  df$direction_class <- factor(
    df$direction_class,
    levels = c(
      "Higher in GFP+ (p < 0.01)",
      "Higher in GFP+ (p < 0.05)",
      "Higher in GFP- (p < 0.01)",
      "Higher in GFP- (p < 0.05)",
      "NS"
    )
  )

  label_pos <- df %>%
    filter(estimate > 0) %>%
    arrange(p_value, desc(abs(estimate))) %>%
    distinct(gene_short_name, .keep_all = TRUE) %>%
    slice_head(n = label_n_each_side)

  label_neg <- df %>%
    filter(estimate < 0) %>%
    arrange(p_value, desc(abs(estimate))) %>%
    distinct(gene_short_name, .keep_all = TRUE) %>%
    slice_head(n = label_n_each_side)

  label_df <- bind_rows(label_pos, label_neg)

  xmax <- max(abs(df$estimate), na.rm = TRUE) * 1.05

  ggplot(df, aes(x = estimate, y = neglog10_p)) +
    geom_point(aes(color = direction_class), alpha = 0.8, size = 2) +
    geom_hline(yintercept = -log10(p_cutoff_1), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cutoff_2), linetype = "dotted", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4) +
    geom_text_repel(
      data = label_df,
      aes(label = gene_short_name),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2
    ) +
    scale_x_continuous(
      limits = c(-xmax, xmax),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(
      values = c(
        "Higher in GFP+ (p < 0.01)" = "#1B7837",
        "Higher in GFP+ (p < 0.05)" = "#A6DBA0",
        "Higher in GFP- (p < 0.01)" = "#B2182B",
        "Higher in GFP- (p < 0.05)" = "#F4A6A6",
        "NS" = "grey75"
      )
    ) +
    labs(
      title = title,
      x = "Model coefficient (left = higher in GFP-, right = higher in GFP+)",
      y = expression(-log[10](p_value)),
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
}

# Adventitial fibroblasts
coef_table_adv <- read.csv(
  file.path(out_dir, "tables/Adventitial_Fibroblasts_coefficient_table_full.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

p_adv <- plot_volcano_directional(
  coef_table = coef_table_adv,
  term_keep = "GFP_statusGFP+",
  #x_var = "normalized_effect",
  label_n = 12,
  title = "Adventitial Fibroblasts: GFP+ vs GFP-"
) 

ggsave(
  filename = file.path(out_dir, "Adventitial_Fibroblasts_volcano.pdf"),
  plot = p_adv,
  width = 7,
  height = 5
)

# Alveolar fibroblasts
coef_table_alv <- read.csv(
  file.path(out_dir, "tables/Alveolar_Fibroblasts_coefficient_table_full.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
)

p_alv <- plot_volcano_directional(
  coef_table = coef_table_alv,
  term_keep = "GFP_statusGFP+",
  #x_var = "normalized_effect",
  label_n = 12,
  title = "Alveolar Fibroblasts: GFP+ vs GFP-"
)

ggsave(
  filename = file.path(out_dir, "Alveolar_Fibroblasts_volcano.pdf"),
  plot = p_alv,
  width = 7,
  height = 5
)

