################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 09-Jul-2025
# This script is to generate and combine precision plots for EWAS results.

################################################################################

rm(list = ls())

# ── 1. Libraries ──────────────────────────────────────────────────────────────
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(glue)
library(ggplot2)
library(cowplot)

# ── 2. Constants ──────────────────────────────────────────────────────────────
root_res  <- "/user/home/zd20208/MatVegDiet_PACE_EWAS/results"
models    <- c("NoCellModel", "MinModel", "FullModel", "AddModel")
exp_order <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")

model_labels <- c(
  NoCellModel = "No-cell model",
  MinModel    = "Minimally-adjusted model",
  FullModel   = "Partially-adjusted model",
  AddModel    = "Fully-adjusted model"
)

exp_labels <- c(
  veggie1 = "Full vs. non-vegetarian",
  veggie2 = "Pesco-/full vs. non-vegetarian",
  PDI  = "Overall PDI",
  hPDI = "Healthful PDI",
  uPDI = "Unhealthful PDI"
)

cohort_labels <- c(
  ALSPAC   = "AP",
  GenR     = "GR",
  INMA     = "IM",
  MoBa1    = "MB1",
  MoBa2    = "MB2",
  MoBa4    = "MB4",
  MoBa8    = "MB8",
  NORTHPOP = "NP",
  Viva     = "PV"
)

out_dir <- file.path(root_res, "meta/Precision")
dir.create(out_dir, showWarnings = FALSE)

# ── 3. Helper functions ───────────────────────────────────────────────────────
precision_plot <- function(df, title_txt) {
  ggplot(df, aes(x = sqrt_n, y = precision, label = lab)) +
    geom_point(size = 2) +
    geom_text(nudge_y = 0.05, family = "sans") +
    scale_x_continuous("sqrt(sample size)") +
    scale_y_continuous("Precision (1 / median(SE))") +
    labs(title = title_txt) +
    theme_bw(10) +
    theme(plot.title = element_text(size = 8, face = "bold"))
}

compute_precision <- function(path) {
  dt <- fread(path, select = c("n", "se"))
  dt <- dt[is.finite(se) & !is.na(se) & se > 0]
  list(n         = unique(dt$n)[1],
       precision = 1 / median(dt$se, na.rm = TRUE))
}

# ── 4. Build master precision table ───────────────────────────────────────────
message("\nCollecting precision statistics from all cohorts…")

data_list <- list()
cohorts   <- list.dirs(root_res, full.names = TRUE, recursive = FALSE)

for (coh_dir in cohorts) {
  cohort <- basename(coh_dir)
  if (cohort %in% c("meta", "BiB", "DCHS"))
    next
  
  csvs <- list.files(coh_dir, pattern = "\\.ewas\\.res\\..+\\.csv$", full.names = TRUE)
  
  for (f in csvs) {
    parts <- str_split(basename(f), "\\.")[[1]]
    if (length(parts) < 6)
      next
    model    <- parts[4]
    exposure <- parts[5]
    if (!model %in% models || !exposure %in% exp_order)
      next
    
    out <- compute_precision(f)
    
    data_list[[length(data_list) + 1]] <- data.frame(
      exposure  = exposure,
      model     = model,
      sqrt_n    = sqrt(out$n),
      precision = out$precision,
      lab       = cohort_labels[[cohort]],
      stringsAsFactors = FALSE
    )
  }
}

precisions <- bind_rows(data_list)

# ── 5. Generate single precision plots ────────────────────────────────────────
plot_list <- list()

for (exp in exp_order) {
  for (model in models) {
    df_sub <- filter(precisions, exposure == exp & model == model)
    if (nrow(df_sub) == 0)
      next
    
    key <- glue("{exp}_{model}")
    title_txt <- glue("{exp_labels[[exp]]} – {model_labels[[model]]}")
    
    p <- precision_plot(df_sub, title_txt)
    
    ggsave(
      file.path(out_dir, glue("Precision_{exp}_{model}.png")),
      p,
      width = 4,
      height = 4,
      dpi = 300
    )
    plot_list[[key]] <- p
  }
}

# ── 6. Combine into 5 × 4 grid ────────────────────────────────────────────────
grid <- plot_grid(
  plotlist = plot_list[order(names(plot_list))],
  nrow = length(exp_order),
  ncol = length(models),
  align = "hv"
)

ggsave(
  file.path(out_dir, "Comb_Precision.png"),
  grid,
  width  = 20,
  # cm
  height = 4 * length(exp_order),
  # cm
  units  = "cm",
  dpi    = 300
)

message("\nAll precision plots finished.")

################################################################################
