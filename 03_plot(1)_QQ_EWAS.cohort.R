################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 29-Jun-2025
# This script is to generate and combine cohort-specific QQ plots for EWAS results.

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
library(png)
library(grid)

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
  ALSPAC   = "ALSPAC (450K)",
  GenR     = "Generation R (450K)",
  INMA     = "INMA (450K)",
  MoBa1    = "MoBa-1 (450K)",
  MoBa2    = "MoBa-2 (450K)",
  MoBa4    = "MoBa-4 (EPIC)",
  MoBa8    = "MoBa-8 (EPIC)",
  NORTHPOP = "NorthPop (EPIC)",
  Viva     = "Project Viva (450K)"
)

lambda_gc <- function(p) {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  median(qchisq(1 - p, 1), na.rm = TRUE) / qchisq(0.5, 1)
}

qq_plot <- function(pvals, title_txt, subtitle_txt) {
  pvals <- pvals[is.finite(pvals) &
                   !is.na(pvals) & pvals > 0 & pvals <= 1]
  if (length(pvals) < 10)
    return(NULL)
  
  observed <- -log10(sort(pvals))
  expected <- -log10(ppoints(length(pvals)))
  
  ggplot(data.frame(Expected = expected, Observed = observed),
         aes(Expected, Observed)) +
    geom_point(size = .4) +
    geom_abline(slope = 1,
                intercept = 0,
                linetype = "dashed") +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Expected –log10(P)",
      y = "Observed –log10(P)"
    ) +
    theme_bw(10) +
    theme(plot.title    = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 15))
}

# ── 3. Loop over cohorts ──────────────────────────────────────────────────────
cohorts <- list.dirs(root_res, full.names = TRUE, recursive = FALSE)

for (coh_dir in cohorts) {
  cohort <- basename(coh_dir)
  if (cohort %in% c("meta", "BiB", "DCHS"))
    next
  
  message("\n=== Cohort: ", cohort, " ===")
  qq_dir <- file.path(coh_dir, "QQ")
  dir.create(qq_dir, showWarnings = FALSE)
  
  # 3a. Generate single QQ plots
  csvs <- list.files(coh_dir, pattern = "\\.ewas\\.res\\..+\\.csv$", full.names = TRUE)
  for (f in csvs) {
    parts <- str_split(basename(f), "\\.")[[1]]
    if (length(parts) < 6)
      next
    model    <- parts[4]
    exposure <- parts[5]
    if (!model %in% models || !exposure %in% exp_order)
      next
    
    dt   <- fread(f, select = c("P", "p"))
    pcol <- intersect(c("P", "p"), names(dt))[1]
    if (is.na(pcol))
      next
    
    lam_txt <- format(round(lambda_gc(dt[[pcol]]), 2), nsmall = 2)
    out_png <- file.path(qq_dir, glue("QQ_{exposure}_{model}.png"))
    
    ggsave(
      out_png,
      qq_plot(
        dt[[pcol]],
        glue("{exp_labels[[exposure]]} - {model_labels[[model]]}"),
        glue("lambda = {lam_txt}")
      ),
      width = 4,
      height = 4,
      dpi = 300
    )
  }
  
  # 3b. Combine QQs into grid
  plot_list <- list()
  for (exp in exp_order) {
    for (model in models) {
      img_path <- file.path(qq_dir, glue("QQ_{exp}_{model}.png"))
      if (file.exists(img_path)) {
        img  <- readPNG(img_path)
        plot_list[[paste(exp, model, sep = "_")]] <-
          ggdraw() + draw_grob(rasterGrob(img, interpolate = FALSE))
      }
    }
  }
  if (length(plot_list) == 0) {
    message("  • No QQs found, skipping combine.")
    next
  }
  
  present_exp <- unique(str_extract(names(plot_list), "^[^_]+"))
  grid <- plot_grid(
    plotlist = plot_list,
    nrow = length(present_exp),
    ncol = length(models),
    align = "hv"
  )
  
  disp_name <- cohort_labels[cohort]
  combined <- ggdraw() +
    draw_label(
      glue("{disp_name}"),
      x = 0.01,
      y = 0.98,
      hjust = 0,
      vjust = 1,
      size = 10,
      fontface = "bold"
    ) +
    draw_plot(
      grid,
      x = 0,
      y = 0,
      width = 1,
      height = 0.95
    )
  
  ggsave(
    file.path(qq_dir, glue("Comb_QQ_{cohort}.png")),
    combined,
    width = 20,
    height = 4 * length(present_exp) + 5,
    units = "cm",
    dpi = 300
  )
  message("  • Combined QQ saved.")
}

message("\nAll cohorts finished.")

################################################################################
