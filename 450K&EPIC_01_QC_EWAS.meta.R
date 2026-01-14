################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 27-Jun-2025
# This script is to perform QC for EWAS meta-analysis results (for probes available on both 450K and EPIC).

## --------------------------------------------------------------------------- ##
## 0.  Clear workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(dplyr)
library(purrr)
library(glue)
library(ggplot2)
library(qqman)
library(meffil)

## --------------------------------------------------------------------------- ##
## 2.  Directories / exposures / models / prefix                               ##
## --------------------------------------------------------------------------- ##

meta_dir   <- "Z:/working/results/EWAS/meta"
plot_dir   <- file.path(meta_dir, "")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(plot_dir, "QQ"), showWarnings = FALSE)
dir.create(file.path(plot_dir, "Manhattan"), showWarnings = FALSE)
dir.create(file.path(plot_dir, "Volcano"), showWarnings = FALSE)

prefix     <- "450K&EPIC"
exp_list   <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
exp_labels <- c(
  veggie1 = "Full vs. non-vegetarian",
  veggie2 = "Pesco-/full vs. non-vegetarian",
  PDI     = "Overall PDI",
  hPDI    = "Healthful PDI",
  uPDI    = "Unhealthful PDI"
)

model_list   <- c("NoCellModel", "MinModel", "FullModel", "AddModel")
model_labels <- c(
  NoCellModel = "No-cell model",
  MinModel    = "Minimally-adjusted model",
  FullModel   = "Partially-adjusted model",
  AddModel    = "Fully-adjusted model"
)

## --------------------------------------------------------------------------- ##
## 3.  Helper functions & annotation (via meffil)                              ##
## --------------------------------------------------------------------------- ##

harmonise_cols <- function(dt) {
  if (!"CpG" %in% names(dt)) {
    probe_candidates <- c("MarkerName", "probeid", "ProbeID", "cpg_id")
    found <- intersect(probe_candidates, names(dt))
    setnames(dt, found[1], "CpG")
  }
  if (!"P" %in% names(dt)) {
    if ("Pvalue" %in% names(dt))
      setnames(dt, "Pvalue", "P")
    if ("p"      %in% names(dt))
      setnames(dt, "p", "P")
  }
  if (!"Effect" %in% names(dt) && "Beta"   %in% names(dt))
    setnames(dt, "Beta", "Effect")
  if (!"SE"     %in% names(dt) && "StdErr" %in% names(dt))
    setnames(dt, "StdErr", "SE")
  dt
}

anno_dt <- as.data.table(meffil.get.features("epic"))[, .(CpG = name, Chr = chromosome, Pos = position)]
anno_dt[, Chr := sub("^chr", "", Chr)]
anno_dt[Chr == "X", Chr := "23"]
anno_dt[Chr == "Y", Chr := "24"]
anno_dt[, Chr := as.integer(Chr)]

lambda_gc <- function(p)
  median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

qq_plot <- function(p, title_txt, subtitle_txt) {
  observed <- -log10(sort(p))
  expected <- -log10(ppoints(length(p)))
  ggplot(data.frame(Expected = expected, Observed = observed),
         aes(Expected, Observed)) +
    geom_point(size = .4) +
    geom_abline(slope = 1,
                intercept = 0,
                linetype = "dashed") +
    labs(
      title    = title_txt,
      subtitle = subtitle_txt,
      x = "Expected –log10(P)",
      y = "Observed –log10(P)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 15))
}

## --------------------------------------------------------------------------- ##
## 4.  Loop over models: QC pipeline                                            ##
## --------------------------------------------------------------------------- ##

for (model in model_list) {
  message("Processing model: ", model, " ...")
  
  # 4.1  Build file list and check existence
  files <- glue("{meta_dir}/{prefix}_ewas.res.{model}.{exp_list}.txt")
  names(files) <- exp_list
  exists   <- file.exists(files)
  
  if (!any(exists)) {
    message("  WARNING: no files found for model '", model, "'. Skipping.")
    next
  }
  if (any(!exists)) {
    missing_exps <- names(files)[!exists]
    message(
      "  WARNING: missing exposures for: ",
      paste(missing_exps, collapse = ", "),
      ". Skipping those."
    )
    files <- files[exists]
  }
  
  # 4.2  Read & harmonise meta-analysis outputs
  meta_dt <- map(files, ~ fread(.x) %>% harmonise_cols())
  meta_dt <- map(meta_dt, ~ merge(.x, anno_dt, by = "CpG", all.x = TRUE))
  
  # 4.3  λGC table
  lambda_tbl <- map_dfr(names(meta_dt), function(e) {
    tibble(
      Exposure  = e,
      Lambda_GC = lambda_gc(meta_dt[[e]]$P),
      n_CpG     = nrow(meta_dt[[e]])
    )
  })
  fwrite(lambda_tbl, file.path(plot_dir, "QQ", glue("{prefix}_lambda.{model}.csv")))
  
  # 4.4  QQ plots
  iwalk(meta_dt, function(dt, exp) {
    lam_txt <- format(round(lambda_gc(dt$P), 2), nsmall = 2)
    qq_file <- file.path(plot_dir, "QQ", glue("{prefix}_QQ_{exp}_{model}.png"))
    ggsave(
      filename = qq_file,
      plot     = qq_plot(
        p            = dt$P,
        title_txt    = glue("{exp_labels[[exp]]} – {model_labels[[model]]}"),
        subtitle_txt = glue("lambda = {lam_txt}")
      ),
      width  = 4,
      height = 4,
      dpi = 300
    )
  })
  
  # 4.5  Manhattan & Volcano plots
  message("  Creating Manhattan & Volcano plots ...")
  suggest_y <- -log10(1e-5)
  
  for (exp in names(meta_dt)) {
    dt      <- meta_dt[[exp]]
    n_tests <- nrow(dt)
    bonf_y  <- -log10(0.05 / n_tests)
    ylim_top <- 8
    
    # Manhattan plot
    png(
      filename = file.path(
        plot_dir,
        "Manhattan",
        glue("{prefix}_Manhattan_{exp}_{model}.png")
      ),
      width  = 2400,
      height = 1500,
      res = 300
    )
    qqman::manhattan(
      dt,
      chr            = "Chr",
      bp             = "Pos",
      snp            = "CpG",
      p              = "P",
      col            = c("grey20", "grey60"),
      genomewideline = FALSE,
      suggestiveline = FALSE,
      ylim           = c(0, ylim_top),
      cex            = 0.5,
      main           = glue("{exp_labels[[exp]]} – {model_labels[[model]]}")
    )
    abline(h = suggest_y, lty = 2)  # p = 1e-5
    abline(h = bonf_y, lty = 3)  # Bonferroni
    dev.off()
    
    # Volcano plot
    dt       <- mutate(dt, log10P = -log10(P))
    vol_file <- file.path(plot_dir,
                          "Volcano",
                          glue("{prefix}_Volcano_{exp}_{model}.png"))
    g_vol <- ggplot(dt, aes(Effect, log10P)) +
      geom_point(size = .7, alpha = .3) +
      geom_hline(yintercept = suggest_y, linetype = "dashed") +
      geom_hline(yintercept = bonf_y, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_y_continuous(limits = c(0, ylim_top), expand = c(0, 0)) +
      labs(
        title = glue("{exp_labels[[exp]]} – {model_labels[[model]]}"),
        x     = "Effect size (β)",
        y     = "–log10(P)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        panel.grid.major = element_line(colour = "grey90", linewidth = .2),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, face = "bold")
      )
    ggsave(vol_file,
           g_vol,
           width = 5,
           height = 7,
           dpi = 300)
  }
  
}

message("QC pipeline finished successfully for all models on 450K&EPIC data!")

## --------------------------------------------------------------------------- ##
## 5.  Combine QQ plots                                                         ##
## --------------------------------------------------------------------------- ##

library(cowplot)
library(magick)

qq_dir   <- file.path(plot_dir, "QQ")
out_png  <- file.path(plot_dir, "QQ", glue("{prefix}_Comb_QQ_meta.png"))

exposures <- exp_list   # row order
models    <- model_list # col order

plot_list <- list()
for (exp in exposures) {
  for (mod in models) {
    img_path <- file.path(qq_dir, glue("{prefix}_QQ_{exp}_{mod}.png"))
    if (!file.exists(img_path))
      stop("Missing file: ", img_path)
    g <- ggdraw() +
      draw_image(image_read(img_path))
    plot_list[[paste(exp, mod, sep = "_")]] <- g
  }
}

# assemble 5×4 grid
grid_5x4 <- plot_grid(
  plotlist = plot_list,
  nrow     = length(exposures),
  ncol     = length(models),
  align    = "hv"
)

# add global title
combined <- ggdraw() +
  draw_label(
    "EWAS meta-analysis (probes on both 450K and EPIC)",
    fontface = "bold",
    x        = 0.01,
    y        = 0.98,
    hjust    = 0,
    vjust    = 1,
    size     = 10
  ) +
  draw_plot(
    grid_5x4,
    x      = 0,
    y      = 0,
    width  = 1,
    height = 0.95
  )

# save
ggsave2(
  out_png,
  plot  = combined,
  width = 20,
  height = 25,
  units = "cm",
  dpi   = 300
)

message("Combined QQ plot saved to: ", out_png)

################################################################################
