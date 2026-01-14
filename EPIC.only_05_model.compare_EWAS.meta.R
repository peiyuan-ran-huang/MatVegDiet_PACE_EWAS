################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 29-Jun-2025
# This script is to compare EWAS meta-analysis results across different models (for probes unique to EPIC).

## --------------------------------------------------------------------------- ##
## 0.  House-keeping                                                           ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggforestplot)
  library(ggforce)
  library(glue)
})

## --------------------------------------------------------------------------- ##
## 2.  Paths & parameters                                                      ##
## --------------------------------------------------------------------------- ##

meta_dir <- "Z:/working/results/EWAS/meta"
prefix   <- "EPIC.only"

models <- c("NoCellModel", "MinModel", "FullModel", "AddModel")
model_labels <- c(
  NoCellModel = "No-cell model",
  MinModel    = "Minimally-adjusted model",
  FullModel   = "Partially-adjusted model",
  AddModel    = "Fully-adjusted model"
)

exp_order <- c("PDI", "hPDI", "uPDI")
exp_labels <- c(PDI  = "Overall plant-based diet index (PDI)", hPDI = "Healthful plant-based diet index (hPDI)", uPDI = "Unhealthful plant-based diet index (uPDI)")

plot_dir <- file.path(meta_dir, "ModelCompare")
dir.create(plot_dir, showWarnings = FALSE)

## --------------------------------------------------------------------------- ##
## 3.  Helper – load results for selected CpGs                                 ##
## --------------------------------------------------------------------------- ##

load_model_results <- function(exp, model, cpg_list) {
  fn <- file.path(meta_dir, glue("{prefix}_ewas.res.{model}.{exp}.txt"))
  if (!file.exists(fn))
    stop("File not found: ", fn)
  dt <- fread(fn)
  if (!all(c("MarkerName", "Direction") %in% names(dt)))
    stop("Missing required columns in: ", fn)
  
  dt <- dt[MarkerName %in% cpg_list]
  dt <- dt[!Direction %in% c("???+", "???-")]
  if (nrow(dt) == 0L)
    return(NULL)
  
  out <- dt[, .(
    CpG  = MarkerName,
    beta = Effect,
    se   = StdErr,
    pval = Pvalue
  )]
  out[, Exposure := exp_labels[[exp]]]
  out[, Model    := model_labels[[model]]]
  return(out)
}

## --------------------------------------------------------------------------- ##
## 4.  Loop through FullModel/AddModel to create forest plots                  ##
## --------------------------------------------------------------------------- ##

for (main_model in c("AddModel", "FullModel")) {
  message("\n======= Forest plot for CpGs with p < 1e-5 in ",
          main_model,
          " =======")
  
  all_cpgs <- list()
  for (exp in exp_order) {
    cpg_file <- file.path(meta_dir,
                          glue("{prefix}_{exp}_{main_model}_TopList_1e-5.txt"))
    if (!file.exists(cpg_file))
      next
    cpgs <- fread(cpg_file, header = FALSE)[[1]]
    if (length(cpgs) > 0)
      all_cpgs[[exp]] <- unique(cpgs)
  }
  
  if (length(all_cpgs) == 0L) {
    message("No CpGs found – skipping plot.")
    next
  }
  
  plot_data <- rbindlist(lapply(names(all_cpgs), function(exp) {
    cpgs <- all_cpgs[[exp]]
    rbindlist(
      lapply(models, function(m)
        load_model_results(exp, m, cpgs)),
      use.names = TRUE,
      fill = TRUE
    )
  }),
  use.names = TRUE,
  fill = TRUE)
  
  if (nrow(plot_data) == 0L) {
    message("No eligible CpG after filtering – skipping.")
    next
  }
  
  plot_data[, Model := factor(Model, levels = model_labels)]
  plot_data[, Exposure := factor(Exposure, levels = exp_labels[exp_order])]
  plot_data[, CpG := factor(CpG)]
  
  p <- ggforestplot::forestplot(
    df       = plot_data,
    name     = CpG,
    estimate = beta,
    se       = se,
    pvalue   = pval,
    psignif  = 0.05,
    shape    = Model,
    colour   = Model,
    xlab     = "Beta and 95% CI",
    logodds  = FALSE,
    title    = glue(
      "CpGs (unique to EPIC) with p < 1e-5 in {model_labels[[main_model]]}"
    )
  ) +
    ggforce::facet_col(facets = ~ Exposure,
                       scales = "free_y",
                       space  = "free") +
    scale_colour_manual(values = c("darkgreen", "darkorange", "navyblue", "darkred")) +
    scale_shape_manual(values = c(21, 21, 21, 21)) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(size = 12))
  
  ggsave(
    filename = file.path(
      plot_dir,
      glue(
        "{prefix}_ModelCompare_{main_model}_TopHits_forestplot.png"
      )
    ),
    plot   = p,
    width  = 14,
    height = 4 * length(unique(plot_data$Exposure)),
    dpi    = 300
  )
  
  message("→ Forest plot saved for ", main_model)
}

message("\nAll forest plots completed.")

################################################################################
