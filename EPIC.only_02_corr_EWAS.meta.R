################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 29-Jun-2025
# This script is to generate correlation matrices across dietary exposures and models for EWAS meta-analysis (for probes available on both 450K and EPIC).

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
library(stringr)
library(glue)
library(ggplot2)
library(tidyr)

## --------------------------------------------------------------------------- ##
## 2.  Directories / constants                                                ##
## --------------------------------------------------------------------------- ##

meta_dir   <- "Z:/working/results/EWAS/meta"
prefix     <- "EPIC.only"

models       <- c("NoCellModel", "MinModel", "FullModel", "AddModel")
model_labels <- c(
  NoCellModel = "No-cell model",
  MinModel    = "Minimally-adjusted model",
  FullModel   = "Partially-adjusted model",
  AddModel    = "Fully-adjusted model"
)

exp_order  <- c("PDI", "hPDI", "uPDI")
exp_labels <- c(PDI  = "Overall PDI", hPDI = "Healthful PDI", uPDI = "Unhealthful PDI")

## --------------------------------------------------------------------------- ##
## 3.  Read meta files & store β vectors                                       ##
## --------------------------------------------------------------------------- ##

get_beta <- function(f) {
  dt <- fread(f)
  # detect CpG column
  id_cands <- c("CpG", "MarkerName", "probeid", "ProbeID", "name")
  id_col   <- id_cands[id_cands %in% names(dt)][1]
  if (is.na(id_col))
    stop("CpG column not found in ", f)
  # detect beta column
  beta_cands <- c("Effect", "Beta", "coef")
  beta_col   <- beta_cands[beta_cands %in% names(dt)][1]
  if (is.na(beta_col))
    stop("Beta column not found in ", f)
  # return CpG and beta
  dt[, .(CpG = get(id_col), beta = get(beta_col))]
}

beta_list <- list()
for (m in models) {
  for (e in exp_order) {
    f <- file.path(meta_dir, glue("{prefix}_ewas.res.{m}.{e}.txt"))
    if (file.exists(f)) {
      beta_list[[paste(e, m, sep = "_")]] <- get_beta(f)
    }
  }
}
if (length(beta_list) < 2)
  stop("Not enough EPIC-only meta files found!")

## --------------------------------------------------------------------------- ##
## 4.  Compute pairwise correlations                                           ##
## --------------------------------------------------------------------------- ##

set_names <- names(beta_list)
n_sets    <- length(set_names)

cors <- matrix(
  NA,
  nrow = n_sets,
  ncol = n_sets,
  dimnames = list(set_names, set_names)
)

for (i in seq_len(n_sets)) {
  for (j in seq_len(n_sets)) {
    dt1 <- beta_list[[set_names[i]]]
    dt2 <- beta_list[[set_names[j]]]
    common <- intersect(dt1$CpG, dt2$CpG)
    v1 <- dt1$beta[match(common, dt1$CpG)]
    v2 <- dt2$beta[match(common, dt2$CpG)]
    cors[i, j] <- cor(v1, v2, use = "complete.obs")
  }
}

## --------------------------------------------------------------------------- ##
## 5.  Prepare data for plotting                                               ##
## --------------------------------------------------------------------------- ##

pretty_names <- imap_chr(beta_list, ~ {
  parts <- str_split(.y, "_")[[1]]
  glue("{exp_labels[[parts[1]]]} – {model_labels[[parts[2]]]}")
})
rownames(cors) <- colnames(cors) <- pretty_names

df_plot <- as.data.frame(as.table(cors), stringsAsFactors = FALSE)
colnames(df_plot) <- c("Var1", "Var2", "Corr")

df_plot$Var1 <- factor(df_plot$Var1, levels = pretty_names)
df_plot$Var2 <- factor(df_plot$Var2, levels = rev(pretty_names))
df_plot <- df_plot %>%
  mutate(i = as.integer(Var1), j = as.integer(Var2)) %>%
  filter(i + j <= n_sets + 1) %>%
  dplyr::select(-i, -j)

## --------------------------------------------------------------------------- ##
## 6.  Plot correlation heat-map                                               ##
## --------------------------------------------------------------------------- ##

p <- ggplot(df_plot, aes(Var1, Var2, fill = Corr)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = sprintf("%.2f", Corr)), size = 4) +
  scale_fill_gradient2(
    limits   = c(-1, 1),
    midpoint = 0,
    low      = "#2166ac",
    mid      = "white",
    high     = "#b2182b",
    name     = "Pearson's r"
  ) +
  scale_x_discrete(position = "bottom", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 315,
      hjust = 0,
      vjust = 0.5,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )

ggsave(
  file.path(meta_dir, glue("{prefix}_CorrMatrix_meta.png")),
  p,
  width  = 30,
  height = 30,
  units  = "cm",
  dpi    = 300
)

message("Correlation heat-map saved to: ", file.path(meta_dir, glue("{prefix}_CorrMatrix_meta.png")))

################################################################################
