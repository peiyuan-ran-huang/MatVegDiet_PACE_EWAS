################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 29-Jun-2025
# This script is to select top hits of EWAS meta-analysis results (for probes available on both 450K and EPIC).
## (1) For fully-adjusted models (i.e., AddModel)

## --------------------------------------------------------------------------- ##
## 0.  Clear workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(dplyr)
library(glue)
library(openxlsx)
library(meffil)

## --------------------------------------------------------------------------- ##
## 2.  Directories & parameters                                                ##
## --------------------------------------------------------------------------- ##

meta_dir  <- "Z:/working/results/EWAS/meta"
out_dir   <- file.path(meta_dir, "")
dir.create(out_dir, showWarnings = FALSE)

prefix    <- "450K&EPIC"
model     <- "AddModel"
exp_list  <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")

## --------------------------------------------------------------------------- ##
## 3.  Load CpG annotation via meffil                                          ##
## --------------------------------------------------------------------------- ##

feat450k <- meffil.get.features("450k")
featEpic <- meffil.get.features("epic")

# pull in annotation from meffil for EPIC array
anno_dt <- as.data.table(meffil.get.features("epic"))[, .(
  CpG                = name,
  Chr                = sub("^chr", "", chromosome),
  Pos                = position,
  UCSC_Gene          = gene.symbol,
  UCSC_Group         = gene.region,
  Relation_to_Island = relation.to.island
)]

# normalize Chr to character, drop any "chr" prefix
anno_dt[, Chr := as.character(sub("^chr", "", Chr))]

## --------------------------------------------------------------------------- ##
## 4.  Extract and save top lists                                              ##
## --------------------------------------------------------------------------- ##

naeem_path   <- "Z:/working/data/EWAS/Naeem_list.csv"
naeem_dt     <- fread(naeem_path, select = c("probe", "Flag(discard/keep)"))
naeem_exclude <- naeem_dt[`Flag(discard/keep)` == "discard", probe]

chen_path   <- "Z:/working/data/EWAS/48639-non-specific-probes-Illumina450k.csv"
chen_dt     <- fread(chen_path, header = TRUE)
chen_exclude <- if ("TargetID" %in% names(chen_dt)) {
  trimws(chen_dt[["TargetID"]])
} else {
  trimws(chen_dt[[1]])
}

snp_exclude <- featEpic$name[featEpic$snp.exclude]
snp_exclude <- union(snp_exclude, chen_exclude)

# --------------------------------------------------------------------------- ##

for (exp in exp_list) {
  message("Processing exposure: ", exp, " ...")
  
  # 4.1  Read AddModel results
  fn <- file.path(meta_dir, glue("{prefix}_ewas.res.{model}.{exp}.txt"))
  dt <- fread(fn)
  
  # 4.2  Harmonise column names
  setnames(dt, "MarkerName", "CpG")
  setnames(dt, "StdErr", "SE")
  setnames(dt, "Pvalue", "P")
  setnames(dt, "HetISq", "I2_orig")
  # recalculate I2 from HetChiSq and HetDf if present
  if (all(c("HetChiSq", "HetDf") %in% names(dt))) {
    dt[, I2 := pmax(0, (HetChiSq - HetDf) / HetChiSq) * 100]
  } else {
    dt[, I2 := I2_orig]  # fallback
  }
  
  # 4.3  Merge in annotation
  dt <- merge(dt, anno_dt, by = "CpG", all.x = TRUE)
  
  # 4.4  Compute FDR
  dt[, FDR := p.adjust(P, "fdr")]
  
  # 4.5  p < 1e-5 CpG list
  list1_fn <- file.path(out_dir, glue("{prefix}_{exp}_{model}_TopList_1e-5.txt"))
  writeLines(dt[P < 1e-5, CpG], list1_fn)
  
  # 4.6  FDR < 0.05 CpG list
  list2_fn <- file.path(out_dir, glue("{prefix}_{exp}_{model}_TopList_FDR.txt"))
  writeLines(dt[FDR < 0.05, CpG], list2_fn)
  
  # 4.7  Detailed table for p < 1e-5
  tab <- dt[P < 1e-5][order(P), .(
    CpG,
    Chr,
    Pos,
    UCSC_Gene,
    UCSC_Group,
    Relation_to_Island,
    Effect = sprintf("%.7f", Effect),
    SE     = sprintf("%.7f", SE),
    Direction,
    I2     = paste0(sprintf("%.1f", I2), "%"),
    P      = format(P, scientific = TRUE, digits = 3),
    FDR    = sprintf("%.3f", FDR),
    Flag_Naeem        = CpG %in% naeem_exclude,
    Flag_Chen         = CpG %in% chen_exclude,
    Flag_SNP          = CpG %in% snp_exclude
  )]
  
  # 4.8  Write Excel
  tab_fn <- file.path(out_dir, glue("{prefix}_{exp}_{model}_TopTab.xlsx"))
  write.xlsx(as.data.frame(tab), tab_fn, rowNames = FALSE)
}

################################################################################
