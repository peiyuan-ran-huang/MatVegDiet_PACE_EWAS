################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 29-Jun-2025
# This script is to select top hits of EWAS meta-analysis results (for probes available on both 450K and EPIC).
## (2) For partially-adjusted models (i.e., FullModel)

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
model     <- "FullModel"
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

for (exp in exp_list) {
  message("Processing exposure: ", exp, " ...")
  
  # 4.1  Read FullModel results
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
  
  # 4.5  Write CpG lists
  writeLines(dt[P < 1e-5 , CpG], file.path(out_dir, glue(
    "{prefix}_{exp}_{model}_TopList_1e-5.txt"
  )))
  
  writeLines(dt[FDR < 0.05, CpG], file.path(out_dir, glue(
    "{prefix}_{exp}_{model}_TopList_FDR.txt"
  )))
}

################################################################################
