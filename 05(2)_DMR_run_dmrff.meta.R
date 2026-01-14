################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 12-Jul-2025
# This script is to perform differentially methylated region (DMR) meta-analysis with dmrff.meta.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(dmrff)
library(data.table)
library(meffil)

## --------------------------------------------------------------------------- ##
## 2.  Parameters                                                              ##
## --------------------------------------------------------------------------- ##

exposures <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
cohorts   <- c("ALSPAC",
               "GenR",
               "INMA",
               "MoBa1",
               "MoBa2",
               "MoBa4",
               "MoBa8",
               "NORTHPOP",
               "Viva")
platforms <- list("450K" = meffil.get.features("450k")$name,
                  "EPIC" = meffil.get.features("epic")$name)

dmr.dir   <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/results"
out.dir   <- file.path(dmr.dir, "meta/DMR_dmrff")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

## --------------------------------------------------------------------------- ##
## 3.  Main loop: exposure × AddModel                                          ##
## --------------------------------------------------------------------------- ##

for (exp in exposures) {
  tag <- paste(exp, "AddModel", sep = ".")
  cat("\nProcessing:", tag, "...\n")
  
  cohort.files <- list()
  cohort.names <- c()
  
  for (cohort in cohorts) {
    path <- file.path(dmr.dir, cohort, "DMR", paste0(tag, ".dmrff.pre.rds"))
    if (!file.exists(path))
      next
    
    obj <- readRDS(path)
    cohort.files[[cohort]] <- obj
    cohort.names <- c(cohort.names, cohort)
  }
  
  if (length(cohort.files) < 2) {
    warning("Not enough cohorts for meta-analysis of ", tag)
    next
  }
  
  meta.res <- dmrff.meta(cohort.files)
  saveRDS(meta.res, file = file.path(out.dir, paste0(tag, "_dmrff.meta.rds")))
  
  cat(
    "  → Done:",
    length(cohort.files),
    "cohorts included; ",
    nrow(meta.res),
    "DMRs found.\n"
  )
}

message("\n=== All meta-analyses complete ===")

################################################################################
