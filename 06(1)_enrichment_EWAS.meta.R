################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 16-Jul-2025
# This script is to perform GO and KEGG enrichment analyses for: (1) top-hit CpGs.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(meffil)
library(missMethyl)
library(openxlsx)

## --------------------------------------------------------------------------- ##
## 2.  Paths & exposure settings                                               ##
## --------------------------------------------------------------------------- ##

meta.dir <- "Z:/working/results/EWAS/meta"
out.dir  <- file.path(meta.dir, "Enrichment")
dir.create(out.dir, showWarnings = FALSE)

exp.list <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
model.tag  <- "AddModel"

## --------------------------------------------------------------------------- ##
## 3.  Helper functions                                                        ##
## --------------------------------------------------------------------------- ##

load_top_cpgs <- function(exp) {
  f1 <- file.path(meta.dir,
                  paste0("450K&EPIC_", exp, "_AddModel_TopList_1e-5.txt"))
  f2 <- file.path(meta.dir,
                  paste0("EPIC.only_", exp, "_AddModel_TopList_1e-5.txt"))
  ids <- character()
  for (f in c(f1, f2)) {
    if (file.exists(f) && file.info(f)$size > 0) {
      ids <- union(ids, fread(f, header = FALSE, showProgress = FALSE)[[1]])
    }
  }
  return(ids)
}

load_all_cpgs <- function(exp) {
  f1 <- file.path(meta.dir,
                  paste0("450K&EPIC_ewas.res.AddModel.", exp, ".txt"))
  f2 <- file.path(meta.dir,
                  paste0("EPIC.only_ewas.res.AddModel.", exp, ".txt"))
  ids <- character()
  for (f in c(f1, f2)) {
    if (file.exists(f)) {
      dt <- fread(f, select = "MarkerName", showProgress = FALSE)
      ids <- union(ids, dt[["MarkerName"]])
    }
  }
  return(ids)
}

## --------------------------------------------------------------------------- ##
## 4.  Main loop: GO and KEGG enrichment for top‐hit CpGs                     ##
## --------------------------------------------------------------------------- ##

for (exp in exp.list) {
  message("\n=== Processing exposure:", exp, "===\n")
  
  # load sig and background CpGs (as before) …
  sig.cpg <- load_top_cpgs(exp)
  if (length(sig.cpg) == 0) {
    message("  No significant CpGs for ", exp, "; skipped.")
    next
  }
  all.cpg <- load_all_cpgs(exp)
  
  # run GO
  go.res <- tryCatch(
    gometh(
      sig.cpg    = sig.cpg,
      all.cpg    = all.cpg,
      collection = "GO",
      array.type = "EPIC"
    ),
    error = function(e) {
      warning("  GO enrichment error for ", exp, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(go.res)) {
    fwrite(go.res, file = file.path(out.dir, paste0(
      exp, ".", model.tag, ".CpG.GO.full.csv"
    )))
    message("  GO full results: ", exp, ".", model.tag, ".CpG.GO.full.csv")
  }
  
  # run KEGG
  kegg.res <- tryCatch(
    gometh(
      sig.cpg    = sig.cpg,
      all.cpg    = all.cpg,
      collection = "KEGG",
      array.type = "EPIC"
    ),
    error = function(e) {
      warning("  KEGG enrichment error for ", exp, ": ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(kegg.res)) {
    fwrite(kegg.res, file = file.path(out.dir, paste0(
      exp, ".", model.tag, ".CpG.KEGG.full.csv"
    )))
    message("  KEGG full results: ",
            exp,
            ".",
            model.tag,
            ".CpG.KEGG.full.csv")
  }
}

message("\n=== CpG‐based enrichment complete for all exposures ===")

################################################################################
