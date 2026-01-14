################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 21-Jul-2025
# This script is to: (1) look-up top CpGs or CpGs in top DMRs on EWAS Catalog (https://ewascatalog.org/).

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(openxlsx)
library(ewascatalog)

## --------------------------------------------------------------------------- ##
## 2.  Paths & exposure settings                                               ##
## --------------------------------------------------------------------------- ##

meta.dir   <- "Z:/working/results/EWAS/meta"
dmr.dir    <- file.path(meta.dir, "DMR_dmrff")
out.dir    <- file.path(meta.dir, "")
dir.create(out.dir, showWarnings = FALSE)

exp.list   <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
model.tag  <- "AddModel"

## --------------------------------------------------------------------------- ##
## 3.  Helper functions to load CpG IDs                                        ##
## --------------------------------------------------------------------------- ##

# Pretty name converter for exposure
pretty.exposure <- function(tag) {
  switch(
    tag,
    "veggie2" = "Pesco-/full vs. non-vegetarian",
    "veggie1" = "Full vs. non-vegetarian",
    "PDI"     = "Overall plant-based diet index (PDI)",
    "hPDI"    = "Healthful plant-based diet index (hPDI)",
    "uPDI"    = "Unhealthful plant-based diet index (uPDI)",
    tag
  )
}

# Load top CpG-level hits
load_top_cpgs <- function(exp) {
  f1 <- file.path(meta.dir,
                  paste0("450K&EPIC_", exp, "_", model.tag, "_TopList_1e-5.txt"))
  f2 <- file.path(meta.dir,
                  paste0("EPIC.only_", exp, "_", model.tag, "_TopList_1e-5.txt"))
  ids <- character()
  for (f in c(f1, f2)) {
    if (file.exists(f) && file.info(f)$size > 0) {
      ids <- union(ids, fread(f, header = FALSE, showProgress = FALSE)[[1]])
    }
  }
  return(data.table(
    CpG_ID   = ids,
    Exposure = pretty.exposure(exp),
    Type     = "Top CpG"
  ))
}

# Load CpGs located in top DMRs
load_dmr_cpgs <- function(exp) {
  f <- file.path(dmr.dir,
                 paste0(exp, "_", model.tag, "_CpG.in.DMR_TopList_FDR.txt"))
  if (!file.exists(f) || file.info(f)$size == 0) {
    return(data.table())
  }
  ids <- fread(f, header = FALSE, showProgress = FALSE)[[1]]
  return(data.table(
    CpG_ID   = ids,
    Exposure = pretty.exposure(exp),
    Type     = "CpG in top DMR"
  ))
}

## --------------------------------------------------------------------------- ##
## 4.  Run EWAS Catalog lookup                                                ##
## --------------------------------------------------------------------------- ##

message("\n=== Performing EWAS Catalog lookup ===")

# Combine all CpGs to be queried
cpg.table <- rbindlist(lapply(exp.list, function(e) {
  rbind(load_top_cpgs(e), load_dmr_cpgs(e), fill = TRUE)
}), use.names = TRUE)

# Keep unique CpG–Exposure–Type combinations
cpg.table <- unique(cpg.table, by = c("CpG_ID", "Exposure", "Type"))

# Lookup via EWAS Catalog
res.list <- list()
for (i in seq_len(nrow(cpg.table))) {
  cg <- cpg.table$CpG_ID[i]
  msg <- paste0("[", i, "/", nrow(cpg.table), "] Querying: ", cg)
  message(msg)
  hits <- tryCatch(
    ewascatalog::ewascatalog(cg),
    error = function(e)
      NULL
  )
  if (!is.null(hits) && nrow(hits) > 0) {
    temp <- cbind(cpg.table[i], hits)
    res.list[[length(res.list) + 1]] <- temp
  }
}

# Merge and export results if any matched
if (length(res.list) == 0) {
  message("No CpGs found in EWAS Catalog.")
} else {
  result.dt <- rbindlist(res.list, use.names = TRUE, fill = TRUE)
  
  # Ensure column names are unique and Excel-compatible
  names(result.dt) <- tolower(names(result.dt))
  names(result.dt) <- make.unique(names(result.dt), sep = "_")
  names(result.dt) <- make.names(names(result.dt), unique = TRUE)
  
  if (any(duplicated(names(result.dt)))) {
    dup.cols <- names(result.dt)[duplicated(names(result.dt))]
    stop("Duplicate column names still exist: ",
         paste(dup.cols, collapse = ", "))
  }
  
  # Write output
  out.xlsx <- file.path(out.dir, "EWASCatalog_Lookup_TopCpG_&_TopDMR.xlsx")
  write.xlsx(result.dt,
             file     = out.xlsx,
             asTable  = TRUE,
             keepNA   = TRUE)
  message("\n=== Lookup complete. Results saved to: ", out.xlsx, " ===")
}

################################################################################
