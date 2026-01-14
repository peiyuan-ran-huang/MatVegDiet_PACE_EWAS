################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 27-Jun-2025
# This script is to screen EWAS meta-analysis results for subsequent analyses.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(meffil)
library(maxprobes)

## --------------------------------------------------------------------------- ##
## 2.  Paths & basic parameters                                                ##
## --------------------------------------------------------------------------- ##

meta.dir  <- "Z:/working/results/EWAS/meta"
out.dir   <- file.path(meta.dir, "")
dir.create(out.dir, showWarnings = FALSE)

zhou.hm450 <- "Z:/working/data/EWAS/HM450.hg38.mask.tsv.gz"
zhou.epic   <- "Z:/working/data/EWAS/EPIC.hg38.mask.tsv.gz"

file.pattern <- "^ewas\\.res\\..*\\.txt$"

## --------------------------------------------------------------------------- ##
## 3.  Define probe universes                                                  ##
## --------------------------------------------------------------------------- ##

cpg.450k <- meffil.get.features("450k")$name
cpg.epic <- meffil.get.features("epic")$name

cpg.overlap  <- intersect(cpg.450k, cpg.epic)
cpg.epiconly <- setdiff(cpg.epic, cpg.450k)

## --------------------------------------------------------------------------- ##
## 4.  Build BAD-probe list                                                   ##
## --------------------------------------------------------------------------- ##

### 4.1  Zhou SNP masking  -------------------------------------------------- ###

read_zhou_mask <- function(fn) {
  if (!file.exists(fn))
    stop("File not found: ", fn)
  hdr <- names(fread(fn, nrows = 0))
  id.col   <- hdr[match(TRUE, tolower(hdr) %in% c("probe_id", "probeid"))]
  if (is.na(id.col))
    stop("Probe-ID column not found in ", basename(fn))
  mask.col <- hdr[match(TRUE, tolower(hdr) %in% c("mask_general", "mask.general"))]
  if (is.na(mask.col))
    stop("MASK_general column not found in ", basename(fn))
  dt <- fread(fn,
              select = c(id.col, mask.col),
              showProgress = FALSE)
  dt[get(mask.col) != 0 & !is.na(get(mask.col)), get(id.col)]
}

cat("Reading Zhou SNP mask …\n")
bad.snp <- union(read_zhou_mask(zhou.hm450), read_zhou_mask(zhou.epic))
cat("  SNP-affected probes (Zhou)    :", length(bad.snp), "\n")

### 4.2  Cross-reactive probes (maxprobes)  --------------------------------- ###

xr_450k       <- maxprobes::xreactive_probes("450K")
xr_epic       <- maxprobes::xreactive_probes("EPIC")
bad.xreactive <- union(xr_450k, xr_epic)
if (!length(bad.xreactive))
  stop("Cross-reactive list is empty! Please check platform args in xreactive_probes().")
cat("  Cross-reactive probes          :", length(bad.xreactive), "\n")

### 4.3  (Double-check) Sex chromosomes & control probes  ------------------------ ###

anno_450k <- meffil.get.features("450k")
anno_epic <- meffil.get.features("epic")

bad.sex.450k <- anno_450k$name[anno_450k$chromosome %in% c("chrX", "chrY")]
bad.sex.epic <- anno_epic$name[anno_epic$chromosome %in% c("chrX", "chrY")]
bad.sex      <- union(bad.sex.450k, bad.sex.epic)

bad.ctrl.450k <- anno_450k$name[!grepl("^cg|^ch", anno_450k$name)]
bad.ctrl.epic <- anno_epic$name[!grepl("^cg|^ch", anno_epic$name)]
bad.ctrl      <- union(bad.ctrl.450k, bad.ctrl.epic)

cat("  Sex-chromosome probes         :", length(bad.sex), "\n")
cat("  Control probes                :", length(bad.ctrl), "\n\n")

# ### 4.4  Naeem & Chen exclusion lists --------------------------------------- ###
#
# naeem_path   <- "Z:/working/data/EWAS/Naeem_list.csv"
# naeem_dt     <- fread(naeem_path, select = c("probe", "Flag(discard/keep)"))
# bad.naeem    <- naeem_dt[`Flag(discard/keep)` == "discard", probe]
#
# chen_path    <- "Z:/working/data/EWAS/48639-non-specific-probes-Illumina450k.csv"
# chen_dt      <- fread(chen_path, header = TRUE)
# bad.chen     <- if ("TargetID" %in% names(chen_dt)) {
#   trimws(chen_dt[["TargetID"]])
# } else {
#   trimws(chen_dt[[1]])
# }
#
# cat("  Naeem discard probes           :", length(bad.naeem), "\n")
# cat("  Chen non-specific probes       :", length(bad.chen), "\n\n")

### 4.5  Combined BAD list --------------------------------------------------- ###

bad.probes <- Reduce(union, list(# bad.naeem, bad.chen,
  bad.snp, bad.xreactive, bad.sex, bad.ctrl))
cat("  TOTAL probes to exclude        :", length(bad.probes), "\n\n")

## --------------------------------------------------------------------------- ##
## 5.  Helper: detect CpG column name                                         ##
## --------------------------------------------------------------------------- ##

find_cpg_col <- function(dt) {
  if ("MarkerName" %in% names(dt))
    return("MarkerName")
  if ("name"       %in% names(dt))
    return("name")
  stop("No CpG-ID column (MarkerName / name) found.")
}

## --------------------------------------------------------------------------- ##
## 6.  Main loop: split every meta-analysis result file                        ##
## --------------------------------------------------------------------------- ##

files <- list.files(meta.dir, pattern = file.pattern, full.names = TRUE)
if (!length(files))
  stop("No *.txt files found in ", meta.dir)

for (f in files) {
  message("Processing: ", basename(f))
  dt     <- fread(f, showProgress = FALSE)
  id.col <- find_cpg_col(dt)
  
  ## --- split first, then clean (per PI request) --------------------------- ##
  both_dt <- dt[get(id.col) %in% cpg.overlap]
  epic_dt <- dt[get(id.col) %in% cpg.epiconly]
  
  ## --- remove BAD probes in each subset ---------------------------------- ##
  r1 <- nrow(both_dt)
  r2 <- nrow(epic_dt)
  both_dt <- both_dt[!get(id.col) %in% bad.probes]
  epic_dt <- epic_dt[!get(id.col) %in% bad.probes]
  rem1 <- r1 - nrow(both_dt)
  rem2 <- r2 - nrow(epic_dt)
  
  ## --- write out ---------------------------------------------------------- ##
  fwrite(both_dt, file = file.path(out.dir, paste0("450K&EPIC_", basename(f))), sep  = "\t")
  fwrite(epic_dt, file = file.path(out.dir, paste0("EPIC.only_", basename(f))), sep  = "\t")
  
  message(
    "  → 450K ∩ EPIC: ",
    nrow(both_dt),
    " rows (",
    rem1,
    " removed) |  ",
    "EPIC-only: ",
    nrow(epic_dt),
    " rows (",
    rem2,
    " removed)."
  )
}

message("\n=== Finished: ", length(files), " files processed ===")

################################################################################
