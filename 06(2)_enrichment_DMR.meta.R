################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 20-Jul-2025
# This script is to perform GO and KEGG enrichment analyses for: (2) top-hit DMRs.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)    # fast I/O
library(readxl)        # read_excel()
library(GenomicRanges) # findOverlaps()
library(missMethyl)    # gometh()
library(ggplot2)       # plotting

## --------------------------------------------------------------------------- ##
## 2.  Paths & parameters                                                      ##
## --------------------------------------------------------------------------- ##

meta.dir   <- "Z:/working/results/EWAS/meta"
base.dir   <- file.path(meta.dir, "DMR_dmrff")
enrich.dir <- file.path(meta.dir, "Enrichment")
dir.create(enrich.dir, recursive = TRUE, showWarnings = FALSE)

exposures <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
model.tag <- "AddModel"

## --------------------------------------------------------------------------- ##
## 3.  Helper: gometh wrapper                                                  ##
## --------------------------------------------------------------------------- ##

run_enrichment <- function(sig, bg, coll) {
  if (length(sig) == 0 || !any(sig %in% bg)) {
    warning("Skipping ", coll, ": no overlap between sig- and bg-CpGs.")
    return(data.frame())
  }
  gometh(
    sig.cpg   = sig,
    all.cpg   = bg,
    array.type = "EPIC",
    collection = coll
  )
}

## --------------------------------------------------------------------------- ##
## 4.  Main loop: iterate exposures                                            ##
## --------------------------------------------------------------------------- ##

for (exp in exposures) {
  cat("\n=== Processing", exp, "===\n")
  
  ## ------------------------------ 4.1  TopTab ------------------------------ ##
  top.file <- file.path(base.dir, paste0(exp, ".", model.tag, "_DMR_TopTab.xlsx"))
  if (!file.exists(top.file)) {
    warning("TopTab missing for ", exp)
    next
  }
  top.dt <- as.data.table(read_excel(top.file))
  if (nrow(top.dt) == 0) {
    warning("No significant DMRs for ", exp)
    next
  }
  setnames(top.dt, c("Chr", "Start", "End"), c("chr", "start", "end"))
  top.dt[, chr := sub("^chr", "", tolower(chr))]
  
  ## ----------------------------- 4.2  meta.rds ----------------------------- ##
  meta.rds <- file.path(base.dir, paste0(exp, ".", model.tag, "_dmrff.meta.rds"))
  if (!file.exists(meta.rds)) {
    warning("Meta RDS missing for ", exp)
    next
  }
  meta.lst <- readRDS(meta.rds)
  
  ## ----------------------------- 4.3  EWAS CpGs ---------------------------- ##
  ewas.dt <- as.data.table(meta.lst$ewas, keep.rownames = "cgID")
  setnames(ewas.dt, c("chr", "pos"), c("chr", "pos"))
  ewas.dt[, chr := sub("^chr", "", tolower(chr))]
  
  ## ----------------------------- 4.4  All DMRs ----------------------------- ##
  dmrs.dt <- as.data.table(meta.lst$dmrs)
  setnames(dmrs.dt, c("chr", "start", "end"), c("chr", "start", "end"))
  dmrs.dt[, chr := sub("^chr", "", tolower(chr))]
  
  ## ----------------------------- 4.5  GRanges ------------------------------ ##
  ## significant DMRs (signal)
  gr.dmr.sig <- GRanges(top.dt$chr, IRanges(top.dt$start, top.dt$end))
  ## all DMRs (background)
  gr.dmr.all <- GRanges(dmrs.dt$chr, IRanges(dmrs.dt$start, dmrs.dt$end))
  ## all CpGs tested in EWAS
  gr.cpgs <- GRanges(ewas.dt$chr, IRanges(ewas.dt$pos, ewas.dt$pos), cgID = ewas.dt$cgID)
  
  ## ----------------------------- 4.6  Overlaps ----------------------------- ##
  ## signal CpGs  = unique CpGs inside TopTab DMRs
  hits.sig <- findOverlaps(gr.dmr.sig, gr.cpgs)
  if (length(hits.sig) == 0) {
    warning("No CpG in significant DMRs for ", exp)
    next
  }
  sigCpGs <- unique(gr.cpgs$cgID[subjectHits(hits.sig)])
  
  ## background CpGs = union of CpGs inside ANY meta-analysis DMR
  hits.bg <- findOverlaps(gr.dmr.all, gr.cpgs)
  if (length(hits.bg) == 0) {
    warning("No CpG in background DMRs for ", exp)
    next
  }
  bgCpGs  <- unique(gr.cpgs$cgID[subjectHits(hits.bg)])
  
  ## ----------------------------- 4.7  Save list ----------------------------- ##
  sigfile <- file.path(base.dir,
                       paste0(exp, "_", model.tag, "_CpG.in.DMR_TopList_FDR.txt"))
  fwrite(data.table(sigCpGs),
         sigfile,
         col.names = FALSE,
         quote     = FALSE)
  message("  Saved DMR CpG list to: ", sigfile)
  
  ## ----------------------------- 4.8  Enrichment --------------------------- ##
  go.res   <- run_enrichment(sigCpGs, bgCpGs, "GO")
  kegg.res <- run_enrichment(sigCpGs, bgCpGs, "KEGG")
  
  fwrite(go.res, file.path(enrich.dir, paste0(exp, ".", model.tag, ".DMR.GO.full.csv")))
  fwrite(kegg.res, file.path(
    enrich.dir,
    paste0(exp, ".", model.tag, ".DMR.KEGG.full.csv")
  ))
}

message("\n=== DMR enrichment complete for all exposures ===")

################################################################################
