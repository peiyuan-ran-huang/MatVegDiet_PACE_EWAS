################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 21-Jul-2025
# This script is to: (2) look-up top CpGs or CpGs in top DMRs on the HELIX database for eQTMs.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())
cat("\014")

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(readxl)
  library(stringr)
  library(dplyr)
})

## --------------------------------------------------------------------------- ##
## 2.  Paths and settings                                                      ##
## --------------------------------------------------------------------------- ##

meta.dir       <- "Z:/working/results/EWAS/meta"
dmr.dir        <- file.path(meta.dir, "DMR_dmrff")
helix.file     <- "Z:/working/data/EWAS/eQTM_autosome_adj.cells.txt.gz"
out.full       <- file.path(meta.dir, "HELIX_eQTM_lookup_ALL.xlsx")
out.sig        <- file.path(meta.dir, "HELIX_eQTM_lookup_TOP.xlsx")

exp.list       <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
model.tag      <- "AddModel"

output.cols    <- c(
  "Exposure",
  "Type",
  "CpG",
  "CpG_chr",
  "CpG_pos",
  "TC",
  "Gene",
  "Beta",
  "SE",
  "P",
  "sigPair"
)

## --------------------------------------------------------------------------- ##
## 3.  Read HELIX eQTM catalogue                                               ##
## --------------------------------------------------------------------------- ##

message("Reading HELIX eQTM catalogue …")

helix.dt <- fread(helix.file, showProgress = TRUE) %>%
  rename(
    CpG     = CpG,
    TC      = TC,
    Gene    = TC_gene,
    Beta    = log2FC,
    SE      = SE,
    P       = p.value,
    sigPair = sigPair
  ) %>%
  select(CpG, TC, Gene, Beta, SE, P, CpG_chr, CpG_pos, sigPair)

## --------------------------------------------------------------------------- ##
## 4.  Helper: lookup function                                                 ##
## --------------------------------------------------------------------------- ##

lookup_helix <- function(cpg.vec) {
  data.table(CpG = unique(cpg.vec)) %>%
    left_join(helix.dt, by = "CpG")
}

## --------------------------------------------------------------------------- ##
## 5.  Main loop over exposures                                                ##
## --------------------------------------------------------------------------- ##

wb.all <- createWorkbook()
wb.sig <- createWorkbook()

for (exp in exp.list) {
  message("Processing exposure: ", exp)
  
  ## --- Try 450K&EPIC first, fallback to EPIC.only ---
  f.cpg.1 <- file.path(meta.dir,
                       paste0("450K&EPIC_", exp, "_", model.tag, "_TopTab.xlsx"))
  f.cpg.2 <- file.path(meta.dir,
                       paste0("EPIC.only_", exp, "_", model.tag, "_TopTab.xlsx"))
  f.cpg   <- if (file.exists(f.cpg.1))
    f.cpg.1
  else if (file.exists(f.cpg.2))
    f.cpg.2
  else
    NA
  
  if (is.na(f.cpg)) {
    message("  No CpG TopTab file found for ", exp)
    next
  }
  
  ## --- DMR file: only one version ---
  f.dmr <- file.path(dmr.dir,
                     paste0(exp, "_", model.tag, "_CpG.in.DMR_TopList_FDR.txt"))
  has.dmr <- file.exists(f.dmr) && file.info(f.dmr)$size > 0
  
  ## --- read CpG-level top hits ---
  cpg.vec <- read_excel(f.cpg) %>% pull(CpG) %>% unique()
  dt.cpg <- lookup_helix(cpg.vec) %>%
    mutate(Exposure = exp, Type = "CpG") %>%
    select(all_of(output.cols))
  
  ## --- read CpGs in top DMRs (if any) ---
  if (has.dmr) {
    cpg.dmr <- fread(f.dmr, header = FALSE)[[1]] %>% unique()
    dt.dmr <- lookup_helix(cpg.dmr) %>%
      mutate(Exposure = exp, Type = "DMR") %>%
      select(all_of(output.cols))
    dt.out <- rbind(dt.cpg, dt.dmr)
  } else {
    dt.out <- dt.cpg
  }
  
  ## --- write to both workbooks ---
  addWorksheet(wb.all, sheetName = exp, gridLines = FALSE)
  writeDataTable(wb.all, exp, dt.out)
  
  dt.sig <- dt.out[sigPair == TRUE]
  if (nrow(dt.sig) > 0) {
    addWorksheet(wb.sig, sheetName = exp, gridLines = FALSE)
    writeDataTable(wb.sig, exp, dt.sig)
  }
}

## --------------------------------------------------------------------------- ##
## 6.  Write results to Excel                                                  ##
## --------------------------------------------------------------------------- ##

saveWorkbook(wb.all, out.full, overwrite = TRUE)
saveWorkbook(wb.sig, out.sig, overwrite = TRUE)

message("\nAll done! Results saved to:\n  → ", out.full, "\n  → ", out.sig)


################################################################################
