################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 21-Jul-2025
# This script is to: (3) select top GO and KEGG enrichment results for top CpGs and DMRs.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)    # fread/fwrite, rbindlist
library(openxlsx)      # createWorkbook, writeDataTable
library(stringr)       # str_split_fixed

## --------------------------------------------------------------------------- ##
## 2.  Paths & settings                                                        ##
## --------------------------------------------------------------------------- ##

meta.dir  <- "Z:/working/results/EWAS/meta"
enrich.dir <- file.path(meta.dir, "Enrichment")
out.file  <- file.path(enrich.dir, "ALL.exp_AddModel_enrichment.xlsx")

exposures <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
types     <- c("CpG", "DMR")
colls     <- c("GO", "KEGG")

## --------------------------------------------------------------------------- ##
## 3.  Read, tag, select top‐5                                                 ##
## --------------------------------------------------------------------------- ##

results <- list()

for (exp in exposures) {
  for (t in types) {
    for (coll in colls) {
      # construct file name
      fname <- file.path(enrich.dir,
                         paste0(exp, ".AddModel.", t, ".", coll, ".full.csv"))
      
      if (!file.exists(fname)) {
        warning("File not found, skipping: ", fname)
        next
      }
      
      # read the full enrichment table
      dt <- fread(fname, showProgress = FALSE)
      if (nrow(dt) == 0) {
        warning("Empty file, skipping: ", fname)
        next
      }
      
      # add Exposure and Type
      dt[, Exposure := exp]
      dt[, Type     := t]
      
      # select/rename columns
      if (coll == "GO") {
        sel <- dt[, .(
          Ontology         = ONTOLOGY,
          Pathway          = TERM,
          N_DMC            = DE,
          N_all_CpGs       = N,
          Enrichment_P     = `P.DE`,
          FDR_P            = FDR
        )]
      } else {
        sel <- dt[, .(
          Ontology         = "KEGG",
          Pathway          = Description,
          N_DMC            = DE,
          N_all_CpGs       = N,
          Enrichment_P     = `P.DE`,
          FDR_P            = FDR
        )]
      }
      
      # attach exposure and type
      sel[, `:=`(Exposure = exp,
                 Type = t,
                 Collection = coll)]
      
      # order by Enrichment_P and take top 5
      setorder(sel, Enrichment_P)
      top5 <- sel[1:min(5, .N)]
      
      # collect
      results[[length(results) + 1]] <- top5
    }
  }
}

# bind all
all.dt <- rbindlist(results, use.names = TRUE, fill = TRUE)

## --------------------------------------------------------------------------- ##
## 4.  Write to Excel                                                          ##
## --------------------------------------------------------------------------- ##

wb <- createWorkbook()
addWorksheet(wb, "Top5_GO_KEGG")
writeDataTable(wb, sheet = 1, x = all.dt)
saveWorkbook(wb, out.file, overwrite = TRUE)

message("Compiled enrichment results saved to:\n  ", out.file)

################################################################################
