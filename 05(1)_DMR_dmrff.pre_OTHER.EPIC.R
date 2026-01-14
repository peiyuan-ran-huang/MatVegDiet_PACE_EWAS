################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 11-Jul-2025
# This script is to perform differentially methylated region (DMR) analysis in each cohort.

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
library(dmrff)

## --------------------------------------------------------------------------- ##
## 2.  Parameters                                                              ##
## --------------------------------------------------------------------------- ##

studies        <- c("MoBa4", "MoBa8", "NORTHPOP")
exposures      <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
models         <- c("NoCellModel", "MinModel", "FullModel", "AddModel")

meth.file      <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/data/GSE154915/beta_matrix_GSE154915.rds"
ewas.base.dir  <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/results"
ewas.prefix    <- "ewas.res"

## --------------------------------------------------------------------------- ##
## 3.  Load GSE154915 reference methylation matrix                                ##
## --------------------------------------------------------------------------- ##

cat("Reading GSE154915 reference methylation matrix …\n")
norm.beta.random <- readRDS(meth.file)
meth <- as.matrix(norm.beta.random)
rm(norm.beta.random)

## --------------------------------------------------------------------------- ##
## 4.  Identify BAD probes                                                     ##
## --------------------------------------------------------------------------- ##

cat("Identifying bad probes …\n")

zhou.hm450 <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/data/resources/HM450.hg38.mask.tsv.gz"
zhou.epic  <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/data/resources/EPIC.hg38.mask.tsv.gz"

read_zhou_mask <- function(fn) {
  hdr <- names(fread(fn, nrows = 0))
  id.col   <- hdr[match(TRUE, tolower(hdr) %in% c("probe_id", "probe.id", "probeid"))]
  mask.col <- hdr[match(TRUE, tolower(hdr) %in% c("mask_general", "mask.general"))]
  dt <- fread(fn,
              select = c(id.col, mask.col),
              showProgress = FALSE)
  dt[get(mask.col) != 0 & !is.na(get(mask.col)), get(id.col)]
}

bad.snp <- union(read_zhou_mask(zhou.hm450), read_zhou_mask(zhou.epic))
xr_450k <- maxprobes::xreactive_probes("450K")
xr_epic <- maxprobes::xreactive_probes("EPIC")
bad.xreactive <- union(xr_450k, xr_epic)

anno_450k <- meffil.get.features("450k")
anno_epic <- meffil.get.features("epic")

bad.sex <- union(anno_450k$name[anno_450k$chromosome %in% c("chrX", "chrY")], anno_epic$name[anno_epic$chromosome %in% c("chrX", "chrY")])

bad.ctrl <- union(anno_450k$name[!grepl("^cg|^ch", anno_450k$name)], anno_epic$name[!grepl("^cg|^ch", anno_epic$name)])

bad.probes <- Reduce(union, list(bad.snp, bad.xreactive, bad.sex, bad.ctrl))

cat("  → A total of ", length(bad.probes), " bad probes identified.\n")

## --------------------------------------------------------------------------- ##
## 5.  Loop through studies × exposures × models                               ##
## --------------------------------------------------------------------------- ##

for (study in studies) {
  ewas.dir <- file.path(ewas.base.dir, study)
  out.dir  <- file.path(ewas.base.dir, study, "DMR")
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  
  for (exp in exposures) {
    for (model in models) {
      tag <- paste(exp, model, sep = ".")
      cat("\n[", study, "] Processing:", tag, "…\n")
      
      ewas.file <- file.path(ewas.dir,
                             sprintf("%s.%s.%s.%s.csv", study, ewas.prefix, model, exp))
      if (!file.exists(ewas.file)) {
        warning("EWAS file not found: ", ewas.file)
        next
      }
      
      ewas <- fread(ewas.file)
      ewas <- ewas[!(probeid %in% bad.probes)]
      
      anno.subset <- rbind(anno_450k, anno_epic)
      anno.sub <- anno.subset[match(ewas$probeid, anno.subset$name), ]
      ewas$chr <- anno.sub$chromosome
      ewas$pos <- anno.sub$position
      
      common.cpg <- intersect(rownames(meth), ewas$probeid)
      meth.sub   <- meth[common.cpg, , drop = FALSE]
      ewas.sub   <- ewas[match(common.cpg, ewas$probeid)]
      
      stopifnot(identical(rownames(meth.sub), ewas.sub$probeid))
      
      cat("  → Matched ", nrow(ewas.sub), " probes.\n")
      
      pre.obj <- dmrff.pre(
        estimate    = ewas.sub$coef,
        se          = ewas.sub$se,
        methylation = meth.sub,
        chr         = ewas.sub$chr,
        pos         = ewas.sub$pos
      )
      
      saveRDS(pre.obj, file = file.path(out.dir, sprintf("%s.dmrff.pre.rds", tag)))
      
      cat("  → Saved:", file.path(out.dir, sprintf("%s.dmrff.pre.rds", tag)), "\n")
    }
  }
}

################################################################################
