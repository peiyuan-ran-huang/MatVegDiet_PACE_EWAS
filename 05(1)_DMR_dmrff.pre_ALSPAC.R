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

study    <- "ALSPAC"
exposures      <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
models         <- c("NoCellModel", "MinModel", "FullModel", "AddModel")

pheno.file     <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/results/ALSPAC/MatVegDiet.ALSPAC.pheno.SV.birth.20250709.rds"
meth.file      <- "/mnt/storage/private/alspac/1/alspac/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj"
ewas.dir       <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/results/ALSPAC/"
ewas.prefix    <- "ewas.res"
platform       <- "450k"

out.dir        <- file.path("/user/work/zd20208/MatVegDiet_PACE_EWAS/results/ALSPAC/DMR")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

## --------------------------------------------------------------------------- ##
## 3.  Convert EWAS Rdata to CSV (one file per exposure × model)               ##
## --------------------------------------------------------------------------- ##

rdata_files <- list.files(path = ewas.dir,
                          pattern = "\\.EWASres\\.birth\\..+\\.Rdata$",
                          full.names = TRUE)

for (rfile in rdata_files) {
  fname     <- basename(rfile)
  parts     <- strsplit(fname, "\\.")[[1]]
  study     <- parts[2]
  exposure  <- parts[3]
  env       <- new.env()
  load(rfile, envir = env)
  
  for (obj_name in ls(env)) {
    obj <- env[[obj_name]]
    if (!is.data.frame(obj))
      next
    
    obj_parts <- strsplit(obj_name, "\\.")[[1]]
    model     <- obj_parts[3]
    exp2      <- obj_parts[4]
    if (exp2 != exposure)
      next
    
    out_name  <- sprintf("%s.%s.%s.%s.csv", study, ewas.prefix, model, exp2)
    obj_out   <- data.frame(
      probeid = rownames(obj),
      obj,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    write.csv(
      obj_out,
      file = file.path(ewas.dir, out_name),
      row.names = FALSE,
      quote = FALSE
    )
  }
}

## --------------------------------------------------------------------------- ##
## 4.  Load and prepare methylation matrix                                     ##
## --------------------------------------------------------------------------- ##

cat("Reading phenotype file …\n")
pheno <- readRDS(pheno.file)

cat("Reading methylation beta matrix …\n")
load(meth.file)
meth <- as.matrix(norm.beta.random[, pheno$sample.id])
rm(norm.beta.random)

## --------------------------------------------------------------------------- ##
## 5.  Identify BAD probes                                                     ##
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

cat("  → A total of ",
    length(bad.probes),
    " bad probes (for either 450K or EPIC) identified.\n")

## --------------------------------------------------------------------------- ##
## 6.  Loop through exposures × models                                         ##
## --------------------------------------------------------------------------- ##

for (exp in exposures) {
  for (model in models) {
    tag <- paste(exp, model, sep = ".")
    cat("\nProcessing:", tag, "…\n")
    
    ewas.file <- file.path(ewas.dir,
                           sprintf("%s.%s.%s.%s.csv", study, ewas.prefix, model, exp))
    if (!file.exists(ewas.file)) {
      warning("EWAS file not found: ", ewas.file)
      next
    }
    ewas <- fread(ewas.file)
    ewas <- ewas[!(probeid %in% bad.probes)]
    
    anno <- meffil.get.features(platform)
    anno.sub <- anno[match(ewas$probeid, anno$name), ]
    ewas$chr <- anno.sub$chromosome
    ewas$pos <- anno.sub$position
    
    common.cpg <- intersect(rownames(meth), ewas$probeid)
    meth.sub   <- meth[common.cpg, , drop = FALSE]
    ewas.sub   <- ewas[match(common.cpg, ewas$probeid)]
    
    cat("  → Done: ", nrow(ewas.sub), "probes.\n")
    
    stopifnot(identical(rownames(meth.sub), ewas.sub$probeid))
    
    pre.obj <- dmrff.pre(
      estimate    = ewas.sub$coef,
      se          = ewas.sub$se,
      methylation = meth.sub,
      chr         = ewas.sub$chr,
      pos         = ewas.sub$pos
    )
    
    saveRDS(pre.obj, file = file.path(out.dir, sprintf("%s.dmrff.pre.rds", tag)))
    
    cat("  → dmrff pre-processing object saved to:",
        file.path(out.dir, sprintf("%s.dmrff.pre.rds", tag)),
        "\n")
  }
}

################################################################################
