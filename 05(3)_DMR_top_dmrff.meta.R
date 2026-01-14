################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 13-Jul-2025
# This script is to select top hits of differentially methylated region (DMR) meta-analysis results.

## --------------------------------------------------------------------------- ##
## 0.  Clean workspace                                                         ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(meffil)
library(openxlsx)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

## --------------------------------------------------------------------------- ##
## 2.  Parameters                                                              ##
## --------------------------------------------------------------------------- ##

exposures <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")

platforms <- list("450K" = meffil.get.features("450k")$name,
                  "EPIC"  = meffil.get.features("epic")$name)

in.dir   <- "Z:/working/results/EWAS/meta/DMR_dmrff/"
out.dir  <- in.dir
p.cutoff <- 0.05

## --------------------------------------------------------------------------- ##
## 3.  Helper functions                                                        ##
## --------------------------------------------------------------------------- ##

annotate_platform <- function(cpg.list) {
  sapply(cpg.list, function(cpgs) {
    in.450k <- cpgs %in% platforms[["450K"]]
    in.epic <- cpgs %in% platforms[["EPIC"]]
    
    if (all(in.450k & in.epic)) {
      "Both"
    } else if (all(in.450k & !in.epic)) {
      "450K only"
    } else if (all(!in.450k & in.epic)) {
      "EPIC only"
    } else {
      "Mixed"
    }
  })
}

get_cpgs_in_region <- function(region_chr,
                               region_start,
                               region_end,
                               cpg.lookup) {
  cpg.lookup[chr == region_chr &
               pos >= region_start &
               pos <= region_end, cpg]
}

## --------------------------------------------------------------------------- ##
## 4.  Main loop                                                               ##
## --------------------------------------------------------------------------- ##

all.results <- list()  # collect per-exposure results for final merge

for (exposure in exposures) {
  tag <- paste0(exposure, ".AddModel")
  cat("\nProcessing:", tag, "...\n")
  
  rds.file <- file.path(in.dir, paste0(tag, "_dmrff.meta.rds"))
  if (!file.exists(rds.file)) {
    warning("RDS file not found: ", rds.file)
    next
  }
  
  dt <- readRDS(rds.file)
  if (!("dmrs" %in% names(dt)) || !("ewas" %in% names(dt))) {
    warning("Invalid RDS structure: ", tag)
    next
  }
  
  dmr <- as.data.table(dt$dmrs)
  ewas <- as.data.table(dt$ewas, keep.rownames = "cpg")
  cpg.lookup <- ewas[, .(cpg, chr, pos)]
  
  if (!nrow(dmr)) {
    warning("No DMRs in result: ", tag)
    next
  }
  
  ## --- Annotate CpGs per region ------------------------------------------- ##
  chr.vec   <- as.character(dmr$chr)
  start.vec <- as.integer(dmr$start)
  end.vec   <- as.integer(dmr$end)
  
  cpg.list <- vector("list", length = nrow(dmr))
  for (i in seq_len(nrow(dmr))) {
    cpg.list[[i]] <- get_cpgs_in_region(chr.vec[i], start.vec[i], end.vec[i], cpg.lookup)
  }
  
  dmr$cpgs <- cpg.list
  dmr <- as.data.table(dmr)
  
  ## --- Annotate platform and compute FDR and Bonferroni -------------------- ##
  dmr[, platform := annotate_platform(cpgs)]
  dmr[, `:=`(FDR.by.platform = NA_real_, Bonf.by.platform = NA_real_)]
  
  for (p in unique(dmr$platform)) {
    ix <- which(dmr$platform == p)
    dmr$FDR.by.platform[ix]  <- p.adjust(dmr$p.value[ix], method = "fdr")
    dmr$Bonf.by.platform[ix] <- p.adjust(dmr$p.value[ix], method = "bonferroni")
  }
  
  ## --- Order and filter --------------------------------------------------- ##
  platform.order <- c("Both", "450K only", "EPIC only", "Mixed")
  dmr[, platform_rank := match(platform, platform.order)]
  setorder(dmr, platform_rank, p.value)
  dmr[, platform_rank := NULL]
  
  top.dmr <- dmr[FDR.by.platform < p.cutoff & n >= 2]
  if (!nrow(top.dmr)) {
    warning("No top DMRs found: ", tag)
    next
  }
  
  ## --- Gene annotation ---------------------------------------------------- ##
  gr <- GRanges(
    seqnames  = ifelse(
      startsWith(top.dmr$chr, "chr"),
      top.dmr$chr,
      paste0("chr", top.dmr$chr)
    ),
    ranges    = IRanges(start = top.dmr$start, end = top.dmr$end),
    region_id = seq_len(nrow(top.dmr))
  )
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  peakAnno <- annotatePeak(
    peak      = gr,
    TxDb      = txdb,
    tssRegion = c(-3000, 3000),
    annoDb    = "org.Hs.eg.db",
    addFlankGeneInfo = FALSE,
    flankDistance    = 0,
    verbose   = FALSE
  )
  
  anno.df <- as.data.table(as.data.frame(peakAnno))
  setnames(anno.df, c("geneId", "SYMBOL"), c("gene_id", "gene_symbol"))
  
  anno.out <- cbind(top.dmr, anno.df[, .(annotation, gene_id, gene_symbol, distanceToTSS)])
  
  ## --- Clean and reformat columns ----------------------------------------- ##
  # Format Bonf_Pvalue as required
  bonf_fmt <- fifelse(
    anno.out$Bonf.by.platform < 0.001,
    "<0.001",
    fifelse(
      anno.out$Bonf.by.platform > 0.999,
      ">0.999",
      sprintf("%.3f", anno.out$Bonf.by.platform)
    )
  )
  
  # Blank out gene_symbol & annotation if Distal Intergenic
  anno.out[annotation == "Distal Intergenic", `:=`(gene_symbol = "", annotation = "")]
  
  # Final tidy output
  anno.out.pretty <- anno.out[, .(
    Exposure     = exposure,
    Chr          = sub("^chr", "", chr),
    Start        = start,
    End          = end,
    Gene_symbol  = gene_symbol,
    Annotation   = annotation,
    N_CpGs       = n,
    Estimate     = sprintf("%.7f", estimate),
    SE           = sprintf("%.7f", se),
    Pvalue       = formatC(p.value, format = "e", digits = 2),
    FDR_Pvalue   = formatC(FDR.by.platform, format = "e", digits = 2),
    Bonf_Pvalue  = bonf_fmt,
    Platform     = platform
  )]
  
  ## --- Write per-exposure output ------------------------------------------ ##
  out.xlsx <- file.path(out.dir, paste0(tag, "_DMR_TopTab.xlsx"))
  write.xlsx(as.data.frame(anno.out.pretty),
             file = out.xlsx,
             rowNames = FALSE)
  cat("Top DMRs for", tag, "written to:\n", out.xlsx, "\n")
  
  ## --- Save to list for combined output ----------------------------------- ##
  all.results[[exposure]] <- anno.out.pretty
}

## --- Merge all exposure tables and export ---------------------------------- ##
if (length(all.results)) {
  combined <- rbindlist(all.results, use.names = TRUE)
  out.file <- file.path(out.dir, "ALL.exp.AddModel_DMR_TopTab.xlsx")
  write.xlsx(as.data.frame(combined),
             file = out.file,
             rowNames = FALSE)
  cat("\nCombined DMR results written to:\n", out.file, "\n")
}

message("\n=== All top DMR tables generated ===")

################################################################################
