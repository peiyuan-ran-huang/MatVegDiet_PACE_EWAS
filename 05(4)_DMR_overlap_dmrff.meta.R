################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 13-Jul-2025
# This script is to identify overlapped differentially methylated regions (DMRs) across exposures.

## --------------------------------------------------------------------------- ##
## 0. Clean workspace                                                          ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1. Libraries                                                                ##
## --------------------------------------------------------------------------- ##

library(data.table)
library(GenomicRanges)
library(VennDiagram)
library(openxlsx)
library(RColorBrewer)

## --------------------------------------------------------------------------- ##
## 2. Parameters                                                               ##
## --------------------------------------------------------------------------- ##

exposures <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")
indir     <- "Z:/working/results/EWAS/meta/DMR_dmrff"
outdir    <- file.path(indir, "")
dir.create(outdir, showWarnings = FALSE)

## --------------------------------------------------------------------------- ##
## 3. Load top DMRs and construct GRanges + CpG sets + Genes                   ##
## --------------------------------------------------------------------------- ##

gr.list   <- list()
cpg.list  <- list()
gene.list <- list()

for (exp in exposures) {
  file <- file.path(indir, paste0(exp, ".AddModel_DMR_TopTab.xlsx"))
  if (!file.exists(file))
    next
  
  dt <- as.data.table(read.xlsx(file))
  if (!nrow(dt))
    next
  
  # Genomic ranges
  gr <- GRanges(
    seqnames = paste0("chr", dt$Chr),
    ranges   = IRanges(start = dt$Start, end = dt$End),
    region_id = paste0(exp, "_", seq_len(nrow(dt)))
  )
  gr.list[[exp]] <- gr
  
  # CpG list: from meta rds
  rds.file <- file.path(indir, paste0(exp, ".AddModel_dmrff.meta.rds"))
  if (file.exists(rds.file)) {
    rds <- readRDS(rds.file)
    dmrs <- as.data.table(rds$dmrs)
    sel  <- dmrs$p.value < 1e-5 & dmrs$n >= 2
    cpgs <- unique(unlist(dmrs$cpgs[sel]))
    cpg.list[[exp]] <- cpgs
  }
  
  # Gene symbols
  gene.list[[exp]] <- unique(na.omit(dt$Gene_symbol))
}

## --------------------------------------------------------------------------- ##
## 4. GenomicRanges overlap matrix (binary)                                    ##
## --------------------------------------------------------------------------- ##

overlap.table <- matrix(
  0,
  nrow = length(exposures),
  ncol = length(exposures),
  dimnames = list(exposures, exposures)
)

for (i in exposures) {
  for (j in exposures) {
    if (i %in% names(gr.list) && j %in% names(gr.list)) {
      olap <- findOverlaps(gr.list[[i]], gr.list[[j]], minoverlap = 1)
      overlap.table[i, j] <- length(unique(queryHits(olap)))
    }
  }
}

write.xlsx(as.data.frame(overlap.table),
           file = file.path(outdir, "DMR_RegionOverlap_Matrix.xlsx"))

## --------------------------------------------------------------------------- ##
## 5. CpG overlap summary                                                      ##
## --------------------------------------------------------------------------- ##

cpg.sets <- cpg.list[intersect(names(cpg.list), exposures)]

# Save per-exposure CpG counts
if (length(cpg.sets) > 0) {
  cpg.counts.df <- data.table(Exposure = names(cpg.sets),
                              N_CpGs   = sapply(cpg.sets, length))
  fwrite(cpg.counts.df, file = file.path(outdir, "DMR_CpG_Counts.csv"))
} else {
  warning("No CpG sets found. Skipping CpG overlap summary.")
}

# Pairwise CpG overlaps
if (length(cpg.sets) >= 2) {
  cpg.intersection <- combn(
    names(cpg.sets),
    2,
    simplify = FALSE,
    FUN = function(x) {
      list(
        Exposure1 = x[1],
        Exposure2 = x[2],
        N_overlap = length(intersect(cpg.sets[[x[1]]], cpg.sets[[x[2]]]))
      )
    }
  )
  
  cpg.inter.df <- rbindlist(cpg.intersection)
  fwrite(cpg.inter.df, file = file.path(outdir, "DMR_CpG_Overlap_Summary.csv"))
  
} else {
  message("Too few exposures with valid CpG sets for pairwise comparison.")
}

## --------------------------------------------------------------------------- ##
## 6. Gene-level overlap summary                                               ##
## --------------------------------------------------------------------------- ##

gene.sets <- gene.list[intersect(names(gene.list), exposures)]
gene.counts <- sapply(gene.sets, length)

# Pairwise overlaps
gene.overlap <- combn(
  names(gene.sets),
  2,
  simplify = FALSE,
  FUN = function(x) {
    list(
      Exposure1 = x[1],
      Exposure2 = x[2],
      N_overlap = length(intersect(gene.sets[[x[1]]], gene.sets[[x[2]]]))
    )
  }
)
gene.overlap.df <- rbindlist(gene.overlap)
write.xlsx(gene.overlap.df,
           file = file.path(outdir, "DMR_GeneOverlap_Summary.xlsx"))

## --------------------------------------------------------------------------- ##
## 7. Optional: Venn diagram for gene symbols                                  ##
## --------------------------------------------------------------------------- ##

library(grid)

pretty.names <- c(
  veggie2 = "Pesco-/full vs. non-vegetarian",
  veggie1 = "Full\nvs.\nnon-\nvegetarian",
  PDI     = "PDI",
  hPDI    = "hPDI",
  uPDI    = "uPDI"
)

# Replace names with pretty labels (only for those in gene.sets)
names(gene.sets) <- pretty.names[names(gene.sets)]

# Generate venn diagram object only (no file/log)
venn.plot <- venn.diagram(
  gene.sets,
  filename       = NULL,
  imagetype      = "png",
  height         = 4000,
  width          = 4000,
  resolution     = 600,
  fill           = c("red", "blue", "purple", "green", "orange"),
  cex            = 1.5,
  cat.cex        = 1.2,
  cat.dist       = c(0.13, 0.12, 0.23, 0.25, 0.18),
  cat.fontface   = "bold",
  fontface       = "plain",
  main           = NULL
)

# Save the image manually (no .log will be produced)
png(
  filename = file.path(outdir, "DMR_GeneSymbol_Venn.png"),
  height = 5000,
  width = 5000,
  res = 600
)
grid.draw(venn.plot)
dev.off()

## --------------------------------------------------------------------------- ##
## 8. Gene symbol overlap: list of genes per combination (multi-sheet Excel)   ##
## --------------------------------------------------------------------------- ##

library(limma)

# Build presence/absence matrix for all gene symbols
all.genes <- unique(unlist(gene.sets))
gene.mat <- sapply(gene.sets, function(x)
  all.genes %in% x)
rownames(gene.mat) <- all.genes

# Generate all overlap combinations with ≥1 gene
venn.table <- vennCounts(gene.mat)

# Mapping for shorter labels to avoid Excel sheet name length overflow
sheet.label.map <- c(
  "Pesco-/full vs. non-vegetarian" = "V2",
  "Full\nvs.\nnon-\nvegetarian" = "V1",
  "PDI"     = "PDI",
  "hPDI"    = "hPDI",
  "uPDI"    = "uPDI"
)

# Create list to store per-combination gene symbols
comb.list <- list()

for (i in seq_len(nrow(venn.table))) {
  pattern <- venn.table[i, 1:length(exposures)] == 1
  
  combo.genes <- all.genes[apply(gene.mat, 1, function(row)
    all(row[names(pattern)] == pattern))]
  
  if (length(combo.genes)) {
    # Construct safe sheet name (short and Excel-compliant)
    combo.key <- paste0(sheet.label.map[names(pattern)[pattern]], collapse = "_")
    sheet.name <- paste0("Combo_", combo.key)
    sheet.name <- substr(sheet.name, 1, 31)  # Excel limit
    
    comb.list[[sheet.name]] <- data.table(Gene_symbol = combo.genes)
  }
}

# Write to Excel (each overlap combination in separate sheet)
out.file <- file.path(outdir, "DMR_GeneSymbol_OverlapDetails.xlsx")
write.xlsx(comb.list, file = out.file, rowNames = FALSE)

## --------------------------------------------------------------------------- ##
## 9. Heatmap for overlapping DMR genes                                        ##
## --------------------------------------------------------------------------- ##

library(ggplot2)
library(data.table)
library(openxlsx)

# Step 1: Keep genes with ≥2 exposures and non-empty names
gene.mat.overlap <- gene.mat[rowSums(gene.mat) >= 2 &
                               !(is.na(rownames(gene.mat)) |
                                   rownames(gene.mat) == ""), , drop = FALSE]
colnames(gene.mat.overlap) <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")

# Step 2: Read and clean top table (with Chr and Start)
top.tab <- read.xlsx(file.path(outdir, "ALL.exp.AddModel_DMR_TopTab.xlsx"))
top.tab <- top.tab[, c("Exposure",
                       "Gene_symbol",
                       "Estimate",
                       "SE",
                       "Bonf_Pvalue",
                       "Chr",
                       "Start")]
top.tab$Bonf_Pvalue <- gsub("<0.001", "0.0001", top.tab$Bonf_Pvalue)
top.tab$Bonf_Pvalue <- gsub(">0.999", "0.9999", top.tab$Bonf_Pvalue)
top.tab$Estimate    <- as.numeric(top.tab$Estimate)
top.tab$SE          <- as.numeric(top.tab$SE)
top.tab$Bonf_Pvalue <- as.numeric(top.tab$Bonf_Pvalue)

# Step 3: Get one DMR per Gene × Exposure (smallest Bonf_Pvalue)
genes.keep <- rownames(gene.mat.overlap)
exps.keep  <- c("veggie2", "veggie1", "PDI", "hPDI", "uPDI")

plot.df <- data.table()
for (g in genes.keep) {
  for (e in exps.keep) {
    dt.sub <- top.tab[top.tab$Gene_symbol == g &
                        top.tab$Exposure == e, ]
    if (nrow(dt.sub) > 0 && any(!is.na(dt.sub$Bonf_Pvalue))) {
      i <- which.min(dt.sub$Bonf_Pvalue)
      est <- dt.sub$Estimate[i]
      se  <- dt.sub$SE[i]
      z   <- if (!is.na(se) && se != 0)
        est / se
      else
        NA
      p   <- dt.sub$Bonf_Pvalue[i]
      chr <- dt.sub$Chr[i]
      pos <- dt.sub$Start[i]
      plot.df <- rbind(plot.df,
                       data.table(
                         Gene = g,
                         Exposure = e,
                         Beta = z,
                         P = p,
                         Chr = chr,
                         Start = pos
                       ))
    }
  }
}

# Step 4: Add all Gene × Exposure combinations (including missing ones)
all.comb <- CJ(Gene = genes.keep, Exposure = exps.keep)
plot.df <- merge(all.comb,
                 plot.df,
                 by = c("Gene", "Exposure"),
                 all.x = TRUE)

# Step 5: Add significance symbols
plot.df[, Sig := NA_character_]
plot.df[!is.na(P), Sig := as.character(cut(
  P,
  breaks = c(0, 0.05, 1),
  labels = c("*", ""),
  right = FALSE
))]
plot.df[is.na(P), Sig := "NA"]

# Step 6: Determine gene order by genomic coordinates
gene.coords <- plot.df[!is.na(Chr), .SD[1], by = Gene]  # one row per gene
gene.coords[, Chr_numeric := as.numeric(gsub("chr", "", Chr))]  # remove "chr" if present
gene.coords[is.na(Chr_numeric), Chr_numeric := 99]  # place unknown at end
gene.coords <- gene.coords[order(Chr_numeric, Start)]
ordered.genes <- gene.coords$Gene

genelist <- file.path(outdir, "DMR_Overlapped_GeneSymbol_List.txt")
fwrite(data.table(ordered.genes),
       genelist,
       col.names = FALSE,
       quote     = FALSE)
message("  Saved overlapped gene list to: ", genelist)

# Step 7: Pretty exposure names and gene factor levels
exp.pretty <- c(
  veggie2 = "Pesco-/full vs.\nnon-vegetarian",
  veggie1 = "Full vs.\nnon-vegetarian",
  PDI     = "PDI",
  hPDI    = "hPDI",
  uPDI    = "uPDI"
)
plot.df[, Exposure := factor(exp.pretty[Exposure], levels = exp.pretty)]
plot.df[, Gene := factor(Gene, levels = rev(ordered.genes))]

# Step 8: Label NA explicitly
plot.df[, Label := ifelse(Sig == "NA", "NA", Sig)]

# Step 9: Draw heatmap
rd_bu_palette <- rev(brewer.pal(n = 11, name = "RdBu"))

p <- ggplot(plot.df, aes(x = Exposure, y = Gene)) +
  geom_tile(aes(fill = Beta), color = "white") +
  geom_text(aes(
    label = Label,
    fontface = ifelse(Sig == "NA", "plain", "bold")
  ), size = 2.5) +
  scale_fill_gradientn(colours = rd_bu_palette,
                       na.value = "grey90",
                       name = "Estimate/SE") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_text(face = "bold.italic"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) + guides(fill = guide_colourbar(
    title.theme = element_text(size = 9, face = "bold"),
    label.theme = element_text(size = 9)
  ))

# Step 10: Save plot
ggsave(
  filename = file.path(outdir, "DMR_GeneSymbol_Heatmap.png"),
  plot = p,
  width = 6.5,
  height = 5,
  dpi = 300
)

## --------------------------------------------------------------------------- ##
## Done                                                                        ##
## --------------------------------------------------------------------------- ##

message("\n=== DMR overlap analysis complete ===\nResults in: ", outdir)

################################################################################
