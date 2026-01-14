################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 13-Jul-2025
# This script is to generate circos plots for differentially methylated regions (DMRs) across exposures.

## --------------------------------------------------------------------------- ##
## 0. Clean workspace and load required libraries                              ##
## --------------------------------------------------------------------------- ##

rm(list = ls())               # remove all existing objects
library(data.table)           # for fast data loading/manipulation
library(readxl)               # to read Excel files
library(circlize)             # for circular plotting
library(scales)               # for graphical scaling functions

## --------------------------------------------------------------------------- ##
## 1. Define input/output paths and exposure labels                            ##
## --------------------------------------------------------------------------- ##

in.file  <- "Z:/working/results/EWAS/meta/DMR_dmrff/ALL.exp.AddModel_DMR_TopTab.xlsx"
out.png  <- "Z:/working/results/EWAS/meta/DMR_dmrff/DMR_CircosPlot.png"

exposure.labels <- c(
  veggie2 = "Pesco-/full vs. non-vegetarian",
  veggie1 = "Full vs. non-vegetarian",
  PDI     = "Overall PDI",
  hPDI    = "Healthful PDI",
  uPDI    = "Unhealthful PDI"
)

bg.palette <- c("#FDE0DD", "#E0ECF4", "#E5F5E0", "#FFF7BC", "#FEE6CE")
# background colors for each exposure ring

## --------------------------------------------------------------------------- ##
## 2. Load data and prepare key columns                                        ##
## --------------------------------------------------------------------------- ##

dt <- as.data.table(read_excel(in.file))  # read the DMR summary statistics

# Determine plotting order by count of DMRs per exposure
exp.dmr.counts <- dt[, .N, by = Exposure][order(N)]
exposure.order <- exp.dmr.counts$Exposure

# Convert columns for plotting
dt[, Exposure     := factor(Exposure, levels = exposure.order)]
dt[, chr          := paste0("chr", Chr)]               # chromosome labels
dt[, mid          := floor((Start + End) / 2)]         # midpoint of each region
dt[, Estimate_num := as.numeric(Estimate)]             # effect estimate (numeric)
dt[, Bonf_num     := as.numeric(gsub("[^0-9eE.-]", "", Bonf_Pvalue))]  # Bonferroni p-value

track.bg.cols <- setNames(bg.palette, exposure.order)

## --------------------------------------------------------------------------- ##
## 3. Select top DMRs for gene labeling                                        ##
## --------------------------------------------------------------------------- ##

# Significant DMRs
gene.labels_sig <- unique(dt[Bonf_num < 0.05 &
                               Gene_symbol != "", .(chr, start = mid, end = mid, Gene_symbol)])

# Extract from dt all DMRs whose gene symbol is in the overlap list
overlap_syms <- scan(
  "Z:/working/results/EWAS/meta/DMR_dmrff/DMR_Overlapped_GeneSymbol_List.txt",
  what = "",
  sep = "\n"
)

gene.labels_overlap <- dt[Gene_symbol %in% overlap_syms &
                            Gene_symbol != "", .(chr = chr[1],
                                                 start = mid[1],
                                                 end = mid[1],
                                                 Gene_symbol), by = Gene_symbol]
gene.labels_overlap <- gene.labels_overlap[, .(chr, start, end, Gene_symbol)]

# Combine and assign colors
gene.labels <- rbind(gene.labels_sig, gene.labels_overlap)
base_col <- ifelse(gene.labels$Gene_symbol %in% overlap_syms,
                   "navyblue",
                   "darkred")

# Assign different colors for repeated gene symbols (if any)
occ <- ave(seq_along(gene.labels$Gene_symbol),
           gene.labels$Gene_symbol,
           FUN = seq_along)

label_cols <- ifelse(occ %% 2 == 1,
                     base_col,
                     ifelse(base_col == "navyblue", "darkred", "navyblue"))

## --------------------------------------------------------------------------- ##
## 4. Initialize PNG device and circos layout                                  ##
## --------------------------------------------------------------------------- ##

png(out.png,
    width = 2400,
    height = 2400,
    res = 300)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after     = c(rep(1, 21), 11),
  track.margin  = c(0.001, 0.001)
)
circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22))

# add gene labels for significant and overlapping DMRs
circos.genomicLabels(
  gene.labels,
  labels.column     = 4,
  side              = "outside",
  cex               = 0.5,
  col               = label_cols,
  line_col          = label_cols,
  connection_height = 0.01,
  font              = 3
)

## --------------------------------------------------------------------------- ##
## 5. Loop over exposures and draw each track                                  ##
## --------------------------------------------------------------------------- ##

exp_list <- rev(exposure.order)
for (i in seq_along(exp_list)) {
  exp <- exp_list[i]
  d   <- dt[Exposure == exp]
  
  # define fixed tick marks and labels for each exposure
  if (exp %in% c("PDI", "hPDI", "uPDI")) {
    ticks <- c(-0.01, -0.005, 0.00, 0.005, 0.01)
    labs  <- formatC(ticks, format = "f", digits = 3)
  } else if (exp == "veggie1") {
    ticks <- c(0.4, 0.2, 0.0, -0.2, -0.4)
    labs  <- formatC(ticks, format = "f", digits = 1)
  } else {
    # veggie2 (innermost ring)
    ticks <- c(0.2, 0.1, 0.0, -0.1, -0.2)
    labs  <- formatC(ticks, format = "f", digits = 1)
  }
  ticks_nz <- ticks[ticks != 0]  # non-zero ticks for dashed grid
  ylim     <- range(ticks)       # y-axis limits for this track
  
  # cache midpoints and annotations for use within panel.fun
  mids_all  <- d$mid
  annos_all <- d$Annotation
  
  circos.genomicTrackPlotRegion(
    data.frame(
      chr   = d$chr,
      start = mids_all,
      end   = mids_all,
      y     = d$Estimate_num,
      bonf  = d$Bonf_num
    ),
    numeric.column = c(4, 5),
    ylim           = ylim,
    stack          = FALSE,
    bg.border      = NA,
    bg.col         = track.bg.cols[exp],
    track.height   = 0.12,
    panel.fun      = function(region, value, ...) {
      # dashed horizontal lines at each non-zero tick
      for (t in ticks_nz) {
        circos.lines(
          CELL_META$xlim,
          c(t, t),
          col = "gray70",
          lty = 2,
          lwd = 0.4
        )
      }
      # solid horizontal line at y = 0
      circos.lines(
        CELL_META$xlim,
        c(0, 0),
        col = "gray30",
        lty = 1,
        lwd = 0.6
      )
      # map midpoints back to annotation vector
      mids_plot <- (region$start + region$end) / 2
      idx       <- match(mids_plot, mids_all)
      anno_vec  <- annos_all[idx]
      # determine filled vs hollow based on Bonferroni p-value
      sig     <- !is.na(value[, 2]) & value[, 2] < 0.05
      pch_vec <- ifelse(sig, 16, 1)
      # color by annotation: red for Promoter, blue otherwise
      is_prom       <- grepl("Promoter", anno_vec, ignore.case = TRUE)
      is_prom[is.na(is_prom)] <- FALSE
      col_vec       <- ifelse(is_prom, "#D73027", "#4575B4")
      # plot points
      circos.points(mids_plot,
                    value[, 1],
                    col = col_vec,
                    pch = pch_vec,
                    cex = 0.6)
      # add y-axis ticks and labels on chr1 only
      if (CELL_META$sector.index == "chr1") {
        circos.yaxis(
          side         = "left",
          at           = ticks,
          labels       = labs,
          tick         = TRUE,
          tick.length  = convert_x(0.3, "mm"),
          labels.cex   = 0.3,
          sector.index = CELL_META$sector.index,
          track.index  = CELL_META$track.index
        )
      }
    }
  )
  
  # add extra y-axis for the innermost ring (veggie2)
  if (exp == "veggie2") {
    current_track <- get.current.track.index()
    circos.yaxis(
      side         = "left",
      at           = ticks,
      labels       = labs,
      tick         = TRUE,
      tick.length  = convert_x(0.3, "mm"),
      labels.cex   = 0.3,
      sector.index = "chr1",
      track.index  = current_track
    )
  }
}

## --------------------------------------------------------------------------- ##
## 6. Add legends                                                              ##
## --------------------------------------------------------------------------- ##

# Legend for DMR points: annotation + significance
legend(
  "topleft",
  inset  = c(0.001, 0.001),
  bty    = "n",
  cex    = 0.7,
  legend = c(
    "Promoter, Bonf-P < 0.05",
    "Promoter, Bonf-P ≥ 0.05",
    "Other annotations, Bonf-P < 0.05",
    "Other annotations, Bonf-P ≥ 0.05"
  ),
  col   = c("#D73027", "#D73027", "#4575B4", "#4575B4"),
  pch   = c(16, 1, 16, 1),
  pt.cex = 1
)

# Legend for exposure rings
dmr_counts <- setNames(exp.dmr.counts$N, exp.dmr.counts$Exposure)
legend_labels <- paste0(exposure.labels[exp_list], " (N DMRs: ", dmr_counts[exp_list], ")")

legend(
  "topright",
  inset  = c(0.001, 0.001),
  bty    = "n",
  cex    = 0.65,
  legend = legend_labels,
  fill   = track.bg.cols[exp_list]
)

dev.off()
message("=== Circos plot saved successfully ===")

################################################################################
