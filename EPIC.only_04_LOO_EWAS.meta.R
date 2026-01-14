################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 29-Jun-2025
# This script is to perform leave-one-out analysis for EWAS meta-analysis results (for probes available on both 450K and EPIC).
## For fully-adjusted models (i.e., AddModel)

## --------------------------------------------------------------------------- ##
## 0.  House-keeping                                                           ##
## --------------------------------------------------------------------------- ##

rm(list = ls())

## --------------------------------------------------------------------------- ##
## 1.  Libraries                                                               ##
## --------------------------------------------------------------------------- ##

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
  library(metafor)
  library(ggplot2)
  library(glue)
})

## --------------------------------------------------------------------------- ##
## 2.  Paths & global parameters                                               ##
## --------------------------------------------------------------------------- ##

ewas_dir <- "/user/work/zd20208/MatVegDiet_PACE_EWAS/results"
out_dir  <- file.path(ewas_dir, "meta/LOO")
dir.create(out_dir, showWarnings = FALSE)

prefix <- "EPIC.only"
model  <- "AddModel"

exp_list <- c("PDI", "hPDI", "uPDI")

exp_labels <- c(PDI  = "Overall plant-based diet index (PDI)", hPDI = "Healthful plant-based diet index (hPDI)", uPDI = "Unhealthful plant-based diet index (uPDI)")

cohorts <- c("MoBa4", "MoBa8", "NORTHPOP")

pretty_tbl <- data.frame(
  study  = cohorts,
  pretty = c("MoBa-4", "MoBa-8", "NorthPop"),
  stringsAsFactors = FALSE
)

## --------------------------------------------------------------------------- ##
## 3.  Helper – read one CpG from one cohort result file                       ##
## --------------------------------------------------------------------------- ##

get_cohort_row <- function(cohort, exposure, cpg) {
  f <- file.path(ewas_dir,
                 cohort,
                 glue("{cohort}.ewas.res.{model}.{exposure}.csv"))
  
  if (!file.exists(f))
    stop("File not found: ", f)
  
  dt <- fread(f, showProgress = FALSE)
  
  required_cols <- c("probeid", "coef", "se", "p")
  if (!all(required_cols %in% names(dt)))
    stop("Missing one or more required columns in ", f)
  
  out <- dt[probeid == cpg, .(study = cohort, beta = coef, se = se)]
  
  return(out)
}

## --------------------------------------------------------------------------- ##
## 4.  Loop over exposures                                                     ##
## --------------------------------------------------------------------------- ##

for (exp in exp_list) {
  message("\n===========  Exposure: ", exp, "  ===========")
  
  exp_pretty <- exp_labels[[exp]]
  
  top_tab <- file.path(ewas_dir, glue("meta/{prefix}_{exp}_{model}_TopTab.xlsx"))
  if (!file.exists(top_tab)) {
    message("  TopTab not found – skipped.")
    next
  }
  
  dt <- data.table(read_xlsx(top_tab))
  sel_cpgs <- dt[as.numeric(FDR) < 0.05, CpG]
  
  if (length(sel_cpgs) == 0L) {
    message("  No CpG with FDR < 0.05.")
    next
  }
  
  loo_collect <- vector("list", length(sel_cpgs))
  names(loo_collect) <- sel_cpgs
  
  for (i in seq_along(sel_cpgs)) {
    cpg <- sel_cpgs[i]
    message(sprintf("  [%d/%d] %s", i, length(sel_cpgs), cpg))
    
    dat <- rbindlist(lapply(
      cohorts,
      get_cohort_row,
      exposure = exp,
      cpg      = cpg
    ),
    use.names = TRUE)
    
    dat <- merge(dat,
                 pretty_tbl,
                 by = "study",
                 all.x = TRUE,
                 sort = FALSE)
    
    if (nrow(dat) < 2L) {
      warning("    < 2 cohorts available – skipped")
      next
    }
    
    fit_all <- rma.uni(
      yi     = beta,
      sei    = se,
      data   = dat,
      method = "REML"
    )
    
    leave_out_list <- vector("list", nrow(dat))
    
    for (j in seq_len(nrow(dat))) {
      dat_j <- dat[-j]
      fit_j <- rma.uni(
        yi = beta,
        sei = se,
        data = dat_j,
        method = "REML"
      )
      
      leave_out_list[[j]] <- data.table(
        Exposure = exp_pretty,
        CpG      = cpg,
        Left_out = dat$pretty[j],
        k        = fit_j$k,
        beta     = coef(fit_j),
        se       = fit_j$se,
        z_val    = fit_j$zval,
        p_val    = fit_j$pval,
        tau2     = fit_j$tau2,
        I2       = fit_j$I2
      )
    }
    
    loo_dt <- rbind(
      data.table(
        Exposure = exp_pretty,
        CpG      = cpg,
        Left_out = "None",
        k        = fit_all$k,
        beta     = coef(fit_all),
        se       = fit_all$se,
        z_val    = fit_all$zval,
        p_val    = fit_all$pval,
        tau2     = fit_all$tau2,
        I2       = fit_all$I2
      ),
      rbindlist(leave_out_list, use.names = TRUE)
    )
    
    loo_collect[[i]] <- loo_dt
    
    if (nrow(loo_dt) > 1) {
      loo_dt[, Left_out := paste0(Left_out, " out")]
      loo_dt[, Left_out := factor(Left_out, levels = Left_out)]  # preserve order
      
      fig_path <- file.path(out_dir, glue("{prefix}_{exp}_{model}_LOO_{cpg}.png"))
      
      png(
        filename = fig_path,
        width    = 1400,
        height   = 600,
        res      = 150
      )
      
      print(
        ggplot(loo_dt, aes(x = Left_out, y = beta)) +
          geom_hline(
            yintercept = 0,
            linetype = "dotted",
            colour = "grey40"
          ) +
          geom_hline(
            yintercept = loo_dt$beta[1],
            linetype = 2,
            colour = "red"
          ) +
          geom_errorbar(
            aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se),
            width = 0.2,
            colour = "steelblue"
          ) +
          geom_point(size = 3, colour = "steelblue") +
          labs(
            title    = glue("{cpg}  |  {exp_pretty}"),
            x        = NULL,
            y        = expression(paste("Beta (95% CI)"))
          ) +
          theme_minimal(base_size = 11) +
          theme(axis.text.x = element_text(
            angle = 45, hjust = 1
          ))
      )
      
      invisible(dev.off())
      message("    → LOO plot saved: ", basename(fig_path))
      
    } else {
      message("    → No LOO plot: only 1 cohort")
    }
  }
  
  all_loo <- rbindlist(loo_collect, use.names = TRUE, fill = TRUE)
  
  fwrite(all_loo, file = file.path(out_dir, glue("{prefix}_{exp}_{model}_LOO.csv")))
  
  write.xlsx(
    as.data.frame(all_loo),
    file      = file.path(out_dir, glue("{prefix}_{exp}_{model}_LOO.xlsx")),
    overwrite = TRUE,
    rowNames  = FALSE
  )
}

message("\nAll exposures finished – LOO outputs saved to:\n", out_dir)

################################################################################
