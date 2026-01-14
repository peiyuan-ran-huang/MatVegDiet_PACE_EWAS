################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETs AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 26-Jun-2025
# This script is to prepare cohort-specific EWAS results for meta-analysis.

################################################################################

# Clear environment
rm(list = ls())

# 0. Define the parent directory that holds one subfolder per cohort
parent_dir <- "/user/home/zd20208/MatVegDiet_PACE_EWAS/results"

# 1. Get all immediate subdirectories (one per cohort or subcohort)
cohort_dirs <- list.dirs(path = parent_dir,
                         full.names = TRUE,
                         recursive = FALSE)

# 2. For each cohort folder...
for (cohort_dir in cohort_dirs) {
  # 2.1 Extract the cohort name for messaging
  cohort_name <- basename(cohort_dir)
  
  # 2.2 Skip BiB and DCHS (empty or intentionally excluded)
  if (cohort_name %in% c("BiB", "DCHS", "meta")) {
    message("Skipping cohort ", cohort_name, " (empty or excluded).")
    next
  }
  
  # 2.3 List the EWAS .Rdata files in that folder
  rdata_files <- list.files(path       = cohort_dir,
                            pattern    = "\\.EWASres\\.birth\\..+\\.Rdata$",
                            full.names = TRUE)
  
  # 2.4 If none found, warn and move on
  if (length(rdata_files) == 0) {
    warning("No EWAS .Rdata files found for ", cohort_name)
    next
  }
  
  # 2.5 Use the same directory for outputs (or define a subfolder)
  output_dir <- cohort_dir
  
  # 3. Loop through each .Rdata file in this cohort
  for (rfile in rdata_files) {
    # 3.1 Parse STUDY and EXPOSURE from the filename
    #     e.g. "MatVegDiet.ALSPAC.PDI.EWASres.birth.20230807.Rdata"
    fname   <- basename(rfile)
    parts   <- strsplit(fname, "\\.")[[1]]
    study   <- parts[2]    # should equal cohort_name
    exposure <- parts[3]
    
    # 3.2 Load the .Rdata into its own environment
    env <- new.env()
    load(rfile, envir = env)
    
    # 3.3 For each object in that environment...
    for (obj_name in ls(env)) {
      obj <- env[[obj_name]]
      
      # Only convert data.frames (your EWAS tables)
      if (!is.data.frame(obj))
        next
      
      # 3.4 Extract model name and check match to exposure
      #     Object names look like "ewas.res.AddModel.PDI"
      obj_parts <- strsplit(obj_name, "\\.")[[1]]
      model     <- obj_parts[3]
      exp2      <- obj_parts[4]
      
      # Skip if it doesn’t match
      if (exp2 != exposure) {
        warning(
          "Skipping ",
          obj_name,
          ": exposure tag '",
          exp2,
          "' != file exposure '",
          exposure,
          "'."
        )
        next
      }
      
      # 3.5 Construct the CSV filename
      out_name <- sprintf("%s.ewas.res.%s.%s.csv", study, model, exp2)
      
      # 3.6 Move rownames into a column called "probeid"
      obj_out <- data.frame(
        probeid        = rownames(obj),
        obj,
        row.names      = NULL,
        stringsAsFactors = FALSE
      )
      
      # 3.7 Write the CSV
      write.csv(
        obj_out,
        file      = file.path(output_dir, out_name),
        row.names = FALSE,
        quote     = FALSE
      )
      
      message("[", cohort_name, "] Saved ", out_name)
    }
  }
}
