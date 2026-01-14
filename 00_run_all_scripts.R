################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIETS AND CORD BLOOD DNA METHYLATION     #
################################################################################

# Last edited date: 27-Jun-2025
# This script is to run all cohort-specific scripts in batch.

################################################################################

# Clear environment
rm(list = ls())

# Collect information about the current R session
sessionInfo()

################################################################################

# # Time the whole process
# start_time <- Sys.time ()
# start_time
# write(
#   start_time,
#   "C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/start_time.txt"
# )

################################################################################

# Run scripts
# 01_prepare_EWAS.cohort.R: Run on HPC
# 02_MatVegDiet_EWAS_meta.sh: Run on HPC
# 03_plot(1)_QQ_EWAS.cohort.R: Run on HPC
# 03_plot(2)_precision_EWAS.cohort.R: Run on HPC

source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/04_screen_EWAS.meta.R", echo = T)
# 05(1)_DMR_dmrff.pre_*.R: Run on HPC
# 05(2)_DMR_run_dmrff.meta.R: Run on HPC
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/05(3)_DMR_top_dmrff.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/05(4)_DMR_overlap_dmrff.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/05(5)_DMR_circos_dmrff.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/06(1)_enrichment_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/06(2)_enrichment_DMR.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/06(3)_enrichment_top.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/07(1)_lookup_EWAS.Catalog.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/07(2)_lookup_HELIX.eQTM.R", echo = T)


source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/450K&EPIC_01_QC_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/450K&EPIC_02_corr_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/450K&EPIC_03(1)_top_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/450K&EPIC_03(2)_top_EWAS.meta.R", echo = T)
# 450K&EPIC_04_LOO_EWAS.meta.R: Run on HPC
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/450K&EPIC_05_model.compare_EWAS.meta.R", echo = T)


source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/EPIC.only_01_QC_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/EPIC.only_02_corr_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/EPIC.only_03(1)_top_EWAS.meta.R", echo = T)
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/EPIC.only_03(2)_top_EWAS.meta.R", echo = T)
# EPIC.only_04_LOO_EWAS.meta.R: Run on HPC
source("C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/EPIC.only_05_model.compare_EWAS.meta.R", echo = T)


################################################################################

# # Calculate the running time
# start_time <- as.POSIXct(as.numeric(
#   read.table(
#     "C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/start_time.txt"
#   )
# ), tz = "GMT", origin = "1970-01-01")
# start_time
#
# end_time <- Sys.time ()
# end_time
# write(
#   end_time,
#   "C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/end_time.txt"
# )
#
# time_used <- end_time - start_time
# time_used
# write(
#   time_used,
#   "C:/Users/zd20208/OneDrive - University of Bristol/Desktop/2020 Bristol PhD/Main Project/Data/Scripts/MainAssocAnalysis/EWAS scripts/MatVegDiet_EWAS_meta/time_used.txt"
# )
