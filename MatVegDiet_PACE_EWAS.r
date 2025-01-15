################################################################################
#     MATERNAL VEGETARIAN/PLANT-BASED DIET AND CORD BLOOD DNA METHYLATION      #
#                           EWAS Analysis Plan Code                            #
#------------------------------------------------------------------------------#
#                             Peiyuan (Ran) Huang                              #
#                                 August 2024                                  #
#         Acknowledgements: Gemma Sharp, Leanne Kupers, Giulia Mancano         #
################################################################################

################################################################################

# The following R code will allow you to complete all the EWAS requested in the vegetarian diet EWAS analysis plan.
# The code also produces .csv files summarising the variables included in the EWASes.
# [NOTE: You should not have to rewrite or add to the following code,
#        UNLESS otherwise stated (see the "[PLEASE CHANGE]" notice).]
#
# There are two input files required for this analysis:
#
# 1. pheno: A dataframe containing all the phenotype data needed for this project. 
#    Each row is a sample (individual) and each column is a different variable. 
#    Necessary variable names are: 
#    - Sample ID (1): "sample.id"
#    - Food groups (18): "wholegrain", "fruit", "vegetable", "nut", "legume", "vegetableoil", "teacoffee", "fruitjuice", "refinedgrain", "potato", "sugarbeverage", "sweetdessert", "animalfat", "dairy", "egg", "fishseafood", "meat", "misc.animal"
#    - Intake frequency of key food groups (4): "meat_freq", "fishseafood_freq", "egg_freq", "dairy_freq"
#    - Covariates (8): "sex", "mat.age", "mat.edu", "mat.parity", "mat.smoking", "mat.bmi", "mat.suppl", "mat.kcal"
#    - Cell types (7): "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"
#    Please use this exact order of variable names and do not include other columns in the pheno file.
#    Please use the Salas reference set for cell type correction
#    If these columns are named differently in your dataset, please rename the columns accordingly
#    Details on how to code these variables are provided in the analysis plan and supplementary documents.
#
# 2. meth: A matrix of methylation Illumina beta values. Each column is a sample and each row is a probe on the array (450k or EPIC). 
#    Column names must correspond to the sample.id column in pheno.
#
# If you have any questions about this script, please contact Peiyuan (Ran) Huang (peiyuan.huang@bristol.ac.uk) and cc Dr Gemma Sharp (g.c.sharp@exeter.ac.uk).

# Update on 2023/08/08: No longer consider semi-vegetarian as a vegetarian subgroup (see changes in Lines 228-250)

################################################################################

#------------------------------------------------------------------------------#
#                             Packages & Functions                             #----
#------------------------------------------------------------------------------#

# Load required packages
# [NOTE: If these are not already installed, you will have to install them as the first step.]
library(tableone)
library(matrixStats)
library(limma)
library(data.table)
library(plyr)
library(sva)
library(qqman)
library(dmrff)

# Setup the necessary functions

## Function to remove outliers using the Tukey IQR*3 method
IQR.removal <- function(meth) {
  rowIQR <- rowIQRs(meth, na.rm = T)
  row2575 <- rowQuantiles(meth, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth < row2575[, 1] - 3 * rowIQR
  maskU <- meth > row2575[, 2] + 3 * rowIQR
  meth[maskL] <- NA
  meth[maskU] <- NA
  return(meth)
}

## Function to run EWAS
## [NOTE: Here we use least squares instead of robust regression lmFit, as we will manually remove outliers before running EWAS.]
ewas.function <- function(meth, pheno, variable.of.interest) {
  model.covariates <- colnames(pheno)[-which(colnames(pheno) %in% c(variable.of.interest, "sample.id"))]
  des <- model.matrix(reformulate(paste0("pheno$", c(variable.of.interest, model.covariates))))
  fit <- lmFit(meth, des, method = "ls")  # Use least squares instead of robust regression (method = "robust")
  fit.ebayes <- eBayes(fit)
  n <- rowSums(!is.na(meth))
  se <- sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[, grep(paste0(variable.of.interest, collapse = "|"), colnames(fit.ebayes$stdev.unscaled))]
  res <- data.frame(n = n,
                    coef = fit.ebayes$coefficient[, grep(paste0(variable.of.interest, collapse = "|"), colnames(fit.ebayes$coefficient))],
                    se = se,
                    p = fit.ebayes$p.value[, grep(paste0(variable.of.interest, collapse = "|"), colnames(fit.ebayes$p.value))])
  return(res)
}

#------------------------------------------------------------------------------#
#                        Load & Prepare Phenotype Data                         #----
#------------------------------------------------------------------------------#

# [PLEASE CHANGE] Set to your working directory where you want to save the results from this script
setwd("...")

# Set initial parameters
# [NOTE: We will use 20 surrogate variables (SVs; named as sv1-sv20) to adjust for batch in this study.
#        See details in the "Generate Surrogate Variables" section below.]
study <- "ALSPAC"  # [PLEASE CHANGE] Set your study identifier
analysis.date <- format(Sys.Date(), "%Y%m%d") # Automatically set to the date on which you perform the analyses
timepoint <- "birth"  # Should be "birth" for this study
cell.names <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")  # Should be these 7 cell types of according to the Salas reference
traits.and.covariates <- c(
  ## 5 exposure variables of interest
  "veggie1", "veggie2", "PDI", "hPDI", "uPDI",
  
  ## 18 food groups
  "wholegrain", "fruit", "vegetable", "nut", "legume", "vegetableoil", "teacoffee", "fruitjuice", "refinedgrain", "potato", "sugarbeverage", "sweetdessert", "animalfat", "dairy", "egg", "fishseafood", "meat", "misc.animal",
  
  ## Covariates (including SVs)
  "sex", "mat.age", "mat.edu", "mat.parity", "mat.smoking", "mat.bmi", "mat.kcal", "mat.suppl", paste0("sv", 1:20)
  )  # You can also add selection factors relevant to your cohort

# Set up 4 models for each of the 5 exposure variables
# [NOTE: The 5 exposure variables include:
#        - "veggie1": Full vegetarian vs. non-vegetarian
#        - "veggie2": Full vegetarian + pesco-vegetarian vs. non-vegetarian
#        - PDI: The overall plant-based diet index
#        - hPDI: The healthful plant-based diet index
#        - uPDI: The unhealthful plant-based diet index]
for (my_exp in c("veggie1", "veggie2", "PDI", "hPDI", "uPDI")) {
  ## 1. The "no-cells" model (NoCellModel): exposure + child sex + batch (20 SVs)
  assign(paste0("covs.NoCellModel.", my_exp), c(my_exp, "sex", paste0("sv", 1:20)))
  
  ## 2. The minimally adjusted model (MinModel): above + cell types
  assign(paste0("covs.MinModel.", my_exp), c(my_exp, "sex", paste0("sv", 1:20), cell.names))
  
  ## 3. The fully adjusted model (FullModel): above + maternal age + education + parity + smoking
  assign(paste0("covs.FullModel.", my_exp), c(my_exp, "sex", paste0("sv", 1:20), cell.names, "mat.age", "mat.edu", "mat.parity", "mat.smoking"))
  
  ## 4. The additional model (AddModel): above + maternal BMI + energy intake + dietary supplementation
  assign(paste0("covs.AddModel.", my_exp), c(my_exp, "sex", paste0("sv", 1:20), cell.names, "mat.age", "mat.edu", "mat.parity", "mat.smoking", "mat.bmi", "mat.kcal", "mat.suppl"))
}

# Load phenotype data including cell types

## [PLEASE CHANGE] Load phenotype data (named as "pheno") from the path where your phenotype data is stored
pheno <- readRDS("...")
dim(pheno)  # Should be 38 variables; in ALSPAC, N = 835 (with both phenotype and methylation data)
# [NOTE: Variable/column names of the 38 variables:
#        - Sample ID (1): "sample.id"
#        - Food groups (18): "wholegrain", "fruit", "vegetable", "nut", "legume", "vegetableoil", "teacoffee", "fruitjuice", "refinedgrain", "potato", "sugarbeverage", "sweetdessert", "animalfat", "dairy", "egg", "fishseafood", "meat", "misc.animal"
#        - Frequency of key food groups (4): "meat_freq", "fishseafood_freq", "egg_freq", "dairy_freq"
#        - Covariates (8): "sex", "mat.age", "mat.edu", "mat.parity", "mat.smoking", "mat.bmi", "mat.suppl", "mat.kcal"
#        - Cell types (7): "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"

## Check phenotype data
for (i in 1:length(c("sample.id", traits.and.covariates, cell.names))) {
  print(ifelse(c("sample.id", traits.and.covariates, cell.names)[i] %in% colnames(pheno) == F,
               paste("Warning: The variable called", c("sample.id", traits.and.covariates, cell.names)[i], "is missing from pheno."),
               paste("OK: The variable called", c("sample.id", traits.and.covariates, cell.names)[i], "is present in pheno.")))
}

#------------------------------------------------------------------------------#
#                       Load & Prepare Methylation Data                        #----
#------------------------------------------------------------------------------#

# [PLEASE CHANGE] Load methylation data (named as "meth") from the path where your methylation data is stored
load("...")
meth <- as.matrix(norm.beta.random[, pheno$sample.id])  # Rename the methylation data, only keep the samples in pheno, and make sure it is a matrix
rm(norm.beta.random)

# [PLEASE CHANGE] Filter probes
# [NOTE: You can use your preferred methods for probe filtering.]
################################################################################
# Below shows as an example of how ALSPAC remove: 1) probes with high detection p-values, 2) probes on sex chromosomes, and 3) probes used as controls:

## 1) Probes with high detection p-values

### Load detection p-values
load("...")
pvals <- as.matrix(detp[, pheno$sample.id])  # Rename the detection p-value data, only keep the samples in pheno, and make sure it is a matrix (similar to meth)
rm(detp)

### Identify probes with high detection p-values
pvalue_over_0.05 <- pvals > 0.05
count_over_0.05 <- rowSums(sign(pvalue_over_0.05))
probes.to.exclude.p <- rownames(pvals)[which(count_over_0.05 > ncol(pvals) * 0.05)]

## 2) Probes on sex chromosomes

### Extract annotation data using the "meffil" package (see package details at: https://github.com/perishky/meffil)
### [PLEASE CHANGE] Please uncomment the code below if the required package has not been installed
# source("http://bioconductor.org/biocLite.R")
# install.packages("devtools")  # If not already installed
# library(devtools)
# install_github("perishky/meffil")
library(meffil)
annotation <- meffil.get.features("450k")

### Identify probes on sex chromosomes
XY <- as.character(annotation$name[which(annotation$chromosome %in% c("chrX", "chrY"))])

## 3) Probes used as controls
SNPs.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])

## Remove all above probes
print(paste("Before filtering, there were ", nrow(meth), " probes"))

annotation <- annotation[-which(annotation$name %in% c(probes.to.exclude.p, XY, SNPs.and.controls)),]
meth <- subset(meth, row.names(meth) %in% annotation$name)

print(paste("After filtering, there are now ", nrow(meth), " probes."))
print(paste(length(probes.to.exclude.p), " removed because they have a high detection p-value."))
print(paste(length(XY), " removed because they are on sex chromosomes."))
print(paste(length(SNPs.and.controls), " removed because they are controls."))

rm(pvalue_over_0.05, count_over_0.05, probes.to.exclude.p, XY, SNPs.and.controls, pvals)
################################################################################

#------------------------------------------------------------------------------#
#                               Remove Outliers                                #----
#------------------------------------------------------------------------------#

# Use the Tukey (IQR*3) method
# [NOTE: Please skip this step if it has already been applied in your cohort data.]
log.iqr <- data.frame(cpgs = row.names(meth), NAs.before.IQR3 = rowSums(is.na(meth)))
meth <- IQR.removal(meth)
log.iqr$NAs.after.IQR3 <- rowSums(is.na(meth))
save(log.iqr, file = paste0("MatVegDiet.", study, ".logIQR.", timepoint, ".", analysis.date, ".Rdata"))

#------------------------------------------------------------------------------#
#                        Generate Variables of Interest                        #----
#------------------------------------------------------------------------------#

# "veggie1" & "veggie2" - Diet-based vegetarianism (vegetarian subgroups)

## Data preparation: Generate weekly intake frequency (times/week) of meat & poultry, fish & seafood, egg, and dairy
## [NOTE: Please do NOT add "na.rm = T" in rowSums below, otherwise all NAs will be classified as vegans.]
pheno$meatfish_freq <- rowSums(pheno[, c("meat_freq", "fishseafood_freq")])  # Meat & poultry + fish & seafood
pheno$eggdairy_freq <- rowSums(pheno[, c("egg_freq", "dairy_freq")])  # Egg + dairy
pheno$animal_freq <- rowSums(pheno[, c("meat_freq", "fishseafood_freq", "egg_freq", "dairy_freq")])  # Meat & poultry + fish & seafood + egg + dairy (all animal-based foods)

## Vegetarian subgroups (detailed) - Coded from the least to the most restricted subgroup

### 0 - Non-vegetarian - meat & poultry: no restriction; fish & seafood: no restriction; egg & dairy: no restriction
pheno$VegDiet_subgroup <- 0
pheno$VegDiet_subgroup[is.na(pheno$meat_freq) | is.na(pheno$fishseafood_freq) | is.na(pheno$egg_freq) | is.na(pheno$dairy_freq)] <- NA

### 1 - Pesco-vegetarian - meat & poultry < 1 time/month (i.e., 0.25 time/week); fish & seafood >= 1 time/month; egg & dairy: no restriction
pheno$VegDiet_subgroup[pheno$meat_freq < 0.25 & pheno$fishseafood_freq >= 0.25] <- 1

### 2 - Lacto-ovo-vegetarian - meat & poultry + fish & seafood < 1 time/month; egg + dairy >= 1 time/month
pheno$VegDiet_subgroup[pheno$meatfish_freq < 0.25 & pheno$eggdairy_freq >= 0.25] <- 2

### 3 - Vegan - meat & poultry + fish & seafood + egg + dairy < 1 time/month
pheno$VegDiet_subgroup[pheno$animal_freq < 0.25] <- 3

## Vegetarian subgroups (combined into 3 categories) - Used for surrogate variable analysis below
### 0 - Non-vegetarian
### 1 - Pesco-vegetarian
### 2 - Full vegetarian (including lacto-ovo-vegetarian & vegan)
pheno$VegDiet_3cat <- NA
pheno$VegDiet_3cat[pheno$VegDiet_subgroup == 0] <- 0
pheno$VegDiet_3cat[pheno$VegDiet_subgroup == 1] <- 1
pheno$VegDiet_3cat[pheno$VegDiet_subgroup %in% c(2, 3)] <- 2

## Vegetarianism (binary)

### "veggie1": Full vegetarian vs. non-vegetarian
#### 0 - Non-vegetarian
#### 1 - Full vegetarian
#### [NOTE: In this exposure variable, all pesco-vegetarians are set as NAs.]
pheno$veggie1 <- NA
pheno$veggie1[pheno$VegDiet_3cat == 0] <- 0
pheno$veggie1[pheno$VegDiet_3cat == 2] <- 1

### "veggie2": Full vegetarian + pesco-vegetarian vs. non-vegetarian
#### 0 - Non-vegetarian
#### 1 - Full vegetarian + pesco-vegetarian
pheno$veggie2 <- NA
pheno$veggie2[pheno$VegDiet_3cat == 0] <- 0
pheno$veggie2[pheno$VegDiet_3cat %in% c(1, 2)] <- 1

################################################################################

# PDI / hPDI / uPDI - Plant-based diet indices
# [NOTE: Please skip the whole section if your study is in tier 2 and DO NOT FORGET to remove the parts related to PDI/hPDI/uPDI in subsequent analyses.
#        Tier 1 studies refer to those with data for BOTH "veggie1"/"veggie2" AND PDI/hPDI/uPDI analyses;
#        Tier 2 studies refer to those with data for "veggie1"/"veggie2" analyses ONLY.]

## First keep complete cases only
## [NOTE: This step is VERY IMPORTANT, as we need to ensure that all tertiles are calculated in the final EWAS samples (i.e., complete cases).
##        "veggie1" has more missing data than "veggie2", as it has dropped pesco-vegetarians.
##        Here, we do NOT consider the missing data in "veggie1", but will address them later at the EWAS stage.
##        In ALSPAC, N dropped from 835 (those with both phenotype and methylation data) to 687 (complete cases).]
pheno <- pheno[complete.cases(subset(pheno, select = -c(veggie1))), ]

## Check distribution of 18 food groups
food_group_name <- c("wholegrain", "fruit", "vegetable", "nut", "legume", "vegetableoil", "teacoffee", "fruitjuice", "refinedgrain", "potato", "sugarbeverage", "sweetdessert", "animalfat", "dairy", "egg", "fishseafood", "meat", "misc.animal")

show_quantile_values <- function(var_name) {
  my_var <- pheno[, var_name]
  x <- quantile(my_var, probs = c(0, 0.1, 0.2, 0.25, 0.3, 0.333, 0.4, 0.5, 0.6, 0.667, 0.7, 0.75, 0.8, 0.9, 1), na.rm = T)
  y <- data.frame("Food_group" = var_name,
                  "N" = length(na.omit(my_var)),
                  `Min` = x[1], `P10` = x[2], `P20` = x[3], `P25` = x[4], `P30`  = x[5],
                  `P33.3` = x[6], `P40` = x[7], `P50` = x[8], `P60` = x[9],
                  `P66.7` = x[10], `P70` = x[11], `P75` = x[12], `P80` = x[13], `P90` = x[14], `Max` = x[15])
  return(y)
}

food_group_distribute <- data.frame(matrix(ncol = 17, nrow = 0))

for (var_name in food_group_name) {
  x <- show_quantile_values(var_name)
  food_group_distribute <- rbind(food_group_distribute, x)
}

write.csv(food_group_distribute, paste0("MatVegDiet.", study, ".food.group.distribute.", timepoint, ".", analysis.date, ".csv"), row.names = F)

## Set up functions for positive and negative scoring based on tertiles of each food group
## [NOTE: You may find some food groups in your data have very skewed distribution or do not have enough variation/granularity to calculate tertiles.
##        For example, in ALSPAC, >33% of participants had no intake of nuts, vegetable oil, sugar-sweetened beverages, or miscellaneous animal food.
##        The functions below automatically address this issue by:
##        - If >33% of participants had no intake - 
##          1) Set the first "tertile" to zero (i.e., the lowest intake in your data, which includes >33% of participants);
##          2) Then dividing the second and third "tertiles" using the median of the rest participants (i.e., those with non-zero intakes).
##        - Otherwise - Use the standard method for tertile calculation.
##        If you still fail to get tertiles (i.e., 3 relatively balanced categories) using the method provided, please get in touch.]

### Positive scoring
p_score_3 <- function(var_name) {
  my_var <- pheno[, var_name]
  
  #### First obtain the cutoffs for tertiles
  P0 <- min(my_var, na.rm = T)
  P33.3 <- quantile(my_var, probs = 0.333, na.rm = T)
  P66.7 <- quantile(my_var, probs = 0.667, na.rm = T)
  P100 <- max(my_var, na.rm = T)
  
  #### When >33.3% with the same lowest intake (most likely to be 0)
  if (P0 == P33.3) {
    ##### Create a warning message and introduce the solution
    writeLines(paste0("Warning for ", var_name, ": >33.3% with the same lowest intake (", P0, ");\nSolution: Set all of the lowest as Q1, then take the median of the rest to split Q2 and Q3."))
    
    ##### Modify the range of tertiles and assign positive scores (1-3)
    pheno[, "my_var_3"] <- NA  # Create an empty column
    pheno[which(my_var == P0), "my_var_3"] <- 1  # Set all of those with the lowest intake as Q1
    pheno[which(my_var > P0 & my_var <= median(my_var[my_var > P0], na.rm = T)), "my_var_3"] <- 2  # For the rest, set those with intake <=median as Q2
    pheno[which(my_var > median(my_var[my_var > P0], na.rm = T)), "my_var_3"] <- 3  # And then set those with intake >median as Q3
    my_var_3 <- c(pheno[, "my_var_3"])  # Convert column into vector for subsequent combination
    pheno <- subset(pheno, select = -c(my_var_3))  # Remove the newly created column from the main dataset
  }
  
  #### In other situations - Take tertiles and assign positive scores
  else {
    my_var_3 <- findInterval(my_var, c(-Inf, quantile(my_var, probs = c(0.333, 0.667), na.rm = T), Inf))
  }
  
  #### Return outputs
  return(my_var_3)
}

### Negative scoring
r_score_3 <- function(var_name) {
  my_var_3 <- p_score_3(var_name)  # Same as above for setting tertiles
  my_var_3 <- dplyr::recode(my_var_3, "1" = 3, "2" = 2, "3" = 1)  # Recoded as reverse scores
  return(my_var_3)
}

## Apply function to create tertiles for each of the 18 food groups
for (var_name in food_group_name) {
  assign(paste0(var_name, "_3"), factor(p_score_3(var_name), label = c("T1", "T2", "T3")))
  pheno <- cbind(pheno, get(paste0(var_name, "_3")))
}
for (i in 1:18) {
  colnames(pheno)[ncol(pheno) - 18 + i] <- paste0(food_group_name[i], "_3")  # Add 18 food group tertile variables (columns) to the end of pheno
}
head(pheno)  # Tertile variables should have been added to the last 18 columns

## Set up functions for showing tertile cutoffs and generating descriptive table for the 18 food groups
tertile_cutoff <- function(var_name) {
  my_var <- pheno[, var_name]
  my_var_3 <- p_score_3(var_name)
  T1_lo <- round(min(my_var[my_var_3 == 1], na.rm = T), digits = 2)
  T1_up <- round(max(my_var[my_var_3 == 1], na.rm = T), digits = 2)
  T2_lo <- round(min(my_var[my_var_3 == 2], na.rm = T), digits = 2)
  T2_up <- round(max(my_var[my_var_3 == 2], na.rm = T), digits = 2)
  T3_lo <- round(min(my_var[my_var_3 == 3], na.rm = T), digits = 2)
  T3_up <- round(max(my_var[my_var_3 == 3], na.rm = T), digits = 2)
  x <- paste0("T1: ", T1_lo, "-", T1_up, " | T2: ", T2_lo, "-", T2_up, " | T3: ", T3_lo, "-", T3_up)
  return(x)
}

describe_food_group <- function(var_name) {
  x <- c(`Food group` = var_name,
         "Cutoffs" = tertile_cutoff(var_name),
         table(factor(p_score_3(var_name), label = c("T1", "T2", "T3")), useNA = "always"))
  return(x)
}

## Apply function to generate descriptive table for the 18 food groups
food_group_tertile <- data.frame(matrix(ncol = 6, nrow = 0))
for (var_name in food_group_name) {
  x <- describe_food_group(var_name)
  food_group_tertile <- rbind(food_group_tertile, x)
}
colnames(food_group_tertile) <- c("Food_group", "Cutoffs", "N_T1", "N_T2", "N_T3", "N_missing")
write.csv(food_group_tertile, paste0("MatVegDiet.", study, ".food.group.tertile.", timepoint, ".", analysis.date, ".csv"), row.names = F)

## Apply functions to calculate PDIs

### Overall plant-based diet index (PDI)
pheno$PDI <-
  #### 7 healthy plant food groups assigned positive scores
  p_score_3("wholegrain") + p_score_3("fruit") + p_score_3("vegetable") + p_score_3("nut") + p_score_3("legume") + p_score_3("vegetableoil") + p_score_3("teacoffee") +
  #### 5 less healthy plant food groups assigned positive scores
  p_score_3("fruitjuice") + p_score_3("refinedgrain") + p_score_3("potato") + p_score_3("sugarbeverage") + p_score_3("sweetdessert") +
  #### 6 animal food groups assigned reverse scores
  r_score_3("animalfat") + r_score_3("dairy") + r_score_3("egg") + r_score_3("fishseafood") + r_score_3("meat") + r_score_3("misc.animal")

### Healthful plant-based diet index (hPDI)
pheno$hPDI <-
  #### 7 healthy plant food groups assigned positive scores
  p_score_3("wholegrain") + p_score_3("fruit") + p_score_3("vegetable") + p_score_3("nut") + p_score_3("legume") + p_score_3("vegetableoil") + p_score_3("teacoffee") +
  #### 5 less healthy plant food groups assigned reverse scores
  r_score_3("fruitjuice") + r_score_3("refinedgrain") + r_score_3("potato") + r_score_3("sugarbeverage") + r_score_3("sweetdessert") +
  #### 6 animal food groups assigned reverse scores
  r_score_3("animalfat") + r_score_3("dairy") + r_score_3("egg") + r_score_3("fishseafood") + r_score_3("meat") + r_score_3("misc.animal")

### Unhealthful plant-based diet index (uPDI)
pheno$uPDI <-
  #### 7 healthy plant food groups assigned reverse scores
  r_score_3("wholegrain") + r_score_3("fruit") + r_score_3("vegetable") + r_score_3("nut") + r_score_3("legume") + r_score_3("vegetableoil") + r_score_3("teacoffee") +
  #### 5 less healthy plant food groups assigned positive scores
  p_score_3("fruitjuice") + p_score_3("refinedgrain") + p_score_3("potato") + p_score_3("sugarbeverage") + p_score_3("sweetdessert") +
  #### 6 animal food groups assigned reverse scores
  r_score_3("animalfat") + r_score_3("dairy") + r_score_3("egg") + r_score_3("fishseafood") + r_score_3("meat") + r_score_3("misc.animal")

## Histograms of PDI/hPDI/uPDI
for (my_exp in c("PDI", "hPDI", "uPDI")) {
  jpeg(paste0("MatVegDiet.", study, ".", my_exp, ".histogram.", timepoint, ".", analysis.date, ".jpg"))
  hist(pheno[, my_exp], main = paste0("Histogram of ", my_exp, " (N = ", sum(is.na(pheno[, my_exp]) == F), ")"), xlab = my_exp)
  dev.off()
}

#------------------------------------------------------------------------------#
#                         Generate Surrogate Variables                         #----
#------------------------------------------------------------------------------#

# [NOTE: For consistency, we highly recommend all cohorts use surrogate variable analysis (SVA) to adjust for batch.
#        The code below will generate 20 surrogate variables (SVs) reflecting systematic variations in methylation data caused by unknown/unmodeled/latent sources of noise.
#        These SVs will be included as covariates in EWAS models to remove the impact from batch effects.
#        Please note that SVA is a COMPLETE CASE ANALYSIS and does NOT allow any missing data in the input files.
#        See details about SVA at https://rdrr.io/bioc/sva/f/inst/doc/sva.pdf.
#        If you do not think SVA is applicable to your data, or if you prefer using methods other than/in addition to SVA (e.g., ComBat), please get in touch.]

# Prepare data

## Phenotype data - Keep complete cases and key variables
sva.pheno <- na.omit(pheno[, colnames(pheno) %in% c("sample.id", "VegDiet_3cat", covs.AddModel.PDI, covs.AddModel.hPDI, covs.AddModel.uPDI)])
dim(sva.pheno)  # Sample size will not change if you have already subset to complete cases at the beginning of PDI calculation above

## Methylation data - Recode missing values
sva.meth <- meth[, match(sva.pheno$sample.id, colnames(meth))]  # Match with sample IDs in pheno
k <- which(is.na(sva.meth), arr.ind = T)  # Find indices of NAs
sva.meth[k] <- rowMedians(sva.meth, na.rm = T)[k[, 1]]  # Replace NAs with the row median

## Generate the full model matrix (includes all adjustment variables and exposure variables of interest)
## [NOTE: Here we generate 1 set of 20 SVs for all 5 exposure variables.
##        To avoid including missing data in veggie1 (pesco-vegetarians set as NAs), we use the 3-category vegetarian subgroup variable ("VegDiet_3cat") to combine the information in veggie1 and veggie2.]
mod <- model.matrix(reformulate(paste0("sva.pheno$", colnames(sva.pheno)[which(colnames(sva.pheno) %in% c("VegDiet_3cat", covs.AddModel.PDI, covs.AddModel.hPDI, covs.AddModel.uPDI))])))

## Generate the null model matrix (contains only the adjustment variables)
mod0 <- mod[, -which(colnames(mod) %in% c("sva.pheno$VegDiet_3cat", "sva.pheno$PDI", "sva.pheno$hPDI", "sva.pheno$uPDI"))]

# Run SVA
# [NOTE: Troubleshooting tips:
#        - If you see the error message "Error in eigen(t(resid) %*% resid) : infinite or missing values in 'x'",
#          please make sure that all columns in sva.meth have been matched with sample.id in sva.pheno with the same order.
#        - If you see the error message "Error in solve.default(t(mod) %*% mod) : Lapack routine dgesv: system is exactly singular: U[19,19] = 0",
#          please make sure that your input files (sva.meth, mod, mod0) do not have any missing data or duplicated columns.
#        Please get in touch if you have any difficulties running SVA.
sva.res <- sva(sva.meth, mod = mod, mod0 = mod0, n.sv = 20)

# Combine data
SVs <- as.data.frame(sva.res$sv)
colnames(SVs) <- paste0("sv", 1:ncol(SVs))  # 20 SVs named as sv1-sv20
SVs$sample.id <- sva.pheno$sample.id
pheno <- merge(pheno, SVs, by = "sample.id")  # [NOTE: This step should NOT reduce the sample size, as we have already subset the data to complete cases.]
saveRDS(pheno, paste0("MatVegDiet.", study, ".pheno.SV.", timepoint, ".", analysis.date, ".rds"))

#------------------------------------------------------------------------------#
#                                Match & Check                                 #----
#------------------------------------------------------------------------------#

# Include identical individuals in pheno and meth files
common <- intersect(colnames(meth), pheno$sample.id)
pheno2 <- pheno[pheno$sample.id %in% common, ]
sum(is.na(pheno2[, colnames(pheno2) != "veggie1"]))  # Should be 0 (i.e., no missing data in all variables except for "veggie1")
meth2 <- meth[, (match(pheno2$sample.id, colnames(meth)))]

# Double check

## Question 1: Are the number of samples identical between pheno2 and meth2 - Should be YES
if (!dim(pheno2)[1] == dim(meth2)[2]) {
  print("Warning: The number of individuals NOT identical between pheno2 and meth2 - PLEASE FIX")
} else {
  print("OK: The SAME number of individuals in pheno2 and meth2 - PLEASE CONTINUE")
}

## Question 2: Are pheno2 and meth2 ordered in the same way? - Should be YES
if (!identical(as.character(pheno2$sample.id), as.character(colnames(meth2))))  {
  print("Warning: Ordering of meth2 and pheno2 NOT identical, or variable type (factor vs, character?) NOT matched - PLEASE CHECK AND FIX")
} else {
  print("OK: Individuals in the SAME order in pheno2 and meth2 - PLEASE CONTINUE")
}

## Question 3: Is meth2 a matrix? - Should be YES
if (!is.matrix(meth2)) {
  print("Warning: meth2 is NOT a matrix - PLEASE FIX VIA: meth2 <- as.matrix(meth2)")
} else {
  print("OK: meth2 is a matrix - PLEASE CONTINUE")
}

## Question 4: Is pheno2 a data frame? - Should be YES
if (!is.data.frame(pheno2)) {
  print("Warning: pheno2 is NOT a dataframe - PLEASE FIX VIA: pheno2 <- as.data.frame(pheno2)")
} else {
  print("OK: pheno2 is a dataframe - PLEASE CONTINUE")
}

# Rename for subsequent analyses
pheno <- pheno2
meth <- meth2
rm(pheno2, meth2)

#------------------------------------------------------------------------------#
#                                EWAS Analyses                                 #----
#------------------------------------------------------------------------------#

# Run EWASes for each exposure variable

## "veggie1" (full vegetarian vs. non-vegetarian)
## [NOTE: We need to generate separate datasets for this exposure variable to remove missing data (pesco-vegetarians).]
pheno2 <- pheno[-which(is.na(pheno$veggie1)), ]  # Keep complete cases for "veggie1" (in ALSPAC, N further dropped from 687 to 666), otherwise an error will be reported due to missing data in pesco-vegetarians
meth2 <- meth[, (match(pheno2$sample.id, colnames(meth)))]
identical(pheno2$sample.id, colnames(meth2))  # Should be TRUE

ewas.res.NoCellModel.veggie1 <- ewas.function(meth2, pheno2[, colnames(pheno) %in% covs.NoCellModel.veggie1], variable.of.interest = "veggie1")
ewas.res.MinModel.veggie1 <- ewas.function(meth2, pheno2[, colnames(pheno) %in% covs.MinModel.veggie1], variable.of.interest = "veggie1")
ewas.res.FullModel.veggie1 <- ewas.function(meth2, pheno2[, colnames(pheno) %in% covs.FullModel.veggie1], variable.of.interest = "veggie1")
ewas.res.AddModel.veggie1 <- ewas.function(meth2, pheno2[, colnames(pheno) %in% covs.AddModel.veggie1], variable.of.interest = "veggie1")

save(list = intersect(ls(), c("ewas.res.NoCellModel.veggie1", "ewas.res.MinModel.veggie1", "ewas.res.FullModel.veggie1", "ewas.res.AddModel.veggie1")),
     file = paste0("MatVegDiet.", study, ".veggie1.EWASres.", timepoint, ".", analysis.date, ".Rdata"))

## "veggie2" (full vegetarian + pesco-vegetarian vs. non-vegetarian)
ewas.res.NoCellModel.veggie2 <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.NoCellModel.veggie2], variable.of.interest = "veggie2")
ewas.res.MinModel.veggie2 <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.MinModel.veggie2], variable.of.interest = "veggie2")
ewas.res.FullModel.veggie2 <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.FullModel.veggie2], variable.of.interest = "veggie2")
ewas.res.AddModel.veggie2 <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.AddModel.veggie2], variable.of.interest = "veggie2")

save(list = intersect(ls(), c("ewas.res.NoCellModel.veggie2", "ewas.res.MinModel.veggie2", "ewas.res.FullModel.veggie2", "ewas.res.AddModel.veggie2")),
     file = paste0("MatVegDiet.", study, ".veggie2.EWASres.", timepoint, ".", analysis.date, ".Rdata"))

## PDI
ewas.res.NoCellModel.PDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.NoCellModel.PDI], variable.of.interest = "PDI")
ewas.res.MinModel.PDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.MinModel.PDI], variable.of.interest = "PDI")
ewas.res.FullModel.PDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.FullModel.PDI], variable.of.interest = "PDI")
ewas.res.AddModel.PDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.AddModel.PDI], variable.of.interest = "PDI")

save(list = intersect(ls(), c("ewas.res.NoCellModel.PDI", "ewas.res.MinModel.PDI", "ewas.res.FullModel.PDI", "ewas.res.AddModel.PDI")),
     file = paste0("MatVegDiet.", study, ".PDI.EWASres.", timepoint, ".", analysis.date, ".Rdata"))

## hPDI
ewas.res.NoCellModel.hPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.NoCellModel.hPDI], variable.of.interest = "hPDI")
ewas.res.MinModel.hPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.MinModel.hPDI], variable.of.interest = "hPDI")
ewas.res.FullModel.hPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.FullModel.hPDI], variable.of.interest = "hPDI")
ewas.res.AddModel.hPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.AddModel.hPDI], variable.of.interest = "hPDI")

save(list = intersect(ls(), c("ewas.res.NoCellModel.hPDI", "ewas.res.MinModel.hPDI", "ewas.res.FullModel.hPDI", "ewas.res.AddModel.hPDI")),
     file = paste0("MatVegDiet.", study, ".hPDI.EWASres.", timepoint, ".", analysis.date, ".Rdata"))

## uPDI
ewas.res.NoCellModel.uPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.NoCellModel.uPDI], variable.of.interest = "uPDI")
ewas.res.MinModel.uPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.MinModel.uPDI], variable.of.interest = "uPDI")
ewas.res.FullModel.uPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.FullModel.uPDI], variable.of.interest = "uPDI")
ewas.res.AddModel.uPDI <- ewas.function(meth, pheno[, colnames(pheno) %in% covs.AddModel.uPDI], variable.of.interest = "uPDI")

save(list = intersect(ls(), c("ewas.res.NoCellModel.uPDI", "ewas.res.MinModel.uPDI", "ewas.res.FullModel.uPDI", "ewas.res.AddModel.uPDI")),
     file = paste0("MatVegDiet.", study, ".uPDI.EWASres.", timepoint, ".", analysis.date, ".Rdata"))

#------------------------------------------------------------------------------#
#                             Descriptive Results                              #----
#------------------------------------------------------------------------------#

# Table 1 - Summary statistics of phenotypes
tableone <- as.data.frame(print(CreateTableOne(
  data = pheno[, !colnames(pheno) %in% c("sample.id", "veggie1", "veggie2")],  # Keep "VegDiet_3cat" only for vegetarian subgroups
  factorVars = c("VegDiet_subgroup", "VegDiet_3cat", "sex", "mat.edu", "mat.parity", "mat.smoking", "mat.suppl"),
  test = F)))

write.csv(tableone, paste0("MatVegDiet.", study, ".tableone.", timepoint, ".", analysis.date, ".csv"), row.names = T)

# Intakes of the 18 food groups and PDIs by vegetarian subgroups
# [NOTE: You still need to generate this table even if your study is in tier 2 (please remove the variables not available in your study from below).]
table_by_VegDiet <- as.data.frame(print(CreateTableOne(
  data = subset(pheno, select = c("wholegrain", "fruit", "vegetable", "nut", "legume",
                                  "vegetableoil", "teacoffee", "fruitjuice", "refinedgrain", "potato",
                                  "sugarbeverage", "sweetdessert", "animalfat", "dairy", "egg",
                                  "fishseafood", "meat", "misc.animal", "PDI", "hPDI", "uPDI", "VegDiet_subgroup")),
  vars = c("wholegrain", "fruit", "vegetable", "nut", "legume",
           "vegetableoil", "teacoffee", "fruitjuice", "refinedgrain", "potato",
           "sugarbeverage", "sweetdessert", "animalfat", "dairy", "egg",
           "fishseafood", "meat", "misc.animal", "PDI", "hPDI", "uPDI"),
  strata = "VegDiet_subgroup"),
  test = F))

write.csv(table_by_VegDiet, paste0("MatVegDiet.", study, ".ByVegDiet.", timepoint, ".", analysis.date, ".csv"), row.names = T)

# Correlation matrix between PDIs and intakes of the 18 food groups
correlogram <- GGally::ggcorr(data = subset(pheno, select = c("misc.animal", "meat", "fishseafood", "egg", "dairy",
                                                              "animalfat", "sweetdessert", "sugarbeverage", "potato", "refinedgrain",
                                                              "fruitjuice", "teacoffee", "vegetableoil", "legume", "nut",
                                                              "vegetable", "fruit", "wholegrain", "uPDI", "hPDI", "PDI")), 
                              method = c("pairwise", "pearson"), label = T, label_alpha = T, label_round = 2)

ggplot2::ggsave(correlogram, file = paste0("MatVegDiet.", study, ".correlogram.", timepoint, ".", analysis.date, ".png"), height = 16, width = 18)

#------------------------------------------------------------------------------#
#                             Q-Q Plots & Lambdas                              #----
#------------------------------------------------------------------------------#

# Generate Q-Q plots
model_list <- c("NoCellModel.veggie1", "MinModel.veggie1", "FullModel.veggie1", "AddModel.veggie1",
                "NoCellModel.veggie2", "MinModel.veggie2", "FullModel.veggie2", "AddModel.veggie2",
                "NoCellModel.PDI", "MinModel.PDI", "FullModel.PDI", "AddModel.PDI",
                "NoCellModel.hPDI", "MinModel.hPDI", "FullModel.hPDI", "AddModel.hPDI",
                "NoCellModel.uPDI", "MinModel.uPDI", "FullModel.uPDI", "AddModel.uPDI")

for (my_mod in model_list) {
  my_res <- get(paste0("ewas.res.", my_mod))
  jpeg(file = paste0("MatVegDiet.", study, ".", my_mod, ".QQ.", timepoint, ".", analysis.date, ".jpg"))
  qq(my_res$p)
  dev.off()
}

# Calculate lambdas
lambda_cal <- function(my_mod) {
  my_res <- get(paste0("ewas.res.", my_mod))
  lambda <- median(qchisq(my_res$p, df = 1, lower.tail = F), na.rm = T) / qchisq(0.5, 1)
  x <- c(my_mod, lambda)
  return(x)
}

lambda_res <- data.frame(matrix(ncol = 2, nrow = 0))
for (my_mod in model_list) {
  x <- lambda_cal(my_mod)
  lambda_res <- rbind(lambda_res, x)
}
colnames(lambda_res) <- c("Model", "Lambda")
write.csv(lambda_res, paste0("MatVegDiet.", study, ".lambda.", timepoint, ".", analysis.date, ".csv"), row.names = F)

#------------------------------------------------------------------------------#
#                   Differentially Methylated Regions (DMR)                    #----
#------------------------------------------------------------------------------#

# Loop through all models and exposure variables to run DMR analysis
for (my_mod in c("NoCellModel", "MinModel", "FullModel", "AddModel")) {
  for (my_exp in c("veggie1", "veggie2", "PDI", "hPDI", "uPDI")) {
    ## Set up input data
    stats <- get(paste0("ewas.res.", my_mod, ".", my_exp))
    ifelse(my_exp == "veggie1", betas <- meth2, betas <- meth)  # meth2 (i.e., meth excluding pesco-vegetarians) was used for "veggie1" EWAS above
    
    ## Match genome annotation
    annotation <- meffil::meffil.get.features("450k")
    annotation <- annotation[match(rownames(betas), annotation$name), ]
    order.idx <- order(annotation$chromosome, annotation$position)
    betas <- betas[order.idx, ]
    annotation <- annotation[order.idx, ]
    stats <- stats[order.idx, ]
    
    ## Run DMR analysis
    dmr.res <- dmrff(estimate = stats$coef, se = stats$se, p.value = stats$p, methylation = betas, chr = annotation$chromosome, pos = annotation$position, maxgap = 500, verbose = T)
    
    ## Select and save top findings
    dmr.res <- dmr.res[which(dmr.res$p.adjust < 0.05 & dmr.res$n > 1), ]
    print(dmr.res)
    write.csv(dmr.res, paste0("MatVegDiet.", study, ".", my_mod, ".", my_exp, ".dmr.", timepoint, ".", analysis.date, ".csv"), row.names = F)
  }
}

#------------------------------------------------------------------------------#
#                                 Output List                                  #----
#------------------------------------------------------------------------------#

# Lastly, please check if you have produced the following output files in the specified formats:
# - Rdata files (5 for tier 1 studies, 2 for tier 2 studies) for EWAS results
#   - Each contains outputs from all four EWAS models for one exposure variable
#   - Named as: “MatVegDiet.[study].[my_exp].EWASres.[timepoint].[analysis.date].Rdata”
# - One Rdata files for the log of the IQR*3 (Tukey) method
#   - Named as: “MatVegDiet.[study].logIQR.[timepoint].[analysis.date].Rdata”
# - One csv file showing the distribution of the 18 food groups
#   - Named as: “MatVegDiet.[study].food.group.distribute.[timepoint].[analysis.date].csv”
#   - Not applicable for tier 2 studies
# - One csv file showing the tertile cutoffs of the 18 food groups
#   - Named as: “MatVegDiet.[study].food.group.tertile.[timepoint].[analysis.date].csv”
#   - Not applicable for tier 2 studies
# - One csv file showing summary statistics of phenotype data
#   - Named as: “MatVegDiet.[study].tableone.[timepoint].[analysis.date].csv”
# - One csv file showing the intake of each available food group by vegetarian subgroups
#   - Named as: “MatVegDiet.[study].ByVegDiet.[timepoint].[analysis.date].csv”
# - One correlogram showing the correlation matrix between PDIs and intakes of the 18 food groups
#   - Named as: “MatVegDiet.[study].correlogram.[timepoint].[analysis.date].png”
#   - Not applicable for tier 2 studies
# - Three histograms for PDI, hPDI, and uPDI
#   - Named as: “MatVegDiet.[study].[my_exp].histogram.[timepoint].[analysis.date].jpg”
#   - Not applicable for tier 2 studies
# - Q-Q plots (20 for tier 1 studies, 8 for tier 2 studies)
#   - Separate plots will be created for each exposure variable and model
#   - Named as: “MatVegDiet.[study].[MODEL].[my_exp].QQ.[timepoint].[analysis.date].jpg”
# - One csv file showing lambdas for all EWAS models
#   - Named as: “MatVegDiet.[study].lambda.[timepoint].[analysis.date].csv”
# - csv files (20 for tier 1 studies, 8 for tier 2 studies) for DMR results
#   - Based on the fully adjusted model for each exposure variable
#   - Named as: “MatVegDiet.[study].[my_exp].dmr.[timepoint].[analysis.date].csv”

# Please do not forget to complete and return the Excel file named "MatVegDiet.[STUDY].cohortinfo.xlsx" (including 2 tabs).

# Massive thanks for your contribution!
