#!/bin/bash -l
#SBATCH --job-name=DMR_meta
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=dmr_analysis_%j.out
#SBATCH --error=dmr_analysis_%j.err

# Load R module
module load languages/R/4.3.3

# Run R scripts
Rscript /user/work/zd20208/MatVegDiet_PACE_EWAS/scripts/05\(2\)_DMR_run_dmrff.meta.R

