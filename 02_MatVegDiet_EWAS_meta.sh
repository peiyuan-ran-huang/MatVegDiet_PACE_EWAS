#!/bin/bash
#SBATCH --job-name=MatVegDiet_Meta
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=50000M

#### 0. Define your base directories and arrays ####

BASE_DIR=/user/home/zd20208/MatVegDiet_PACE_EWAS
RESULTS_DIR=${BASE_DIR}/results
META_DIR=${RESULTS_DIR}/meta

# cohorts to include (skip BiB & DCHS since they’re empty)
COHORTS=(ALSPAC GenR INMA MoBa1 MoBa2 MoBa4 MoBa8 NORTHPOP Viva)

# all exposures and models you want to meta-analyze
EXPOSURES=(PDI hPDI uPDI veggie1 veggie2)
MODELS=(AddModel FullModel MinModel NoCellModel)

#### 1. Prepare the meta output dir and METAL command file ####

mkdir -p "${META_DIR}"
CMD_FILE=${META_DIR}/metal_commands.txt

# Write global METAL settings
cat > "${CMD_FILE}" <<EOF
REMOVEFILTERS
SCHEME STDERR
USESTRAND OFF
GENOMICCONTROL OFF
AVERAGEFREQ OFF
MINMAXFREQ OFF
COLUMNCOUNTING LENIENT
SEPARATOR COMMA
PVALUELABEL p
MARKER probeid
EFFECTLABEL coef
STDERRLABEL se
EOF

#### 2. Loop over each exposure × model ####
for exp in "${EXPOSURES[@]}"; do
  for mdl in "${MODELS[@]}"; do

    # Count how many cohort CSVs actually exist for this combo
    count=0
    for coh in "${COHORTS[@]}"; do
      CSV="${RESULTS_DIR}/${coh}/${coh}.ewas.res.${mdl}.${exp}.csv"
      if [ -f "${CSV}" ]; then
        ((count++))
      fi
    done

    # If none found, skip entire block
    if [ "${count}" -eq 0 ]; then
      echo "## Skipping ${mdl}.${exp}: no files found" >> "${CMD_FILE}"
      continue
    fi

    # Otherwise, write a human-readable header and the OUTFILE line
    echo -e "\n# ===== ${exp}: ${mdl} (${count} cohorts) =====" >> "${CMD_FILE}"
    echo "OUTFILE ${META_DIR}/ewas.res.${mdl}.${exp} .txt" >> "${CMD_FILE}"

    # Add PROCESSFILE for each existing CSV, and note missing ones
    for coh in "${COHORTS[@]}"; do
      CSV="${RESULTS_DIR}/${coh}/${coh}.ewas.res.${mdl}.${exp}.csv"
      if [ -f "${CSV}" ]; then
        echo "PROCESSFILE      ${CSV}" >> "${CMD_FILE}"
      else
        echo "## Missing ${CSV}, skipping" >> "${CMD_FILE}"
      fi
    done

    # Finally run the random‐effects analysis and clear
    echo "ANALYZE RANDOM" >> "${CMD_FILE}"
    echo "CLEAR"         >> "${CMD_FILE}"
    echo "CLEAR"         >> "${CMD_FILE}"

  done
done

#### 3. Compile & run random-metal ####

module add tools/git/2.42.0-a3he       # load git with modules
git clone https://github.com/explodecomputer/random-metal "${BASE_DIR}/random-metal"
cd "${BASE_DIR}/random-metal"
make

# Execute METAL once with our dynamically generated commands
cd executables
./metal "${CMD_FILE}"
