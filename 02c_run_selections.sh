
# Parameters that need to be adjusted each run
# Choose config file and output directory (give OUTPUT_DIR an experiment specific name)
CONFIG_FILE="./02b_selection_configs/01_selection_config_Mad20_basic_methods_n150.yaml"
OUTPUT_DIR="./results/01_selection"

STEP=0 #1, 2, for a given config go through 0, 1, 2 after each of them is finished (step 1 is parallelised)

#--------------------------------------------------------------------------------------

# Script
YAMLS_DIR="${OUTPUT_DIR}/tmp_yamls"
SELECTIONS_DIR="${OUTPUT_DIR}/tmp_select"
JOBS_DIR="${OUTPUT_DIR}/tmp_jobs"

PARTITION="cpu_p"

mkdir $OUTPUT_DIR
mkdir $YAMLS_DIR
mkdir $SELECTIONS_DIR
mkdir $JOBS_DIR

# Process 1, Generate yamls for each selection
sleep 0.1
#STEP=0
if [ $STEP -eq 0 ]
then
job_file="${JOBS_DIR}/selection_step_${STEP}.cmd"
echo "#!/bin/bash
#SBATCH -J selection_step_${STEP}
#SBATCH -o ${JOBS_DIR}/selection_step_${STEP}.out
#SBATCH -e ${JOBS_DIR}/selection_step_${STEP}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python select_script.py ${STEP} ${CONFIG_FILE} ${YAMLS_DIR}
" > ${job_file}
sbatch $job_file
fi

# Processes 2 select for each configuration. (parallelize over yaml files)
#STEP=1
if [ $STEP -eq 1 ]
then
COUNT=1
for yaml in "$YAMLS_DIR"/* # yaml is the full file path
do
    sleep 0.1
    job_file="${JOBS_DIR}/selection_step_${STEP}_${COUNT}.cmd"
    echo "#!/bin/bash
#SBATCH -J selection_step_${STEP}_${yaml}
#SBATCH -o ${JOBS_DIR}/selection_step_${STEP}_${COUNT}.out
#SBATCH -e ${JOBS_DIR}/selection_step_${STEP}_${COUNT}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python select_script.py ${STEP} ${yaml} ${SELECTIONS_DIR}
" > ${job_file}
    sbatch $job_file
    COUNT=$((COUNT+1))
done
fi

# Process 3 Summarize results
sleep 0.1
#STEP=2
if [ $STEP -eq 2 ]
then
job_file="${JOBS_DIR}/selection_step_${STEP}.cmd"
echo "#!/bin/bash
#SBATCH -J selection_step_${STEP}
#SBATCH -o ${JOBS_DIR}/selection_step_${STEP}.out
#SBATCH -e ${JOBS_DIR}/selection_step_${STEP}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 01:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python select_script.py ${STEP} ${SELECTIONS_DIR} ${OUTPUT_DIR}
" > ${job_file}
sbatch $job_file
fi

#rm -rf $SELECTIONS_DIR
#rm -rf $JOBS_DIR
#rm -rf $YAMLS_DIR