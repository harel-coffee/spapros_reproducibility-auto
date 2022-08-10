#!/bin/bash


# Parameters that need to be adjusted each run
# Choose config file and output directory (give OUTPUT_DIR an experiment specific name)
CONFIG_FILE="./04b_evaluation_configs/01_evaluation_config_Mad.yaml"
PROBESETS_CSV="./results/01_selection/selections.csv"
OUTPUT_DIR="./results/01_eval_01_select"

STEP=0 #1, 2, 3 for a given CONFIG_FILE and PROBESETS_CSV go through 0, 1, 2, 3 after each of them is finished


#--------------------------------------------------------------------------------------

CODE_DIR="./" # Abs path needed?
PARTITION="cpu_p"
JOBS_DIR="${OUTPUT_DIR}/tmp_jobs"

# Script

###

# probesets set_ids # TODO: check if the csv file index also has a column name, guess we need to filter that one out in the python script.
PROBESETS=( $(head -n +1 ${PROBESETS_CSV}) )
PROBESETS=${PROBESETS//,/ } # call list via ${PROBESETS[@]}

#
#OUT_PATH=${OUT_PATH_BASE}/grid_searches_gen/${GS_KEY}
#
#rm -rf ${OUT_PATH}/jobs
#rm -rf ${OUT_PATH}/logs
#rm -rf ${OUT_PATH}/results
mkdir -p ${OUTPUT_DIR}
mkdir -p ${JOBS_DIR}
#mkdir -p ${OUT_PATH}/logs
#mkdir -p ${OUT_PATH}/results


## step 0
#STEP=0
if [ $STEP -eq 0 ]
then
for metric in "cluster_similarity" "knn_overlap" "gene_corr-marker_corr" 
do
    sleep 0.1
    job_file="${JOBS_DIR}/eval_step_${STEP}_${metric}.cmd"
    echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${JOBS_DIR}/eval_step_${STEP}_${metric}.out
#SBATCH -e ${JOBS_DIR}/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_DIR}/eval_script.py ${OUTPUT_DIR} ${CONFIG_FILE} ${PROBESETS_CSV} ${metric} ${STEP} ${PROBESETS[@]}
" > ${job_file}
    sbatch $job_file      
#echo ${CODE_DIR}/eval_script.py ${OUTPUT_DIR} ${METRICS_YAML} ${DATA_YAML} ${PROBESETS_CSV} ${metric} 0 ${PROBESETS[@]}
done
fi

## step 1
#STEP=1
if [ $STEP -eq 1 ]
then
for metric in "cluster_similarity" "knn_overlap"
    do
    for set_id in ${PROBESETS[@]}
        do
        sleep 0.1
        job_file="${JOBS_DIR}/eval_step_${STEP}_${metric}_${set_id}.cmd"
        echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${JOBS_DIR}/eval_step_${STEP}_${metric}_${set_id}.out
#SBATCH -e ${JOBS_DIR}/eval_step_${STEP}_${metric}_${set_id}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_DIR}/eval_script.py ${OUTPUT_DIR} ${CONFIG_FILE} ${PROBESETS_CSV} ${metric} ${STEP} ${set_id}
" > ${job_file}
        sbatch $job_file
    done
done
fi

## step 2
#STEP=2
if [ $STEP -eq 2 ]
then
for metric in "forest_clfs"
    do
    for set_id in ${PROBESETS[@]}
        do
        sleep 0.1
        job_file="${JOBS_DIR}/eval_step_${STEP}_${metric}_${set_id}.cmd"
        echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${JOBS_DIR}/eval_step_${STEP}_${metric}_${set_id}.out
#SBATCH -e ${JOBS_DIR}/eval_step_${STEP}_${metric}_${set_id}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_DIR}/eval_script.py ${OUTPUT_DIR} ${CONFIG_FILE} ${PROBESETS_CSV} ${metric} ${STEP} ${set_id}
" > ${job_file}
        sbatch $job_file
    done
done
fi

#STEP=2
if [ $STEP -eq 2 ]
then
for metric in "cluster_similarity" "knn_overlap" "gene_corr-marker_corr"
    do
    sleep 0.1
    job_file="${JOBS_DIR}/eval_step_${STEP}_${metric}.cmd"
    echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${JOBS_DIR}/eval_step_${STEP}_${metric}.out
#SBATCH -e ${JOBS_DIR}/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_DIR}/eval_script.py ${OUTPUT_DIR} ${CONFIG_FILE} ${PROBESETS_CSV} ${metric} ${STEP} ${PROBESETS[@]}
" > ${job_file}
    sbatch $job_file
done
fi

#step 3
#STEP=3
if [ $STEP -eq 3 ]
then
for metric in "cluster_similarity-knn_overlap-forest_clfs-gene_corr-marker_corr"
    do
    sleep 0.1
    job_file="${JOBS_DIR}/eval_step_${STEP}_${metric}.cmd"
    echo "#!/bin/bash
#SBATCH -J eval_step_${STEP}_${metric}
#SBATCH -o ${JOBS_DIR}/eval_step_${STEP}_${metric}.out
#SBATCH -e ${JOBS_DIR}/eval_step_${STEP}_${metric}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000
source ${HOME}/.bash_profile
conda activate spapros
python3 ${CODE_DIR}/eval_script.py ${OUTPUT_DIR} ${CONFIG_FILE} ${PROBESETS_CSV} ${metric} ${STEP} ${PROBESETS[@]}
" > ${job_file}
    sbatch $job_file
done
fi
