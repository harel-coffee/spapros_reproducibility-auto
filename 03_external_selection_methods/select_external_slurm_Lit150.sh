#!/bin/bash

###############
# Description #
###############
# The script runs all selections that are defined with config yaml files in $YAMLS_DIR.
# Selection results are saved in $SELECTIONS_DIR.
# 1. Read all yaml files
# 2. for each yaml file:
#    2.1 Read the method
#    2.2 Activate the according conda env
#    2.3 Run the according selection script
###############

# Parameters

YAMLS_DIR="./config_yamls_Lit150"
ENV_DIR="./environments"
SELECTIONS_DIR="./external_selection/selections_Lit150"
JOBS_DIR="./external_selection/jobs_Lit150"
PARTITION="cpu_p"

################
# Define your subset of methods

METHODS=(asfs cosg nsforest scgenefit scmer smash triku)
R_METHODS=(selfe)  # scpnmf genebasis

################

CWD_var=$(pwd)

source ${HOME}/.bash_profile
conda activate spapros

mkdir -p $SELECTIONS_DIR
mkdir -p $JOBS_DIR


# Function to read method from yaml
read_yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

COUNT=1
for yaml in "$YAMLS_DIR"/* # yaml is the full file path
do
    sleep 0.1
    conda activate "${ENV_DIR}""/venv_yaml"
    METHOD_ENV=$(read_yaml "${yaml}" "['venv']")
    METHOD=$(read_yaml "${yaml}" "['method']")
    if [[ " ${METHODS[*]} " =~ " ${METHOD} " ]]; then
    job_file="${JOBS_DIR}/selection_external_${COUNT}.cmd"
    echo "#!/bin/bash
#SBATCH -J ${COUNT}_selection_external_${yaml}
#SBATCH -o ${JOBS_DIR}/selection_external_${COUNT}.out
#SBATCH -e ${JOBS_DIR}/selection_external_${COUNT}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000   
source ${HOME}/.bash_profile
conda activate "${ENV_DIR}""/${METHOD_ENV}"
python3 external_select_script.py ${yaml} ${SELECTIONS_DIR}
" > ${job_file}
    sbatch $job_file
fi
    if [[ " ${R_METHODS[*]} " =~ " ${METHOD} " ]]; then
    job_file="${JOBS_DIR}/selection_external_${COUNT}.cmd"
    echo "#!/bin/bash
#SBATCH -J ${COUNT}_selection_external_${yaml}
#SBATCH -o ${JOBS_DIR}/selection_external_${COUNT}.out
#SBATCH -e ${JOBS_DIR}/selection_external_${COUNT}.err
#SBATCH -p ${PARTITION}
#SBATCH -t 2-00:00:00
#SBATCH -c 12
#SBATCH --mem=50G
#SBATCH --nice=10000 
echo start script
source ${HOME}/.bash_profile
source selection_script_charliecloud.sh ${METHOD} ${yaml} ${SELECTIONS_DIR} ${CWD_var}
echo script finished
" > ${job_file}
    sbatch $job_file    
fi
    COUNT=$((COUNT+1))
done


