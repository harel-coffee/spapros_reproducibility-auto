#! /bin/bash

# Set Variables
METHOD=$1
CONT_NAME=sc_analysis_r_py_210114
yaml=$2
SELECTIONS_DIR=$3
CWD_var=$4
COMMAND="cd ${CWD_var}; python3 ./external_select_script.py ${yaml} ${SELECTIONS_DIR}"

CH_PARENT=/localscratch/"$USER"/ch/"$SLURM_JOB_ID"
CH_ROOTDIR="$CH_PARENT"/"$CONT_NAME"
CH_HOMEDIR="$CH_PARENT"/ch_home

# Cleanup: Delete all charliecloud containers whose corresponding job is no longer running
mkdir -p /localscratch/$USER/ch/dummy
touch /home/icb/${USER}/slurm_jupyter_dummy.job
for d in /localscratch/$USER/ch/*/ ; do
  CHID=$(echo "$d" | awk -F'/' '{print $5}')
  if ! squeue -j "$CHID"&>/dev/null; then
      rm -rf "${d::-1}"
      rm /home/icb/${USER}/slurm_jupyter_$CHID.job
  fi
done

mkdir -p $CH_PARENT
ch-tar2dir /storage/groups/ml01/workspace/ml_charliecloud_containers/"$(echo "$CONT_NAME" | sed 's/\(.*\)_/\1\//')"/$CONT_NAME.tar.gz $CH_PARENT

mkdir -p $CH_ROOTDIR/home/$USER
mkdir -p $CH_ROOTDIR/home/icb/$USER
mkdir -p $CH_ROOTDIR/storage/groups
mkdir -p $CH_ROOTDIR/storage/scratch

echo "
USER=${USER}
LOGNAME=${USER}
HOME=/home/${USER}
" >> $CH_PARENT/env

mkdir -p "$HOME"/charliecloud_homedirs/

#if [ -e "$HOME"/charliecloud_homedirs/charliehome_"$CONT_NAME"_"$USER".tar.gz ]
if [ -e "$HOME"/charliecloud_homedirs/charliehome_"$METHOD"_"$USER".tar.gz ]
then
echo Existing home directory found. Unpacking...
#tar -xzf "$HOME"/charliecloud_homedirs/charliehome_"$CONT_NAME"_"$USER".tar.gz -C "$CH_PARENT"/
tar -xzf "$HOME"/charliecloud_homedirs/charliehome_"$METHOD"_"$USER".tar.gz -C "$CH_PARENT"/
else
echo No existing charliecloud home directory found. Creating a new one...
mkdir -p $CH_HOMEDIR/.ch_tmp
mkdir -p $CH_HOMEDIR/.local/bin
mkdir -p $CH_HOMEDIR/bin
ln -s /home/icb/$USER $CH_HOMEDIR/server_home
ln -s /storage/scratch $CH_HOMEDIR/storage_scratch
ln -s /storage/groups $CH_HOMEDIR/storage_groups

echo "
# ~/.profile: executed by the command interpreter for login shells.
# This file is not read by bash(1), if ~/.bash_profile or ~/.bash_login exists.

# if running bash, source ~/.bashrc
if [ -n '$BASH_VERSION' ]; then
    # include .bashrc if it exists
    if [ -f '/home/${USER}/.bashrc' ]; then
        . '/home/${USER}/.bashrc'
    fi
fi
" >> $CH_HOMEDIR/.profile


echo "
# .bashrc
# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# Setting environment variables
export USER=${USER}
export LOGNAME=${USER}
export HOME=/home/${USER}
export PATH=/home/${USER}/.local/bin:/home/${USER}/bin:/opt/python/bin:/opt/R/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export LD_LIBRARY_PATH=/opt/R/lib/R/lib
export SHELL=/bin/bash
export TERM=xterm-256color
export LANG=en_US.UTF-8
export TMP=/home/${USER}/.ch_tmp
export TMPDIR=/home/${USER}/.ch_tmp
export TEMP=/home/${USER}/.ch_tmp

# Activate some jupyter notebook extensions
jupyter nbextension enable collapsible_headings/main &> /dev/null
jupyter nbextension enable notify/notify &> /dev/null
jupyter nbextension enable scroll_down/main &> /dev/null
jupyter nbextension enable toc2/main &> /dev/null
jupyter nbextension enable varInspector/main &> /dev/null
jupyter nbextension enable execute_time/ExecuteTime &> /dev/null

# Setting aliases
alias ll='ls -lha --color=auto'
alias jn='jupyter-notebook --no-browser --ip=0.0.0.0'
alias jl='jupyter-lab --no-browser --ip=0.0.0.0'
alias cx='cellxgene launch --host 0.0.0.0 --port 8888'
" >> $CH_HOMEDIR/.bashrc
fi

chmod 600 /home/icb/${USER}/slurm_batch_$SLURM_JOB_ID.job

ch-run $CH_ROOTDIR --no-home --unset-env='*' --set-env=$CH_PARENT/env -b /storage/groups/:/storage/groups -b /storage/scratch/:/storage/scratch -b $HOME:/home/icb/$USER -b $CH_HOMEDIR:/home/$USER --cd=/home/$USER -- /bin/bash -c "source /home/${USER}/.bashrc; ${COMMAND}"

#echo Saving your charliecloud home directory to "$HOME"/charliecloud_homedirs/charliehome_"$CONT_NAME"_"$USER".tar.gz
echo Saving your charliecloud home directory to "$HOME"/charliecloud_homedirs/charliehome_"$METHOD"_"$USER".tar.gz
rm -rf $CH_HOMEDIR/.ch_tmp/*
#tar -czf "$HOME"/charliecloud_homedirs/charliehome_"$CONT_NAME"_"$USER".tar.gz -C "$CH_PARENT"/ ch_home
tar -czf "$HOME"/charliecloud_homedirs/charliehome_"$METHOD"_"$USER".tar.gz -C "$CH_PARENT"/ ch_home
echo Saved your charliecloud home directory. It will be loaded and reused in your next charliecloud session unless you delete it in the meantime.

rm -rf $CH_PARENT
# Do not remove job output files for batch jobs
#rm /home/icb/${USER}/slurm_batch_$SLURM_JOB_ID.job
