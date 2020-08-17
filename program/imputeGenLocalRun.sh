#!/bin/bash

##########################################################
############## EXPORT PARAMETER FILE VARIABLES ###########
##########################################################
export REF_GENOME=$(readlink -f $(cat $PARAMETERS | grep -w "|REF_GENOME" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export REF_PANEL_DIR=$(readlink -f $(cat $PARAMETERS | grep -w "|REF_PANEL_DIR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export REF_PANEL_PREFIX=$(cat $PARAMETERS | grep -w "|REF_PANEL_PREFIX" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export VAL_PANEL_DIR=$(readlink -f $(cat $PARAMETERS | grep -w "|VAL_PANEL_DIR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export VAL_PANEL_PREFIX=$(cat $PARAMETERS | grep -w "|VAL_PANEL_PREFIX" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export GENO_PANEL=$(readlink -f $(cat $PARAMETERS | grep -w "|GENO_PANEL" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export BAM_DIR=$(readlink -f $(cat $PARAMETERS | grep -w "|BAM_DIR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export PEDIGREE=$(readlink -f $(cat $PARAMETERS | grep -w "|PEDIGREE" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export FINDHAPF90=$(readlink -f $(cat $PARAMETERS | grep -w "|FINDHAPF90" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export GENETIC_MAP=$(cat $PARAMETERS | grep -w "|GENETIC_MAP" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export OUTPUT_DIR=$(readlink -f $(cat $PARAMETERS | grep -w "|OUTPUT_DIR" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1) 2> /dev/null)
export OUTPUT_NAME=$(cat $PARAMETERS | grep -w "|OUTPUT_NAME" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)

export PHASING_PROGRAM=$(cat $PARAMETERS | grep -w "|PHASING_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export IMPUTE_PROGRAM=$(cat $PARAMETERS | grep -w "|IMPUTE_PROGRAM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)

export CHROM_SET=$(cat $PARAMETERS | grep -w "|CHROM_SET" | cut -d "=" -f 2- | cut -d "#" -f 1 | sed -e "s/[[:space:]]\+//g" | tr "," "\n")
export INCLUDE_CHROM_X=$(cat $PARAMETERS | grep -w "|INCLUDE_CHROM_X" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export SELECT_CHROM=$(cat $PARAMETERS | grep "|SELECT_CHROM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1 | tr "," "\n")
export MAF=$(cat $PARAMETERS | grep -w "|MAF" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MISSING=$(cat $PARAMETERS | grep -w "|MISSING" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)

export DOWN_SCALING=$(cat $PARAMETERS | grep -w "|DOWN_SCALING" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export THREADS=$(cat $PARAMETERS | grep -w "|THREADS" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export MEM=$(cat $PARAMETERS | grep -w "|MEM" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
export BATCH=$(cat $PARAMETERS | grep -w "|BATCH" | sed -e "s/[[:space:]]\+//g" | cut -d "=" -f 2- | cut -d "#" -f 1)
##########################################################


### CHECK CONDA ENVIRONMENTS
eval "$(conda shell.bash hook)"
conda activate imputegen
conda activate --stack r-env4.0
export CONDA_DEFAULT_ENV=$CONDA_DEFAULT_ENV
export CONDA=$(conda env list | grep "envs/imputegen" | awk '{print $2}')

if [[ "$CONDA_DEFAULT_ENV" != "imputegen" && "$CONDA_DEFAULT_ENV" != "r-env4.0" ]]; then
   echo
   echo "ERROR: Activate conda imputegen environment before run ImputeGen "
   echo "Usage: conda activate imputegen"
   echo "Please, check if the ImputeGen dependencies are installed (install_imputegen_dependencies.sh file)"
   echo
   exit 1
fi


### CHECK IMPUTEGEN PARAMETERS
source $IMPUTEGEN_DIR/program/runCheckParameters.sh
wait


### LOAD OUTPUT FOLDERS
export SWITCH=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/data/switchRefAlt
export DOWN_SCALED=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/data/downScaled
export FILTER=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/data/filtered
export PHASE=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/data/phased
export IMPUTE=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/data/imputed
export ACC=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/accuracy

mkdir -p $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log 2>&1 > /dev/null
mkdir -p $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/params 2>&1 > /dev/null
export ANALYSIS_ID=$(echo $RANDOM)


##### Run ImputeGen
cp $PARAMETERS $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/params/"$ANALYSIS_ID"_"$(basename $PARAMETERS)"
# https://unix.stackexchange.com/questions/49042/is-it-possible-for-nohup-to-write-the-output-both-to-the-file-nohup-out-and-to-t
# https://stackoverflow.com/questions/8363519/how-do-i-terminate-all-the-subshell-processes
set -m
nohup $IMPUTEGEN_DIR/program/runStepsLocal.sh >/dev/null 2>&1 &
echo $! > $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/pid.txt
echo "Running ImputeGen"
echo "ANALYSIS NAME: $OUTPUT_NAME"
echo "General log: tail -f -n +1 $OUTPUT_DIR/imputeGen_${OUTPUT_NAME}/log/log_${OUTPUT_NAME}.txt"
echo "Check the log files per chromosome in $OUTPUT_DIR/imputeGen_${OUTPUT_NAME}/log"
echo "To open the log file of a chromosome, type on terminal: tail -f -n +1 $OUTPUT_DIR/imputeGen_${OUTPUT_NAME}/log/${OUTPUT_NAME}_CHR.txt"
echo "Just replace CHR by the chromosome name"
echo "To kill this ImputeGen analysis, run the command line: kill -- -$(cat $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/pid.txt)"


