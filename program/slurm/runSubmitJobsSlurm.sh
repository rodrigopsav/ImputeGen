#!/bin/bash
#IMPORTANT: NEVER LEAVE A SPACE AFTER ANY \ WHEN USING SBATCH


##### DOWNSCALING, FILTERING, PHASING, IMPUTATION
for REGION in $SELECT_CHROM; do
   export REGION=$REGION
   export CHR=$(echo $REGION | awk -F":" '{print $1}')
   export BPSTART=$(echo $REGION | awk -F":" '{print $2}' | awk -F"-" '{print $1}')
   export BPEND=$(echo $REGION | awk -F":" '{print $2}' | awk -F"-" '{print $2}')
   
   export REF_PANEL=$REF_PANEL_DIR/${REF_PANEL_PREFIX}${CHR}${EXTENSION}
   if [[ $VAL_PANEL_DIR != "none" ]]; then export VAL_PANEL=$VAL_PANEL_DIR/${VAL_PANEL_PREFIX}${CHR}${EXTENSION_VAL}; fi
   
   if [[ -z "$BPSTART" ]]; then export BPSTART=0; fi
   if [[ -z "$BPEND" ]]; then export BPEND=1000000000; fi
   
   if [[ "$BPSTART" == 0 && "$BPEND" == 1000000000 ]]; then
      export FileExtension=$CHR
   else
      export FileExtension=$CHR_"$BPSTART"_"$BPEND"
   fi   
   
   JOBID1=$(sbatch --job-name="${ANALYSIS_ID}"_"$OUTPUT_NAME"_"$FileExtension" --time=$IMPUTE_TIME --cpus-per-task=$IMPUTE_CPU --mem=$IMPUTE_MEM \
   --output $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/logSlurm/"$FileExtension".txt \
   --export=CONDA=$CONDA,IMPUTEGEN_DIR=$IMPUTEGEN_DIR,REF_GENOME=$REF_GENOME,REF_PANEL=$REF_PANEL,VAL_PANEL_DIR=$VAL_PANEL_DIR,VAL_PANEL=$VAL_PANEL,GENO_PANEL=$GENO_PANEL,GENETIC_MAP=$GENETIC_MAP,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,PHASING_PROGRAM=$PHASING_PROGRAM,IMPUTE_PROGRAM=$IMPUTE_PROGRAM,CHROM_SET=$CHROM_SET,INCLUDE_CHROM_X=$INCLUDE_CHROM_X,MAF=$MAF,MAF=$MAF,MISSING=$MISSING,SWITCH=$SWITCH,DOWN_SCALING=$DOWN_SCALING,FileExtension=$FileExtension,REGION=$REGION,CHR=$CHR,BPSTART=$BPSTART,BPEND=$BPEND,DOWN_SCALED=$DOWN_SCALED,FILTER=$FILTER,PHASE=$PHASE,IMPUTE=$IMPUTE,ACC=$ACC \
   $IMPUTEGEN_DIR/program/slurm/runStepsSlurm.sb)
   echo $JOBID1 >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/logSlurm/jobid1.txt
done


WAIT1=$(cat $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/logSlurm/jobid1.txt | awk '{print $4}' | tr "\n" ":" | sed 's/.$//')
wait

##### IMPUTATION ACCURACY
sbatch --dependency=afterok:$WAIT1 \
--job-name="${ANALYSIS_ID}"_"$OUTPUT_NAME"_accR --time=$ACC_TIME --cpus-per-task=$ACC_CPU --mem=$ACC_MEM \
--output=$OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/logSlurm/accR.txt \
--export=CONDA=$CONDA,IMPUTEGEN_DIR=$IMPUTEGEN_DIR,REF_GENOME=$REF_GENOME,VAL_PANEL_DIR=$VAL_PANEL_DIR,GENO_PANEL=$GENO_PANEL,DOWN_SCALING=$DOWN_SCALING,PHASE=$PHASE,ACC=$ACC,OUTPUT_DIR=$OUTPUT_DIR,OUTPUT_NAME=$OUTPUT_NAME,PARAMETERS=$PARAMETERS \
$IMPUTEGEN_DIR/program/slurm/run_accR_Slurm.sb
wait

