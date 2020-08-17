#!/bin/bash

echo
echo
echo "#@############################################################"
echo "#@                 CHECK ImputeGen PARAMETERS"
#echo


#############################################
########## LOAD DEFAULT PARAMETERS ##########
#############################################
COUNT_ERROR=0
COUNT_WARNING=0

### REF_GENOME
if [[ -z "$REF_GENOME" || "$REF_GENOME" == "none" ]]; then
   export REF_GENOME=none
   #echo "ERROR: Missing reference genome file (REF_GENOME variable) "
   #COUNT_ERROR=$((COUNT_ERROR + 1))
   #echo 
    
else

   if [[ ! -f "$REF_GENOME" ]]; then
      echo "ERROR: Invalid reference genome file (REF_GENOME variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   else
      export REF_GENOME=$REF_GENOME
      if [[ ! -f ${REF_GENOME}.fai ]]; then
         samtools faidx $REF_GENOME
      fi
   fi
fi


### REF_PANEL_DIR
if [[ -z "$REF_PANEL_DIR" || "$REF_PANEL_DIR" == "none" ]]; then
   export REF_PANEL_DIR=none
   #echo "ERROR: Missing reference panel directory in parameter file (REF_PANEL_DIR variable) "
   #COUNT_ERROR=$((COUNT_ERROR + 1))
   #echo
else
   if [[ ! -d "$REF_PANEL_DIR" ]]; then
      echo "ERROR: Invalid reference panel directory in parameter file (REF_PANEL_DIR variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   else
      export REF_PANEL_DIR=$REF_PANEL_DIR
   fi
fi


if [[ $(ls $REF_PANEL_DIR/*vcf.gz | wc -l) != 0 ]]; then
   export EXTENSION=.vcf.gz
elif [[ $(ls $REF_PANEL_DIR/*vcf | wc -l) != 0 ]]; then
   export EXTENSION=.vcf
else
   export EXTENSION=.bcf
fi


if [[ -z "$REF_PANEL_PREFIX" || "$REF_PANEL_PREFIX" == "none" ]]; then
   export REF_PANEL_PREFIX=none
   #echo "ERROR: Missing reference panel prefix file in parameter file (REF_PANEL_PREFIX variable) "
   #COUNT_ERROR=$((COUNT_ERROR + 1))
   #echo
else
   if [[ $(ls $REF_PANEL_DIR/"$REF_PANEL_PREFIX"* | wc -l) == 0 ]]; then
      echo "ERROR: Invalid reference panel prefix file in parameter file (REF_PANEL_PREFIX variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   else
      export REF_PANEL_PREFIX=$REF_PANEL_PREFIX
   fi
fi


### VAL_PANEL_DIR
if [[ -z "$VAL_PANEL_DIR" || "$VAL_PANEL_DIR" == "none" ]]; then
   export VAL_PANEL_DIR=none
else
   if [[ ! -d "$VAL_PANEL_DIR" ]]; then
      echo "ERROR: Invalid validation panel directory in parameter file (VAL_PANEL_DIR variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   else
      export VAL_PANEL_DIR=$VAL_PANEL_DIR
   fi
fi


if [[ $VAL_PANEL_DIR != "none" ]]; then
   if [[ $(ls $VAL_PANEL_DIR/*vcf.gz | wc -l) != 0 ]]; then
      export EXTENSION_VAL=.vcf.gz
   elif [[ $(ls $VAL_PANEL_DIR/*vcf | wc -l) != 0 ]]; then
      export EXTENSION_VAL=.vcf
   else
      export EXTENSION_VAL=.bcf
   fi
fi


if [[ -z "$VAL_PANEL_PREFIX" || "$VAL_PANEL_PREFIX" == "none" ]]; then
   export VAL_PANEL_PREFIX=none
else
   if [[ $(ls $VAL_PANEL_DIR/"$VAL_PANEL_PREFIX"* | wc -l) == 0 ]]; then
      echo "ERROR: Invalid validation panel prefix file in parameter file (VAL_PANEL_PREFIX variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   else
      export VAL_PANEL_PREFIX=$VAL_PANEL_PREFIX
   fi
fi


### GENO_PANEL
if [[ -z "$GENO_PANEL" ]]; then
   echo "ERROR: Missing low density file directory in parameter file (GENO_PANEL variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
else
   if [[ ! -f "$GENO_PANEL" ]]; then
      echo "ERROR: Invalid low density file directory in parameter file (GENO_PANEL variable) "
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   else
      export GENO_PANEL=$GENO_PANEL
   fi
fi


if echo $GENO_PANEL | grep -i -q ".gz"; then
   if [[ ! -f "$GENO_PANEL.tbi" ]]; then
      echo "Indexing LD panel file. Please wait."
      bcftools index -f -t --threads $THREADS $GENO_PANEL
   fi

else
   echo "Compressing and indexing LD panel file. Please wait."
   bgzip -f -@ $THREADS $GENO_PANEL
   bcftools index -f -t --threads $THREADS $GENO_PANEL.gz
   export GENO_PANEL=$GENO_PANEL.gz
fi
wait


### OUTPUT_DIR
if [[ -z "$OUTPUT_DIR" ]]; then
   echo "WARNING: Missing output directory in parameter file (OUTPUT_DIR variable) "
   echo "The output folder will be created in $"HOME" directory"
   export OUTPUT_DIR=$HOME
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export OUTPUT_DIR=$OUTPUT_DIR
fi


### OUTPUT_NAME
if [[ -z "$OUTPUT_NAME" ]]; then
   echo "WARNING: Missing analysis name in parameter file (OUTPUT_NAME variable) "
   echo "The analysis name will be $"USER" "
   export OUTPUT_NAME=$USER
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export OUTPUT_NAME=$OUTPUT_NAME
fi


### PHASING_PROGRAM
if [[ -z "$PHASING_PROGRAM" ]]; then
   echo "ERROR: Missing phasing program (PHASING_PROGRAM variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
else
   export PHASING_PROGRAM=${PHASING_PROGRAM,,}
fi


### IMPUTE_PROGRAM
if [[ -z "$IMPUTE_PROGRAM" ]]; then
   echo "ERROR: Missing imputation program (IMPUTE_PROGRAM variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
else
   export IMPUTE_PROGRAM=${IMPUTE_PROGRAM,,}
fi


### GENETIC_MAP
if [[ -z "$GENETIC_MAP" || "$GENETIC_MAP" == "1cMperMb" ]]; then

   if [[ "PHASING_PROGRAM" == "eagle" ]]; then
      export GENETIC_MAP=$IMPUTEGEN_DIR/program/genetic_map_1cMperMb.txt
   else
      export GENETIC_MAP="1cMperMb"
   fi
else
   export GENETIC_MAP=$(readlink -f $GENETIC_MAP)
fi


### CHROM_SET
if [[ -z "$CHROM_SET" ]]; then
   echo "ERROR: Missing chromosome set parameter file (CHROM_SET variable) "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
else
   export CHROM_SET=$CHROM_SET
fi


### INCLUDE_CHROM_X
if [[ -z "$INCLUDE_CHROM_X" ]]; then 
   export INCLUDE_CHROM_X=no
fi

if [[ "${INCLUDE_CHROM_X,,}" != "no" && "${INCLUDE_CHROM_X,,}" != "yes" ]]; then
   echo "#@ WARNING: Invalid option in parameter file (INCLUDE_CHROM_X variable) "
   echo "#@ IVDP will consider INCLUDE_CHROM_X=no"
   export INCLUDE_CHROM_X=no
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export INCLUDE_CHROM_X=${INCLUDE_CHROM_X,,}
fi


### SELECT_CHROM
# Defaults
if [[ -z "$SELECT_CHROM" ]]; then 
   export SELECT_CHROM=$(seq 1 1 $CHROM_SET)
else
   if [[ "$SELECT_CHROM" == "all" ]]; then
      export SELECT_CHROM=$(seq 1 1 $CHROM_SET)
   else
      export SELECT_CHROM=$SELECT_CHROM
   fi
fi


if [[ "$(echo $SELECT_CHROM | tr " " "\n" | wc -l)" -gt "$CHROM_SET" ]]; then
   echo "#@ ERROR: Number of selected autosomes (SELECT_CHROM variable) is greater than CHROM_SET."
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi


for chr in $SELECT_CHROM; do 
   if [[ "$chr" -gt "$CHROM_SET" ]]; then 
      echo "#@ ERROR: Wrong SELECT_CHROM value ($chr). Some values are greater than CHROM_SET variable."
      COUNT_ERROR=$((COUNT_ERROR + 1))
      echo
   fi
done
wait


### MAF
if [[ -z "$MAF" ]]; then
   echo "WARNING: Missing variable in parameter file (MAF variable) "
   echo "ImputeGen will skip maf for filtering"
   export MAF=0
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export MAF=$MAF
fi
 

### MISSING
if [[ -z "$MISSING" ]]; then
   echo "WARNING: Missing variable in parameter file (MISSING variable) "
   echo "ImputeGen will skip missingness snps for filtering"
   export MISSING=0
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export MISSING=$MISSING
fi


### DOWN_SCALING
if [[ -z "$DOWN_SCALING" ]]; then
   export DOWN_SCALING=0
else
   export DOWN_SCALING=$DOWN_SCALING
fi


### NUMBER OF THREADS
#https://stackoverflow.com/questions/394230/how-to-detect-the-os-from-a-bash-script
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
   if [[ -z "$THREADS" ]]; then
      if [[ $(nproc) -le 8 ]]; then
         export THREADS=$(nproc)
      elif [[ $(nproc) -le 12 ]]; then
         export THREADS=$(nproc)
      elif [[ $(nproc) -le 16 ]]; then
         export THREADS=$(nproc)
      else
         export THREADS=16
      fi
      
   else
      export THREADS=$(echo $THREADS | tr -d [:alpha:])
   fi

elif [[ "$OSTYPE" == "darwin"* ]]; then
   if [[ -z "$THREADS" ]]; then
      if [[ $(sysctl -n hw.ncpu) -le 8 ]]; then
         export THREADS=$(sysctl -n hw.ncpu)
      elif [[ $(sysctl -n hw.ncpu) -le 12 ]]; then
         export THREADS=$(sysctl -n hw.ncpu)
      elif [[ $(sysctl -n hw.ncpu) -le 16 ]]; then
         export THREADS=$(sysctl -n hw.ncpu)
      else
         export THREADS=16
      fi
      
   else
      export THREADS=$(echo $THREADS | tr -d [:alpha:])
   fi

else
   echo "ERROR: This is not a linux or mac machine. The analysis will stop "
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo

fi


### MEM
if [[ -z "$MEM" ]]; then
   echo "WARNING: Missing variable in parameter file (MEM variable) "
   echo "ImputeGen will consider MEM=16"
   export MEM=16
   COUNT_WARNING=$((COUNT_WARNING + 1))
   echo
else
   export MEM=$MEM
fi


### NUMBER OF BATCH SAMPLES
if [[ -z "$BATCH" ]]; then
   export BATCH=5
else
   export BATCH=$(echo $BATCH | tr -d [:alpha:])
fi


### CHECK PROGRAMS

###
if ! command -v vcftools &> /dev/null; then
   echo "ERROR: vcftools could not be found "
   echo "Please install vcftools using install_imputegen_dependencies.sh"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi
wait


if ! command -v bcftools &> /dev/null; then
   echo "ERROR: bcftools could not be found ${deafault}"
   echo "Please install bcftools using install_imputegen_dependencies.sh"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi
wait


##
echo $CONDA
if ! command -v $CONDA/bin/eagle &> /dev/null; then
   echo "ERROR: eagle could not be found "
   echo "Please install eagle using install_imputegen_dependencies.sh"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi
wait


##
if ! command -v java -jar $CONDA/bin/beagle &> /dev/null; then
   echo "ERROR: beagle could not be found "
   echo "Please install beagle using install_imputegen_dependencies.sh"
   COUNT_ERROR=$((COUNT_ERROR + 1))
   echo
fi
wait


###############################################
########## CHECK ERRORS AND WARNINGS ##########
###############################################

echo "***** TOTAL OF ERRORS: $COUNT_ERROR *****"
echo "***** TOTAL OF WARNINGS: $COUNT_WARNING *****"

if [[ "$COUNT_ERROR" -gt "0" ]]; then
   echo
   echo "Solve all the ERRORS before run ImputeGen again "
   echo "Don't worry about WARNINGS: ImputeGen will use default parameters "
   echo "Aborting analysis" 
   exit 1
else
   if [[ "$COUNT_WARNING" -gt "0" ]]; then
      echo
      echo "Don't worry about WARNINGS: ImputeGen will use default parameters "
   fi
fi

#echo
#echo "#@#############################################################"
