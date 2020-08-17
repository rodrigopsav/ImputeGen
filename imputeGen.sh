#!/bin/bash


##### SET PARAMETERS imputeGen_run.sh #####
usage() { echo "Usage: $0 -p params.txt -c configSlurm.txt
       -p: parameter file
       -c: config slurm parameter file (required only for the HPCC with slurm)
       
Parameter file variables

REF_GENOME         # File directory for reference genome (*.fa). The reference genome must be uncompressed file and not (*.fa.gz).
REF_PANEL_DIR      # Directory path with high density panel files (Reference panel).
REF_PANEL_PREFIX   # Prefix name of reference panel files (reference files must be divided by chromosome and named as: (prefix name + chromosome name + .vcf) or (prefix name + chromosome name + vcf.gz).
VAL_PANEL_DIR      # Directory path with validaiton panel files (Validation panel).
VAL_PANEL_PREFIX   # Prefix name of validaiton panel files (validation files must be divided by chromosome and named as: (prefix name + chromosome name + .vcf) or (prefix name + chromosome name + vcf.gz).
GENO_PANEL         # File path for the low density panel to be imputed.
BAM_DIR            # Directory path with the bam files (for whatshap and stitch)
PEDIGREE           # File path with the pedigree (for phasing with whatshap - when the pedigree is available)
FINDHAPF90         # Directory path with the files to be used with findhap90 (not implemented yet)
GENETIC_MAP        # File path for the recombination genetic map (leave it empty to consider 1cMperMb)
OUTPUT_DIR         # Directory path for the outputs (Default is the home directory if left blank: OUTPUT_DIR= )
OUTPUT_NAME        # Analysis name. Choose a single name without spaces (Default is user name if left blank: OUTPUT_NAME= )

PHASING_PROGRAM    # Choose the software for phasing (Options: eagle, beagle, shapeit4, glimpse)
IMPUTE_PROGRAM     # Choose the software for imputation (leave it empty if you won't use minimac3)

CHROM_SET         # Number of autosome pairs (Eg: CHROM_SET=22 for humans, CHROM_SET=29 for cow, CHROM_SET=18 for pig)
INCLUDE_CHROM_X   # Do the variants call for chromosome X? INCLUDE_CHROM_X=no or INCLUDE_CHROM_X=yes (Default: INCLUDE_CHROM_X=no)
SELECT_CHROM      # Comma-separated list of chromosomes that will be used for variant calling (Eg. SELECT_CHROM=1,2,3). If omitted, IVDP will consider all the autosomes. 
MAF               # Mimimum allele frequency (Use MAF=0 or leave it empty if vcf was already filtered for maf).
MISSING           # Threshold for SNP missingness (Use MISSING=0 or leave it empty if low density panel was already filtered for call rate). Eg. MISSING=0.1 exclude positions with more than 10% of missing calls.

DOWN_SCALING    # Down-scaling HD panel to LD panel level to test imputation method. Choose the proportion of of samples to be used in the test set. Eg. DOWN_SCALING=0.3 will split the reference panel in 70% of samples for reference population and 30% for imputation.
THREADS         # Number of threads for analysis
MEM             # Amount of memory for analysis
BATCH           # Number of analysis at the same time, i.e. number of imputations to be done at the same time.

       " 1>&2; exit 1; }

while getopts :p:c: option; do
   case "${option}" in
   p) PARAMETERS=${OPTARG};;
   c) CONFIGSLURM=${OPTARG};;
   *) usage;;
   esac
done
shift $((OPTIND -1))


##### Run imputeGenLocal or imputeGenSlurm #####
export IMPUTEGEN_DIR=$(dirname $(readlink -f $0))


#
if [[ -z "$PARAMETERS" ]]; then
   echo "ERROR: -p flag is empty "
   echo "Aborting analysis"
   usage
   exit 1

else
   export PARAMETERS=$(readlink -f $PARAMETERS)
   
   if [[ ! -f "$PARAMETERS" ]]; then
      echo "ERROR: parameter file does not exist "
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait


#
if [[ "$CONFIGSLURM" ]]; then
   export CONFIGSLURM=$(readlink -f $CONFIGSLURM)

   if [[ ! -f "$CONFIGSLURM" ]]; then
      echo "ERROR: slurm parameter file does not exist"
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait


#
if [[ ! -f "$CONFIGSLURM" ]]; then
   source $IMPUTEGEN_DIR/program/imputeGenLocalRun.sh
else
   for par in $(cat $CONFIGSLURM); do export $par; done
   source $IMPUTEGEN_DIR/program/imputeGenSlurmRun.sh
fi
