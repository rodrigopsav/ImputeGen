#!/bin/bash

echo
echo "#@#################################################################"
echo "#@                         INPUT PARAMETERS"
echo "#@   High Density Panek: $REF_PANEL"
echo "#@   Low Density Panel: $GENO_PANEL"
echo "#@   Output Folder: $OUTPUT_NAME"
echo "#@"
echo "#@   Phasing Program: $PHASING_PROGRAM"
echo "#@   Imputation Program: $IMPUTE_PROGRAM"
echo "#@"
echo "#@   Chromosome: $REGION"
echo "#@   Filtering: MAF=$MAF  |  MISSING=$MISSING"
echo "#@"
D1=$(date "+%D    %T")
echo "#@   Date and Time: $D1"
echo "#@#################################################################"
echo