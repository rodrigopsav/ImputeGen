#!/bin/bash

mkdir -p $IMPUTE 2>&1 > /dev/null
wait

if [[ $PHASING_PROGRAM == "eagle" || $PHASING_PROGRAM == "beagle" || $PHASING_PROGRAM == "shapeit4" || $PHASING_PROGRAM == "glimpse" || $PHASING_PROGRAM == "whatshap" || $PHASING_PROGRAM == "stitch" ]]; then
   VCF=$PHASE/phased_"$FileExtension".vcf.gz
   
else

   if [[ $MAF != 0 || $MISSING != 0 ]]; then
      VCF=$FILTER/filtered_"$FileExtension".vcf.gz
   else

      if [[ $DOWN_SCALING != 0 ]]; then
         VCF=$DOWN_SCALED/test/test_"$FileExtension".vcf.gz
      else
      
         if [[ -f $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz ]]; then
            VCF=$SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz
         else
            VCF=$GENO_PANEL
         fi

      fi
   fi
fi
wait

if [[ $DOWN_SCALING != 0 ]]; then
   REF_PANEL=$DOWN_SCALED/train/train_"$FileExtension".vcf.gz
fi
wait


echo
echo "#@################################################"
echo "#@ IMPUTATION : Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_NAME"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

if [[ "$BPSTART" == "0" ]]; then
   BPSTART_1=1
else
   BPSTART_1=$BPSTART
fi

cd $IMPUTE
$CONDA/bin/minimac3 \
--refHaps $REF_PANEL --rsid \
--haps $VCF \
--prefix imputed_"$FileExtension" \
--chr $CHR \
--start $BPSTART_1 \
--end $BPEND \
--cpus $THREADS \
--log
wait

bcftools index -f -t --threads $THREADS $IMPUTE/imputed_"$FileExtension".dose.vcf.gz
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ IMPUTATION CHROMOSOME $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"
echo

