#!/bin/bash

mkdir -p $PHASE 2>&1 > /dev/null
wait

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
wait


if [[ $DOWN_SCALING != 0 ]]; then
   REF_PANEL=$DOWN_SCALED/train/train_"$FileExtension".vcf.gz
fi
wait


echo
echo "#@################################################"
echo "#@ PHASING LOW DENSITY PANEL: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_NAME"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

echo
echo "##### PHASING FILTERED LD PANEL"
echo


### Get bam files
#https://science.sciencemag.org/content/sci/suppl/2020/07/15/369.6501.eaba4674.DC1/aba4674-Fuller-SM.pdf
samtools view -bs ${RNUM}.${DOWNSAMPLE} ${IN}.marked_duplicates.bam "$CHR" \
| samtools mpileup - > ${IN}.downsampled_${COVERAGE}x.pileup

samtools view -bs ${RNUM}.${DOWNSAMPLE} ${IN}.marked_duplicates.bam "$CHR" \
| samtools mpileup -f ${REF} -v - > ${IN}.downsampled_${COVERAGE}x.pileup.vcf

gzip ${IN}.downsampled_${COVERAGE}x.pileup

REF_PANEL="chr1.HaplotypeData.vcf.gz"

loimpute \
-i ${IN}.downsampled_${COVERAGE}x.pileup.gz \
-h ${REF_PANEL} \
-o ${IN}.${COVERAGE}_imputed \
-k 80
-id ${IN}

gunzip ${IN}.${COVERAGE}_imputed.vcf.gz

bcftools call -c ${IN}.downsampled_${COVERAGE}x.pileup.vcf > ${IN}.${COVERAGE}x_sites.vcf


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ FILTER AND PHASING LOW DENSITY PANEL CHROMOSOME $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"

echo
echo





