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


if [[ "$GENETIC_MAP" == "1cMperMb" ]]; then
   echo $'pos chr cM' > $PHASE/genMap_1cMperMb_"$FileExtension".txt
   echo 0 $CHR 0 >> $PHASE/genMap_1cMperMb_"$FileExtension".txt
   echo 1000000000 $CHR 1000 >> $PHASE/genMap_1cMperMb_"$FileExtension".txt
   GENETIC_MAP=$PHASE/genMap_1cMperMb_"$FileExtension".txt
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


##### Split the reference panel into chunks
cd $PHASE

$CONDA/bin/GLIMPSE_chunk \
--thread $THREADS \
--input $REF_PANEL \
--region "$CHR" \
--window-size 2000000 \
--buffer-size 20000 \
--output chunks."$CHR".txt
wait


##### Impute and phase a whole chromosome
while IFS="" read -r LINE || [ -n "$LINE" ]; do   
	printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
	IRG=$(echo $LINE | cut -d" " -f3)
	ORG=$(echo $LINE | cut -d" " -f4)
	OUT=$PHASE/imputed_"${ID}"_"$IRG".bcf

	$CONDA/bin/GLIMPSE_phase \
	--thread $THREADS \
	--input $VCF \
	--reference $REF_PANEL \
	--input-region $IRG \
	--map $GENETIC_MAP \
	--output-region $ORG \
	--output $OUT
	wait
	
	bcftools index -f -t --threads $THREADS $OUT
done < chunks."$CHR".txt
wait

##### Ligate multiple chunks together
LIST=$PHASE/list_$CHR.txt
ls $PHASE/imputed_*_$CHR:*.bcf > $LIST
$CONDA/bin/GLIMPSE_ligate \
--thread $THREADS \
--input $LIST \
--output merged_"$FileExtension".bcf
wait

bcftools index -f -t --threads $THREADS merged_"$FileExtension".bcf
wait

##### Sample haplotypes
$CONDA/bin/GLIMPSE_sample \
--thread $THREADS \
--input merged_"$FileExtension".bcf \
--solve \
--output phased_"$FileExtension".vcf
wait

bgzip -@ $THREADS phased_"$FileExtension".vcf
bcftools index -f -t --threads $THREADS phased_"$FileExtension".vcf.gz
wait

rm chunks."$CHR".txt
rm imputed_*_$CHR:*
rm merged_"$FileExtension".bcf*
rm $LIST
rm $IMPUTE/genMap_1cMperMb_"$FileExtension".txt 2>&1 > /dev/null
wait

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ PHASING LOW DENSITY PANEL CHROMOSOME $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"
echo






