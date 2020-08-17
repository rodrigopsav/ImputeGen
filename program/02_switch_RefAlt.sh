#!/bin/bash

mkdir -p $SWITCH 2>&1 > /dev/null
wait

VCF=$GENO_PANEL
wait


echo
echo "#@################################################"
echo "#@ SWITCH REF/ALT LOW DENSITY PANEL: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_NAME"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

# Keep only biallelic snps
bcftools view --threads $THREADS -m2 -M2 -v snps -Oz -o $SWITCH/biallelic_"$FileExtension".vcf.gz $VCF -r "$CHR":"$BPSTART"-"$BPEND"
bcftools index -f -t --threads $THREADS $SWITCH/biallelic_"$FileExtension".vcf.gz
wait

# Match / Mismatch stats
bcftools +fixref --threads $THREADS $SWITCH/biallelic_"$FileExtension".vcf.gz -- \
-f $REF_GENOME  > $SWITCH/mismatchStats_"$FileExtension".txt 2>&1
wait
mismatch=$(echo $(grep "ref mismatch" $SWITCH/mismatchStats_"$FileExtension".txt | cut -f 4- | awk -F"%" '{printf "%.0f\n", $1}') | bc -l)
wait

# Switch Ref/Alt if necessary
if [[ $mismatch -gt 3 ]]; then
   bcftools +fixref --threads $THREADS $SWITCH/biallelic_"$FileExtension".vcf.gz \
   -Ob -o $SWITCH/biallelic_"$FileExtension"_fixref.bcf -- \
   -f $REF_GENOME \
   -m flip -d
   wait
   
   bcftools sort -Oz -o $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz \
   $SWITCH/biallelic_"$FileExtension"_fixref.bcf
   wait
   
   bcftools index -f -t --threads $THREADS $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz
   wait
   rm $SWITCH/biallelic_"$FileExtension".vcf.gz 2>&1 > /dev/null
   wait
   rm $SWITCH/biallelic_"$FileExtension"_fixref.bcf 2>&1 > /dev/null
   wait

else

   mv $SWITCH/biallelic_"$FileExtension".vcf.gz $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz  2>&1 > /dev/null
   mv $SWITCH/biallelic_"$FileExtension".vcf.gz.tbi $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz.tbi  2>&1 > /dev/null
   wait
   
fi
wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ SWITCH REF/ALT LOW DENSITY PANEL CHROMOSOME $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"
echo

