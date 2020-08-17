#!/bin/bash

if [[ -f $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz ]]; then
   VCF=$SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz
else
   VCF=GENO_PANEL
fi
wait


echo
echo "#@##################################################################"
echo "#@ DOWN-SCALING FROM HIGH TO LOW DENSITY: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_NAME"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)


##### CREATE FOLDERS
mkdir -p $DOWN_SCALED/subset 2>&1 > /dev/null
mkdir -p $DOWN_SCALED/train 2>&1 > /dev/null
mkdir -p $DOWN_SCALED/test 2>&1 > /dev/null
mkdir -p $DOWN_SCALED/validation 2>&1 > /dev/null


##### SELECT TRAIN, AND VALIDATION SETS PER CHROMOSOME
# SELECT SAMPLES FROM VCF samples from vcf
# https://www.biostars.org/p/355360/
# RANDOM SAMPLING
# https://unix.stackexchange.com/questions/108581/how-to-randomly-sample-a-subset-of-a-file

# GET ALL SAMPLE NAMES FROM REFERENCE PANEL
bcftools query -l $REF_PANEL > $DOWN_SCALED/subset/samplesAll_"$FileExtension".txt

# SUBSETTING SAMPLE NAMES FOR TEST AND VALIDATION
N=$(cat $DOWN_SCALED/subset/samplesAll_"$FileExtension".txt | wc -l)
SUBSET=$( printf %.0f $(echo "scale=3; $DOWN_SCALING*$N" | bc))
shuf -n $SUBSET $DOWN_SCALED/subset/samplesAll_"$FileExtension".txt > $DOWN_SCALED/subset/samplesTest_"$FileExtension".txt

# SUBSETTING SAMPLE NAMES FOR TRAIN
comm -23 <(sort $DOWN_SCALED/subset/samplesAll_"$FileExtension".txt) <(sort DOWN_SCALED/subset/samplesTest_"$FileExtension".txt) > $DOWN_SCALED/subset/samplesTrain_"$FileExtension".txt


##### TRAINING FILE
bcftools view --threads $THREADS -S $DOWN_SCALED/subset/samplesTrain_"$FileExtension".txt \
-Oz -o $DOWN_SCALED/train/train_"$FileExtension".vcf.gz \
$REF_PANEL
wait
bcftools index -f -t --threads $THREADS $DOWN_SCALED/train/train_"$FileExtension".vcf.gz
wait

##### VALIDATION FILE
bcftools view --threads $THREADS -S $DOWN_SCALED/subset/samplesTest_"$FileExtension".txt \
-Oz -o $DOWN_SCALED/validation/validation_"$FileExtension".vcf.gz \
$REF_PANEL
wait
bcftools index -f -t --threads $THREADS $DOWN_SCALED/validation/validation_"$FileExtension".vcf.gz
wait

##### TEST FILE
# GET POSITIONS LOW DENSITY PANEL
bcftools query -f '%CHROM %POS %POS\n' $VCF | awk -v OFS='\t' '{$1=$1};1' > $DOWN_SCALED/subset/lociLowPanel_"$FileExtension".txt

# SUBSETTING POSITIONS FROM VALIDATION SET
bcftools view --threads $THREADS -R $DOWN_SCALED/subset/lociLowPanel_"$FileExtension".txt \
-Oz -o $DOWN_SCALED/test/test_"$FileExtension".vcf.gz \
$DOWN_SCALED/validation/validation_"$FileExtension".vcf.gz
wait
bcftools index -f -t --threads $THREADS $DOWN_SCALED/test/test_"$FileExtension".vcf.gz
wait

# or using bedtools
#bedtools intersect -header -wa -a $DOWN_SCALED/validation/validation_"$FileExtension".vcf \
#-b $DOWN_SCALED/subset/lociLowPanel_"$FileExtension".txt > DOWN_SCALED/test/test_"$FileExtension".vcf
#bgzip -@ $THREADS DOWN_SCALED/test/test_"$FileExtension".vcf
#bcftools index -f -t --threads $THREADS $DOWN_SCALED/test/test_"$FileExtension".vcf.gz
#wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ DOWN-SCALED FROM HIGH TO LOW DENSITY CHROMOSOME $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"

echo
echo

