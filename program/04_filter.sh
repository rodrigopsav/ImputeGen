#!/bin/bash

mkdir -p $FILTER 2>&1 > /dev/null
wait

if [[ $DOWN_SCALING != 0 ]]; then
   VCF=$DOWN_SCALED/test/test_"$FileExtension".vcf.gz
else
   if [[ -f $SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz ]]; then
      VCF=$SWITCH/biallelic_"$FileExtension"_fixref_sorted.vcf.gz
   else
      VCF=$GENO_PANEL
   fi
fi
wait


echo
echo "#@################################################"
echo "#@ FILTERING LOW DENSITY PANEL: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_NAME"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

plink2 \
--vcf $VCF \
--chr-set $CHROM_SET no-xy no-mt \
--chr $CHR \
--from-bp $BPSTART \
--to-bp $BPEND \
--snps-only \
--min-alleles 2 \
--max-alleles 2 \
--set-missing-var-ids @:# \
--maf $MAF \
--geno $MISSING \
--keep-allele-order \
--const-fid 0 \
--recode vcf id-paste=iid \
--out $FILTER/filtered_"$FileExtension"
wait


# Grep header of former vcf file to avoid problems with bcftools isec and vcf-compare snf beagle
zcat $VCF | grep "^#" > $FILTER/filtered_"$FileExtension".header
grep -v "^#" $FILTER/filtered_"$FileExtension".vcf > $FILTER/filtered_"$FileExtension".tmp
cat $FILTER/filtered_"$FileExtension".header $FILTER/filtered_"$FileExtension".tmp > $FILTER/filtered_"$FileExtension".vcf
wait

rm -r $FILTER/filtered_"$FileExtension".log 2> /dev/null
rm -r $FILTER/filtered_"$FileExtension".header 2> /dev/null
rm -r $FILTER/filtered_"$FileExtension".tmp 2> /dev/null
rm -r $FILTER/filtered_"$FileExtension".nosex 2> /dev/null
wait

bgzip -f -@ $THREADS $FILTER/filtered_"$FileExtension".vcf
bcftools index -f -t --threads $THREADS $FILTER/filtered_"$FileExtension".vcf.gz

END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ FILTERING LOW DENSITY PANEL CHROMOSOME $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"
echo


##########################################################
###################### Using vcftools #####################
##########################################################

#vcftools --gzvcf $VCF \
#--remove-indels \
#--min-alleles 2 --max-alleles 2 \
#--chr $CHR \
#--from-bp $BPSTART \
#--to-bp $BPEND \
#--maf $MAF \
#--max-missing $MISSING \
#--recode --recode-INFO-all \
#--out $FILTER/filtered_"$FileExtension"
#wait
##minQ 30 \

#mv $FILTER/filtered_"$FileExtension".recode.vcf \
#   $FILTER/filtered_"$FileExtension".vcf

#bgzip -f -@ $THREADS $FILTER/filtered_"$FileExtension".vcf
#bcftools index -f -t --threads $THREADS $FILTER/filtered_"$FileExtension".vcf.gz
#wait

#########################################################
###################### Using Plink1 #####################
#########################################################

#plink \
#--vcf $VCF \
#--chr-set $CHROM_SET no-xy no-mt \
#--chr $CHR \
#--from-bp $BPSTART \
#--to-bp $BPEND \
#--snps-only \
#--biallelic-only strict \
#--set-missing-var-ids @:# \
#--geno $MISSING \
#--maf $MAF \
#--keep-allele-order \
#--id-delim " " \
#--const-fid 0 \
#--recode vcf-iid \
#--out $FILTER/filtered_"$FileExtension"
#wait
##--vcf-min-dp 3 \
##--mind $MISSING \
##--vcf-min-qual 30 \

#rm $FILTER/*nosex 2>&1 > /dev/null
#bgzip -f -@ $THREADS $FILTER/filtered_"$FileExtension".vcf
#bcftools index -f -t --threads $THREADS $FILTER/filtered_"$FileExtension".vcf.gz


#########################################################
###################### Using Plink2 #####################
# I had problem with ploid with Eagle. I don't know why #
#########################################################

### Convert vcf to plink format
#plink2 \
#--vcf $VCF \
#--chr-set $CHROM_SET no-xy no-mt \
#--chr $CHR \
#--from-bp $BPSTART \
#--to-bp $BPEND \
#--snps-only \
#--min-alleles 2 \
#--max-alleles 2 \
#--keep-allele-order --make-bed \
#--out $FILTER/plink_"$FileExtension"
#wait
##--vcf-min-dp 3 \
##--vcf-min-qual 30 \


### Filtering
#plink2 \
#--bfile $FILTER/plink_"$FileExtension" \
#--chr-set $CHROM_SET no-xy no-mt \
#--chr $CHR \
#--from-bp $BPSTART \
#--to-bp $BPEND \
#--set-missing-var-ids @:# \
#--geno $MISSING \
#--maf $MAF
#--keep-allele-order --make-bed \
#--out $FILTER/plinkFiltered_"$FileExtension"
#wait
##--mind $MISSING \


# Convert plink format to vcf
#plink2 \
#--bfile $FILTER/plinkFiltered_"$FileExtension" \
#--chr-set $CHROM_SET no-xy no-mt \
#--chr $CHR \
#--from-bp $BPSTART \
#--to-bp $BPEND \
#--keep-allele-order \
#--id-delim " " \
#--const-fid 0 \
#--recode vcf \
#--out $FILTER/filtered_"$FileExtension"
#wait

#rm -r $FILTER/plink* 2>&1 > /dev/null
#rm -r $FILTER/filtered*.log 2>&1 > /dev/null
#bgzip -f -@ $THREADS $FILTER/filtered_"$FileExtension".vcf
#bcftools index -f -t --threads $THREADS $FILTER/filtered_"$FileExtension".vcf.gz
#wait



