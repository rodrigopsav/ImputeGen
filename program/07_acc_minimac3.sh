#!/bin/bash

#mkdir -p $ACC/genotypes 2>&1 > /dev/null
mkdir -p $ACC/012 2>&1 > /dev/null
mkdir -p $ACC/r2_maf 2>&1 > /dev/null
mkdir -p $ACC/stats 2>&1 > /dev/null

IMP_PANEL=$IMPUTE/imputed_"$FileExtension".dose.vcf.gz


echo
echo "#@##################################################################"
echo "#@ IMPUTATION ACCURACY: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_PREFIX"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

##### IMPUTED FILE
### MAP
bcftools query -i 'ER2>0' -f '%CHROM\t%POS\tGenotyped\n' $IMPUTE/imputed_"$FileExtension".dose.vcf.gz > $ACC/posGenotyped_"$FileExtension".txt
bcftools query -e 'ER2>0' -f '%CHROM\t%POS\tImputed\n' $IMPUTE/imputed_"$FileExtension".dose.vcf.gz > $ACC/posImputed_"$FileExtension".txt
cat $ACC/posGenotyped_"$FileExtension".txt $ACC/posImputed_"$FileExtension".txt | sort -t$'\t' -k1,1 -k2,2 --numeric-sort  > $ACC/012/map_"$FileExtension".txt
wait

rm $ACC/posGenotyped_"$FileExtension".txt 2>&1 > /dev/null
rm $ACC/posImputed_"$FileExtension".txt 2>&1 > /dev/null
wait

### EXTRACT R2 and MAF
cat $IMPUTE/imputed_"$FileExtension".info | awk -F"\t" '{print $7"\t"$5}' | tail -n +2 > $ACC/r2_maf/r2_maf_"$FileExtension".txt
paste -d "\t" $ACC/012/map_"$FileExtension".txt $ACC/r2_maf/r2_maf_"$FileExtension".txt > $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp && mv $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp $ACC/r2_maf/r2_maf_"$FileExtension".txt
echo -e "chrom\tpos\tclass\tr2\tmaf" > $ACC/r2_maf/header_"$FileExtension".txt
cat $ACC/r2_maf/header_"$FileExtension".txt $ACC/r2_maf/r2_maf_"$FileExtension".txt > $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp && mv $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp $ACC/r2_maf/r2_maf_"$FileExtension".txt
rm $ACC/r2_maf/header_"$FileExtension".txt 2>&1 > /dev/null
wait

### SAMPLE NAMES
bcftools query -l $IMP_PANEL | sort | awk '{print $1"\t"$1}' > $ACC/012/sampleNames_"$FileExtension".txt

### CONVERT 012
# Imputed set
echo "PLINK: Converting vcf to 012 format (imputation set)"
plink --vcf $IMP_PANEL \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --indiv-sort f $ACC/012/sampleNames_"$FileExtension".txt \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --make-bed \
   --out $ACC/012/plink_"$FileExtension" > /dev/null 2>&1 
wait
   
plink --bfile $ACC/012/plink_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $ACC/012/geno012_"$FileExtension" > /dev/null 2>&1
wait

### PCA
echo "GCTA: PCA (validation set)"
gcta64 \
   --grm-bin $ACC/012/plink_"$FileExtension" \
   --pca 2 \
   --threads $THREADS \
   --out $ACC/012/genoGRM_"$FileExtension" > /dev/null 2>&1
wait

### STATS
echo "PLINK: Calculating allele frequencies (imputation set)"
plink2 --bfile $ACC/012/plink_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --freq --keep-allele-order \
   --out $ACC/stats/pred_"$FileExtension" > /dev/null 2>&1 
wait

#echo "PLINK: Calculating linkage disequilibrium (imputation set)"
#plink --bfile $ACC/012/plink_"$FileExtension" \
#   --chr-set $CHROM_SET no-xy no-mt \
#   --chr $CHR \
#   --from-bp $BPSTART \
#   --to-bp $BPEND \
#   --threads $THREADS \
#   --set-missing-var-ids @:# \
#   --r2 --ld-window-kb 500 \
#   --ld-window 99999 \
#   --ld-window-r2 0 \
#   --out $ACC/stats/pred_"$FileExtension" > /dev/null 2>&1 
#wait


#if [[ $(wc -l $ACC/stats/pred_"$FileExtension".ld | awk '{print $1}') -gt 500000000 ]]; then
#   head -n 500000000 $ACC/stats/pred_"$FileExtension".ld > $ACC/stats/pred_"$FileExtension".ld.tmp && mv $ACC/stats/pred_"$FileExtension".ld.tmp $ACC/stats/pred_"$FileExtension".ld
#fi
#wait

echo "LDkit: Calculating linkage disequilibrium (imputation set)"
java -jar $IMPUTEGEN_DIR/program/LDkit.jar \
--infile $IMP_PANEL --out $ACC/stats/LDkit_pred_"$FileExtension" \
--ws 1000 --chr $CHR --maf 0.001 --threads $THREADS --type 1 --Intermediate yes
wait
bcftools index -f -t --threads $THREADS $IMP_PANEL
wait
zcat $ACC/stats/LDkit_pred_"$FileExtension"/result.txt.gz > $ACC/stats/pred_"$FileExtension".LDkit
wait

echo "PLINK: Calculating heterozygosity (imputation set)"
plink --bfile $ACC/012/plink_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --het --keep-allele-order \
   --out $ACC/stats/pred_"$FileExtension" > /dev/null 2>&1 
wait

rm $ACC/012/plink_"$FileExtension".* > /dev/null 2>&1 
rm $ACC/012/geno012_"$FileExtension".log > /dev/null 2>&1 
rm $ACC/012/geno012_"$FileExtension".nosex > /dev/null 2>&1 
rm $ACC/012/genoGRM_"$FileExtension".log > /dev/null 2>&1 
rm $ACC/stats/pred_"$FileExtension".log > /dev/null 2>&1
rm $ACC/stats/pred_"$FileExtension".nosex > /dev/null 2>&1
wait

##### REFERENCE POPULATION #####

# ### STATS
# echo
# echo "PLINK: Calculating allele frequencies (reference set)"
# plink2 --vcf $REF_PANEL \
#    --chr-set $CHROM_SET no-xy no-mt \
#    --chr $CHR \
#    --from-bp $BPSTART \
#    --to-bp $BPEND \
#    --threads $THREADS \
#    --set-missing-var-ids @:# \
#    --double-id \
#    --freq --keep-allele-order \
#    --out $ACC/stats/ref_"$FileExtension" > /dev/null 2>&1 
# wait
# 
# #echo "PLINK: Calculating linkage disequilibrium (reference set)"
# #plink --vcf $REF_PANEL \
# #   --chr-set $CHROM_SET no-xy no-mt \
# #   --chr $CHR \
# #   --from-bp $BPSTART \
# #   --to-bp $BPEND \
# #   --threads $THREADS \
# #   --set-missing-var-ids @:# \
# #   --double-id \
# #   --r2 --ld-window-kb 500 \
# #   --ld-window 99999 \
# #   --ld-window-r2 0 \
# #   --out $ACC/stats/ref_"$FileExtension" > /dev/null 2>&1 
# #wait
# 
# 
# #if [[ $(wc -l $ACC/stats/ref_"$FileExtension".ld | awk '{print $1}') -gt 500000000 ]]; then
# #   head -n 500000000 $ACC/stats/ref_"$FileExtension".ld > $ACC/stats/ref_"$FileExtension".ld.tmp && mv $ACC/stats/ref_"$FileExtension".ld.tmp $ACC/stats/ref_"$FileExtension".ld
# #fi
# #wait
# 
# echo "LDkit: Calculating linkage disequilibrium (reference set)"
# java -jar $IMPUTEGEN_DIR/program/LDkit.jar \
# --infile $REF_PANEL --out $ACC/stats/LDkit_ref_"$FileExtension" \
# --ws 1000 --chr $CHR --maf 0.001 --threads $THREADS --type 1 --Intermediate yes
# wait
# bcftools index -f -t --threads $THREADS $REF_PANEL
# wait
# zcat $ACC/stats/LDkit_ref_"$FileExtension"/result.txt.gz > $ACC/stats/ref_"$FileExtension".LDkit
# wait
# 
# echo "PLINK: Calculating heterozygosity (reference set)"
# plink --vcf $REF_PANEL \
#    --chr-set $CHROM_SET no-xy no-mt \
#    --chr $CHR \
#    --from-bp $BPSTART \
#    --to-bp $BPEND \
#    --threads $THREADS \
#    --set-missing-var-ids @:# \
#    --double-id \
#    --het --keep-allele-order \
#    --out $ACC/stats/ref_"$FileExtension" > /dev/null 2>&1 
# wait
# 
# 
rm -r $ACC/stats/*.log 2> /dev/null
rm -r $ACC/stats/*.nosex 2> /dev/null


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

echo
echo "#@ IMPUTATION ACCURACY Chromosome $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@##################################################################"

echo
echo
