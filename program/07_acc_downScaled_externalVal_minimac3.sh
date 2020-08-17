#!/bin/bash

mkdir -p $ACC/012 2>&1 > /dev/null
mkdir -p $ACC/r2_maf 2>&1 > /dev/null
mkdir -p $ACC/stats 2>&1 > /dev/null


if [[ $DOWN_SCALING != 0 ]]; then
   VAL_PANEL=$DOWN_SCALED/validation/validation_"$FileExtension".vcf.gz
   IMP_PANEL=$IMPUTE/imputed_"$FileExtension".dose.vcf.gz

else
   VAL_PANEL=$VAL_PANEL
   IMP_PANEL=$IMPUTE/imputed_"$FileExtension".dose.vcf.gz
   
fi
wait


echo
echo "#@##################################################################"
echo "#@ IMPUTATION ACCURACY: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_PREFIX"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START1=$(date +%s)

##### VALIDATION POPULATION #####
### MAP
bcftools query -f '%CHROM\t%POS\n' $VAL_PANEL > $ACC/012/mapObs_"$FileExtension".txt
wait

### SAMPLE NAMES
bcftools query -l $VAL_PANEL | sort | awk '{print $1"\t"$1}' > $ACC/012/sampleNamesObs_"$FileExtension".txt

### CONVERT 012
echo "PLINK: Converting vcf to 012 format (validation set)"
plink --vcf $VAL_PANEL \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --indiv-sort f $ACC/012/sampleNamesObs_"$FileExtension".txt \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --make-bed \
   --out $ACC/012/plinkObs_"$FileExtension" > /dev/null 2>&1 
wait

plink --bfile $ACC/012/plinkObs_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $ACC/012/obs012_"$FileExtension" > /dev/null 2>&1
wait

### Building GRM
echo "GCTA: Building GRM (validation set)"
gcta64 \
   --bfile $ACC/012/plinkObs_"$FileExtension" \
   --autosome-num $CHROM_SET \
   --chr $CHR \
   --threads $THREADS \
   --make-grm \
   --out $ACC/012/obsGRM_"$FileExtension" > /dev/null 2>&1
wait

### STATS
echo "PLINK: Calculating allele frequencies (validation set)"
plink2 --bfile $ACC/012/plinkObs_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --freq --keep-allele-order \
   --out $ACC/stats/obs_"$FileExtension" > /dev/null 2>&1 
wait

#echo "PLINK: Calculating linkage disequilibrium (validation set)"
#plink --bfile $ACC/012/plinkObs_"$FileExtension" \
#   --chr-set $CHROM_SET no-xy no-mt \
#   --chr $CHR \
#   --from-bp $BPSTART \
#   --to-bp $BPEND \
#   --threads $THREADS \
#   --set-missing-var-ids @:# \
#   --r2 --ld-window-kb 500 \
#   --ld-window 99999 \
#   --ld-window-r2 0 \
#   --out $ACC/stats/obs_"$FileExtension" > /dev/null 2>&1 
#wait


#if [[ $(wc -l $ACC/stats/obs_"$FileExtension".ld | awk '{print $1}') -gt 500000000 ]]; then
#   head -n 500000000 $ACC/stats/obs_"$FileExtension".ld > $ACC/stats/obs_"$FileExtension".ld.tmp && mv $ACC/stats/obs_"$FileExtension".ld.tmp $ACC/stats/obs_"$FileExtension".ld
#fi
#wait

echo "LDkit: Calculating linkage disequilibrium (validation set)"
java -jar $IMPUTEGEN_DIR/program/LDkit.jar \
--infile $VAL_PANEL --out $ACC/stats/LDkit_obs_"$FileExtension" \
--ws 1000 --chr $CHR --maf 0.001 --threads $THREADS --type 1 --Intermediate yes
wait
bcftools index -f -t --threads $THREADS $VAL_PANEL
wait
zcat $ACC/stats/LDkit_obs_"$FileExtension"/result.txt.gz > $ACC/stats/obs_"$FileExtension".LDkit
wait

echo "PLINK: Calculating heterozygosity (validation set)"
plink --bfile $ACC/012/plinkObs_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --het --keep-allele-order \
   --out $ACC/stats/obs_"$FileExtension" > /dev/null 2>&1 
wait

rm $ACC/012/obs012_"$FileExtension".log > /dev/null 2>&1 
rm $ACC/012/obs012_"$FileExtension".nosex > /dev/null 2>&1 
rm $ACC/stats/obs_"$FileExtension".log > /dev/null 2>&1
rm $ACC/stats/obs_"$FileExtension".nosex > /dev/null 2>&1
wait

##### IMPUTED POPULATION #####
### MAP
bcftools query -i 'ER2>0' -f '%CHROM\t%POS\tGenotyped\n' $IMPUTE/imputed_"$FileExtension".dose.vcf.gz > $ACC/posGenotyped_"$FileExtension".txt
bcftools query -e 'ER2>0' -f '%CHROM\t%POS\tImputed\n' $IMPUTE/imputed_"$FileExtension".dose.vcf.gz > $ACC/posImputed_"$FileExtension".txt
cat $ACC/posGenotyped_"$FileExtension".txt $ACC/posImputed_"$FileExtension".txt | sort -t$'\t' -k1,1 -k2,2 --numeric-sort > $ACC/012/mapPred_"$FileExtension".txt
wait

rm $ACC/posGenotyped_"$FileExtension".txt 2>&1 > /dev/null
rm $ACC/posImputed_"$FileExtension".txt 2>&1 > /dev/null
wait

### EXTRACT R2 and MAF
cat $IMPUTE/imputed_"$FileExtension".info | awk -F"\t" '{print $7"\t"$5}' | tail -n +2 > $ACC/r2_maf/r2_maf_"$FileExtension".txt
paste -d "\t" $ACC/012/mapPred_"$FileExtension".txt $ACC/r2_maf/r2_maf_"$FileExtension".txt > $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp && mv $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp $ACC/r2_maf/r2_maf_"$FileExtension".txt
echo -e "chrom\tpos\tclass\tr2\tmaf" > $ACC/r2_maf/header_"$FileExtension".txt
cat $ACC/r2_maf/header_"$FileExtension".txt $ACC/r2_maf/r2_maf_"$FileExtension".txt > $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp && mv $ACC/r2_maf/r2_maf_"$FileExtension".txt.tmp $ACC/r2_maf/r2_maf_"$FileExtension".txt
rm $ACC/r2_maf/header_"$FileExtension".txt 2>&1 > /dev/null
wait

### SAMPLE NAMES
bcftools query -l $IMP_PANEL | sort | awk '{print $1"\t"$1}' > $ACC/012/sampleNamesPred_"$FileExtension".txt

### CONVERT 012
# Imputed set
echo
echo "PLINK: Converting vcf to 012 format (imputation set)"
plink --vcf $IMP_PANEL \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --indiv-sort f $ACC/012/sampleNamesPred_"$FileExtension".txt \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --make-bed \
   --out $ACC/012/plinkPred_"$FileExtension" > /dev/null 2>&1 
wait
   
plink --bfile $ACC/012/plinkPred_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $ACC/012/pred012_"$FileExtension" > /dev/null 2>&1 
wait

### Building GRM
echo "GCTA: Building GRM (imputation set)"
gcta64 \
   --bfile $ACC/012/plinkPred_"$FileExtension" \
   --autosome-num $CHROM_SET \
   --chr $CHR \
   --threads $THREADS \
   --make-grm \
   --out $ACC/012/predGRM_"$FileExtension" > /dev/null 2>&1
wait

### STATS
echo "PLINK: Calculating allele frequencies (imputation set)"
plink2 --bfile $ACC/012/plinkPred_"$FileExtension" \
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
#plink --bfile $ACC/012/plinkPred_"$FileExtension" \
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
plink --bfile $ACC/012/plinkPred_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --het --keep-allele-order \
   --out $ACC/stats/pred_"$FileExtension" > /dev/null 2>&1 
wait

rm $ACC/012/pred012_"$FileExtension".log > /dev/null 2>&1 
rm $ACC/012/pred012_"$FileExtension".nosex > /dev/null 2>&1 
rm $ACC/stats/pred_"$FileExtension".log > /dev/null 2>&1
rm $ACC/stats/pred_"$FileExtension".nosex > /dev/null 2>&1
wait


##### PCA OF VALIDATION AND IMPUTED FILES
awk -v OFS=' ' '{print $1"_obs", $2"_obs", $3, $4, $5, $6}' $ACC/012/plinkObs_"$FileExtension".fam > $ACC/012/plinkObs_"$FileExtension".fam.tmp && mv $ACC/012/plinkObs_"$FileExtension".fam.tmp $ACC/012/plinkObs_"$FileExtension".fam
awk -v OFS=' ' '{print $1"_pred", $2"_pred", $3, $4, $5, $6}' $ACC/012/plinkPred_"$FileExtension".fam > $ACC/012/plinkPred_"$FileExtension".fam.tmp && mv $ACC/012/plinkPred_"$FileExtension".fam.tmp $ACC/012/plinkPred_"$FileExtension".fam

plink --bfile $ACC/012/plinkObs_"$FileExtension" \
   --bmerge $ACC/012/plinkPred_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --pca 8000 \
   --out $ACC/012/pca_"$FileExtension"> /dev/null 2>&1 
wait

rm $ACC/012/plinkObs_"$FileExtension".* > /dev/null 2>&1
rm $ACC/012/plinkPred_"$FileExtension".* > /dev/null 2>&1 
rm $ACC/012/pca_"$FileExtension".log > /dev/null 2>&1 
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

echo
echo "#@##################################################################"
echo "#@ PHASING IMPUTATION ACCURACY: Chromosome $REGION"
echo "#@ NAME OF ANALYSIS: $OUTPUT_PREFIX"
D1=`date "+%D    %T"`
echo "#@ Date and Time: $D1"
echo
       #################################

START2=$(date +%s)

##### GET LOW DENSITY MAP FROM PHASING STEP
plink --vcf $PHASE/phased_"$FileExtension".vcf.gz \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --make-bed \
   --out $PHASE/plink_"$FileExtension" > /dev/null 2>&1 
wait

awk '{print $2}' $PHASE/plink_"$FileExtension".bim > $PHASE/pos_"$FileExtension".txt
wait

##### IMPUTATION FILE
echo "PLINK: Converting vcf to 012 format (imputation set)"
plink --vcf $IMP_PANEL \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --extract $PHASE/pos_"$FileExtension".txt \
   --indiv-sort f $ACC/012/sampleNamesPred_"$FileExtension".txt \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --make-bed \
   --out $ACC/012/plinkPhasedPred_"$FileExtension" > /dev/null 2>&1 
wait
   
plink --bfile $ACC/012/plinkPhasedPred_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $ACC/012/phasedPred012_"$FileExtension" > /dev/null 2>&1 
wait

##### VALIDATION FILE
echo "PLINK: Converting vcf to 012 format (validation set)"
plink --vcf $VAL_PANEL \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --extract $PHASE/pos_"$FileExtension".txt \
   --indiv-sort f $ACC/012/sampleNamesObs_"$FileExtension".txt \
   --set-missing-var-ids @:# \
   --double-id \
   --keep-allele-order --make-bed \
   --out $ACC/012/plinkPhasedObs_"$FileExtension" > /dev/null 2>&1 
wait
   
plink --bfile $ACC/012/plinkPhasedObs_"$FileExtension" \
   --chr-set $CHROM_SET no-xy no-mt \
   --chr $CHR \
   --from-bp $BPSTART \
   --to-bp $BPEND \
   --threads $THREADS \
   --keep-allele-order --recode A \
   --out $ACC/012/phasedObs012_"$FileExtension" > /dev/null 2>&1 
wait

END2=$(date +%s)
DIFF2=$(( $END2 - $START2 ))

echo
echo "#@ PHASING IMPUTATION ACCURACY Chromosome $REGION TOOK $(printf '%dh:%dm:%ds\n' $(($DIFF2/3600)) $(($DIFF2%3600/60)) $(($DIFF2%60)))"
echo "#@##################################################################"
echo
echo

