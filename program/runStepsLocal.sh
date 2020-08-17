#!/bin/bash

(
echo "#@#############################################################"
echo "#@                    ImputeGen: $OUTPUT_NAME"
echo
D1=$(date "+%D    %T")
echo "#@ Date and Time: $D1"
echo
) > $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/log_"$OUTPUT_NAME".txt

START1=$(date +%s)


for REGION in $SELECT_CHROM; do
(
   export REGION=$REGION
   export CHR=$(echo $REGION | awk -F":" '{print $1}')
   export BPSTART=$(echo $REGION | awk -F":" '{print $2}' | awk -F"-" '{print $1}')
   export BPEND=$(echo $REGION | awk -F":" '{print $2}' | awk -F"-" '{print $2}')
   
   export REF_PANEL=$REF_PANEL_DIR/${REF_PANEL_PREFIX}${CHR}${EXTENSION}
   if [[ $VAL_PANEL_DIR != "none" ]]; then export VAL_PANEL=$VAL_PANEL_DIR/${VAL_PANEL_PREFIX}${CHR}${EXTENSION_VAL}; fi
   
   if [[ -z "$BPSTART" ]]; then export BPSTART=0; fi
   if [[ -z "$BPEND" ]]; then export BPEND=1000000000; fi
   
   if [[ "$BPSTART" == 0 && "$BPEND" == 1000000000 ]]; then
      export FileExtension=$CHR
   else
      export FileExtension=$CHR_"$BPSTART"_"$BPEND"
   fi
   
   echo "Running chromosome $FileExtension" >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/log_"$OUTPUT_NAME".txt
   wait
   
   
   ##################
   ##### Header #####
   ##################
   source $IMPUTEGEN_DIR/program/01_header.sh >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
   ################################
   ##### Check Switch Ref/Alt #####
   ################################
   (
   if [[ $REF_GENOME != "none" ]]; then
      source $IMPUTEGEN_DIR/program/02_switch_RefAlt.sh
   fi
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
   #######################
   ##### Downscaling #####
   #######################
   (
   if [[ $DOWN_SCALING != 0 ]]; then
      source $IMPUTEGEN_DIR/program/03_downScaling.sh
   fi
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
   ##################
   ##### Filter #####
   ##################
   (
   if [[ $MAF != 0 || $MISSING != 0 ]]; then
      source $IMPUTEGEN_DIR/program/04_filter.sh
   fi
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
   ###################
   ##### Phasing #####
   ###################
   (
   if [[ $PHASING_PROGRAM == "eagle" || $PHASING_PROGRAM == "beagle" || $PHASING_PROGRAM == "shapeit4" || $PHASING_PROGRAM == "glimpse" || $PHASING_PROGRAM == "whatsHap" || $PHASING_PROGRAM == "stitch" ]]; then
      source $IMPUTEGEN_DIR/program/05_phasing_"$PHASING_PROGRAM".sh
   fi
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
   ######################
   ##### Imputation #####
   ######################
   (
   if [[ "$IMPUTE_PROGRAM" != 'none' ]]; then
      if [[ "$REF_PANEL_DIR" != "none" ]]; then
         source $IMPUTEGEN_DIR/program/06_imputation_"$IMPUTE_PROGRAM".sh
      
      else
         
         echo
         echo "SKIP IMPUTATION: there is no Reference panel for imputation analysis"
         echo "Phased low density panel files are in $PHASE"
         exit 1
         
      fi
   fi
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
   if [[ $IMPUTE_PROGRAM == "beagle" || $IMPUTE_PROGRAM == "impute5" ]]; then
      if [[ $(bcftools view -H "$IMPUTE/imputed_"$FileExtension".vcf.gz" | wc -l) == 0 ]]; then
         echo "#@ ERROR: No imputed_"$FileExtension".vcf.gz in $IMPUTE. Please check any possible error on Imputation Step"
         exit 1
      fi
   else
      if [[ $(bcftools view -H "$IMPUTE/imputed_"$FileExtension".dose.vcf.gz" | wc -l) == 0 ]]; then
         echo "#@ ERROR: No imputed_"$FileExtension".vcf.gz in $IMPUTE. Please check any possible error on Imputation Step"
         exit 1
      fi
   fi
   wait

   #######################################
   ##### Accuracy - prepare the data #####
   #######################################
   (
   if [[ $DOWN_SCALING != 0 || $VAL_PANEL_DIR != "none" ]]; then
   
      if [[ $IMPUTE_PROGRAM == "beagle" ]]; then
         source $IMPUTEGEN_DIR/program/07_acc_downScaled_externalVal_beagle.sh

      elif [[ $IMPUTE_PROGRAM == "impute5" ]]; then
         source $IMPUTEGEN_DIR/program/07_acc_downScaled_externalVal_impute5.sh
         
      else
         source $IMPUTEGEN_DIR/program/07_acc_downScaled_externalVal_minimac3.sh

      fi
      wait
      
   else
   
      if [[ $IMPUTE_PROGRAM == "beagle" ]]; then
         source $IMPUTEGEN_DIR/program/acc_beagle.sh

      elif [[ $IMPUTE_PROGRAM == "impute5" ]]; then
         source $IMPUTEGEN_DIR/program/acc_impute5.sh

      else
         source $IMPUTEGEN_DIR/program/acc_minimac3.sh

      fi
      wait
   
   fi
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   wait
   
) &
   
   # limit jobs
   if (( $(($((++n)) % $BATCH)) == 0 )) ; then
   wait # wait until all have finished (not optimal, but most times good enough)
   echo $n Files completed wait
   fi
   
done
wait


#############################################
##### Create initial files for accuracy #####
#############################################
export NPOS=$(echo $SELECT_CHROM | tr " " "\n" | wc -l)
COUNT=0

if [[ $DOWN_SCALING != 0 || $VAL_PANEL_DIR != "none" ]]; then
   ### CREATE AN EMPTY ARRAY HERE
   Rscript -e "accPath=\"$ACC\";npos=\"$NPOS\";source(\"$IMPUTEGEN_DIR/program/08_acc1_downScaled_externalVal.R\")"
else
   ### CREATE AN EMPTY ARRAY HERE
   Rscript -e "accPath=\"$ACC\";npos=\"$NPOS\";source(\"$IMPUTEGEN_DIR/program/08_acc1.R\")"
fi
wait


###
for REGION in $SELECT_CHROM; do
   export COUNT=$(( $COUNT + 1 ))
   export REGION=$REGION
   export CHR=$(echo $REGION | awk -F":" '{print $1}')
   export BPSTART=$(echo $REGION | awk -F":" '{print $2}' | awk -F"-" '{print $1}')
   export BPEND=$(echo $REGION | awk -F":" '{print $2}' | awk -F"-" '{print $2}')
   
   export REF_PANEL=$REF_PANEL_DIR/${REF_PANEL_PREFIX}${CHR}${EXTENSION}
   if [[ $VAL_PANEL_DIR != "none" ]]; then export VAL_PANEL=$VAL_PANEL_DIR/${VAL_PANEL_PREFIX}${CHR}${EXTENSION_VAL}; fi
   
   if [[ -z "$BPSTART" ]]; then export BPSTART=0; fi
   if [[ -z "$BPEND" ]]; then export BPEND=1000000000; fi
   
   if [[ "$BPSTART" == 0 && "$BPEND" == 1000000000 ]]; then
      export FileExtension=$CHR
   else
      export FileExtension=$CHR_"$BPSTART"_"$BPEND"
   fi
   wait

   #########################
   ##### Accuracy in R #####
   #########################
   
   (
   echo "#@##################################################################"
   echo "#@ Rscript (plot and summary): Chromosome $REGION"
   echo "#@" 
   D1=$(date "+%D    %T")
   echo "Running Rscript: chromosome $FileExtension  $D1"
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1
   
   START2=$(date +%s)
   
   if [[ $DOWN_SCALING != 0 || $VAL_PANEL_DIR != "none" ]]; then
      NROW_OBS=$(cat $ACC/012/sampleNamesObs_"$FileExtension".txt | wc -l)
      NCOL_OBS=$(cat $ACC/012/mapObs_"$FileExtension".txt | wc -l)
      NROW_PRED=$(cat $ACC/012/sampleNamesPred_$FileExtension.txt | wc -l)
      NCOL_PRED=$(cat $ACC/012/mapPred_"$FileExtension".txt | wc -l)

      Rscript -e "accPath=\"$ACC\";referenceGenome=\"$REF_GENOME\";genoPanel=\"$GENO_PANEL\";count=\"$COUNT\";numColObs=\"$NCOL_OBS\";numRowObs=\"$NROW_OBS\";numColPred=\"$NCOL_PRED\";numRowPred=\"$NROW_PRED\";pos=\"$FileExtension\";source(\"$IMPUTEGEN_DIR/program/08_acc2_downScaled_externalVal.R\")"
      wait
      Rscript -e "accPath=\"$ACC\";pos=\"$FileExtension\";source(\"$IMPUTEGEN_DIR/program/09.alleleFreq_LD_downScaled_externalVal.R\")"
      wait
      gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=$ACC/plots/plots_"$FileExtension".pdf $ACC/plots/page1_"$FileExtension".pdf $ACC/plots/page2_"$FileExtension".pdf $ACC/plots/page3_"$FileExtension".pdf $ACC/plots/page4_"$FileExtension".pdf $ACC/plots/page5_"$FileExtension".pdf
      wait
      rm $ACC/plots/page1_"$FileExtension".pdf $ACC/plots/page2_"$FileExtension".pdf $ACC/plots/page3_"$FileExtension".pdf $ACC/plots/page4_"$FileExtension".pdf $ACC/plots/page5_"$FileExtension".pdf $ACC/plots/Rplots.pdf $ACC/stats/Rplots.pdf 2> /dev/null
      wait
      rm $ACC/012/obsGRM_"$FileExtension".* $ACC/012/predGRM_"$FileExtension".* $ACC/012/pca_"$FileExtension".* 2> /dev/null
      wait

      NROW_OBS=$(cat $ACC/012/plinkPhasedObs_"$FileExtension".fam | wc -l)
      NCOL_OBS=$(cat $ACC/012/plinkPhasedObs_"$FileExtension".bim | wc -l)
      NROW_PRED=$(cat $ACC/012/plinkPhasedPred_"$FileExtension".fam | wc -l)
      NCOL_PRED=$(cat $ACC/012/plinkPhasedPred_"$FileExtension".bim | wc -l)      
      Rscript -e "accPath=\"$ACC\";numColObs=\"$NCOL_OBS\";numRowObs=\"$NROW_OBS\";numColPred=\"$NCOL_PRED\";numRowPred=\"$NROW_PRED\";pos=\"$FileExtension\";source(\"$IMPUTEGEN_DIR/program/08_acc3_downScaled_externalVal_phasing.R\")"
      wait
      rm $PHASE/plink* > /dev/null 2>&1
      rm $PHASE/pos_"$FileExtension".txt > /dev/null 2>&1
      rm $ACC/012/plinkPhasedPred_"$FileExtension".* > /dev/null 2>&1 
      rm $ACC/012/phasedPred012_"$FileExtension".* > /dev/null 2>&1 
      rm $ACC/012/plinkPhasedObs_"$FileExtension".* > /dev/null 2>&1
      rm $ACC/012/phasedObs012_"$FileExtension".* > /dev/null 2>&1 
      wait
      

   else
      NROW=$(cat $ACC/012/sampleNames_"$FileExtension".txt | wc -l)
      NCOL=$(cat $ACC/012/map_"$FileExtension".txt | wc -l)

      Rscript -e "accPath=\"$ACC\";referenceGenome=\"$REF_GENOME\";genoPanel=\"$GENO_PANEL\";count=\"$COUNT\";numCol=\"$NCOL\";numRow=\"$NROW\";pos=\"$FileExtension\";source(\"$IMPUTEGEN_DIR/program/08_acc2.R\")"
      wait
      Rscript -e "accPath=\"$ACC\";pos=\"$FileExtension\";source(\"$IMPUTEGEN_DIR/program/09.alleleFreq_LD.R\")"
      wait
      gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=$ACC/plots/plots_"$FileExtension".pdf $ACC/plots/page1_"$FileExtension".pdf $ACC/plots/page2_"$FileExtension".pdf $ACC/plots/page3_"$FileExtension".pdf
      wait
      rm $ACC/plots/page1_"$FileExtension".pdf $ACC/plots/page2_"$FileExtension".pdf $ACC/plots/page3_"$FileExtension".pdf $ACC/plots/Rplots.pdf $ACC/stats/Rplots.pdf 2> /dev/null
      wait
      rm $ACC/012/genoGRM_"$FileExtension".* 2> /dev/null
      wait
      
   fi
   wait


   END2=$(date +%s)
   DIFF2=$(( $END2 - $START2 ))


   (
   echo
   echo "#@ Rscript (plot and summary) Chromosome $REGION took $(printf '%dh:%dm:%ds\n' $(($DIFF2/3600)) $(($DIFF2%3600/60)) $(($DIFF2%60)))"
   echo "#@#############################################################"
   ) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/"$OUTPUT_NAME"_"$FileExtension".txt 2>&1

   
done
wait


cd "$OUTPUT_DIR"/imputeGen_"$OUTPUT_NAME"
mkdir -p imputeGen_"$OUTPUT_NAME"_report 2> /dev/null
cp -r accuracy/plots imputeGen_"$OUTPUT_NAME"_report 2> /dev/null
cp -r accuracy/summary imputeGen_"$OUTPUT_NAME"_report 2> /dev/null
zip -q -u -r imputeGen_"$OUTPUT_NAME"_report.zip imputeGen_"$OUTPUT_NAME"_report 2> /dev/null
rm -r imputeGen_"$OUTPUT_NAME"_report 2> /dev/null
wait


#while ! grep -q "took" $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/log_"$OUTPUT_NAME".txt; do
#   sleep 5m
#   echo >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/log_"$OUTPUT_NAME".txt
#   echo "ImputeGen still running ... $(date +'%r')  -  $(date +'%m/%d/%Y')" >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/log_"$OUTPUT_NAME".txt
#done
#wait


END1=$(date +%s)
DIFF1=$(( $END1 - $START1 ))

(
echo
echo "#@ ImputeGen $OUTPUT_NAME took $(printf '%dh:%dm:%ds\n' $(($DIFF1/3600)) $(($DIFF1%3600/60)) $(($DIFF1%60)))"
echo "#@#############################################################"
) >> $OUTPUT_DIR/imputeGen_"$OUTPUT_NAME"/log/log_"$OUTPUT_NAME".txt
