#!/bin/bash


##### SET PARAMETERS ivdp_run.sh #####
usage() { echo "Usage: $0 -a <analysis id> -o <ivdp output path>
       " 1>&2; exit 1; }


while getopts ":a:o:" option; do
   case "${option}" in
   a) ANALYSIS_ID=${OPTARG};;
   o) OUTPUT_PATH=${OPTARG};;
   *) usage;;
   esac
done
shift $((OPTIND -1))

##### TEST PARAMETERS #####
if [[ -z "OUTPUT_PATH" ]]; then
   echo "ERROR: missing parameter "
   usage
   exit 1

else

   export OUTPUT_PATH=$(readlink -f $OUTPUT_PATH)
   
   if [[ ! -d "$OUTPUT_PATH" ]]; then
      echo "ERROR: wrong directory path. Please check -o flag"
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait

###########################

OUTPUT_NAME=$(basename $(echo $OUTPUT_PATH) | cut -f 2- -d"_")

echo
echo
echo
echo "##### Check if ANALYSIS:$ANALYSIS_ID still exist #####"
if `sacct --user pelicion --format="JobID,JobName%50,AllocCPUS,CPUTime,Elapsed,State,ExitCode,DerivedExitCode" | grep  -q $ANALYSIS_ID`; then
   echo ANALYSIS NAME: $OUTPUT_NAME
   echo JOBS COMPLETED: `sacct --user pelicion --format="JobID,JobName%50,AllocCPUS,CPUTime,Elapsed,State,ExitCode,DerivedExitCode" | grep  $ANALYSIS_ID | grep "COMPLETED" | wc -l`
   echo JOBS RUNNING: `sacct --user pelicion --format="JobID,JobName%50,AllocCPUS,CPUTime,Elapsed,State,ExitCode,DerivedExitCode" | grep  $ANALYSIS_ID | grep "RUNNING" | wc -l`
   #echo JOBS PENDING: `sacct --user pelicion --format="JobID,JobName%50,AllocCPUS,CPUTime,Elapsed,State,ExitCode,DerivedExitCode" | grep  $ANALYSIS_ID | grep "PENDING" | wc -l`
   echo JOBS PENDING: `squeue -u $USER | grep  $ANALYSIS_ID | grep "(Dependency)\|(Priority)\|(QOSMaxCpuPerUserLimit)" | wc -l`
   echo JOBS FAILED: `sacct --user pelicion --format="JobID,JobName%50,AllocCPUS,CPUTime,Elapsed,State,ExitCode,DerivedExitCode" | grep  $ANALYSIS_ID | grep "DependencyNeverSatisfied" | wc -l`
   echo JOBS OUT OF MEMORY: `sacct --user pelicion --format="JobID,JobName%50,AllocCPUS,CPUTime,Elapsed,State,ExitCode,DerivedExitCode" | grep  $ANALYSIS_ID | grep "OUT_OF_ME" | wc -l`
   echo JOBS DEPENDENCY NEVER SATISTIED: `squeue -u $USER | grep  $ANALYSIS_ID | grep "(DependencyNeverSatisfied)" | wc -l`
else
   echo "INFORMATION ABOUT ANALYSIS:$OUTPUT_NAME IS NOT AVAILABLE ANYMORE"
   exit 1
fi


echo
echo
echo "##### All jobs in queue list #####"
echo ANALYSIS NAME: $OUTPUT_NAME
squeue -u $USER --format="%.18i %.6P %.50j %.2t %.10M %.6D %R" | grep $ANALYSIS_ID


echo
echo
echo "##### Jobs with status: DependencyNeverSatisfied #####"
echo ANALYSIS NAME: $OUTPUT_NAME
squeue -u $USER --format="%.18i %.6P %.50j %.2t %.10M %.6D %R" | grep $ANALYSIS_ID | grep "DependencyNeverSatisfied"


echo
echo
echo "##### Jobs with status: running, completed, and failed #####"
sacct --user pelicion --format="JobID%10,JobName%50,AllocCPUS%3,CPUTime%12,Elapsed%12,State,ExitCode%4,DerivedExitCode%4" | grep $ANALYSIS_ID > $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info
wait

rm $OUTPUT_PATH/logSlurm/memo_"$OUTPUT_NAME".info > /dev/null 2>&1
for job in $(sacct --user pelicion --format="JobID%10,JobName%50,AllocCPUS%3,CPUTime%12,Elapsed%12,State,ExitCode%4,DerivedExitCode%4" | grep $ANALYSIS_ID | awk '{print $1}'); do 
   seff $job | grep "Memory Utilized\|Memory Efficiency" | awk -F: '{print $2}' | sed -e "s/[[:space:]]\+/_/g" | cut -d '_' -f 2- > $OUTPUT_PATH/logSlurm/tmp.txt

   awk '
   { 
       for (i=1; i<=NF; i++)  {
           a[NR,i] = $i
       }
   }
   NF>p { p = NF }
   END {    
       for(j=1; j<=p; j++) {
           str=a[1,j]
           for(i=2; i<=NR; i++){
               str=str" "a[i,j];
           }
           print str
       }
   }' $OUTPUT_PATH/logSlurm/tmp.txt >> $OUTPUT_PATH/logSlurm/memo_"$OUTPUT_NAME".info
   wait
   
   rm $OUTPUT_PATH/logSlurm/tmp.txt > /dev/null 2>&1

done
wait

paste $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info $OUTPUT_PATH/logSlurm/memo_"$OUTPUT_NAME".info -d '\t' > $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info.tmp && mv $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info.tmp $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info
echo ANALYSIS NAME: $OUTPUT_NAME > $OUTPUT_PATH/logSlurm/header1.txt
echo >> $OUTPUT_PATH/logSlurm/header1.txt
echo >> $OUTPUT_PATH/logSlurm/header1.txt
echo $'JobID\tJobName\tAllocCPUS\tCPUTime\tElapsed(wall-time)\tState\tExitCode\tDerivedExitCode\tMemUsed\tMemEfficiency' > $OUTPUT_PATH/logSlurm/header2.txt
cat $OUTPUT_PATH/logSlurm/header1.txt $OUTPUT_PATH/logSlurm/header2.txt $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info > $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info.tmp && mv $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info.tmp $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info
cat $OUTPUT_PATH/logSlurm/analysis_"$OUTPUT_NAME".info
wait

rm $OUTPUT_PATH/logSlurm/header1.txt $OUTPUT_PATH/logSlurm/header2.txt > /dev/null 2>&1
rm $OUTPUT_PATH/logSlurm/memo_"$OUTPUT_NAME".info > /dev/null 2>&1
