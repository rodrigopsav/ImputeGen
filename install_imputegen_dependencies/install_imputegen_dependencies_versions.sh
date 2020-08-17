#!/bin/bash

################################
##### INSTALL DEPENDENCIES #####
################################

##### SET PARAMETERS
usage() { echo "Usage: $0 -d path/to/install/dependencies/programs
       -d: directory path to install dependencies programs
       " 1>&2; exit 1; }

while getopts :s:d: option; do
   case "${option}" in
   d) INSTALL_FOLDER=${OPTARG};;
   *) usage;;
   esac
done
#shift "$((OPTIND-1))"

##### TEST INSTALL_FOLDER VARIABLE
if [[ -z "$INSTALL_FOLDER" ]]; then
   echo "Error: -d flag is empty"
   usage
   exit 1

else

   export INSTALL_FOLDER=$(readlink -f $INSTALL_FOLDER)
   
   if [[ ! -d "$INSTALL_FOLDER" ]]; then
      echo "ERROR: wrong directory path. Please check -d flag"
      echo "Aborting analysis"
      usage
      exit 1
   fi
fi
wait
   
echo "Install programs in: "$INSTALL_FOLDER


##### Make Installation folder
#export INSTALL_FOLDER=$INSTALL_FOLDER
#mkdir -p $INSTALL_FOLDER


echo "
#--------------------------------#
##### INSTALLING MINICONDA 3 #####
#--------------------------------#
"

if ! command -v conda &> /dev/null; then
   echo "conda CANNOT BE FOUND"
   
   # Make Installation folder
   export INSTALL_FOLDER=$INSTALL_FOLDER
   mkdir -p $INSTALL_FOLDER

   cd $INSTALL_FOLDER
   if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      chmod +x Miniconda3-latest-Linux-x86_64.sh
      # https://docs.anaconda.com/anaconda/install/silent-mode/
      ./Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALL_FOLDER/miniconda3 -f
      source $INSTALL_FOLDER/miniconda3/bin/activate
      conda init bash
      #https://www.kangzhiq.com/2020/05/02/how-to-activate-a-conda-environment-in-a-remote-machine/
      eval "$(conda shell.bash hook)"
      rm Miniconda3-latest-Linux-x86_64.sh
   elif [[ "$OSTYPE" == "darwin"* ]]; then
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
      chmod +x Miniconda3-latest-MacOSX-x86_64.sh
      ./Miniconda3-latest-MacOSX-x86_64.sh -b -p $INSTALL_FOLDER/miniconda3 -f
      source $INSTALL_FOLDER/miniconda3/bin/activate
      conda init bash
      #https://www.kangzhiq.com/2020/05/02/how-to-activate-a-conda-environment-in-a-remote-machine/
      eval "$(conda shell.bash hook)"
      rm Miniconda3-latest-MacOSX-x86_64.sh
   else
      "Error: This is neither linux nor MacOS system"
       exit 1
   fi

else
   echo "conda is installed" 
fi
   
### Configure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels anaconda
conda config --add channels conda-forge

### Create conda environment to save Bioinformatic programs
#https://www.kangzhiq.com/2020/05/02/how-to-activate-a-conda-environment-in-a-remote-machine/
conda init bash
eval "$(conda shell.bash hook)"
conda install mamba -y -n base -c conda-forge
mamba create -y --name imputegen python=3.8
mamba create -y --name r-env4.0 -c conda-forge r-base=4.0.3 python=3.8


### Activate conda imputegen environment
conda activate imputegen


echo "
#-----------------------------#
##### INSTALLING parallel #####
#-----------------------------#
"
mamba install -y -n imputegen -c conda-forge parallel=20210222


echo "
#-----------------------------#
##### INSTALLING dos2unix #####
#-----------------------------#
"
mamba install -y -n imputegen -c trent dos2unix=6.0.3


echo "
#---------------------------#
##### INSTALLING rename #####
#---------------------------#
"
mamba install -y -n imputegen -c bioconda rename


echo "
#------------------------#
##### INSTALLING sed #####
#------------------------#
"
mamba install -y -n imputegen -c conda-forge sed=4.8 


echo "
#------------------------#
##### INSTALLING git #####
#------------------------#
"
mamba install -y -n imputegen -c anaconda git=2.28.0


echo "
#-----------------------------------------------#
##### INSTALLING htslib, samtools, bcftools #####
#-----------------------------------------------#
"
mamba install -y -n imputegen -c bioconda htslib=1.9
mamba install -y -n imputegen -c bioconda samtools=1.9 --force-reinstall
mamba install -y -n imputegen -c bioconda bcftools=1.9


echo "
#-----------------------------#
##### INSTALLING BEDTOOLS #####
#-----------------------------#
"
mamba install -y -n imputegen -c bioconda bedtools=2.29.2


echo "
#-----------------------------#
##### INSTALLING vcftools #####
#-----------------------------#
"
mamba install -y -n imputegen -c bioconda vcftools=0.1.16
mamba install -y -n imputegen -c bioconda perl-vcftools-vcf


echo "
#---------------------------#
##### INSTALLING VCFLIB #####
#---------------------------#
"
mamba install -y -n imputegen -c bioconda vcflib=1.0.2-3


echo "
#--------------------------#
##### INSTALLING plink #####
#--------------------------#
"
mamba install -y -n imputegen -c bioconda plink=1.90b6.18
mamba install -y -n imputegen -c bioconda plink2=v2.00a2.3LM

wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip -O $CONDA_PREFIX/bin/plink.zip
unzip -o $CONDA_PREFIX/bin/plink.zip -d $CONDA_PREFIX/bin
rm $CONDA_PREFIX/bin/plink.zip

wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20220129.zip -O $CONDA_PREFIX/bin/plink2.zip
unzip -o $CONDA_PREFIX/bin/plink2.zip -d $CONDA_PREFIX/bin
rm $CONDA_PREFIX/bin/plink2.zip


echo "
#-------------------------#
##### INSTALLING GCTA #####
#-------------------------#
"
mamba install -y -n imputegen -c bioconda gcta


echo "
#--------------------------#
##### INSTALLING eagle #####
#--------------------------#
"
#conda install -y -n imputegen -c soil eagle-phase=v2.3.5
wget -P $INSTALL_FOLDER https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
wait
tar -xf $INSTALL_FOLDER/Eagle_v2.4.1.tar.gz
wait
CONDA=$(dirname $CONDA_PREFIX)
wait
mv $INSTALL_FOLDER/Eagle_v2.4.1/eagle $CONDA/imputegen/bin/eagle
wait
rm -rf $INSTALL_FOLDER/Eagle_v2* 2> /dev/null
wait


echo "
#----------------------------#
##### INSTALLING glimpse #####
#----------------------------#
"
git clone https://github.com/odelaneau/GLIMPSE.git $INSTALL_FOLDER/glimpse
wait
CONDA=$(dirname $CONDA_PREFIX)
wait
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_chunk_static $CONDA/imputegen/bin/GLIMPSE_chunk
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_concordance_static $CONDA/imputegen/bin/GLIMPSE_concordance
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_ligate_static $CONDA/imputegen/bin/GLIMPSE_ligate
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_phase_static $CONDA/imputegen/bin/GLIMPSE_phase
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_sample_static $CONDA/imputegen/bin/GLIMPSE_sample
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_snparray_static $CONDA/imputegen/bin/GLIMPSE_snparray
mv $INSTALL_FOLDER/glimpse/static_bins/GLIMPSE_stats_static $CONDA/imputegen/bin/GLIMPSE_stats
wait
rm -rf $INSTALL_FOLDER/glimpse 2> /dev/null
wait

echo "
#-----------------------------#
##### INSTALLING shapeit4 #####
#-----------------------------#
"
mamba install -y -n imputegen -c bioconda shapeit4=4.1.3


echo "
#-----------------------------#
##### INSTALLING whatshap #####
#-----------------------------#
"
pip install whatshap


echo "
#---------------------------#
##### INSTALLING stitch #####
#---------------------------#
"
mamba install  -y -n r-env4.0 -c bioconda r-stitch
cp ../program/STITCH.R $(dirname $CONDA_PREFIX)/imputegen/bin/STITCH.R


echo "
#---------------------------#
##### INSTALLING quilt #####
#---------------------------#
"
mamba install  -y -n r-env4.0 -c defaults -c bioconda -c conda-forge r-quilt
cp ../program/QUILT.R $(dirname $CONDA_PREFIX)/imputeme/bin/QUILT.R


echo "
#----------------------------#
##### INSTALLING impute2 #####
#----------------------------#
"
mamba install -y -n imputegen -c bioconda impute2=2.3.2


echo "
#---------------------------#
##### INSTALLING beagle #####
#---------------------------#
"
#conda install -y -n imputegen -c bioconda beagle=5.1_24Aug19.3e8
wget -P $INSTALL_FOLDER https://faculty.washington.edu/browning/beagle/beagle.21Apr21.304.jar
wait
CONDA=$(dirname $CONDA_PREFIX)
wait
mv $INSTALL_FOLDER/beagle*jar $CONDA/imputegen/bin/beagle
wait


echo "
#-----------------------------#
##### INSTALLING liftover #####
#-----------------------------#
"
mamba install -y -n imputegen -c bioconda ucsc-liftover


echo "
#---------------------------#
##### INSTALLING hapbin #####
#---------------------------#
"
# hapbin: Extended Haplotype Homozygosity (EHH), the Integrated Haplotype Score (iHS) and the Cross Population Extended Haplotype Homozogysity (XP-EHH) statistic.
# https://github.com/evotools/hapbin
mamba install -y -n imputegen -c bioconda/label/cf201901 hapbin


##### Population structure software

echo "
#------------------------------#
##### INSTALLING admixture #####
#------------------------------#
"
mamba install -y -n imputegen -c bioconda admixture=1.3


echo "
#------------------------------#
##### INSTALLING structure #####
#------------------------------#
"
mamba install -y -n imputegen -c bioconda structure=2.3.4


echo "
#--------------------------------#
##### INSTALLING ghostscript #####
#--------------------------------#
"
mamba install -y -n imputegen -c conda-forge ghostscript


echo "
#----------------------------#
##### INSTALLING OPENJDK #####
#----------------------------#
#"
mamba install -y -n imputegen -c anaconda openjdk=8.0.282


echo "
#------------------------------#
##### INSTALLING R Program #####
#------------------------------#
"
mamba install -y -n r-env4.0 "libblas=*=*openblas"
mamba install -y -n r-env4.0 -c conda-forge r-dplyr
mamba install -y -n r-env4.0 -c conda-forge r-data.table
mamba install -y -n r-env4.0 -c conda-forge r-ggplot2
mamba install -y -n r-env4.0 -c conda-forge r-ggrepel
mamba install -y -n r-env4.0 -c conda-forge r-gridextra
mamba install -y -n r-env4.0 -c conda-forge r-vcfr
mamba install -y -n r-env4.0 -c conda-forge r-wesanderson
mamba install -y -n r-env4.0 -c conda-forge r-fields
mamba install -y -n r-env4.0 -c conda-forge r-mass
mamba install -y -n r-env4.0 -c conda-forge r-optparse


##### Deactivate ivdp conda environment
conda deactivate


echo "
#----------------------------------#
##### CHECK INSTALLED PROGRAMS #####
#----------------------------------#
"

echo "#############################################################################"
echo "############## CHECK IF ALL THE PROGRAMS ARE INSTALLED ######################"
echo "If any installed program failed, re-install a previous version. To do that:"
echo "Example using beagle package"
echo "Search for versions: conda search beagle --info"
echo "Choose a previous version: eg. 5.1_24Aug19.3e8"
echo "Install: conda install --force-reinstall -y -n imputegen -c bioconda beagle=5.1_24Aug19.3e8"
echo "#############################################################################"
echo

### Activate conda ivdp environment 
eval "$(conda shell.bash hook)"
conda activate imputegen
conda activate --stack r-env4.0

###
if ! command -v parallel &> /dev/null; then
   echo "parallel COULD NOT BE FOUND"
else
   echo "parallel is installed" 
fi

###
if ! command -v dos2unix &> /dev/null; then
   echo "dos2unix COULD NOT BE FOUND"
else
   echo "dos2unix is installed" 
fi

###
if ! command -v rename &> /dev/null; then
   echo "rename COULD NOT BE FOUND"
else
   echo "rename is installed" 
fi

###
if ! command -v sed &> /dev/null; then
   echo "sed COULD NOT BE FOUND"
else
   echo "sed is installed" 
fi

###
if ! command -v git &> /dev/null; then
   echo "git COULD NOT BE FOUND"
else
   echo "git is installed" 
fi

###
if ! command -v bgzip &> /dev/null; then
   echo "htslib COULD NOT BE FOUND"
else
   echo "htslib is installed" 
fi

###
if ! command -v samtools &> /dev/null; then
   echo "samtools COULD NOT BE FOUND"
else
   echo "samtools is installed" 
fi

###
if ! command -v bcftools &> /dev/null; then
   echo "bcftools COULD NOT BE FOUND"
else
   echo "bcftools is installed" 
fi

###
if ! command -v bedtools &> /dev/null; then
   echo "bedtools CANNOT BE FOUND"
else
   echo "bedtools is installed" 
fi

###
if ! command -v vcftools &> /dev/null; then
   echo "vcftools COULD NOT BE FOUND"
else
   echo "vcftools is installed" 
fi

###
if ! command -v vcffilter &> /dev/null; then
   echo "vcflib CANNOT BE FOUND"
else
   echo "vcflib is installed" 
fi

###
if ! command -v plink &> /dev/null; then
   echo "plink COULD NOT BE FOUND"
else
   echo "plink is installed" 
fi

###
if ! command -v plink2 &> /dev/null; then
   echo "plink2 COULD NOT BE FOUND"
else
   echo "plink2 is installed" 
fi

###
if ! command -v $CONDA/imputegen/bin/eagle &> /dev/null; then
   echo "eagle COULD NOT BE FOUND"
else
   echo "eagle is installed" 
fi

###
if ! command -v $CONDA/imputegen/bin/GLIMPSE_phase &> /dev/null; then
   echo "glimpse COULD NOT BE FOUND"
else
   echo "glimpse is installed" 
fi

###
if ! command -v shapeit4 &> /dev/null; then
   echo "shapeit4 COULD NOT BE FOUND"
else
   echo "shapeit4 is installed" 
fi

###
if ! command -v whatshap &> /dev/null; then
   echo "whatshap COULD NOT BE FOUND"
else
   echo "whatshap is installed" 
fi

###
if ! command -v STITCH.R &> /dev/null; then
   echo "stitch COULD NOT BE FOUND"
else
   echo "stitch is installed" 
fi

###
if ! command -v QUILT.R &> /dev/null; then
   echo "quilt COULD NOT BE FOUND"
else
   echo "quilt is installed" 
fi

###
if ! command -v impute2 &> /dev/null; then
   echo "impute2 COULD NOT BE FOUND"
else
   echo "impute2 is installed" 
fi

###
if ! command -v java -jar $CONDA/imputegen/bin/beagle &> /dev/null; then
   echo "beagle COULD NOT BE FOUND"
else
   echo "beagle is installed" 
fi

###
if ! command -v liftOver &> /dev/null; then
   echo "liftOver COULD NOT BE FOUND"
else
   echo "liftOver is installed" 
fi

###
if ! command -v ehhbin &> /dev/null; then
   echo "hapbin COULD NOT BE FOUND"
else
   echo "hapbin is installed" 
fi

###
if ! command -v admixture &> /dev/null; then
   echo "admixture COULD NOT BE FOUND"
else
   echo "admixture is installed" 
fi

###
if ! command -v structure &> /dev/null; then
   echo "structure COULD NOT BE FOUND"
else
   echo "structure is installed" 
fi

###
if ! command -v gs &> /dev/null; then
   echo "ghostscript COULD NOT BE FOUND"
else
   echo "ghostscript is installed" 
fi

###
if ! command -v R &> /dev/null; then
   echo "R Program COULD NOT BE FOUND"
else
   echo "R Program is installed" 
fi


##### Deactivate imputegen conda environment
conda deactivate
conda deactivate
echo


echo "#####################################################################################"
echo "# TO RUN THESE PROGRAMS FOR OTHER APPLICATION, ACTIVATE IMPUTEGEN CONDA ENVIRONMENT #"
echo "############ TO ACTIVATE, TYPE ON TERMINAL: conda activate imputegen ################"
echo "#####################################################################################"

echo "
#--------------------------------#
##### INSTALLATION COMPLETED #####
#--------------------------------#
"

