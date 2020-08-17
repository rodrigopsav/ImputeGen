#---------- INITIAL PARAMETERS ----------#
my_packages <- c("data.table", "dplyr", "ggplot2", "ggrepel", "gridExtra", "wesanderson", "fields", "MASS") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

accPath=accPath
pos=pos
numRowObs=as.numeric(numRowObs)
numColObs=as.numeric(numColObs)
numRowPred=as.numeric(numRowPred)
numColPred=as.numeric(numColPred)

#---------- Observed genotypes ----------#
setwd(paste0(accPath, "/012"))
print("loading observed genotypes (012): Validation population")
geno_obs <- matrix(scan(paste0("phasedObs012_", pos, ".raw"),what = "character"), nrow = (numRowObs + 1), ncol = (numColObs + 6), byrow=T)
geno_obs <- geno_obs[-c(1), -c(1:6)]
geno_obs <- t(geno_obs)
dim(geno_obs)
geno_obs <- as.data.frame(geno_obs)


#---------- Predicted genotypes ----------#
setwd(paste0(accPath, "/012"))
print("loading predicted genotypes (012): Imputed population")
geno_pred <- matrix(scan(paste0("phasedPred012_", pos, ".raw"),what = "character"), nrow = (numRowPred + 1), ncol = (numColPred + 6), byrow=T)
geno_pred <- geno_pred[-c(1), -c(1:6)]
geno_pred <- t(geno_pred)
dim(geno_pred)
geno_pred <- as.data.frame(geno_pred)


#---------- Concordance Rate ----------#
print("Calculating concordance rate")
accuracy_CR <- list(3)

compare <- geno_obs == geno_pred
dim(compare)
CR <- sum(compare)/(as.numeric(dim(compare)[1]) * as.numeric(dim(compare)[2]))

CR_ind <- colSums(compare, na.rm = TRUE)/dim(compare)[1]
CR_ind <- matrix(CR_ind, ncol=1)
colnames(CR_ind) <- "concRate" 

CR_snp <- rowSums(compare, na.rm = TRUE)/dim(compare)[2]
CR_snp <- matrix(CR_snp, ncol=1)
colnames(CR_snp) <- c("concRate")

accuracy_CR[[1]] <- CR
accuracy_CR[[2]] <- CR_ind
accuracy_CR[[3]] <- CR_snp


#---------- Correlation ----------#
print("Calculating correlation")

# add small number to avoid error in correlation
small <- matrix(runif(nrow(geno_pred) * ncol(geno_pred),
                      1e-15,1e-12),
                nrow = nrow(geno_pred), ncol = ncol(geno_pred))
geno_obs <- as.matrix(geno_obs)
mode(geno_obs) <- "numeric"
geno_obs <- geno_obs + small

geno_pred <- as.matrix(geno_pred)
mode(geno_pred) <- "numeric"
geno_pred <- geno_pred + small

accuracy_COR <- list(3)
COR <- cor(c(geno_obs), c(geno_pred))

COR_ind <- sapply(1:ncol(geno_obs), function(i) cor(geno_obs[,i], geno_pred[,i]))
COR_ind <- matrix(COR_ind, ncol=1)
colnames(COR_ind) <- "r2_calc" 

COR_snp <- sapply(1:nrow(geno_obs), function(i) cor(geno_obs[i,], geno_pred[i,]))
COR_snp <- matrix(COR_snp, ncol=1)
colnames(COR_snp) <- c("r2_calc")

accuracy_COR[[1]] <- COR
accuracy_COR[[2]] <- COR_ind
accuracy_COR[[3]] <- COR_snp


#---------- Summary Statistics ----------#
dir.create(paste0(accPath, "/summary"), recursive = T)

sink(paste0(accPath, "/summary/summaryStatsPhasing_", pos, ".txt"))
cat("#----- General Concordance Rate -----#\n")
cat(accuracy_CR[[1]])
cat("\n")
cat("\n")
cat("#----- Concordance Rate by Subject -----#\n")
print(summary(accuracy_CR[[2]][ ,1]))
cat("\n")
cat("\n")
cat("#----- Concordance Rate by SNPs -----#\n")
print(summary(accuracy_CR[[3]][ ,1]))
cat("\n")
cat("\n")

cat("#----- General Correlation -----#\n")
cat(accuracy_COR[[1]])
cat("\n")
cat("\n")
cat("#----- Correlation Rate by Subject -----#\n")
print(summary(accuracy_COR[[2]][ ,1]))
cat("\n")
cat("\n")
cat("#----- Correlation Rate by SNPs -----#\n")
print(summary(accuracy_COR[[3]][ ,1]))
sink()