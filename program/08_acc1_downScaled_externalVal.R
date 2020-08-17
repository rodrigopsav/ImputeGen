#---------- INITIAL PARAMETERS ----------#
accPath=accPath
npos=as.numeric(npos)
dir.create(accPath, recursive = TRUE)

#---------- Create Empty Array for Imputation Accuracy ----------#
# https://stat.ethz.ch/pipermail/r-help/2008-August/170305.html
accuracy <- array(list(NULL), c(npos, 3))
row.names(accuracy) <- 1:npos
colnames(accuracy) <- c("concRate", "r2_calc", "r2_maf")
saveRDS(accuracy, file=paste0(accPath, "/accImputation.rds"))
