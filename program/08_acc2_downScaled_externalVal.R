#---------- INITIAL PARAMETERS ----------#
my_packages <- c("data.table", "dplyr", "reshape2", "ggplot2", "ggrepel", "gridExtra", "wesanderson", "fields", "MASS") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# GCTA: R script to read the GRM binary file (https://cnsgenomics.com/software/gcta/#MakingaGRM)
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

accPath=accPath
referenceGenome=referenceGenome
genoPanel=genoPanel
pos=pos
count=as.numeric(count)
numRowObs=as.numeric(numRowObs)
numColObs=as.numeric(numColObs)
numRowPred=as.numeric(numRowPred)
numColPred=as.numeric(numColPred)

#---------- Load accImputation.rds ----------#
print(paste0("Chromosome: ", pos))
print("loading accImputation.rds")
accuracy <- readRDS(file = paste0(accPath, "/accImputation.rds"))

#---------- R2 and MAF  ----------#
setwd(paste0(accPath, "/r2_maf"))
r2_maf <- read.table(paste0("r2_maf_", pos, ".txt"), header = T, sep = "\t")

#---------- Observed genotypes ----------#
setwd(paste0(accPath, "/012"))
print("loading individuals and position of SNPs: Validation population")
ind_obs <- read.table(paste0("sampleNamesObs_", pos, ".txt"), header=F)
ind_obs <- ind_obs[,1]
ind_obs <- t(ind_obs)
map_obs <- read.table(paste0("mapObs_", pos, ".txt"), header=F)

print("loading observed genotypes (012): Validation population")
geno_obs <- matrix(scan(paste0("obs012_", pos, ".raw"),what = "character"), nrow = (numRowObs + 1), ncol = (numColObs + 6), byrow=T)
geno_obs <- geno_obs[-c(1), -c(1:6)]
geno_obs <- t(geno_obs)
dim(geno_obs)
geno_obs <- as.data.frame(geno_obs)
colnames(geno_obs) <- ind_obs
#geno_obs <- geno_obs[ , order(names(geno_obs))]


#---------- Predicted genotypes ----------#
setwd(paste0(accPath, "/012"))
print("loading individuals and position of SNPs: Imputed population")
ind_pred <- read.table(paste0("sampleNamesPred_", pos, ".txt"), header=F)
ind_pred <- ind_pred[ ,1]
ind_pred <- t(ind_pred)
map_pred <- read.table(paste0("mapPred_", pos, ".txt"), header=F)

print("loading predicted genotypes (012): Imputed population")
geno_pred <- matrix(scan(paste0("pred012_", pos, ".raw"),what = "character"), nrow = (numRowPred + 1), ncol = (numColPred + 6), byrow=T)
geno_pred <- geno_pred[-c(1), -c(1:6)]
geno_pred <- t(geno_pred)
dim(geno_pred)
geno_pred <- as.data.frame(geno_pred)
colnames(geno_pred) <- ind_pred
#geno_pred <- geno_pred[ , order(names(geno_pred))]


#---------- Subset observed and predicted ----------#
index_obs <- paste(map_obs[,1], map_obs[,2], sep=":")
index_pred <- paste(map_pred[,1], map_pred[,2], sep=":")
common <- intersect(index_obs, index_pred)
pos_obs <- index_obs %in% common
pos_pred <- index_pred %in% common

map_obs <- map_obs[pos_obs, ]
map_pred <- map_pred[pos_pred, ]
geno_obs <- geno_obs[pos_obs, ]
geno_pred <- geno_pred[pos_pred, ]


#---------- Concordance Rate ----------#
print("Calculating concordance rate")
compare <- geno_obs == geno_pred
dim(compare)
CR <- sum(compare)/(as.numeric(dim(compare)[1]) * as.numeric(dim(compare)[2]))

CR_ind <- colSums(compare, na.rm = TRUE)/dim(compare)[1]
CR_ind <- matrix(CR_ind, ncol=1)
colnames(CR_ind) <- "concRate" 
rownames(CR_ind) <- ind_pred[order(ind_pred)]

CR_snp <- rowSums(compare, na.rm = TRUE)/dim(compare)[2]
CR_snp <- matrix(CR_snp, ncol=1)
CR_snp <- cbind(map_pred, CR_snp)
colnames(CR_snp) <- c("chrom", "pos", "class", "concRate")

row.names(accuracy)[count] <- pos
accuracy[[count, 1]] <- vector(mode="list", length=3)
accuracy[[count, 1]][[1]] <- CR
accuracy[[count, 1]][[2]] <- CR_ind
accuracy[[count, 1]][[3]] <- CR_snp

#accuracy[[count, 1]][[1]][[1]] # general concordance rate
#accuracy[[count, 1]][[2]][1:5, ] # concordance rate by individual
#accuracy[[count, 1]][[3]][1:5, ] #concordance rate by snp
#accuracy[[count, 2]][1:5,] #R2 and MAF from beagle or minimac3


#------------- Correlation ------------#
# https://stackoverflow.com/questions/9136116/correlation-between-two-dataframes-by-row
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

COR <- cor(c(geno_obs), c(geno_pred))

COR_ind <- sapply(1:ncol(geno_obs), function(i) cor(geno_obs[,i], geno_pred[,i]))
COR_ind <- matrix(COR_ind, ncol=1)
colnames(COR_ind) <- "r2_calc" 
rownames(COR_ind) <- ind_pred[order(ind_pred)]

COR_snp <- sapply(1:nrow(geno_obs), function(i) cor(geno_obs[i,], geno_pred[i,]))
COR_snp <- matrix(COR_snp, ncol=1)
COR_snp <- cbind(map_pred, COR_snp)
colnames(COR_snp) <- c("chrom", "pos", "class", "r2_calc")

accuracy[[count, 2]] <- vector(mode="list", length=3)
accuracy[[count, 2]][[1]] <- COR
accuracy[[count, 2]][[2]] <- COR_ind
accuracy[[count, 2]][[3]] <- COR_snp


#--- R2 and MAF from imputation programs ---#
accuracy[[count, 3]] <- r2_maf


#---------- PLOTS WITH GGPLOT2 ----------#
print("Creating plots")
dir.create(paste0(accPath, "/plots"), recursive = T)

##### Concordance rate plots by markers and individuals #####
#https://www.datamentor.io/r-programming/box-plot/

concRate.tmp <- rbind(accuracy[[count, 1]][[2]], matrix(accuracy[[count, 1]][[3]][,4], ncol=1))
concRate.tmp <- as.data.frame(concRate.tmp)
type <- c(rep("Subject", length(accuracy[[count, 1]][[2]][,1])), rep("SNPs", length(accuracy[[count, 1]][[3]][,4])))
concRate.tmp <-cbind(concRate.tmp, type)
rm(type)

p1 <- ggplot(concRate.tmp, aes(x=type, y=concRate, color=type)) + 
  geom_boxplot(color=c("#999933", "#44AA99"), fill=c("cornsilk1", "azure")) +
  ggtitle(paste0("Imputation concordance rate: chr", pos)) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        legend.position = "right") +
  xlab("") + 
  ylab("Concordance rate")


##### Accuracy vs MAF #####
ImpAcc_maf <- data.frame(concRate = accuracy[[count, 1]][[3]][,4],
                         correl = accuracy[[count, 2]][[3]][,4],
                         maf = accuracy[[count, 3]][ ,5])

ImpAcc_maf$mafBin <- cut(ImpAcc_maf$maf, 
                         seq(0, 0.5, 0.01), 
                         include.lowest=TRUE, right=T,
                         labels=seq(0.01, 0.5, 0.01)
)

ImpAcc_maf_summary <- ImpAcc_maf %>% group_by(mafBin) %>% summarise(concRateSummary=median(concRate), correlSummary=median(correl)) %>% as.data.frame()
ImpAcc_maf_summary <- na.omit(ImpAcc_maf_summary)
ImpAcc_maf_summary$mafBin <- as.numeric.factor(ImpAcc_maf_summary$mafBin)

lm1 <- lm(concRateSummary ~ poly(mafBin, 3), data = ImpAcc_maf_summary)
lm2 <- lm(correlSummary ~ poly(mafBin, 3), data = ImpAcc_maf_summary)

ImpAcc_maf_summary_long <- melt(ImpAcc_maf_summary, id.var = "mafBin")
colnames(ImpAcc_maf_summary_long) <- c("mafBin", "category", "value")
ImpAcc_maf_summary_long$pred <- c(predict(lm1, data=seq(0.01, 0.5, 0.01)),
                                  predict(lm2, data=seq(0.01, 0.5, 0.01))
)

p2 <- ggplot(ImpAcc_maf_summary_long, aes(x = mafBin)) +
  geom_point(aes(y = pred, colour = category), shape=18, size=2.5, alpha=1) +
  scale_color_manual(values = c("darkred", "steelblue"), labels = c("Conc Rate", "Correl")) +
  stat_smooth(aes(x = mafBin, y = value, colour = category, linetype= category), 
              method = "lm", formula = y ~ poly(x, 3), se = FALSE, size=0.8, alpha=0.6) +
  scale_fill_manual(values = c("darkred", "steelblue"), labels = c("Conc Rate", "Correl")) +
  scale_linetype_manual(values = c("solid","dashed"), labels = c("Conc Rate", "Correl")) +
  scale_x_continuous(limits=c(0, 0.5)) + 
  scale_y_continuous(limits=c(0,1)) +
  ggtitle(paste("Imputation Accuracy vs MAF: chr", pos)) +
  theme_light() +
  theme(panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        legend.position = "right",
        legend.text = element_text(size=12, face="bold")) +
  xlab("Minor allele frequency (MAF)") + 
  ylab("Accuracy") +
  labs(color  = "", linetype = "", shape = "")


##### Concordance Rate density plot #####
concRate.tmp2 <- accuracy[[count, 1]][[3]]
# Divide chromosome in segments length = window #####
window=100000
commandAWK <- paste0('awk \'$1 == ', pos, ' {print $2}\' ', referenceGenome, '.fai')
chromLenght <- as.numeric(system(commandAWK, intern = T))

seqBreaks <- seq(0, chromLenght, by=window)
seqBreaks[which.max(seqBreaks)] <- chromLenght
concRate.tmp2$group <- cut(concRate.tmp2$pos,
                           breaks = seqBreaks,
                           include.lowest=TRUE,
                           labels = seqBreaks[-1])

# Calculate average Concordance Rate by SNP by Chromosome segment
conRate_average <- aggregate(concRate.tmp2$concRate,
                             list(concRate.tmp2$group), 
                             mean)
colnames(conRate_average) <- c("group", "average")
conRate_average$height <- 1

# Convert group from factor to numeric
#https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
conRate_average$group <- as.numeric.factor(conRate_average$group)
conRate_average$group <- conRate_average$group/1000000
min(conRate_average$average)
max(conRate_average$average)

library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
p3 <- ggplot(conRate_average, aes(x = group , y = height, fill = average)) +
  geom_line(size=0.01, alpha=0.2) +
  geom_tile() + 
  ggtitle(paste0("Average Concordance Rate (Window of ", window/1000, "kb): Chromosome ", pos)) +
  scale_fill_gradientn(colours = pal, limits=c(0,1)) +
  scale_x_discrete(name = "Position (Mb)", expand = c(0, 0),
                   limits = seq(0, round(max(conRate_average$group), 0), 5)) +
  scale_y_discrete(name = paste0("Chrom", pos), expand = c(0, 0)) +
  coord_fixed(ratio = 5) + 
  theme(axis.title.y=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12,face="bold",angle = 0, vjust=0.5),
        axis.text.x=element_text(size=12),
        plot.title = element_text(size=12,face="bold", hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "cm"), #top, right, bottom, left
        legend.position = "bottom", 
        legend.key.size = unit(1.5, "cm")) 


##### PCA plot of SDV of G matrix #####
eigenvalues <- fread(paste0("pca_", pos, ".eigenval"))
eigenvalues <- as.data.frame(eigenvalues)
varPCA <- round( (eigenvalues$V1^2) / sum(eigenvalues$V1^2) *100, 3)

eigenvectors <- fread(paste0("pca_", pos, ".eigenvec"))
eigenvectors <- as.data.frame(eigenvectors)
eigenvectors <- eigenvectors[,c(3:4, 2)]
colnames(eigenvectors) <- c("pc1", "pc2", "id")
eigenvectors$pop <- rep(c("Observed", "Imputed"), dim(eigenvectors)[1]/2)

# plot PC1 vs PC2 using scores
pc1=1
pc2=2

p4 <- ggplot(eigenvectors, aes_string(x=colnames(eigenvectors)[pc1], y=colnames(eigenvectors)[pc2], col=colnames(eigenvectors)[which(colnames(eigenvectors) == "pop")])) + 
  geom_point(shape=18, size=5) +
  scale_color_manual(values=c("#1F77B4","#FF9896")) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray40", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "gray40", size=1) +
  #geom_text_repel(aes(label=id), size=4, segment.size=0.10, nudge_x=0.2, direction="y") +
  ggtitle(paste0("Singular value decomposition of GRM: chr", pos)) +
  theme_light() +
  theme(#panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face="bold"),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = "right") +
  xlab(paste("PC",pc1, " (",varPCA[pc1],"%)",sep="")) + 
  ylab(paste("PC",pc2, " (",varPCA[pc2],"%)",sep=""))


##### Scatterplot between Observed and Imputed GRM #####
# Get relationships from gcta
G_obs <- ReadGRMBin(paste0("obsGRM_", pos))
G_pred <- ReadGRMBin(paste0("predGRM_", pos))

G_obs_triangular <- c(G_obs[[1]], G_obs[[2]])
G_pred_triangular <- c(G_pred[[1]], G_pred[[2]])

# Correlation between Observed and Imputed GRM
G_correl <- round(cor(G_obs_triangular, G_pred_triangular), 3)
lm_G_correl <- lm(G_pred_triangular ~ G_obs_triangular)
slope <- round(lm_G_correl$coefficients[2], 3)

G_correl_df <- cbind(G_obs_triangular, G_pred_triangular)
G_correl_df <- as.data.frame(G_correl_df)

# I want the color to change in a direction that characterizes the correlation.
# To do this, we can color points by the first principal component
# https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
G_correl_df$pc <- predict(prcomp(~G_obs_triangular+G_pred_triangular, G_correl_df))[,1]

library(MASS)
library(fields)
G_correl_df$density <- fields::interp.surface(
  MASS::kde2d(G_correl_df$G_obs_triangular, 
              G_correl_df$G_pred_triangular),
  G_correl_df[,c("G_obs_triangular", "G_pred_triangular")])


# Plot correlation matrix (G_correl_df)
p5 <- ggplot(G_correl_df, aes(x=G_obs_triangular, y=G_pred_triangular, color = pc, alpha = 1/density)) +
  geom_point(shape = 18, size = 5, show.legend = FALSE) +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  scale_alpha(range = c(0.50, 0.60)) +
  ggtitle(paste0("Imputed vs Observed GRM: chr", pos, "\ncorrelation = ", G_correl, "\nprediction bias (slope) = ", slope)) +
  theme_light() +
  theme(#panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face="bold"),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = "right") +
  xlab("Observed relationships") + 
  ylab("Imputed relationships")

### Plot correlation matrix (with scattermore package)
#https://github.com/exaexa/scattermore
#install.packages("scattermore")
# library("scattermore")
# pdf(paste0(accPath, "/plots/correl_GRM_", pos, ".pdf"))
# scattermoreplot(G_correl_df$G_obs_triangular, G_correl_df$G_pred_triangular, 
#                 main = paste0("Chromosome ", pos, " ==> correl = ", correl),
#                 xlab = "Observed Relationship", 
#                 ylab = "Imputed Relationship",
#                 col = rgb(0.00, 0.50, 0.80, 1))
# dev.off()

##### Scatterplot between Observed and Imputed Heterozygosity #####
setwd(paste0(accPath, "/stats"))

het_obs <- fread(paste0("obs_", pos, ".het"), header=T)
het_obs <- het_obs[ , c(1,6)]
colnames(het_obs) <- c("id", "Fobs")

het_pred <- fread(paste0("pred_", pos, ".het"), header=T)
het_pred <- het_pred[ , c(1,6)]
colnames(het_pred) <- c("id", "Fpred")

het <- inner_join(het_obs, het_pred, by="id")
cor_het <- round(cor(het$Fobs, het$Fpred, use = "complete.obs"), 3)
lm_het <- lm(het$Fpred ~ het$Fobs)
slope <- round(lm_het$coefficients[2], 3)

p6 <- ggplot(het, aes(x=Fobs, y=Fpred)) +
  geom_point(shape = 18, size = 5, alpha = 0.4, color = "orange3", show.legend = FALSE) +
  ggtitle(paste0("Imputed vs Observed Heterozygosity: chr", pos, "\ncorrelation = ", cor_het, "\nprediction bias (slope) = ", slope)) +
  theme_light() +
  theme(#panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5, face="bold"),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12),
    legend.position = "right") +
  xlab("Observed heterozygosity") + 
  ylab("Imputed heterozygosity")


### Save plots
setwd(paste0(accPath, "/plots"))
plots_page1 <- list(p1, p2, p3)
margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))
g1 <- arrangeGrob(grobs=lapply(plots_page1,"+", margin), nrow=2, 
                  layout_matrix = matrix(c(1,2,
                                           3,3), nrow=2, byrow=T))

ggsave(file=paste0("page1_", pos, ".pdf"), 
       g1, dpi=300, width = 16, height = 12)

ggsave(file=paste0("page2_", pos, ".pdf"),
       p4, dpi=300, width = 16, height = 12)

ggsave(file=paste0("page3_", pos, ".pdf"),
       p5, dpi=300, width = 16, height = 12)

ggsave(file=paste0("page4_", pos, ".pdf"),
       p6, dpi=300, width = 16, height = 12)

#---------- Summary Statistics ----------#
dir.create(paste0(accPath, "/summary"), recursive = T)

sink(paste0(accPath, "/summary/summaryStatsImputation_", pos, ".txt"))
cat("#----- General Concordance Rate -----#\n")
cat(accuracy[[count, 1]][[1]][[1]])
cat("\n")
cat("\n")
cat("#----- Concordance Rate by Subject -----#\n")
print(summary(accuracy[[count, 1]][[2]][ ,1]))
cat("\n")
cat("\n")
cat("#----- Concordance Rate by SNPs -----#\n")
print(summary(accuracy[[count, 1]][[3]][ ,4]))
cat("\n")
cat("\n")

cat("#----- General Correlation -----#\n")
cat(accuracy[[count, 2]][[1]][[1]])
cat("\n")
cat("\n")
cat("#----- Correlation by Subject -----#\n")
print(summary(accuracy[[count, 2]][[2]][ ,1]))
cat("\n")
cat("\n")
cat("#----- Correlation by SNPs -----#\n")
print(summary(accuracy[[count, 2]][[3]][ ,4]))
cat("\n")
cat("\n")

cat("#----- Minor Allele Frequency (MAF) -----#\n")
print(summary(accuracy[[count, 3]][, 5]))
sink()

#---------- Save RDS file ----------#
saveRDS(accuracy, file=paste0(accPath, "/accImputation.rds"))




