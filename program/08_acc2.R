#---------- INITIAL PARAMETERS ----------#
my_packages <- c("data.table", "ggplot2", "ggrepel", "gridExtra", "vcfR", "wesanderson") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

library(data.table)
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

accPath=accPath
referenceGenome=referenceGenome
genoPanel=genoPanel
pos=pos
count=as.numeric(count)
numRow=as.numeric(numRow)
numCol=as.numeric(numCol)

#---------- Load accImputation.rds ----------#
print(paste0("Chromosome: ", pos))
print("loading accImputation.rds")
accuracy <- readRDS(file = paste0(accPath, "/accImputation.rds"))

#---------- R2 and MAF  ----------#
setwd(paste0(accPath, "/r2_maf"))
r2_maf <- read.table(paste0("r2_maf_", pos, ".txt"), header = T, sep = "\t")

#---------- Imputed genotypes ----------#
setwd(paste0(accPath, "/012"))
print("loading individuals and position of SNPs")
ind <- read.table(paste0("sampleNames_", pos, ".txt"), header=F)
ind <- ind[,1]
ind <- t(ind)
map <- read.table(paste0("map_", pos, ".txt"), header=F)

print("loading predicted genotypes (012)")
geno <- matrix(scan(paste0("geno012_", pos, ".raw"),what = "character"), nrow = (numRow + 1), ncol = (numCol + 6), byrow=T)
geno <- geno[-c(1), -c(1:6)]
geno <- t(geno)
dim(geno)
geno <- as.data.frame(geno)
colnames(geno) <- ind
geno <- geno[ , order(names(geno))]

#---------- R2 and MAF ----------#
accuracy[[count, 1]] <- r2_maf
#accuracy[[count, 1]][1:5,] #R2 and MAF from beagle or minimac3


#---------- PLOTS WITH GGPLOT2 ----------#
print("Creating plots")
dir.create(paste0(accPath, "/plots"), recursive = T)


### Scatterplot R2 vs MAF
r2_mafDF <- data.frame(r2 = accuracy[[count, 1]][ ,4],
                       maf = accuracy[[count, 1]][ ,5])

r2_mafDF$mafBin <- cut(r2_mafDF$maf, 
                       seq(0, 0.5, 0.01), 
                       include.lowest=TRUE, right=T,
                       labels=seq(0.01, 0.5, 0.01)
)

r2_mafDF_summary <- r2_mafDF %>% group_by(mafBin) %>% summarise(r2=median(r2)) %>% as.data.frame()
r2_mafDF_summary <- na.omit(r2_mafDF_summary)
r2_mafDF_summary$mafBin <- as.numeric.factor(r2_mafDF_summary$mafBin)

lm1 <- lm(r2 ~ poly(mafBin, 3), data = r2_mafDF_summary)
r2_mafDF_summary$pred <- predict(lm1, data=seq(0.01, 0.5, 0.01))

p1 <- ggplot(r2_mafDF_summary, aes(x = mafBin)) +
  geom_point(aes(y = pred), colour="darkred" , label="R2", shape=18, size=2.5, alpha=1) +
  stat_smooth(aes(x = mafBin, y = r2), 
              method = "lm", formula = y ~ poly(x, 3), se = FALSE, 
              colour="darkred", label="R2", size=0.8, alpha=0.6) +
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
        legend.position = "none",
        legend.text = element_text(size=12, face="bold")) +
  xlab("Minor allele frequency (MAF)") + 
  ylab("Correlation") +
  labs(color  = "", linetype = "", shape = "")


##### R2 density plot #####
r2_mafDF <- accuracy[[count, 1]]
# Divide chromosome in segments length = window #####
window=100000
commandAWK <- paste0('awk \'$1 == ', pos, ' {print $2}\' ', referenceGenome, '.fai')
chromLenght <- as.numeric(system(commandAWK, intern = T))

seqBreaks <- seq(0, chromLenght, by=window)
seqBreaks[which.max(seqBreaks)] <- chromLenght
r2_mafDF$group <- cut(r2_mafDF$pos,
                           breaks = seqBreaks,
                           include.lowest=TRUE,
                           labels = seqBreaks[-1])

# Calculate average Concordance Rate by SNP by Chromosome segment
r2_average <- aggregate(r2_mafDF$r2,
                             list(r2_mafDF$group), 
                             mean)
colnames(r2_average) <- c("group", "average")
r2_average$height <- 1

# Convert group from factor to numeric
#https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
r2_average$group <- as.numeric.factor(r2_average$group)
r2_average$group <- r2_average$group/1000000
min(r2_average$average)
max(r2_average$average)

library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
p2 <- ggplot(r2_average, aes(x = group , y = height, fill = average)) +
        geom_line(size=0.01, alpha=0.2) +
        geom_tile() + 
        ggtitle(paste0("Average R2 (Window of ", window/1000, "kb): Chromosome ", pos)) +
        scale_fill_gradientn(colours = pal, limits=c(0,1)) +
        scale_x_discrete(name = "Position (Mb)", expand = c(0, 0),
                         limits = seq(0, round(max(r2_average$group), 0), 5)) +
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
# PCA gcta
eigenvalues <- fread(paste0("genoGRM_", pos, ".eigenval"))
eigenvalues <- as.data.frame(eigenvalues)
varPCA <- round( (eigenvalues$V1^2) / sum(eigenvalues$V1^2) *100, 3)

eigenvectors <- fread(paste0("genoGRM_", pos, ".eigenvec"))
eigenvectors <- as.data.frame(eigenvectors)
eigenvectors <- eigenvec[,c(3:4, 2)]
colnames(eigenvectors) <- c("pc1", "pc2", "id")

# plot PC1 vs PC2 using scores
pc1=1
pc2=2

p3 <- ggplot(eigenvectors, aes_string(x=colnames(eigenvectors)[pc1], y=colnames(eigenvectors)[pc2])) + 
        geom_point(shape=18, size=5) +
        scale_color_manual(values=c("#1F77B4")) +
        geom_hline(yintercept=0, linetype="dashed", color = "gray40", size=1) +
        geom_vline(xintercept=0, linetype="dashed", color = "gray40", size=1) +
        geom_text_repel(aes(label=id), size=4, segment.size=0.10, nudge_x=0.2, direction="y") +
        ggtitle(paste0("Singular value decomposition of GRM: chr", pos)) +
        theme_light() +
        theme(panel.border = element_blank(),
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


### Save plots
setwd(paste0(accPath, "/plots"))
plots_page1 <- list(p1, p2)
margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))
g1 <- arrangeGrob(grobs=lapply(plots_page1,"+", margin), nrow=2, 
                  layout_matrix = matrix(c(1,1,
                                           2,2), nrow=2, byrow=T))

ggsave(file=paste0("page1_", pos, ".pdf"), 
       g1, dpi=300, width = 16, height = 12)

ggsave(file=paste0("page2_", pos, ".pdf"), 
       p3, dpi=300, width = 16, height = 12)


#---------- Summary Statistics ----------#
dir.create(paste0(accPath, "/summary"), recursive = T)

sink(paste0(accPath, "/summary/summaryStats_", pos, ".txt"))
cat("#----- R2 -----#\n")
print(summary(accuracy[[count, 2]][, 4]))
cat("#----- Minor Allele Frequency (MAF) -----#\n")
print(summary(accuracy[[count, 2]][, 5]))
sink()

#---------- Save RDS file ----------#
saveRDS(accuracy, file=paste0(accPath, "/accImputation.rds"))




