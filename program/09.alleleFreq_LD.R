#---------- INITIAL PARAMETERS ----------#
my_packages <- c("data.table", "ggplot2", "gridExtra") # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed, repos = "https://cran.microsoft.com/", dependencies=TRUE)

library(data.table)
library(ggplot2)
library(gridExtra)

accPath=accPath
pos=pos

#---------- ALLELE FREQUENCY ----------#
### Read files
print("Reading files: allele frequency")
setwd(paste0(accPath, "/stats"))

freqPred <- fread(paste0("pred_",pos, ".afreq"), header=T, select=c(5), sep = "\t")
freqPred$pop <- "imputation"

### Concatenate files
freq <- as.data.frame(freqPred)
freq[,1] <- as.numeric(unlist(freq[,1]))
colnames(freq) <- c("alt", "pop")
freq <- na.omit(freq)
freq$maf <- ifelse(freq[,1] <= 0.5, freq[,1], 1 - freq[,1])
freq <- freq[ ,c(3,1,2)]

### Plots
print("Creating plots: allele frequency")
dir.create(paste0(accPath, "/plots"), recursive = T)

p1 <- ggplot(freq, aes(x=maf, fill=pop, color=pop)) + 
  geom_histogram(position="identity", alpha=0.3, na.rm = T) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  ggtitle(paste0("Minor allele frequency: chr", pos)) +
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
  xlab("MAF") + 
  ylab("Density")

p2 <- ggplot(freq, aes(x=alt, fill=pop, color=pop)) + 
  geom_histogram(position="identity", alpha=0.3, na.rm = T) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  ggtitle(paste0("ALT allele frequency: chr", pos)) +
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
  xlab("ALT allele frequency") + 
  ylab("Density")


#---------- LD DECAY: LDkit ----------#
### Read files
print("Reading files: linkage disequilibrium")
setwd(paste0(accPath, "/stats"))

# ldRef <- fread(paste0("ref_",pos, ".LDkit"), header=F, skip = 1, select=c(2,3), sep = "\t")
# ldRef$pop <- "reference"
ldPred <- fread(paste0("pred_",pos, ".LDkit"), header=F, skip = 1, select=c(2,3), sep = "\t")
ldPred$pop <- "imputation"

### Concatenate files
ld <- rbind(ldPred)
ld <- as.data.frame(ld)
colnames(ld) <-c("Position", "Rsquared", "Pop")
ld[,1] <- as.numeric(ld[,1])
ld[,2] <- as.numeric(ld[,2])
index <- seq(0, 1000000, 1000)
index[1] <- 1
ld <- ld[which(ld$Position %in% index), ]
ld[,1] <- ld[,1]/1000 # change bp to kb

### Plots
print("Creating plots: linkage disequilibrium")
p3a <- ggplot(ld, aes(x=Position, y=Rsquared, group=Pop, shape=Pop, colour=Pop)) +
  geom_line() +
  geom_point(size = 2.5, alpha = 0.5) +
  scale_x_continuous(limits=c(0,1000), breaks = seq(0, 1000, by = 50)) +
  scale_y_continuous(limits=c(0,1)) +
  ggtitle(paste0("Linkage Disequilibrium Decay: chr", pos)) +
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
  xlab("Distance (kb)") +
  ylab(expression(LD~(r^{2})))


p3b <- ggplot(ld, aes(x=Position, y=Rsquared, group=Pop, shape=Pop, colour=Pop)) +
  geom_line() +
  geom_point(size = 2.5, alpha = 0.5) +
  scale_x_continuous(limits=c(0,50), breaks = seq(0, 50, by = 10)) +
  scale_y_continuous(limits=c(0,1)) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        legend.position = "none") +
  xlab("Distance (kb)") +
  ylab(expression(LD~(r^{2})))


p3 <- p3a + annotation_custom(ggplotGrob(p3b), xmin = 600, xmax = 1000,
                              ymin = 0.5, ymax = 1)


# #---------- LD DECAY: PLINK ----------#
# ### Read files
# print("Reading files: linkage disequilibrium")
# setwd(paste0(accPath, "/stats"))
# 
# ldRef <- fread(paste0("ref_",pos, ".ld"), header=T, select=c(2,5,7), sep = " ")
# ldRef$pop <- "reference"
# ldPred <- fread(paste0("pred_",pos, ".ld"), header=T, select=c(2,5,7), sep = " ")
# ldPred$pop <- "imputation"
# 
# ### Concatenate files
# ld <- rbind(ldRef, ldPred)
# ld <- as.data.frame(ld)
# 
# # dataFrame: LD along chromosome 
# ld_along_chrom <- ld[ , c(1,3,4)]
# colnames(ld_along_chrom) <- c("Position", "Rsquared", "Pop")
# ld_along_chrom <- ld_along_chrom[order(ld_along_chrom$Position), ]
# 
# # dataFrame: LD between distances
# ld$dist <- ld$BP_B - ld$BP_A # difference between bp
# #ld$dist.kb <- ld$dist/1000 # distance in kb
# ld <- ld[order(ld$dist), ] # sort by distance between bp
# ld <- ld[ ,c(5,3,4)]
# colnames(ld) <- c("Distance", "Rsquared", "Pop")
# 
# ### Create bins of distance
# # https://stackoverflow.com/questions/9946630/colour-points-in-a-plot-differently-depending-on-a-vector-of-values
# # https://rstudio-pubs-static.s3.amazonaws.com/297778_5fce298898d64c81a4127cf811a9d486.html
# ldList <- list()
# population <- c("reference", "imputation")
# 
# for (i in 1:3){
#   bins=seq(1,max(ld[which(ld$Pop == population[i]), "Distance"]), by=100)
#   my.means=rep(0, length(bins))
#   LD.averages=data.frame(bins, my.means, pop=population[i])
#   
#   for (j in 1:length(bins)) {
#     data.interval=subset(ld, (ld[which(ld$Pop == population[i]), "Distance"] >= bins[j] & ld[which(ld$Pop == population[i]), "Distance"] < (bins[j]+500))) 
#     LD.averages$my.means[j]=mean(data.interval$Rsquared) 
#   }
#   
#   LD.averages$bins.kb <- LD.averages$bins/1000
#   ldList[[i]] <- LD.averages
# }
# 
# LD.averages <- do.call(rbind, ldList)
# LD.averages <- LD.averages[ , c(1,2,4,3)]
# ld$Distance.kb <- ld$Distance/1000
# 
# ### Plots
# print("Creating plots: linkage disequilibrium")
# p3 <- ggplot(LD.averages, aes(x=bins.kb, y=my.means, group=pop, shape=pop, colour=pop)) + 
#   geom_line() +
#   geom_point(size = 2.5, alpha = 0.5) +
#   ggtitle(paste0("Linkage Disequilibrium Decay: chr", pos)) +
#   theme_light() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.title.x = element_text(size=12, face="bold"),
#         axis.title.y = element_text(size=12, face="bold"),
#         plot.title = element_text(size=12, hjust = 0.5, face="bold"),
#         legend.position = "right") +
#   xlab("Distance (kb)") + 
#   ylab(expression(LD~(r^{2})))


##### ACTIVATE THIS PART ONLY IF YOU WANT TO ZOOM ANY REGION
# # Plot LD along chromosome (select a region to zoom)
# startbp=100000
# endbp=1000000
# 
# subset_ld_along_chrom <- ld_along_chrom[which(ld_along_chrom$Position >= startbp & ld_along_chrom$Position <= endbp), ]
# ldList <- list()
# population <- c("reference", "imputation")
# 
# for (i in 1:3){
#   average <- aggregate(subset_ld_along_chrom[which(subset_ld_along_chrom$Pop == population[i]), 2], 
#                        by=list(subset_ld_along_chrom[which(subset_ld_along_chrom$Pop == population[i]), 1]), 
#                        mean)
#   average <- data.frame(Position=average[,1],
#                         Rsquared=average[,2],
#                         pop=population[i])
#   ldList[[i]] <- average
# }
# 
# avg_ld_by_position <- do.call(rbind, ldList)
# 
# ggplot(avg_ld_by_position, aes(x=Position, y=Rsquared, group=pop, shape=pop, colour=pop)) + 
#   geom_point(size = 2.5, alpha = 0) + 
#   geom_line() +
#   ggtitle(paste0("Linkage Disequilibrium along chromosome ", pos)) +
#   theme_light() +
#   theme(panel.grid.major.x = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.title.x = element_text(size=12, face="bold"),
#         axis.title.y = element_text(size=12, face="bold"),
#         plot.title = element_text(size=12, hjust = 0.5, face="bold"),
#         legend.position = "right") +
#   xlab("Distance (bp)") + 
#   ylab(expression(LD~(r^{2})))


### Save plots
setwd(paste0(accPath, "/plots"))
plots_page1 <- list(p1, p2, p3)
margin = theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))
g1 <- arrangeGrob(grobs=lapply(plots_page1,"+", margin), nrow=2, 
                  layout_matrix = matrix(c(1,2,
                                           3,3), nrow=2, byrow=T))

ggsave(file=paste0("page3_", pos, ".pdf"), 
       g1, dpi=300, width = 16, height = 12)

