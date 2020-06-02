
library(ggplot2)
library(matrixStats)
library(cowplot)


# thresholds for medium-to-heavy infection
mediumHeavyInfectionThreshold_epg <- 2000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
mediumHeavyInfectionThreshold_counts <- mediumHeavyInfectionThreshold_epg / diagnosticDivisor


# thresholds for heavy infection
heavyInfectionThreshold_epg <- 4000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
heavyInfectionThreshold_counts <- heavyInfectionThreshold_epg / diagnosticDivisor



outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\Figures\\"

# for column names and iterations
path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file <- "fertilisedFemaleWorms_old_20-50_treat1x.RData"
results <- get(load(paste0(path, file)))

################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION FROM ICL EGG OUTPUT
################################################################################################################################ 

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file1 <- "eggCounts_old_20-50_treat1x.RData"
file2 <- "eggCounts_new_20-50_treat1x.RData"
stub1 <- substring(file1, first=10, last=27)
stub2 <- substring(file2, first=10, last=27)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))


##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################


MHI <- eggCounts_old > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- unlist(MHICountsList)

years <- 100
timerep=rep(seq(0, years, 1/12), 100)
time=seq(0, years, 1/12)


temp <- do.call(cbind, MHICountsList)
meanMHICounts <- rowMeans(temp)
sdMHICounts <- sqrt(rowVars(temp))
medianMHICounts <- rowMedians(temp)
p95MHICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))


df.msd.mhi.icl.old <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, sdMHICounts=sdMHICounts/500*100)
df.mp.mhi.icl.old <- data.frame(time=time, medianMHICounts=medianMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)

df.mhi.icl.old <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)


##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- unlist(MHICountsList)

years <- 100
timerep=rep(seq(0, years, 1/12), 100)
time=seq(0, years, 1/12)


temp <- do.call(cbind, MHICountsList)
meanMHICounts <- rowMeans(temp)
sdMHICounts <- sqrt(rowVars(temp))
medianMHICounts <- rowMedians(temp)
p95MHICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))


df.msd.mhi.icl.new <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, sdMHICounts=sdMHICounts/500*100)
df.mp.mhi.icl.new <- data.frame(time=time, medianMHICounts=medianMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)

df.mhi.icl.new <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)


##############################################################################################################

p1 <- ggplot(df.mhi.icl.old)
p1 <- p1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=12))
p1 <- p1 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
p1 <- p1 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p1 <- p1 + scale_x_continuous(name="Age (years)", expand=c(0,0), limits=c(0, 72), breaks=seq(0, 100, 5))
p1 <- p1 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p1)
#ggsave(filename=paste0(outpath, "meanWomenMHI95pc", stub1, "_ICL_20", ".jpeg"), plot=p1, dpi=300)


p2 <- ggplot(df.mhi.icl.new)
p2 <- p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=12))
p2 <- p2 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
p2 <- p2 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p2 <- p2 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p2 <- p2 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p2)
#ggsave(filename=paste0(outpath, "meanWomenMHI95pc", stub2, "_ICL_20", ".jpeg"), plot=p2, dpi=300)



#####################################################################################
# determine number of women in cohort with medium-to-heavy infection over time
#####################################################################################

MHI <- eggCounts_old > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- do.call(cbind, MHICountsList)
MHIFrac <- MHICounts / 500


pp15_50 <- mean(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
pp15_19 <- mean(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
pp20_50 <- mean(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 

pc_95_15_50 <- quantile(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ]), probs=c(0.05, 0.95)) 
pc_95_15_19 <- quantile(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ]), probs=c(0.05, 0.95)) 
pc_95_20_50 <- quantile(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ]), probs=c(0.05, 0.95))  

df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), pc5=c(pc_95_15_50[1], pc_95_15_19[1], pc_95_20_50[1]), pc95=c(pc_95_15_50[2], pc_95_15_19[2], pc_95_20_50[2]))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "pc95MHI", stub1, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


#####################################################################################

MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- do.call(cbind, MHICountsList)
MHIFrac <- MHICounts / 500


pp15_50 <- mean(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
pp15_19 <- mean(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
pp20_50 <- mean(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 

pc_95_15_50 <- quantile(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ]), probs=c(0.05, 0.95)) 
pc_95_15_19 <- quantile(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ]), probs=c(0.05, 0.95)) 
pc_95_20_50 <- quantile(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ]), probs=c(0.05, 0.95)) 

df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), pc5=c(pc_95_15_50[1], pc_95_15_19[1], pc_95_20_50[1]), pc95=c(pc_95_15_50[2], pc_95_15_19[2], pc_95_20_50[2]))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "pc95MHI", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION FROM EMC EGG OUTPUT
################################################################################################################################ 

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "eggCounts_old_20-50_treat1x.RData"
file2 <- "eggCounts_new_20-50_treat1x.RData"
stub1 <- substring(file1, first=10, last=27)
stub2 <- substring(file2, first=10, last=27)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

MHI <- eggCounts_old > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- unlist(MHICountsList)

years <- 100
timerep=rep(seq(0, years, 1/12), 100)
time=seq(0, years, 1/12)


temp <- do.call(cbind, MHICountsList)
meanMHICounts <- rowMeans(temp)
sdMHICounts <- sqrt(rowVars(temp))
medianMHICounts <- rowMedians(temp)
p95MHICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))


df.msd.mhi.emc.old <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, sdMHICounts=sdMHICounts/500*100)
df.mp.mhi.emc.old <- data.frame(time=time, medianMHICounts=medianMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)

df.mhi.emc.old <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- unlist(MHICountsList)

years <- 100
timerep=rep(seq(0, years, 1/12), 100)
time=seq(0, years, 1/12)


temp <- do.call(cbind, MHICountsList)
meanMHICounts <- rowMeans(temp)
sdMHICounts <- sqrt(rowVars(temp))
medianMHICounts <- rowMedians(temp)
p95MHICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))


df.msd.mhi.emc.new <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, sdMHICounts=sdMHICounts/500*100)
df.mp.mhi.emc.new <- data.frame(time=time, medianMHICounts=medianMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)

df.mhi.emc.new <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)


##############################################################################################################

p3 <- ggplot(df.mhi.emc.old)
p3 <- p3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=12))
p3 <- p3 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
p3 <- p3 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p3 <- p3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p3 <- p3 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p3)
#ggsave(filename=paste0(outpath, "meanWomenMHI95pc", stub1, "_EMC_20", ".jpeg"), plot=p3, dpi=300)


p4 <- ggplot(df.mhi.emc.new)
p4 <- p4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=12))
p4 <- p4 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
p4 <- p4 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p4 <- p4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p4 <- p4 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p4)
#ggsave(filename=paste0(outpath, "meanWomenMHI95pc", stub2, "_EMC_20", ".jpeg"), plot=p4, dpi=300)


#####################################################################################
# determine number of women in cohort with medium-to-heavy infection over time
#####################################################################################

MHI <- eggCounts_old > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- do.call(cbind, MHICountsList)
MHIFrac <- MHICounts / 500



pp15_50 <- mean(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
pp15_19 <- mean(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
pp20_50 <- mean(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 

pc_95_15_50 <- quantile(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ]), probs=c(0.05, 0.95)) 
pc_95_15_19 <- quantile(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ]), probs=c(0.05, 0.95)) 
pc_95_20_50 <- quantile(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ]), probs=c(0.05, 0.95))


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), pc5=c(pc_95_15_50[1], pc_95_15_19[1], pc_95_20_50[1]), pc95=c(pc_95_15_50[2], pc_95_15_19[2], pc_95_20_50[2]))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "pc95MHI", stub1, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

#####################################################################################

MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHICounts <- do.call(cbind, MHICountsList)
MHIFrac <- MHICounts / 500


pp15_50 <- mean(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
pp15_19 <- mean(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
pp20_50 <- mean(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ]))

pc_95_15_50 <- quantile(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ]), probs=c(0.05, 0.95)) 
pc_95_15_19 <- quantile(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ]), probs=c(0.05, 0.95)) 
pc_95_20_50 <- quantile(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ]), probs=c(0.05, 0.95)) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), pc5=c(pc_95_15_50[1], pc_95_15_19[1], pc_95_20_50[1]), pc95=c(pc_95_15_50[2], pc_95_15_19[2], pc_95_20_50[2]))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "pc95MHI", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


#####################################################################################
# COMPOSITE FIGURE
#####################################################################################

stub <- substring(stub1, first=5, last=18)

pic_mhi_means <- plot_grid(p1, p2, p3, p4, labels=c('a', 'b', 'c', 'd'), label_size=18, ncol=2, hjust=0, vjust=1.0)
save_plot(file=paste0(outpath, "pc95MHI", stub, ".jpeg"), pic_mhi_means, base_height=8.5, base_width=9.3)




