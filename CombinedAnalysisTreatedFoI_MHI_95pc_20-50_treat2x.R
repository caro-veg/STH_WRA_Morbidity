
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



outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\Tables_OneList\\"



################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION FROM ICL EGG OUTPUT
################################################################################################################################ 

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat2x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat2x.RData"
stub1 <- substring(file1, first=23, last=40)
stub2 <- substring(file2, first=23, last=40)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))



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
file1 <- "redoneEggCountsOneList_old_20-50_treat2x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat2x.RData"
stub1 <- substring(file1, first=23, last=40)
stub2 <- substring(file2, first=23, last=40)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))




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






