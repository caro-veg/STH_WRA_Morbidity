
library(ggplot2)
library(matrixStats)
library(cowplot)


################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING ICL EGG OUTPUT
################################################################################################################################ 
# MODERATE PREVALENCE, ANNUAL TREATMENT
#################################################################################################################################


# thresholds for medium-to-heavy infection
mediumHeavyInfectionThreshold_epg <- 2000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
mediumHeavyInfectionThreshold_counts <- mediumHeavyInfectionThreshold_epg / diagnosticDivisor


# thresholds for heavy infection
heavyInfectionThreshold_epg <- 4000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
heavyInfectionThreshold_counts <- heavyInfectionThreshold_epg / diagnosticDivisor



outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\FiguresAntonio\\"


# for column names and iterations
path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file <- "fertilisedFemaleWorms_old_20-50_treat1x.RData"
results <- get(load(paste0(path, file)))




#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat1x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat1x.RData"
stub1 <- substring(file1, first=23, last=40)
stub2 <- substring(file2, first=23, last=40)
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

###########################################################################################################################

#red.ratio.mhi.icl <- (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts
red.ratio.mhi.icl <- ifelse(df.msd.mhi.icl.old$meanMHICounts != 0, (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts, 0)

red.mean.mhi.icl.mod.1x.15_50 <- mean(red.ratio.mhi.icl[181:601], na.rm=TRUE)
red.sd.mhi.icl.mod.1x.15_50 <-  sd(red.ratio.mhi.icl[181:601], na.rm=TRUE)
red.pc.mhi.icl.mod.1x.15_50 <- quantile(red.ratio.mhi.icl[181:601], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.icl.mod.1x.15_19 <- mean(red.ratio.mhi.icl[181:240], na.rm=TRUE)
red.sd.mhi.icl.mod.1x.15_19 <-  sd(red.ratio.mhi.icl[181:240], na.rm=TRUE)
red.pc.mhi.icl.mod.1x.15_19 <- quantile(red.ratio.mhi.icl[181:240], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.icl.mod.1x.20_50 <- mean(red.ratio.mhi.icl[241:601], na.rm=TRUE)
red.sd.mhi.icl.mod.1x.20_50 <-  sd(red.ratio.mhi.icl[241:601], na.rm=TRUE)
red.pc.mhi.icl.mod.1x.20_50 <- quantile(red.ratio.mhi.icl[241:601], probs=c(0.05, 0.95), na.rm=TRUE)


df <- data.frame(Mean=c(red.mean.mhi.icl.mod.1x.15_50, red.mean.mhi.icl.mod.1x.15_19, red.mean.mhi.icl.mod.1x.20_50),
			pc5=c(red.pc.mhi.icl.mod.1x.15_50[1], red.pc.mhi.icl.mod.1x.15_19[1], red.pc.mhi.icl.mod.1x.20_50[1]),
			pc95=c(red.pc.mhi.icl.mod.1x.15_50[2], red.pc.mhi.icl.mod.1x.15_19[2], red.pc.mhi.icl.mod.1x.20_50[2]))
df <- df*100
df <- signif(df, digits=5)
stub <- substring(stub2, first=5, last=18)
#write.table(df, file=paste0(outpath, "relativeReduction_pc95MHI", stub, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)




################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG OUTPUT
################################################################################################################################ 
################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat1x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat1x.RData"
stub1 <- substring(file1, first=23, last=40)
stub2 <- substring(file2, first=23, last=40)
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



##################################################################################################################################


red.ratio.mhi.emc <- (df.msd.mhi.emc.old$meanMHICounts - df.msd.mhi.emc.new$meanMHICounts) / df.msd.mhi.emc.old$meanMHICounts

red.mean.mhi.emc.mod.1x.15_50 <- mean(red.ratio.mhi.emc[181:601], na.rm=TRUE)
red.sd.mhi.emc.mod.1x.15_50 <-  sd(red.ratio.mhi.emc[181:601], na.rm=TRUE)
red.pc.mhi.emc.mod.1x.15_50 <- quantile(red.ratio.mhi.emc[181:601], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.emc.mod.1x.15_19 <- mean(red.ratio.mhi.emc[181:240], na.rm=TRUE)
red.sd.mhi.emc.mod.1x.15_19 <-  sd(red.ratio.mhi.emc[181:240], na.rm=TRUE)
red.pc.mhi.emc.mod.1x.15_19 <- quantile(red.ratio.mhi.emc[181:240], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.emc.mod.1x.20_50 <- mean(red.ratio.mhi.emc[241:601], na.rm=TRUE)
red.sd.mhi.emc.mod.1x.20_50 <-  sd(red.ratio.mhi.emc[241:601], na.rm=TRUE)
red.pc.mhi.emc.mod.1x.20_50 <- quantile(red.ratio.mhi.emc[241:601], probs=c(0.05, 0.95), na.rm=TRUE)



df <- data.frame(Mean=c(red.mean.mhi.emc.mod.1x.15_50, red.mean.mhi.emc.mod.1x.15_19, red.mean.mhi.emc.mod.1x.20_50),
			pc5=c(red.pc.mhi.emc.mod.1x.15_50[1], red.pc.mhi.emc.mod.1x.15_19[1], red.pc.mhi.emc.mod.1x.20_50[1]),
			pc95=c(red.pc.mhi.emc.mod.1x.15_50[2], red.pc.mhi.emc.mod.1x.15_19[2], red.pc.mhi.emc.mod.1x.20_50[2]))
df <- df*100
df <- signif(df, digits=5)
stub <- substring(stub2, first=5, last=18)
#write.table(df, file=paste0(outpath, "relativeReduction_pc95MHI", stub, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING ICL EGG OUTPUT
################################################################################################################################ 
# MODERATE PREVALENCE, SEMI-ANNUAL TREATMENT
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


#############################################################################################

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


##################################################################################################


#red.ratio.mhi.icl <- (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts
red.ratio.mhi.icl <- ifelse(df.msd.mhi.icl.old$meanMHICounts != 0, (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts, 0)


red.mean.mhi.icl.mod.2x.15_50 <- mean(red.ratio.mhi.icl[181:601], na.rm=TRUE)
red.sd.mhi.icl.mod.2x.15_50 <-  sd(red.ratio.mhi.icl[181:601], na.rm=TRUE)
red.pc.mhi.icl.mod.2x.15_50 <- quantile(red.ratio.mhi.icl[181:601], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.icl.mod.2x.15_19 <- mean(red.ratio.mhi.icl[181:240], na.rm=TRUE)
red.sd.mhi.icl.mod.2x.15_19 <-  sd(red.ratio.mhi.icl[181:240], na.rm=TRUE)
red.pc.mhi.icl.mod.2x.15_19 <- quantile(red.ratio.mhi.icl[181:240], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.icl.mod.2x.20_50 <- mean(red.ratio.mhi.icl[241:601], na.rm=TRUE)
red.sd.mhi.icl.mod.2x.20_50 <-  sd(red.ratio.mhi.icl[241:601], na.rm=TRUE)
red.pc.mhi.icl.mod.2x.20_50 <- quantile(red.ratio.mhi.icl[241:601], probs=c(0.05, 0.95), na.rm=TRUE)



df <- data.frame(Mean=c(red.mean.mhi.icl.mod.2x.15_50, red.mean.mhi.icl.mod.2x.15_19, red.mean.mhi.icl.mod.2x.20_50),
			pc5=c(red.pc.mhi.icl.mod.2x.15_50[1], red.pc.mhi.icl.mod.2x.15_19[1], red.pc.mhi.icl.mod.2x.20_50[1]),
			pc95=c(red.pc.mhi.icl.mod.2x.15_50[2], red.pc.mhi.icl.mod.2x.15_19[2], red.pc.mhi.icl.mod.2x.20_50[2]))
df <- df*100
df <- signif(df, digits=5)
stub <- substring(stub2, first=5, last=18)
#write.table(df, file=paste0(outpath, "relativeReduction_pc95MHI", stub, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)




################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG OUTPUT
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



##################################################################################################################################


red.ratio.mhi.emc <- (df.msd.mhi.emc.old$meanMHICounts - df.msd.mhi.emc.new$meanMHICounts) / df.msd.mhi.emc.old$meanMHICounts

red.mean.mhi.emc.mod.2x.15_50 <- mean(red.ratio.mhi.emc[181:601], na.rm=TRUE)
red.sd.mhi.emc.mod.2x.15_50 <-  sd(red.ratio.mhi.emc[181:601], na.rm=TRUE)
red.pc.mhi.emc.mod.2x.15_50 <- quantile(red.ratio.mhi.emc[181:601], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.emc.mod.2x.15_19 <- mean(red.ratio.mhi.emc[181:240], na.rm=TRUE)
red.sd.mhi.emc.mod.2x.15_19 <-  sd(red.ratio.mhi.emc[181:240], na.rm=TRUE)
red.pc.mhi.emc.mod.2x.15_19 <- quantile(red.ratio.mhi.emc[181:240], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.emc.mod.2x.20_50 <- mean(red.ratio.mhi.emc[241:601], na.rm=TRUE)
red.sd.mhi.emc.mod.2x.20_50 <-  sd(red.ratio.mhi.emc[241:601], na.rm=TRUE)
red.pc.mhi.emc.mod.2x.20_50 <- quantile(red.ratio.mhi.emc[241:601], probs=c(0.05, 0.95), na.rm=TRUE)



df <- data.frame(Mean=c(red.mean.mhi.emc.mod.2x.15_50, red.mean.mhi.emc.mod.2x.15_19, red.mean.mhi.emc.mod.2x.20_50),
			pc5=c(red.pc.mhi.emc.mod.2x.15_50[1], red.pc.mhi.emc.mod.2x.15_19[1], red.pc.mhi.emc.mod.2x.20_50[1]),
			pc95=c(red.pc.mhi.emc.mod.2x.15_50[2], red.pc.mhi.emc.mod.2x.15_19[2], red.pc.mhi.emc.mod.2x.20_50[2]))
df <- df*100
df <- signif(df, digits=5)
stub <- substring(stub2, first=5, last=18)
#write.table(df, file=paste0(outpath, "relativeReduction_pc95MHI", stub, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)





################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING ICL EGG OUTPUT
################################################################################################################################ 
# HIGH PREVALENCE, SEMI-ANNUAL TREATMENT
#################################################################################################################################

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file1 <- "redoneEggCountsOneList_old_50-70_treat2x.RData"
file2 <- "redoneEggCountsOneList_new_50-70_treat2x.RData"
stub1 <- substring(file1, first=23, last=40)
stub2 <- substring(file2, first=23, last=40)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))


#############################################################################################

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


#############################################################################################################

red.ratio.mhi.icl <- (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts

red.mean.mhi.icl.high.2x.15_50 <- mean(red.ratio.mhi.icl[181:601], na.rm=TRUE)
red.sd.mhi.icl.high.2x.15_50 <-  sd(red.ratio.mhi.icl[181:601], na.rm=TRUE)
red.pc.mhi.icl.high.2x.15_50 <- quantile(red.ratio.mhi.icl[181:601], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.icl.high.2x.15_19 <- mean(red.ratio.mhi.icl[181:240], na.rm=TRUE)
red.sd.mhi.icl.high.2x.15_19 <-  sd(red.ratio.mhi.icl[181:240], na.rm=TRUE)
red.pc.mhi.icl.high.2x.15_19 <- quantile(red.ratio.mhi.icl[181:240], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.icl.high.2x.20_50 <- mean(red.ratio.mhi.icl[241:601], na.rm=TRUE)
red.sd.mhi.icl.high.2x.20_50 <-  sd(red.ratio.mhi.icl[241:601], na.rm=TRUE)
red.pc.mhi.icl.high.2x.20_50 <- quantile(red.ratio.mhi.icl[241:601], probs=c(0.05, 0.95), na.rm=TRUE)



df <- data.frame(Mean=c(red.mean.mhi.icl.high.2x.15_50, red.mean.mhi.icl.high.2x.15_19, red.mean.mhi.icl.high.2x.20_50),
			pc5=c(red.pc.mhi.icl.high.2x.15_50[1], red.pc.mhi.icl.high.2x.15_19[1], red.pc.mhi.icl.high.2x.20_50[1]),
			pc95=c(red.pc.mhi.icl.high.2x.15_50[2], red.pc.mhi.icl.high.2x.15_19[2], red.pc.mhi.icl.high.2x.20_50[2]))
df <- df*100
df <- signif(df, digits=5)
stub <- substring(stub2, first=5, last=18)
#write.table(df, file=paste0(outpath, "relativeReduction_pc95MHI", stub, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


##############################################################################################################################################
################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG OUTPUT
################################################################################################################################ 


#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "redoneEggCountsOneList_old_50-100_treat2x.RData"
file2 <- "redoneEggCountsOneList_new_50-100_treat2x.RData"
stub1 <- substring(file1, first=23, last=41)
stub2 <- substring(file2, first=23, last=41)
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



##################################################################################################################################


red.ratio.mhi.emc <- (df.msd.mhi.emc.old$meanMHICounts - df.msd.mhi.emc.new$meanMHICounts) / df.msd.mhi.emc.old$meanMHICounts

red.mean.mhi.emc.high.2x.15_50 <- mean(red.ratio.mhi.emc[181:601], na.rm=TRUE)
red.sd.mhi.emc.high.2x.15_50 <-  sd(red.ratio.mhi.emc[181:601], na.rm=TRUE)
red.pc.mhi.emc.high.2x.15_50 <- quantile(red.ratio.mhi.emc[181:601], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.emc.high.2x.15_19 <- mean(red.ratio.mhi.emc[181:240], na.rm=TRUE)
red.sd.mhi.emc.high.2x.15_19 <-  sd(red.ratio.mhi.emc[181:240], na.rm=TRUE)
red.pc.mhi.emc.high.2x.15_19 <- quantile(red.ratio.mhi.emc[181:240], probs=c(0.05, 0.95), na.rm=TRUE)

red.mean.mhi.emc.high.2x.20_50 <- mean(red.ratio.mhi.emc[241:601], na.rm=TRUE)
red.sd.mhi.emc.high.2x.20_50 <-  sd(red.ratio.mhi.emc[241:601], na.rm=TRUE)
red.pc.mhi.emc.high.2x.20_50 <- quantile(red.ratio.mhi.emc[241:601], probs=c(0.05, 0.95), na.rm=TRUE)

df <- data.frame(Mean=c(red.mean.mhi.emc.high.2x.15_50, red.mean.mhi.emc.high.2x.15_19, red.mean.mhi.emc.high.2x.20_50),
			pc5=c(red.pc.mhi.emc.high.2x.15_50[1], red.pc.mhi.emc.high.2x.15_19[1], red.pc.mhi.emc.high.2x.20_50[1]),
			pc95=c(red.pc.mhi.emc.high.2x.15_50[2], red.pc.mhi.emc.high.2x.15_19[2], red.pc.mhi.emc.high.2x.20_50[2]))
df <- df*100
df <- signif(df, digits=5)
stub <- substring(stub2, first=5, last=19)
#write.table(df, file=paste0(outpath, "relativeReduction_pc95MHI", stub, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


################################################################################

ratios <- c(red.mean.mhi.icl.mod.1x.15_50, red.mean.mhi.emc.mod.1x.15_50, red.mean.mhi.icl.mod.2x.15_50, red.mean.mhi.emc.mod.2x.15_50, red.mean.mhi.icl.high.2x.15_50, red.mean.mhi.emc.high.2x.15_50)*100
#sds <- c(red.sd.mhi.icl.mod.1x.15_50, red.sd.mhi.emc.mod.1x.15_50, red.sd.mhi.icl.mod.2x.15_50, red.sd.mhi.emc.mod.2x.15_50, red.sd.mhi.icl.high.2x.15_50, red.sd.mhi.emc.high.2x.15_50)*100
pc5 <- c(red.pc.mhi.icl.mod.1x.15_50[1], red.pc.mhi.emc.mod.1x.15_50[1], red.pc.mhi.icl.mod.2x.15_50[1], red.pc.mhi.emc.mod.2x.15_50[1], red.pc.mhi.icl.high.2x.15_50[1], red.pc.mhi.emc.high.2x.15_50[1])*100
pc95 <- c(red.pc.mhi.icl.mod.1x.15_50[2], red.pc.mhi.emc.mod.1x.15_50[2], red.pc.mhi.icl.mod.2x.15_50[2], red.pc.mhi.emc.mod.2x.15_50[2], red.pc.mhi.icl.high.2x.15_50[2], red.pc.mhi.emc.high.2x.15_50[2])*100

scenarios <- seq(1, 12, 2)
unis <- c(rep(c("ICL", "EMC"), 3))

df <- data.frame(ratios=ratios, scenarios=scenarios, pc5=pc5, pc95=pc95, unis=unis)

p <- ggplot(data=df, aes(x=scenarios, y=ratios, colour=unis))
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank(), legend.text=element_text(size=16), legend.position="bottom")
p <- p + geom_point(size=4)
p <- p + geom_errorbar(aes(ymin=pc5, ymax=pc95))
p <- p + scale_y_continuous(breaks=seq(0, 110, 20)) + coord_cartesian(ylim=c(0,110)) 
p <- p + xlab("Scenarios") + ylab("Relative reduction in \nM&HI infection prevalence (%)")
p <- p + scale_x_continuous(breaks=c(2, 6, 10), labels=c("Moderate prevalence,\nannual treatment", "Moderate prevalence,\nsemi-annual treatment", "High prevalence,\nsemi-annual treatment"))
print(p)
ggsave(paste0(outpath, "summaryFigureTreatedFoI_95percentiles.pdf"), p, dpi=300)




