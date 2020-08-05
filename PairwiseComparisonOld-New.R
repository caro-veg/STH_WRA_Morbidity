library(matrixStats)

########################################################################################################################################
# WILCOXON RANK SIGN TEST FOR PAIRWISE COMPARISON OF MHI INFECTION PREVALENCE WITH OLD (SB ONLY) VS NEW (SB + WRA) TREATMENT STRATEGIES
################################################################################################################################ #######
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



outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\Tables_OneList\\"


# for column names and iterations
path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file <- "fertilisedFemaleWorms_old_20-50_treat1x.RData"
results <- get(load(paste0(path, file)))


################################################################################################
# FUNCTION TO CALCULATE AVERAGE POINT PREVALENCE FOR EACH AGE GROUP AND ITERATION
################################################################################################

pointPrevAgeClass <- function(x, i1, i2)
{
	return(mean(x[i1:i2]))
}

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat1x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat1x.RData"
stub1 <- substring(file1, first=27, last=40)
stub2 <- substring(file2, first=27, last=40)
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
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)


################################################################################################
# CALCULATE RELATIVE REDUCTION (MEAN AND 95% CREDIBLE INTERVAL)
################################################################################################

rr_15_50_20_50treat1x_ICL = (meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old * 100
rr_15_19_20_50treat1x_ICL = (meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old * 100
rr_20_50_20_50treat1x_ICL = (meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old * 100

mean_rr_15_50_20_50treat1x_ICL = mean(rr_15_50_20_50treat1x_ICL)
mean_rr_15_19_20_50treat1x_ICL = mean(rr_15_19_20_50treat1x_ICL)
mean_rr_20_50_20_50treat1x_ICL = mean(rr_20_50_20_50treat1x_ICL)

pc95_rr_15_50_20_50treat1x_ICL = quantile(rr_15_50_20_50treat1x_ICL, c(0.025, 0.975))
pc95_rr_15_19_20_50treat1x_ICL = quantile(rr_15_19_20_50treat1x_ICL, c(0.025, 0.975))
pc95_rr_20_50_20_50treat1x_ICL = quantile(rr_20_50_20_50treat1x_ICL, c(0.025, 0.975))

df <- data.frame(mean=c(mean_rr_15_50_20_50treat1x_ICL, mean_rr_15_19_20_50treat1x_ICL, mean_rr_20_50_20_50treat1x_ICL),
			pc0.025=c(pc95_rr_15_50_20_50treat1x_ICL[[1]], pc95_rr_15_19_20_50treat1x_ICL[[1]], pc95_rr_20_50_20_50treat1x_ICL[[1]]),
			pc0.975=c(pc95_rr_15_50_20_50treat1x_ICL[[2]], pc95_rr_15_19_20_50treat1x_ICL[[2]], pc95_rr_20_50_20_50treat1x_ICL[[2]]))


write.table(df, file=paste0(outpath, "rrMHI", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



################################################################################################
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")





################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG OUTPUT
################################################################################################################################ 
################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat1x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat1x.RData"
stub1 <- substring(file1, first=27, last=40)
stub2 <- substring(file2, first=27, last=40)
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
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)


################################################################################################
# CALCULATE RELATIVE REDUCTION (MEAN AND 95% CREDIBLE INTERVAL)
################################################################################################

rr_15_50_20_50treat1x_EMC = (meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old * 100
rr_15_19_20_50treat1x_EMC = (meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old * 100
rr_20_50_20_50treat1x_EMC = (meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old * 100

mean_rr_15_50_20_50treat1x_EMC = mean(rr_15_50_20_50treat1x_EMC)
mean_rr_15_19_20_50treat1x_EMC = mean(rr_15_19_20_50treat1x_EMC)
mean_rr_20_50_20_50treat1x_EMC = mean(rr_20_50_20_50treat1x_EMC)

pc95_rr_15_50_20_50treat1x_EMC = quantile(rr_15_50_20_50treat1x_EMC, c(0.025, 0.975))
pc95_rr_15_19_20_50treat1x_EMC = quantile(rr_15_19_20_50treat1x_EMC, c(0.025, 0.975))
pc95_rr_20_50_20_50treat1x_EMC = quantile(rr_20_50_20_50treat1x_EMC, c(0.025, 0.975))

df <- data.frame(mean=c(mean_rr_15_50_20_50treat1x_EMC, mean_rr_15_19_20_50treat1x_EMC, mean_rr_20_50_20_50treat1x_EMC),
			pc0.025=c(pc95_rr_15_50_20_50treat1x_EMC[[1]], pc95_rr_15_19_20_50treat1x_EMC[[1]], pc95_rr_20_50_20_50treat1x_EMC[[1]]),
			pc0.975=c(pc95_rr_15_50_20_50treat1x_EMC[[2]], pc95_rr_15_19_20_50treat1x_EMC[[2]], pc95_rr_20_50_20_50treat1x_EMC[[2]]))


write.table(df, file=paste0(outpath, "rrMHI", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



################################################################################################
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")



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
stub1 <- substring(file1, first=27, last=40)
stub2 <- substring(file2, first=27, last=40)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))


#############################################################################################

MHI <- eggCounts_old > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)


################################################################################################
# CALCULATE RELATIVE REDUCTION (MEAN AND 95% CREDIBLE INTERVAL)
################################################################################################

rr_15_50_20_50treat2x_ICL = ifelse(meanPointPrevs_15_50_old==0, 0, (meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old * 100)
rr_15_19_20_50treat2x_ICL = ifelse(meanPointPrevs_15_19_old==0, 0, (meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old * 100)
rr_20_50_20_50treat2x_ICL = ifelse(meanPointPrevs_20_50_old==0, 0, (meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old * 100)

mean_rr_15_50_20_50treat2x_ICL = mean(rr_15_50_20_50treat2x_ICL, na.rm=TRUE)
mean_rr_15_19_20_50treat2x_ICL = mean(rr_15_19_20_50treat2x_ICL, na.rm=TRUE)
mean_rr_20_50_20_50treat2x_ICL = mean(rr_20_50_20_50treat2x_ICL, na.rm=TRUE)

pc95_rr_15_50_20_50treat2x_ICL = quantile(rr_15_50_20_50treat2x_ICL, c(0.025, 0.975), na.rm=TRUE)
pc95_rr_15_19_20_50treat2x_ICL = quantile(rr_15_19_20_50treat2x_ICL, c(0.025, 0.975), na.rm=TRUE)
pc95_rr_20_50_20_50treat2x_ICL = quantile(rr_20_50_20_50treat2x_ICL, c(0.025, 0.975), na.rm=TRUE)

df <- data.frame(mean=c(mean_rr_15_50_20_50treat2x_ICL, mean_rr_15_19_20_50treat2x_ICL, mean_rr_20_50_20_50treat2x_ICL),
			pc0.025=c(pc95_rr_15_50_20_50treat2x_ICL[[1]], pc95_rr_15_19_20_50treat2x_ICL[[1]], pc95_rr_20_50_20_50treat2x_ICL[[1]]),
			pc0.975=c(pc95_rr_15_50_20_50treat2x_ICL[[2]], pc95_rr_15_19_20_50treat2x_ICL[[2]], pc95_rr_20_50_20_50treat2x_ICL[[2]]))


write.table(df, file=paste0(outpath, "rrMHI", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)




################################################################################################
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")



################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG OUTPUT
################################################################################################################################ 
#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat2x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat2x.RData"
stub1 <- substring(file1, first=27, last=40)
stub2 <- substring(file2, first=27, last=40)
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
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)


##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)

mean((meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old)
quantile((meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old, c(0.025, 0.975))

mean((meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old)
quantile((meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old, c(0.025, 0.975))

mean((meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old)
quantile((meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old, c(0.025, 0.975))
    

################################################################################################
# CALCULATE RELATIVE REDUCTION (MEAN AND 95% CREDIBLE INTERVAL)
################################################################################################

rr_15_50_20_50treat2x_EMC = (meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old * 100
rr_15_19_20_50treat2x_EMC = (meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old * 100
rr_20_50_20_50treat2x_EMC = (meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old * 100

mean_rr_15_50_20_50treat2x_EMC = mean(rr_15_50_20_50treat2x_EMC)
mean_rr_15_19_20_50treat2x_EMC = mean(rr_15_19_20_50treat2x_EMC)
mean_rr_20_50_20_50treat2x_EMC = mean(rr_20_50_20_50treat2x_EMC)

pc95_rr_15_50_20_50treat2x_EMC = quantile(rr_15_50_20_50treat2x_EMC, c(0.025, 0.975))
pc95_rr_15_19_20_50treat2x_EMC = quantile(rr_15_19_20_50treat2x_EMC, c(0.025, 0.975))
pc95_rr_20_50_20_50treat2x_EMC = quantile(rr_20_50_20_50treat2x_EMC, c(0.025, 0.975))

df <- data.frame(mean=c(mean_rr_15_50_20_50treat2x_EMC, mean_rr_15_19_20_50treat2x_EMC, mean_rr_20_50_20_50treat2x_EMC),
			pc0.025=c(pc95_rr_15_50_20_50treat2x_EMC[[1]], pc95_rr_15_19_20_50treat2x_EMC[[1]], pc95_rr_20_50_20_50treat2x_EMC[[1]]),
			pc0.975=c(pc95_rr_15_50_20_50treat2x_EMC[[2]], pc95_rr_15_19_20_50treat2x_EMC[[2]], pc95_rr_20_50_20_50treat2x_EMC[[2]]))


write.table(df, file=paste0(outpath, "rrMHI", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

   

################################################################################################
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")



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
stub1 <- substring(file1, first=27, last=40)
stub2 <- substring(file2, first=27, last=40)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))


#############################################################################################

MHI <- eggCounts_old > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)


################################################################################################
# CALCULATE RELATIVE REDUCTION (MEAN AND 95% CREDIBLE INTERVAL)
################################################################################################

rr_15_50_50_100treat2x_ICL = (meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old * 100
rr_15_19_50_100treat2x_ICL = (meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old * 100
rr_20_50_50_100treat2x_ICL = (meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old * 100

mean_rr_15_50_50_100treat2x_ICL = mean(rr_15_50_50_100treat2x_ICL)
mean_rr_15_19_50_100treat2x_ICL = mean(rr_15_19_50_100treat2x_ICL)
mean_rr_20_50_50_100treat2x_ICL = mean(rr_20_50_50_100treat2x_ICL)

pc95_rr_15_50_50_100treat2x_ICL = quantile(rr_15_50_50_100treat2x_ICL, c(0.025, 0.975))
pc95_rr_15_19_50_100treat2x_ICL = quantile(rr_15_19_50_100treat2x_ICL, c(0.025, 0.975))
pc95_rr_20_50_50_100treat2x_ICL = quantile(rr_20_50_50_100treat2x_ICL, c(0.025, 0.975))

df <- data.frame(mean=c(mean_rr_15_50_50_100treat2x_ICL, mean_rr_15_19_50_100treat2x_ICL, mean_rr_20_50_50_100treat2x_ICL),
			pc0.025=c(pc95_rr_15_50_50_100treat2x_ICL[[1]], pc95_rr_15_19_50_100treat2x_ICL[[1]], pc95_rr_20_50_50_100treat2x_ICL[[1]]),
			pc0.975=c(pc95_rr_15_50_50_100treat2x_ICL[[2]], pc95_rr_15_19_50_100treat2x_ICL[[2]], pc95_rr_20_50_50_100treat2x_ICL[[2]]))


write.table(df, file=paste0(outpath, "rrMHI", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)




################################################################################################
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")


################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG OUTPUT
################################################################################################################################ 


#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "redoneEggCountsOneList_old_50-100_treat2x.RData"
file2 <- "redoneEggCountsOneList_new_50-100_treat2x.RData"
stub1 <- substring(file1, first=27, last=41)
stub2 <- substring(file2, first=27, last=41)
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
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_old <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_old <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- eggCounts_new > mediumHeavyInfectionThreshold_counts
MHI <- as.data.frame(cbind(results[, 1], MHI))
names(MHI) <- colnames(results)
MHIList <- split(MHI, MHI$rep)
MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

MHICountsList <- lapply(MHIList, rowSums)
MHIPrevsList <- lapply(MHICountsList, function(x){x / 500 * 100})


meanPointPrevs_15_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 601)
meanPointPrevs_15_19_new <- sapply(MHIPrevsList, pointPrevAgeClass, 181, 240)
meanPointPrevs_20_50_new <- sapply(MHIPrevsList, pointPrevAgeClass, 241, 601)



################################################################################################
# CALCULATE RELATIVE REDUCTION (MEAN AND 95% CREDIBLE INTERVAL)
################################################################################################

rr_15_50_50_100treat2x_EMC = (meanPointPrevs_15_50_old - meanPointPrevs_15_50_new) / meanPointPrevs_15_50_old * 100
rr_15_19_50_100treat2x_EMC = (meanPointPrevs_15_19_old - meanPointPrevs_15_19_new) / meanPointPrevs_15_19_old * 100
rr_20_50_50_100treat2x_EMC = (meanPointPrevs_20_50_old - meanPointPrevs_20_50_new) / meanPointPrevs_20_50_old * 100

mean_rr_15_50_50_100treat2x_EMC = mean(rr_15_50_50_100treat2x_EMC)
mean_rr_15_19_50_100treat2x_EMC = mean(rr_15_19_50_100treat2x_EMC)
mean_rr_20_50_50_100treat2x_EMC = mean(rr_20_50_50_100treat2x_EMC)

pc95_rr_15_50_50_100treat2x_EMC = quantile(rr_15_50_50_100treat2x_EMC, c(0.025, 0.975))
pc95_rr_15_19_50_100treat2x_EMC = quantile(rr_15_19_50_100treat2x_EMC, c(0.025, 0.975))
pc95_rr_20_50_50_100treat2x_EMC = quantile(rr_20_50_50_100treat2x_EMC, c(0.025, 0.975))

df <- data.frame(mean=c(mean_rr_15_50_50_100treat2x_EMC, mean_rr_15_19_50_100treat2x_EMC, mean_rr_20_50_50_100treat2x_EMC),
			pc0.025=c(pc95_rr_15_50_50_100treat2x_EMC[[1]], pc95_rr_15_19_50_100treat2x_EMC[[1]], pc95_rr_20_50_50_100treat2x_EMC[[1]]),
			pc0.975=c(pc95_rr_15_50_50_100treat2x_EMC[[2]], pc95_rr_15_19_50_100treat2x_EMC[[2]], pc95_rr_20_50_50_100treat2x_EMC[[2]]))


write.table(df, file=paste0(outpath, "rrMHI", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)




################################################################################################
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")



