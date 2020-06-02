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



outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\Figures\\"


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
file1 <- "eggCounts_old_20-50_treat2x.RData"
file2 <- "eggCounts_new_20-50_treat2x.RData"
stub1 <- substring(file1, first=10, last=27)
stub2 <- substring(file2, first=10, last=27)
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
file1 <- "eggCounts_old_20-50_treat2x.RData"
file2 <- "eggCounts_new_20-50_treat2x.RData"
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
file1 <- "eggCounts_old_50-70_treat2x.RData"
file2 <- "eggCounts_new_50-70_treat2x.RData"
stub1 <- substring(file1, first=10, last=27)
stub2 <- substring(file2, first=10, last=27)
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
file1 <- "eggCounts_old_50-100_treat2x.RData"
file2 <- "eggCounts_new_50-100_treat2x.RData"
stub1 <- substring(file1, first=10, last=28)
stub2 <- substring(file2, first=10, last=28)
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
# WILCOXON SIGNED-RANK TESTS TO COMPARE OLD AND NEW TREATMENT STRATEGIES
################################################################################################

res_15_50 <- wilcox.test(meanPointPrevs_15_50_old, meanPointPrevs_15_50_new, paired=TRUE, alternative="greater")
res_15_19 <- wilcox.test(meanPointPrevs_15_19_old, meanPointPrevs_15_19_new, paired=TRUE, alternative="greater")
res_20_50 <- wilcox.test(meanPointPrevs_20_50_old, meanPointPrevs_20_50_new, paired=TRUE, alternative="greater")



