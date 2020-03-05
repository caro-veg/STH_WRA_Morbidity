
library(ggplot2)
library(matrixStats)
library(cowplot)


set.seed(123)

################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING ICL EGG PRODUCTION FUNCTION
################################################################################################################################ 


################################################################################################################################ 
## Returns a random egg count readings from a single sample based on number of fertilised female worms inside a host
################################################################################################################################

getSetOfEggCounts <- function(fertilisedFemales, lambda, gamma, k_epg)
{  
	z <- exp(-gamma)	# Fecundity parameter z
	meanCount <- fertilisedFemales*lambda*z^fertilisedFemales
  
  	readings <- rnbinom(length(meanCount), mu=meanCount, size=k_epg)
  	return(readings)
} 

#################################################################################################################################
# draws egg counts for given number of fertilised females using presampled egg counts
#################################################################################################################################

drawEggCounts <- function(fertilisedFemales, eggList)
{
	if(fertilisedFemales==0)
		return(0)
	index <- sample(1:length(eggList[[1]]), 1)
	return(eggList[[fertilisedFemales]][[index]])
}


#################################################################################################################################

lambda <- 3
gamma <- 0.02
k_epg <- 0.35


# thresholds for medium-to-heavy infection
mediumHeavyInfectionThreshold_epg <- 2000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
mediumHeavyInfectionThreshold_counts <- mediumHeavyInfectionThreshold_epg / diagnosticDivisor


# thresholds for heavy infection
heavyInfectionThreshold_epg <- 4000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
heavyInfectionThreshold_counts <- heavyInfectionThreshold_epg / diagnosticDivisor



#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Untreated\\correctedTimeStep\\"
file1 <- "fertilisedFemaleWorms_old_20-50KK_treat1x_100reps.RData"
file2 <- "fertilisedFemaleWorms_new_20-50KK_treat1x_100reps.RData"
stub1 <- substring(file1, first=22, last=41)
stub2 <- substring(file2, first=22, last=41)
results_old <- get(load(paste0(path, file1)))
results_new <- get(load(paste0(path, file2)))


############################################################################################
# get egg counts for old treatment strategy
############################################################################################

ffwMax <- max(results_old[,-1])


eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_old <- apply(results_old[, -1], c(1,2), drawEggCounts, eggList)


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
# HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

HI <- eggCounts_old > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.icl.old <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.icl.old <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)



############################################################################################
# get egg counts for new treatment strategy
############################################################################################

ffwMax <- max(results_new[,-1])


eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_new <- apply(results_new[, -1], c(1,2), drawEggCounts, eggList)


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
# HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

HI <- eggCounts_new > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.icl.new <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.icl.new <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)


red.ratio.mhi.icl <- (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts

red.mean.mhi.icl.mod.1x.15_50 <- mean(red.ratio.mhi.icl[181:601])
red.sd.mhi.icl.mod.1x.15_50 <-  sd(red.ratio.mhi.icl[181:601])

red.mean.mhi.icl.mod.1x.15_19 <- mean(red.ratio.mhi.icl[181:240])
red.sd.mhi.icl.mod.1x.15_19 <-  sd(red.ratio.mhi.icl[181:240])

red.mean.mhi.icl.mod.1x.20_50 <- mean(red.ratio.mhi.icl[241:601])
red.sd.mhi.icl.mod.1x.20_50 <-  sd(red.ratio.mhi.icl[241:601])


##############################################################################################################################################
################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG PRODUCTION FUNCTION
################################################################################################################################ 


################################################################################################################################ 
## Returns a random epgs from a single sample based on number of fertilised female worms inside a host
################################################################################################################################

getSetOfEpgs <- function(fertilisedFemales, a, k_epg)
{  
	delta <- rgamma(length(fertilisedFemales), shape=50, rate=50)
	b <- 1500 * delta
	xi <- a * fertilisedFemales / (1 + a * fertilisedFemales / b)
  
  	readings <- rnbinom(length(xi), mu=xi, size=k_epg)
  	return(readings)
} 


#################################################################################################################################
# draws epgs for given number of fertilised females using presampled epgs
#################################################################################################################################


drawEpgs <- function(fertilisedFemales, eggList)
{
	if(fertilisedFemales==0)
		return(0)
	index <- sample(1:length(eggList[[1]]), 1)
	return(eggList[[fertilisedFemales]][[index]])
}

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\EMC_data\\correctedTimeStep\\"
file1 <- "fertilisedFemaleWorms_old_20-50KK_treat1x_100reps.RData"
file2 <- "fertilisedFemaleWorms_new_20-50KK_treat1x_100reps.RData"
stub1 <- substring(file1, first=22, last=41)
stub2 <- substring(file2, first=22, last=41)
results_old <- get(load(paste0(path, file1)))
results_new <- get(load(paste0(path, file2)))


a <- 200
k_epg <- 0.35


############################################################################################
# get epgs for old treatment strategy
############################################################################################


ffwMax <- max(results_old[,-1])
eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEpgs(i, a, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
epgs_old <- apply(results_old[, -1], c(1,2), drawEpgs, eggList)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

MHI <- epgs_old > mediumHeavyInfectionThreshold_epg
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
# HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

HI <- epgs_old > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.emc.old <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.emc.old <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)



############################################################################################
# get epgs for new treatment strategy
############################################################################################

ffwMax <- max(results_new[,-1])
eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEpgs(i, a, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
epgs_new <- apply(results_new[, -1], c(1,2), drawEpgs, eggList)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- epgs_new > mediumHeavyInfectionThreshold_epg
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


##################################################################################################
# HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

HI <- epgs_new > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.emc.new <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.emc.new <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)


##################################################################################################################################


red.ratio.mhi.emc <- (df.msd.mhi.emc.old$meanMHICounts - df.msd.mhi.emc.new$meanMHICounts) / df.msd.mhi.emc.old$meanMHICounts

red.mean.mhi.emc.mod.1x.15_50 <- mean(red.ratio.mhi.emc[181:601])
red.sd.mhi.emc.mod.1x.15_50 <-  sd(red.ratio.mhi.emc[181:601])

red.mean.mhi.emc.mod.1x.15_19 <- mean(red.ratio.mhi.emc[181:240])
red.sd.mhi.emc.mod.1x.15_19 <-  sd(red.ratio.mhi.emc[181:240])

red.mean.mhi.emc.mod.1x.20_50 <- mean(red.ratio.mhi.emc[241:601])
red.sd.mhi.emc.mod.1x.20_50 <-  sd(red.ratio.mhi.emc[241:601])


##################################################################################################################################


################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING ICL EGG PRODUCTION FUNCTION
################################################################################################################################ 


################################################################################################################################ 
## Returns a random egg count readings from a single sample based on number of fertilised female worms inside a host
################################################################################################################################

getSetOfEggCounts <- function(fertilisedFemales, lambda, gamma, k_epg)
{  
	z <- exp(-gamma)	# Fecundity parameter z
	meanCount <- fertilisedFemales*lambda*z^fertilisedFemales
  
  	readings <- rnbinom(length(meanCount), mu=meanCount, size=k_epg)
  	return(readings)
} 

#################################################################################################################################
# draws egg counts for given number of fertilised females using presampled egg counts
#################################################################################################################################

drawEggCounts <- function(fertilisedFemales, eggList)
{
	if(fertilisedFemales==0)
		return(0)
	index <- sample(1:length(eggList[[1]]), 1)
	return(eggList[[fertilisedFemales]][[index]])
}


#################################################################################################################################

lambda <- 3
gamma <- 0.02
k_epg <- 0.35


# thresholds for medium-to-heavy infection
mediumHeavyInfectionThreshold_epg <- 2000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
mediumHeavyInfectionThreshold_counts <- mediumHeavyInfectionThreshold_epg / diagnosticDivisor


# thresholds for heavy infection
heavyInfectionThreshold_epg <- 4000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
heavyInfectionThreshold_counts <- heavyInfectionThreshold_epg / diagnosticDivisor



#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Untreated\\correctedTimeStep\\"
file1 <- "fertilisedFemaleWorms_old_20-50KK_treat2x_100reps.RData"
file2 <- "fertilisedFemaleWorms_new_20-50KK_treat2x_100reps.RData"
stub1 <- substring(file1, first=22, last=41)
stub2 <- substring(file2, first=22, last=41)
results_old <- get(load(paste0(path, file1)))
results_new <- get(load(paste0(path, file2)))


############################################################################################
# get egg counts for old treatment strategy
############################################################################################

ffwMax <- max(results_old[,-1])


eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_old <- apply(results_old[, -1], c(1,2), drawEggCounts, eggList)


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
# HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

HI <- eggCounts_old > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.icl.old <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.icl.old <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)



############################################################################################
# get egg counts for new treatment strategy
############################################################################################

ffwMax <- max(results_new[,-1])


eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_new <- apply(results_new[, -1], c(1,2), drawEggCounts, eggList)


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
# HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

HI <- eggCounts_new > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.icl.new <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.icl.new <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)


red.ratio.mhi.icl <- (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts

red.mean.mhi.icl.mod.2x.15_50 <- mean(red.ratio.mhi.icl[181:601])
red.sd.mhi.icl.mod.2x.15_50 <-  sd(red.ratio.mhi.icl[181:601])

red.mean.mhi.icl.mod.2x.15_19 <- mean(red.ratio.mhi.icl[181:240])
red.sd.mhi.icl.mod.2x.15_19 <-  sd(red.ratio.mhi.icl[181:240])

red.mean.mhi.icl.mod.2x.20_50 <- mean(red.ratio.mhi.icl[241:601])
red.sd.mhi.icl.mod.2x.20_50 <-  sd(red.ratio.mhi.icl[241:601])


##############################################################################################################################################
################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG PRODUCTION FUNCTION
################################################################################################################################ 


################################################################################################################################ 
## Returns a random epgs from a single sample based on number of fertilised female worms inside a host
################################################################################################################################

getSetOfEpgs <- function(fertilisedFemales, a, k_epg)
{  
	delta <- rgamma(length(fertilisedFemales), shape=50, rate=50)
	b <- 1500 * delta
	xi <- a * fertilisedFemales / (1 + a * fertilisedFemales / b)
  
  	readings <- rnbinom(length(xi), mu=xi, size=k_epg)
  	return(readings)
} 


#################################################################################################################################
# draws epgs for given number of fertilised females using presampled epgs
#################################################################################################################################


drawEpgs <- function(fertilisedFemales, eggList)
{
	if(fertilisedFemales==0)
		return(0)
	index <- sample(1:length(eggList[[1]]), 1)
	return(eggList[[fertilisedFemales]][[index]])
}

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\EMC_data\\correctedTimeStep\\"
file1 <- "fertilisedFemaleWorms_old_20-50KK_treat2x_100reps.RData"
file2 <- "fertilisedFemaleWorms_new_20-50KK_treat2x_100reps.RData"
stub1 <- substring(file1, first=22, last=41)
stub2 <- substring(file2, first=22, last=41)
results_old <- get(load(paste0(path, file1)))
results_new <- get(load(paste0(path, file2)))


a <- 200
k_epg <- 0.35


############################################################################################
# get epgs for old treatment strategy
############################################################################################


ffwMax <- max(results_old[,-1])
eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEpgs(i, a, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
epgs_old <- apply(results_old[, -1], c(1,2), drawEpgs, eggList)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

MHI <- epgs_old > mediumHeavyInfectionThreshold_epg
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
# HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

HI <- epgs_old > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.emc.old <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.emc.old <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)



############################################################################################
# get epgs for new treatment strategy
############################################################################################

ffwMax <- max(results_new[,-1])
eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEpgs(i, a, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
epgs_new <- apply(results_new[, -1], c(1,2), drawEpgs, eggList)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- epgs_new > mediumHeavyInfectionThreshold_epg
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


##################################################################################################
# HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

HI <- epgs_new > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.emc.new <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.emc.new <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)


##################################################################################################################################


red.ratio.mhi.emc <- (df.msd.mhi.emc.old$meanMHICounts - df.msd.mhi.emc.new$meanMHICounts) / df.msd.mhi.emc.old$meanMHICounts

red.mean.mhi.emc.mod.2x.15_50 <- mean(red.ratio.mhi.emc[181:601])
red.sd.mhi.emc.mod.2x.15_50 <-  sd(red.ratio.mhi.emc[181:601])

red.mean.mhi.emc.mod.2x.15_19 <- mean(red.ratio.mhi.emc[181:240])
red.sd.mhi.emc.mod.2x.15_19 <-  sd(red.ratio.mhi.emc[181:240])

red.mean.mhi.emc.mod.2x.20_50 <- mean(red.ratio.mhi.emc[241:601])
red.sd.mhi.emc.mod.2x.20_50 <-  sd(red.ratio.mhi.emc[241:601])

#######################################################################################################################################


################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING ICL EGG PRODUCTION FUNCTION
################################################################################################################################ 


################################################################################################################################ 
## Returns a random egg count readings from a single sample based on number of fertilised female worms inside a host
################################################################################################################################

getSetOfEggCounts <- function(fertilisedFemales, lambda, gamma, k_epg)
{  
	z <- exp(-gamma)	# Fecundity parameter z
	meanCount <- fertilisedFemales*lambda*z^fertilisedFemales
  
  	readings <- rnbinom(length(meanCount), mu=meanCount, size=k_epg)
  	return(readings)
} 

#################################################################################################################################
# draws egg counts for given number of fertilised females using presampled egg counts
#################################################################################################################################

drawEggCounts <- function(fertilisedFemales, eggList)
{
	if(fertilisedFemales==0)
		return(0)
	index <- sample(1:length(eggList[[1]]), 1)
	return(eggList[[fertilisedFemales]][[index]])
}


#################################################################################################################################

lambda <- 3
gamma <- 0.02
k_epg <- 0.35


# thresholds for medium-to-heavy infection
mediumHeavyInfectionThreshold_epg <- 2000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
mediumHeavyInfectionThreshold_counts <- mediumHeavyInfectionThreshold_epg / diagnosticDivisor


# thresholds for heavy infection
heavyInfectionThreshold_epg <- 4000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
heavyInfectionThreshold_counts <- heavyInfectionThreshold_epg / diagnosticDivisor



#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Untreated\\correctedTimeStep\\"
file1 <- "fertilisedFemaleWorms_old_50-100K_treat2x_100reps.RData"
file2 <- "fertilisedFemaleWorms_new_50-100K_treat2x_100reps.RData"
stub1 <- substring(file1, first=22, last=41)
stub2 <- substring(file2, first=22, last=41)
results_old <- get(load(paste0(path, file1)))
results_new <- get(load(paste0(path, file2)))


############################################################################################
# get egg counts for old treatment strategy
############################################################################################

ffwMax <- max(results_old[,-1])


eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_old <- apply(results_old[, -1], c(1,2), drawEggCounts, eggList)


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
# HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

HI <- eggCounts_old > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.icl.old <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.icl.old <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)



############################################################################################
# get egg counts for new treatment strategy
############################################################################################

ffwMax <- max(results_new[,-1])


eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_new <- apply(results_new[, -1], c(1,2), drawEggCounts, eggList)


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
# HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

HI <- eggCounts_new > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.icl.new <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.icl.new <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)


red.ratio.mhi.icl <- (df.msd.mhi.icl.old$meanMHICounts - df.msd.mhi.icl.new$meanMHICounts) / df.msd.mhi.icl.old$meanMHICounts

red.mean.mhi.icl.high.2x.15_50 <- mean(red.ratio.mhi.icl[181:601])
red.sd.mhi.icl.high.2x.15_50 <-  sd(red.ratio.mhi.icl[181:601])

red.mean.mhi.icl.high.2x.15_19 <- mean(red.ratio.mhi.icl[181:240])
red.sd.mhi.icl.high.2x.15_19 <-  sd(red.ratio.mhi.icl[181:240])

red.mean.mhi.icl.high.2x.20_50 <- mean(red.ratio.mhi.icl[241:601])
red.sd.mhi.icl.high.2x.20_50 <-  sd(red.ratio.mhi.icl[241:601])


##############################################################################################################################################
################################################################################################################################ 
## DETERMINE PERCENTAGE OF WOMEN WITH MEDIUM-HIGH INTENSITY INFECTION USING EMC EGG PRODUCTION FUNCTION
################################################################################################################################ 


################################################################################################################################ 
## Returns a random epgs from a single sample based on number of fertilised female worms inside a host
################################################################################################################################

getSetOfEpgs <- function(fertilisedFemales, a, k_epg)
{  
	delta <- rgamma(length(fertilisedFemales), shape=50, rate=50)
	b <- 1500 * delta
	xi <- a * fertilisedFemales / (1 + a * fertilisedFemales / b)
  
  	readings <- rnbinom(length(xi), mu=xi, size=k_epg)
  	return(readings)
} 


#################################################################################################################################
# draws epgs for given number of fertilised females using presampled epgs
#################################################################################################################################


drawEpgs <- function(fertilisedFemales, eggList)
{
	if(fertilisedFemales==0)
		return(0)
	index <- sample(1:length(eggList[[1]]), 1)
	return(eggList[[fertilisedFemales]][[index]])
}

#################################################################################################################################
# read in data
#################################################################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\EMC_data\\correctedTimeStep\\"
file1 <- "fertilisedFemaleWorms_old_50-100K_treat2x_100reps.RData"
file2 <- "fertilisedFemaleWorms_new_50-100K_treat2x_100reps.RData"
stub1 <- substring(file1, first=22, last=41)
stub2 <- substring(file2, first=22, last=41)
results_old <- get(load(paste0(path, file1)))
results_new <- get(load(paste0(path, file2)))


a <- 200
k_epg <- 0.35


############################################################################################
# get epgs for old treatment strategy
############################################################################################


ffwMax <- max(results_old[,-1])
eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEpgs(i, a, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
epgs_old <- apply(results_old[, -1], c(1,2), drawEpgs, eggList)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

MHI <- epgs_old > mediumHeavyInfectionThreshold_epg
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
# HEAVY INTENSITY INFECTIONS - OLD TREATMENT STRATEGY
##################################################################################################

HI <- epgs_old > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.emc.old <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.emc.old <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)



############################################################################################
# get epgs for new treatment strategy
############################################################################################

ffwMax <- max(results_new[,-1])
eggList <- list()
for(i in 1:ffwMax)
{
	eggSamples <- c()
	for(j in 1:500)
	{
		newSample <- getSetOfEpgs(i, a, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
epgs_new <- apply(results_new[, -1], c(1,2), drawEpgs, eggList)



##################################################################################################
# MEDIUM-TO-HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################


MHI <- epgs_new > mediumHeavyInfectionThreshold_epg
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


##################################################################################################
# HEAVY INTENSITY INFECTIONS - NEW TREATMENT STRATEGY
##################################################################################################

HI <- epgs_new > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- unlist(HICountsList)

temp <- do.call(cbind, HICountsList)
meanHICounts <- rowMeans(temp)
sdHICounts <- sqrt(rowVars(temp))
medianHICounts <- rowMedians(temp)
p95HICounts <- rowQuantiles(temp, probs=c(0.05, 0.95))

df.msd.hi.emc.new <- data.frame(time=time, meanHICounts=meanHICounts/500*100, sdHICounts=sdHICounts/500*100)
df.mp.hi.emc.new <- data.frame(time=time, medianHICounts=medianHICounts/500*100, p5=p95HICounts[,1]/500*100, p95=p95HICounts[,2]/500*100)


##################################################################################################################################


red.ratio.mhi.emc <- (df.msd.mhi.emc.old$meanMHICounts - df.msd.mhi.emc.new$meanMHICounts) / df.msd.mhi.emc.old$meanMHICounts

red.mean.mhi.emc.high.2x.15_50 <- mean(red.ratio.mhi.emc[181:601])
red.sd.mhi.emc.high.2x.15_50 <-  sd(red.ratio.mhi.emc[181:601])

red.mean.mhi.emc.high.2x.15_19 <- mean(red.ratio.mhi.emc[181:240])
red.sd.mhi.emc.high.2x.15_19 <-  sd(red.ratio.mhi.emc[181:240])

red.mean.mhi.emc.high.2x.20_50 <- mean(red.ratio.mhi.emc[241:601])
red.sd.mhi.emc.high.2x.20_50 <-  sd(red.ratio.mhi.emc[241:601])



################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

ratios <- c(red.mean.mhi.icl.mod.1x.15_50, red.mean.mhi.emc.mod.1x.15_50, red.mean.mhi.icl.mod.2x.15_50, red.mean.mhi.emc.mod.2x.15_50, red.mean.mhi.icl.high.2x.15_50, red.mean.mhi.emc.high.2x.15_50)*100
sds <- c(red.sd.mhi.icl.mod.1x.15_50, red.sd.mhi.emc.mod.1x.15_50, red.sd.mhi.icl.mod.2x.15_50, red.sd.mhi.emc.mod.2x.15_50, red.sd.mhi.icl.high.2x.15_50, red.sd.mhi.emc.high.2x.15_50)*100
scenarios <- seq(1, 12, 2)
unis <- c(rep(c("ICL", "EMC"), 3))

df <- data.frame(ratios=ratios, scenarios=scenarios, sds=sds, unis=unis)

p <- ggplot(data=df, aes(x=scenarios, y=ratios, colour=unis))
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank(), legend.text=element_text(size=16), legend.position="bottom")
p <- p + geom_point(size=4)
p <- p + geom_errorbar(aes(ymin=ratios-2*sds, ymax=ratios+2*sds))
p <- p + scale_y_continuous(breaks=seq(0, 110, 20)) + coord_cartesian(ylim=c(0,100)) 
p <- p + xlab("Scenarios") + ylab("Relative reduction in \nM&HI infection prevalence (%)")
p <- p + scale_x_continuous(breaks=c(2, 6, 10), labels=c("Moderate prevalence,\nannual treatment", "Moderate prevalence,\nsemi-annual treatment", "High prevalence,\nsemi-annual treatment"))
print(p)
ggsave(paste0(outpath, "summaryFigureCorrectedMarch2020.jpeg"), p, dpi=300)




