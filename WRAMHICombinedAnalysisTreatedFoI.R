
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

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\L_FoI_Treated\\"
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

##############################################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

p1 <- ggplot(df.msd.mhi.icl.old)
p1 <- p1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p1 <- p1 + geom_ribbon(aes(x=time, ymin=meanMHICounts-(2*sdMHICounts), ymax=meanMHICounts+(2*sdMHICounts)), fill="lightpink", alpha=0.5)
p1 <- p1 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p1 <- p1 + scale_x_continuous(name="Age (years)", expand=c(0,0), limits=c(0, 72), breaks=seq(0, 100, 5))
p1 <- p1 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p1)
ggsave(filename=paste0(outpath, "meanWomenMHI", stub1, "_ICL_20", ".jpeg"), plot=p1, dpi=300)


p2 <- ggplot(df.msd.mhi.icl.new)
p2 <- p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p2 <- p2 + geom_ribbon(aes(x=time, ymin=meanMHICounts-(2*sdMHICounts), ymax=meanMHICounts+(2*sdMHICounts)), fill="lightpink", alpha=0.5)
p2 <- p2 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p2 <- p2 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p2 <- p2 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p2)
ggsave(filename=paste0(outpath, "meanWomenMHI", stub2, "_ICL_20", ".jpeg"), plot=p2, dpi=300)


mp1 <- ggplot(df.mp.mhi.icl.old)
mp1 <- mp1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp1 <- mp1 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp1 <- mp1 + geom_line(aes(x=time, y=medianMHICounts), colour="red", size=1)
mp1 <- mp1 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp1 <- mp1 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(mp1)
ggsave(filename=paste0(outpath, "medianWomenMHI", stub1, "_ICL_20", ".jpeg"), plot=mp1, dpi=300)


mp2 <- ggplot(df.mp.mhi.icl.new)
mp2 <- mp2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp2 <- mp2 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp2 <- mp2 + geom_line(aes(x=time, y=medianMHICounts), colour="red", size=1)
mp2 <- mp2 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp2 <- mp2 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(mp2)
ggsave(filename=paste0(outpath, "medianWomenMHI", stub2, "_ICL_20", ".jpeg"), plot=mp2, dpi=300)


hp1 <- ggplot(df.msd.hi.icl.old)
hp1 <- hp1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hp1 <- hp1 + geom_ribbon(aes(x=time, ymin=meanHICounts-(2*sdHICounts), ymax=meanHICounts+(2*sdHICounts)), fill="lightpink", alpha=0.5)
hp1 <- hp1 + geom_line(aes(x=time, y=meanHICounts), colour="red", size=1)
hp1 <- hp1 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hp1 <- hp1 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hp1)
ggsave(filename=paste0(outpath, "meanWomenHI", stub1, "_ICL_20", ".jpeg"), plot=hp1, dpi=300)


hp2 <- ggplot(df.msd.hi.icl.new)
hp2 <- hp2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hp2 <- hp2 + geom_ribbon(aes(x=time, ymin=meanHICounts-(2*sdHICounts), ymax=meanHICounts+(2*sdHICounts)), fill="lightpink", alpha=0.5)
hp2 <- hp2 + geom_line(aes(x=time, y=meanHICounts), colour="red", size=1)
hp2 <- hp2 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hp2 <- hp2 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hp2)
ggsave(filename=paste0(outpath, "meanWomenHI", stub2, "_ICL_20", ".jpeg"), plot=hp2, dpi=300)


hmp1 <- ggplot(df.mp.hi.icl.old)
hmp1 <- hmp1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hmp1 <- hmp1 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
hmp1 <- hmp1 + geom_line(aes(x=time, y=medianHICounts), colour="red", size=1)
hmp1 <- hmp1 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hmp1 <- hmp1 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hmp1)
ggsave(filename=paste0(outpath, "medianWomenHI", stub1, "_ICL_20", ".jpeg"), plot=hmp1, dpi=300)


hmp2 <- ggplot(df.mp.hi.icl.new)
hmp2 <- hmp2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hmp2 <- hmp2 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
hmp2 <- hmp2 + geom_line(aes(x=time, y=medianHICounts), colour="red", size=1)
hmp2 <- hmp2 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hmp2 <- hmp2 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hmp2)
ggsave(filename=paste0(outpath, "medianWomenHI", stub2, "_ICL_20", ".jpeg"), plot=hmp2, dpi=300)




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



##############################################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

p3 <- ggplot(df.msd.mhi.emc.old)
p3 <- p3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p3 <- p3 + geom_ribbon(aes(x=time, ymin=meanMHICounts-(2*sdMHICounts), ymax=meanMHICounts+(2*sdMHICounts)), fill="lightpink", alpha=0.5)
p3 <- p3 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p3 <- p3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p3 <- p3 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p3)
ggsave(filename=paste0(outpath, "meanWomenMHI", stub1, "_EMC_20", ".jpeg"), plot=p3, dpi=300)


p4 <- ggplot(df.msd.mhi.emc.new)
p4 <- p4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p4 <- p4 + geom_ribbon(aes(x=time, ymin=meanMHICounts-(2*sdMHICounts), ymax=meanMHICounts+(2*sdMHICounts)), fill="lightpink", alpha=0.5)
p4 <- p4 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p4 <- p4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p4 <- p4 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(p4)
ggsave(filename=paste0(outpath, "meanWomenMHI", stub2, "_EMC_20", ".jpeg"), plot=p4, dpi=300)


mp3 <- ggplot(df.mp.mhi.emc.old)
mp3 <- mp3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp3 <- mp3 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp3 <- mp3 + geom_line(aes(x=time, y=medianMHICounts), colour="red", size=1)
mp3 <- mp3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp3 <- mp3 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(mp3)
ggsave(filename=paste0(outpath, "medianWomenMHI", stub1, "_EMC_20", ".jpeg"), plot=mp3, dpi=300)


mp4 <- ggplot(df.mp.mhi.emc.new)
mp4 <- mp4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp4 <- mp4 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp4 <- mp4 + geom_line(aes(x=time, y=medianMHICounts), colour="red", size=1)
mp4 <- mp4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp4 <- mp4 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(mp4)
ggsave(filename=paste0(outpath, "medianWomenMHI", stub2, "_EMC_20", ".jpeg"), plot=mp4, dpi=300)


hp3 <- ggplot(df.msd.hi.emc.old)
hp3 <- hp3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hp3 <- hp3 + geom_ribbon(aes(x=time, ymin=meanHICounts-(2*sdHICounts), ymax=meanHICounts+(2*sdHICounts)), fill="lightpink", alpha=0.5)
hp3 <- hp3 + geom_line(aes(x=time, y=meanHICounts), colour="red", size=1)
hp3 <- hp3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hp3 <- hp3 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hp3)
ggsave(filename=paste0(outpath, "meanWomenHI", stub1, "_EMC_20", ".jpeg"), plot=hp3, dpi=300)


hp4 <- ggplot(df.msd.hi.emc.new)
hp4 <- hp4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hp4 <- hp4 + geom_ribbon(aes(x=time, ymin=meanHICounts-(2*sdHICounts), ymax=meanHICounts+(2*sdHICounts)), fill="lightpink", alpha=0.5)
hp4 <- hp4 + geom_line(aes(x=time, y=meanHICounts), colour="red", size=1)
hp4 <- hp4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hp4 <- hp4 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hp4)
ggsave(filename=paste0(outpath, "meanWomenHI", stub2, "_EMC_20", ".jpeg"), plot=hp4, dpi=300)


hmp3 <- ggplot(df.mp.hi.emc.old)
hmp3 <- hmp3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hmp3 <- hmp3 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
hmp3 <- hmp3 + geom_line(aes(x=time, y=medianHICounts), colour="red", size=1)
hmp3 <- hmp3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
hmp3 <- hmp3 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hmp3)
ggsave(filename=paste0(outpath, "medianWomenHI", stub1, "_EMC_20", ".jpeg"), plot=hmp3, dpi=300)


hmp4 <- ggplot(df.mp.hi.emc.new)
hmp4 <- hmp4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
hmp4 <- hmp4 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
hmp4 <- hmp4 + geom_line(aes(x=time, y=medianHICounts), colour="red", size=1)
hmp4 <- hmp4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=c(0seq(0, 70, 5))
hmp4 <- hmp4 + scale_y_continuous(name="% women with HI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,20))
print(hmp4)
ggsave(filename=paste0(outpath, "medianWomenHI", stub2, "_EMC_20", ".jpeg"), plot=hmp4, dpi=300)


#####################################################################################

stub <- substring(stub1, first=5, last=20)

pic_mhi_means <- plot_grid(p1, p2, p3, p4, labels=c('a', 'b', 'c', 'd'), label_size=18, ncol=2, hjust=0, vjust=1.0)
save_plot(file=paste0(outpath, "MHI", stub1, ".jpeg"), pic_mhi_means, base_height=8.5, base_width=9.3)

pic_hi_means <- plot_grid(hp1, hp2, hp3, hp4, labels=c('a', 'b', 'c', 'd'), label_size=18, ncol=2, hjust=0, vjust=1.0)
save_plot(file=paste0(outpath, "HI", stub1, ".jpeg"), pic_hi_means, base_height=8.5, base_width=9.3)


#####################################################################################
# ICL
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

sd15_50 <- sd(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
sd15_19 <- sd(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
sd20_50 <- sd(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 

df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "MHI", stub1, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


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

sd15_50 <- sd(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
sd15_19 <- sd(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
sd20_50 <- sd(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 

df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "MHI", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



#####################################################################################
# determine number of women in cohort with heavy infection over time
#####################################################################################

HI <- eggCounts_old > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- do.call(cbind, HICountsList)
HIFrac <- HICounts / 500


pp15_50 <- mean(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
pp15_19 <- mean(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
pp20_50 <- mean(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 

sd15_50 <- sd(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
sd15_19 <- sd(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
sd20_50 <- sd(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "HI", stub1, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

#####################################################################################

HI <- eggCounts_new > heavyInfectionThreshold_counts
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- do.call(cbind, HICountsList)
HIFrac <- HICounts / 500


pp15_50 <- mean(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
pp15_19 <- mean(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
pp20_50 <- mean(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 

sd15_50 <- sd(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
sd15_19 <- sd(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
sd20_50 <- sd(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "HI", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



#####################################################################################
# EMC
#####################################################################################
# determine number of women in cohort with medium-to-heavy infection over time
#####################################################################################

MHI <- epgs_old > mediumHeavyInfectionThreshold_epg
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

sd15_50 <- sd(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
sd15_19 <- sd(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
sd20_50 <- sd(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "MHI", stub1, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

#####################################################################################

MHI <- epgs_new > mediumHeavyInfectionThreshold_epg
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

sd15_50 <- sd(colSums(MHIFrac[180:600, ]) / nrow(MHIFrac[180:600, ])) 
sd15_19 <- sd(colSums(MHIFrac[180:239, ]) / nrow(MHIFrac[180:239, ])) 
sd20_50 <- sd(colSums(MHIFrac[240:600, ]) / nrow(MHIFrac[240:600, ])) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "MHI", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


#####################################################################################
# determine number of women in cohort with heavy infection over time
#####################################################################################

HI <- epgs_old > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- do.call(cbind, HICountsList)
HIFrac <- HICounts / 500


pp15_50 <- mean(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
pp15_19 <- mean(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
pp20_50 <- mean(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 

sd15_50 <- sd(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
sd15_19 <- sd(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
sd20_50 <- sd(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "HI", stub1, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


#####################################################################################

HI <- epgs_new > heavyInfectionThreshold_epg
HI <- as.data.frame(cbind(results[, 1], HI))
names(HI) <- colnames(results)
HIList <- split(HI, HI$rep)
HIList <- lapply(HIList, function(x){ x[, 1] <- NULL; x})

HICountsList <- lapply(HIList, rowSums)
HICounts <- do.call(cbind, HICountsList)
HIFrac <- HICounts / 500


pp15_50 <- mean(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
pp15_19 <- mean(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
pp20_50 <- mean(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 

sd15_50 <- sd(colSums(HIFrac[180:600, ]) / nrow(HIFrac[180:600, ])) 
sd15_19 <- sd(colSums(HIFrac[180:239, ]) / nrow(HIFrac[180:239, ])) 
sd20_50 <- sd(colSums(HIFrac[240:600, ]) / nrow(HIFrac[240:600, ])) 


df <- data.frame(Mean=c(pp15_50, pp15_19, pp20_50), SD=c(sd15_50, sd15_19, sd20_50))
df <- df*100
df <- signif(df, digits=5)
write.table(df, file=paste0(outpath, "HI", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)







