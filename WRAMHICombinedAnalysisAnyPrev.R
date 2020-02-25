
library(ggplot2)
library(matrixStats)
library(cowplot)


set.seed(123)

################################################################################################################################ 
## DETERMINE PREVALENCE OF ANY INFECTION USING ICL EGG PRODUCTION FUNCTION
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


#################################################################################################################################
# read in data
#################################################################################################################################


#path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Untreated\\correctedTimeStep\\"
#file1 <- "fertilisedFemaleWorms_old_50-100K_treat2x_100reps.RData"
#file2 <- "fertilisedFemaleWorms_new_50-100K_treat2x_100reps.RData"

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\L_FoI_Treated\\"
file1 <- "fertilisedFemaleWorms_old_50-70_treat2x.RData"
file2 <- "fertilisedFemaleWorms_new_50-70_treat2x.RData"

#stub1 <- substring(file1, first=22, last=41)
#stub2 <- substring(file2, first=22, last=41)

stub1 <- substring(file1, first=22, last=39)
stub2 <- substring(file2, first=22, last=39)

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
# PREVALENCE OF ANY INFECTION - OLD TREATMENT STRATEGY
##################################################################################################


prev <- eggCounts_old >= 1
prev <- as.data.frame(cbind(results[, 1], prev))
names(prev) <- colnames(results)
prevList <- split(prev, prev$rep)
prevList <- lapply(prevList, function(x){ x[, 1] <- NULL; x})

prevCountsList <- lapply(prevList, rowSums)
prevCounts <- unlist(prevCountsList)

years <- 100
timerep=rep(seq(0, years, 1/12), 100)
time=seq(0, years, 1/12)


temp <- do.call(cbind, prevCountsList)
meanprevCounts <- rowMeans(temp)
sdprevCounts <- sqrt(rowVars(temp))
medianprevCounts <- rowMedians(temp)
p95prevCounts <- rowQuantiles(temp, probs=c(0.05, 0.95))


df.msd.icl.old <- data.frame(time=time, meanprevCounts=meanprevCounts/500*100, sdprevCounts=sdprevCounts/500*100)
df.mp.icl.old <- data.frame(time=time, medianprevCounts=medianprevCounts/500*100, p5=p95prevCounts[,1]/500*100, p95=p95prevCounts[,2]/500*100)



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
		newSample <- getSetOfEggCounts(i, lambda, gamma, k_epg) 
		eggSamples <- c(eggSamples, newSample)
	}
	eggList[[i]] <- eggSamples	
}


# get egg output based on number of female worms
eggCounts_new <- apply(results_new[, -1], c(1,2), drawEggCounts, eggList)

##################################################################################################
# PREVALENCE OF ANY INFECTION - NEW TREATMENT STRATEGY
##################################################################################################


prev <- eggCounts_new >= 1
prev <- as.data.frame(cbind(results[, 1], prev))
names(prev) <- colnames(results)
prevList <- split(prev, prev$rep)
prevList <- lapply(prevList, function(x){ x[, 1] <- NULL; x})

prevCountsList <- lapply(prevList, rowSums)
prevCounts <- unlist(prevCountsList)

years <- 100
timerep=rep(seq(0, years, 1/12), 100)
time=seq(0, years, 1/12)


temp <- do.call(cbind, prevCountsList)
meanprevCounts <- rowMeans(temp)
sdprevCounts <- sqrt(rowVars(temp))
medianprevCounts <- rowMedians(temp)
p95prevCounts <- rowQuantiles(temp, probs=c(0.05, 0.95))


df.msd.icl.new <- data.frame(time=time, meanprevCounts=meanprevCounts/500*100, sdprevCounts=sdprevCounts/500*100)
df.mp.icl.new <- data.frame(time=time, medianprevCounts=medianprevCounts/500*100, p5=p95prevCounts[,1]/500*100, p95=p95prevCounts[,2]/500*100)



##############################################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

p1 <- ggplot(df.msd.icl.old)
p1 <- p1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p1 <- p1 + geom_ribbon(aes(x=time, ymin=meanprevCounts-(2*sdprevCounts), ymax=meanprevCounts+(2*sdprevCounts)), fill="lightpink", alpha=0.5)
p1 <- p1 + geom_line(aes(x=time, y=meanprevCounts), colour="red", size=1)
p1 <- p1 + scale_x_continuous(name="Age (years)", expand=c(0,0), limits=c(0, 72), breaks=seq(0, 100, 5))
p1 <- p1 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(p1)
ggsave(filename=paste0(outpath, "meanWomenPrev", stub1, "_ICL", ".jpeg"), plot=p1, dpi=300)


mp1 <- ggplot(df.mp.icl.old)
mp1 <- mp1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp1 <- mp1 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp1 <- mp1 + geom_line(aes(x=time, y=medianprevCounts), colour="red", size=1)
mp1 <- mp1 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp1 <- mp1 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(mp1)
ggsave(filename=paste0(outpath, "medianWomenPrev", stub1, "_ICL", ".jpeg"), plot=mp1, dpi=300)


p2 <- ggplot(df.msd.icl.new)
p2 <- p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p2 <- p2 + geom_ribbon(aes(x=time, ymin=meanprevCounts-(2*sdprevCounts), ymax=meanprevCounts+(2*sdprevCounts)), fill="lightpink", alpha=0.5)
p2 <- p2 + geom_line(aes(x=time, y=meanprevCounts), colour="red", size=1)
p2 <- p2 + scale_x_continuous(name="Age (years)", expand=c(0,0), limits=c(0, 72), breaks=seq(0, 100, 5))
p2 <- p2 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(p2)
ggsave(filename=paste0(outpath, "meanWomenPrev", stub2, "_ICL", ".jpeg"), plot=p2, dpi=300)


mp2 <- ggplot(df.mp.icl.new)
mp2 <- mp2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp2 <- mp2 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp2 <- mp2 + geom_line(aes(x=time, y=medianprevCounts), colour="red", size=1)
mp2 <- mp2 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp2 <- mp2 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(mp2)
ggsave(filename=paste0(outpath, "medianWomenPrev", stub2, "_ICL", ".jpeg"), plot=mp2, dpi=300)




################################################################################################################################ 
## DETERMINE PREVALENCE OF ANY INFECTION USING EMC EGG PRODUCTION FUNCTION
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

MHI <- epgs_old > 0
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


df.msd.emc.old <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, sdMHICounts=sdMHICounts/500*100)
df.mp.emc.old <- data.frame(time=time, medianMHICounts=medianMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)


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


MHI <- epgs_new > 0
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


df.msd.emc.new <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, sdMHICounts=sdMHICounts/500*100)
df.mp.emc.new <- data.frame(time=time, medianMHICounts=medianMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)

##############################################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

p3 <- ggplot(df.msd.emc.old)
p3 <- p3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p3 <- p3 + geom_ribbon(aes(x=time, ymin=meanMHICounts-(2*sdMHICounts), ymax=meanMHICounts+(2*sdMHICounts)), fill="lightpink", alpha=0.5)
p3 <- p3 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p3 <- p3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p3 <- p3 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(p3)
ggsave(filename=paste0(outpath, "meanWomenPrev", stub1, "_EMC", ".jpeg"), plot=p3, dpi=300)


p4 <- ggplot(df.msd.emc.new)
p4 <- p4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
p4 <- p4 + geom_ribbon(aes(x=time, ymin=meanMHICounts-(2*sdMHICounts), ymax=meanMHICounts+(2*sdMHICounts)), fill="lightpink", alpha=0.5)
p4 <- p4 + geom_line(aes(x=time, y=meanMHICounts), colour="red", size=1)
p4 <- p4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
p4 <- p4 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(p4)
ggsave(filename=paste0(outpath, "meanWomenPrev", stub2, "_EMC", ".jpeg"), plot=p4, dpi=300)


mp3 <- ggplot(df.mp.emc.old)
mp3 <- mp3 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp3 <- mp3 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp3 <- mp3 + geom_line(aes(x=time, y=medianMHICounts), colour="red", size=1)
mp3 <- mp3 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp3 <- mp3 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(mp3)
ggsave(filename=paste0(outpath, "medianWomenPrev", stub1, "_EMC", ".jpeg"), plot=mp3, dpi=300)


mp4 <- ggplot(df.mp.emc.new)
mp4 <- mp4 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), legend.title=element_text(size=16), legend.text=element_text(size=16), legend.position="bottom")
mp4 <- mp4 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95), fill="lightpink", alpha=0.5)
mp4 <- mp4 + geom_line(aes(x=time, y=medianMHICounts), colour="red", size=1)
mp4 <- mp4 + scale_x_continuous(name="Age (years)", expand=c(0, 0), limits=c(0, 72), breaks=seq(0, 100, 5))
mp4 <- mp4 + scale_y_continuous(name="Prevalence (%)", expand=c(0, 0)) + coord_cartesian(ylim=c(0,100))
print(mp4)
ggsave(filename=paste0(outpath, "medianWomenPrev", stub2, "_EMC", ".jpeg"), plot=mp4, dpi=300)


#####################################################################################

stub <- substring(stub1, first=5, last=20)

pic_prev_means <- plot_grid(p1, p2, p3, p4, labels=c('a', 'b', 'c', 'd'), label_size=18, ncol=2, hjust=0, vjust=1.0)
save_plot(file=paste0(outpath, "prev", stub1, ".jpeg"), pic_prev_means, base_height=8.5, base_width=9.3)


#####################################################################################
# ICL
#####################################################################################
# determine number of women in cohort with medium-to-heavy infection over time
#####################################################################################

MHI <- eggCounts_old >= 1
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
#write.table(df, file=paste0(outpath, "prev", stub1, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)


#####################################################################################

MHI <- eggCounts_new >= 1
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
#write.table(df, file=paste0(outpath, "prev", stub2, "_ICL.txt"), sep="\t", col.names=TRUE, row.names=FALSE)



#####################################################################################
# EMC
#####################################################################################
# determine number of women in cohort with medium-to-heavy infection over time
#####################################################################################

MHI <- epgs_old > 0
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
#write.table(df, file=paste0(outpath, "prev", stub1, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

#####################################################################################

MHI <- epgs_new > 0
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
#write.table(df, file=paste0(outpath, "prev", stub2, "_EMC.txt"), sep="\t", col.names=TRUE, row.names=FALSE)







