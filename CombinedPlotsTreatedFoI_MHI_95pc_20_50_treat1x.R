
library(ggplot2)
library(matrixStats)
library(cowplot)


# thresholds for medium-to-heavy infection
mediumHeavyInfectionThreshold_epg <- 2000	# WHO value for hookworm
diagnosticDivisor <- 24				# 24 for Kato-Katz. 
mediumHeavyInfectionThreshold_counts <- mediumHeavyInfectionThreshold_epg / diagnosticDivisor



outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\FiguresAntonio\\"

# for column names and iterations
path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file <- "fertilisedFemaleWorms_old_20-50_treat1x.RData"
stub <- substring(file, first=22, last=45)
results <- get(load(paste0(path, file)))



###############################################################################################
# FUNCTION TO DETERMINE MHI infections
###############################################################################################

getMHI <- function(eggCounts, threshold, scenario, refResults)
{
	
	MHI <- eggCounts > threshold
	MHI <- as.data.frame(cbind(refResults[, 1], MHI))
	names(MHI) <- colnames(refResults)
	MHIList <- split(MHI, MHI$rep)
	MHIList <- lapply(MHIList, function(x){ x[, 1] <- NULL; x})

	MHICountsList <- lapply(MHIList, rowSums)
	
	years <- 100
	timerep=rep(seq(0, years, 1/12), 100)
	time=seq(0, years, 1/12)

	temp <- do.call(cbind, MHICountsList)
	meanMHICounts <- rowMeans(temp)
	sdMHICounts <- sqrt(rowVars(temp))
	medianMHICounts <- rowMedians(temp)
	p95MHICounts <- rowQuantiles(temp, probs=c(0.025, 0.975))

	df.mhi <- data.frame(time=time, meanMHICounts=meanMHICounts/500*100, p5=p95MHICounts[,1]/500*100, p95=p95MHICounts[,2]/500*100)
	df.mhi$scenario <- scenario

	return(df.mhi)
}


###############################################################################################
# Figure 1a: 20-50% prevalence treat 1x per year, ICL model
###############################################################################################
# read in data
###############################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\ICLsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat1x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat1x.RData"
file3 <- "eggCounts_20-50_untreated.RData"

eggCounts_untreated <- get(load(paste0(path, file3)))
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))



###############################################################################################
# determine MHI infections in untreated scenario
###############################################################################################

df.mhi.icl.untreat <- getMHI(eggCounts_untreated, mediumHeavyInfectionThreshold_counts, "untreated", results)


###############################################################################################
# determine MHI infections in old treatment scenario
###############################################################################################

df.mhi.icl.old <- getMHI(eggCounts_old, mediumHeavyInfectionThreshold_counts, "old", results)


###############################################################################################
# determine MHI infections in new treatment scenario
###############################################################################################


df.mhi.icl.new <- getMHI(eggCounts_new, mediumHeavyInfectionThreshold_counts, "new", results)


################################################################################################
# combine data frames and plot
################################################################################################

df.mhi.icl.20_50 <- rbind(df.mhi.icl.untreat, df.mhi.icl.old, df.mhi.icl.new)
df.mhi.icl.20_50$scenario <- factor(df.mhi.icl.20_50$scenario, levels=c("new", "old", "untreated"))

#df.mhi.icl.20_50 <- rbind(df.mhi.icl.old, df.mhi.icl.new)
#df.mhi.icl.20_50$scenario <- factor(df.mhi.icl.20_50$scenario, levels=c("new", "old"))


p1 <- ggplot(df.mhi.icl.20_50)
p1 <- p1 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=24), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))
p1 <- p1 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95, fill=scenario), alpha=0.3) + scale_fill_manual(values=c("red", "blue", "seagreen"))
p1 <- p1 + geom_line(aes(x=time, y=meanMHICounts, colour=scenario), size=1, alpha=0.8) + scale_colour_manual(values=c("red", "blue", "seagreen"))
p1 <- p1 + scale_x_continuous(name="Age (years)", expand=c(0,0), limits=c(0, 72), breaks=seq(0, 100, 5))
p1 <- p1 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,21))
print(p1)




###############################################################################################
# Figure 1b: 20-50% prevalence treat 1x per year, EMC model
###############################################################################################
# read in data
###############################################################################################

path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\ForPublication\\EMCsimoutput\\"
file1 <- "redoneEggCountsOneList_old_20-50_treat1x.RData"
file2 <- "redoneEggCountsOneList_new_20-50_treat1x.RData"
file3 <- "eggCounts_20-50_untreated.RData"
stub1 <- substring(file1, first=10, last=27)
stub2 <- substring(file2, first=10, last=27)
eggCounts_old <- get(load(paste0(path, file1)))
eggCounts_new <- get(load(paste0(path, file2)))
eggCounts_untreated <- get(load(paste0(path, file3)))


###############################################################################################
# determine MHI infections in untreated scenario
###############################################################################################

df.mhi.emc.untreat <- getMHI(eggCounts_untreated, mediumHeavyInfectionThreshold_counts, "untreated", results)

###############################################################################################
# determine MHI infections in untreated scenario
###############################################################################################

df.mhi.emc.old <- getMHI(eggCounts_old, mediumHeavyInfectionThreshold_counts, "old", results)

###############################################################################################
# determine MHI infections in untreated scenario
###############################################################################################

df.mhi.emc.new <- getMHI(eggCounts_new, mediumHeavyInfectionThreshold_counts, "new", results)

################################################################################################
# combine data frames and plot
################################################################################################

df.mhi.emc.20_50 <- rbind(df.mhi.emc.untreat, df.mhi.emc.old, df.mhi.emc.new)
df.mhi.emc.20_50$scenario <- factor(df.mhi.emc.20_50$scenario, levels=c("new", "old", "untreated"))

p2 <- ggplot(df.mhi.emc.20_50)
p2 <- p2 + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=24), axis.text=element_text(size=20), legend.title=element_text(size=20), legend.text=element_text(size=20))
p2 <- p2 + geom_ribbon(aes(x=time, ymin=p5, ymax=p95, fill=scenario), alpha=0.3) + scale_fill_manual(values=c("red", "blue", "seagreen"))
p2 <- p2 + geom_line(aes(x=time, y=meanMHICounts, colour=scenario), size=1, alpha=0.8) + scale_colour_manual(values=c("red", "blue", "seagreen"))
p2 <- p2 + scale_x_continuous(name="Age (years)", expand=c(0,0), limits=c(0, 72), breaks=seq(0, 100, 5))
p2 <- p2 + scale_y_continuous(name="% women with MHI", expand=c(0, 0)) + coord_cartesian(ylim=c(0,21))
print(p2)

pic_mhi_means <- plot_grid(p1, p2, labels=c('a', 'b'), label_size=22, ncol=2, hjust=0, vjust=1.0)
#save_plot(file=paste0(outpath, "pc95MHI", stub, ".pdf"), pic_mhi_means, base_height=8.5, base_width=19)




