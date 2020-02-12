
library(ggplot2)

##################################################################################
# Standard deviation of a ratio of independent variables
##################################################################################

sdRatio <- function(mean1, mean2, sd1, sd2)
{
	sdr <- mean1/mean2 * sqrt(sd1^2 / mean1^2 + sd2^2 / mean2^2)
	return(sdr)
}

sdSum <- function(sd1, sd2)
{
	return(sqrt(sd1^2 + sd2^2))
}

##################################################################################
# Moderate prevalence, annual SB
##################################################################################

meanSB_ICL <- 0.0394
meanWRA_ICL <- 0.0234

sdSB_ICL <- 0.00288
sdWRA_ICL <- 0.00184


meanSB_EMC <- 0.0571 
meanWRA_EMC <- 0.0458

sdSB_EMC <- 0.00308
sdWRA_EMC <- 0.00255


#reduction_ICL <- meanSB_ICL - meanWRA_ICL
#sd_reduction_ICL <- sdSum(sdSB_ICL, sdWRA_ICL)

ratio_mod_1_ICL <- meanWRA_ICL / meanSB_ICL
sd_mod_1_ICL <- sdRatio(meanWRA_ICL, meanSB_ICL, sdWRA_ICL, sdSB_ICL)



#reduction_EMC <- meanSB_EMC - meanWRA_EMC
#sd_reduction_EMC <- sdSum(sdSB_EMC, sdWRA_EMC)

ratio_mod_1_EMC <- meanWRA_EMC / meanSB_EMC
sd_mod_1_EMC <- sdRatio(meanWRA_EMC, meanSB_EMC, sdWRA_EMC, sdSB_EMC)



##################################################################################
# Moderate prevalence, semi-annual SB
##################################################################################

meanSB_ICL <- 0.0322
meanWRA_ICL <- 0.0232

sdSB_ICL <- 0.00222
sdWRA_ICL <- 0.00197


meanSB_EMC <- 0.0574 
meanWRA_EMC <- 0.0446

sdSB_EMC <- 0.00314
sdWRA_EMC <- 0.00245



#reduction_ICL <- meanSB_ICL - meanWRA_ICL
#sd_reduction_ICL <- sdSum(sdSB_ICL, sdWRA_ICL)

ratio_mod_2_ICL <- meanWRA_ICL / meanSB_ICL
sd_mod_2_ICL <- sdRatio(meanWRA_ICL, meanSB_ICL, sdWRA_ICL, sdSB_ICL)


#reduction_EMC <- meanSB_EMC - meanWRA_EMC
#sd_reduction_EMC <- sdSum(sdSB_EMC, sdWRA_EMC)

ratio_mod_2_EMC <- meanWRA_EMC / meanSB_EMC
sd_mod_2_EMC <- sdRatio(meanWRA_EMC, meanSB_EMC, sdWRA_EMC, sdSB_EMC)



##################################################################################
# High prevalence, semi-annual SB
##################################################################################

meanSB_ICL <- 0.0798
meanWRA_ICL <- 0.0778

sdSB_ICL <- 0.00363
sdWRA_ICL <- 0.00329


meanSB_EMC <- 0.124 
meanWRA_EMC <- 0.110

sdSB_EMC <- 0.00384
sdWRA_EMC <- 0.00382



#reduction_ICL <- meanSB_ICL - meanWRA_ICL
#sd_reduction_ICL <- sdSum(sdSB_ICL, sdWRA_ICL)

ratio_high_2_ICL <- meanWRA_ICL / meanSB_ICL
sd_high_2_ICL <- sdRatio(meanWRA_ICL, meanSB_ICL, sdWRA_ICL, sdSB_ICL)


#reduction_EMC <- meanSB_EMC - meanWRA_EMC
#sd_reduction_EMC <- sdSum(sdSB_EMC, sdWRA_EMC)

ratio_high_2_EMC <- meanWRA_EMC / meanSB_EMC
sd_high_2_EMC <- sdRatio(meanWRA_EMC, meanSB_EMC, sdWRA_EMC, sdSB_EMC)


################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

#ratios <- c(ratio_mod_1_ICL, ratio_mod_2_ICL, ratio_high_2_ICL, ratio_mod_1_EMC, ratio_mod_2_EMC, ratio_high_2_EMC)*100
ratios <- c(ratio_mod_1_ICL, ratio_mod_1_EMC, ratio_mod_2_ICL, ratio_mod_2_EMC, ratio_high_2_ICL, ratio_high_2_EMC)*100
sds <- c(sd_mod_1_ICL, sd_mod_2_ICL, sd_high_2_ICL, sd_mod_1_EMC, sd_mod_2_EMC, sd_high_2_EMC)*100
scenarios <- seq(1, 12, 2)
unis <- c(rep(c("ICL", "EMC"), 3))

df <- data.frame(ratios=ratios, scenarios=scenarios, sds=sds, unis=unis)

p <- ggplot(data=df, aes(x=scenarios, y=ratios, colour=unis))
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank(), legend.text=element_text(size=16), legend.position="bottom")
p <- p + geom_point(size=4)
p <- p + geom_errorbar(aes(ymin=ratios-2*sds, ymax=ratios+2*sds))
p <- p + scale_y_continuous(breaks=seq(0, 110, 20), limits=c(0, 110)) 
p <- p + xlab("Scenarios") + ylab("Ratio MHI prevalence \nnew treatment strategy : \nold treatment strategy")
p <- p + scale_x_continuous(breaks=c(2, 6, 10), labels=c("Moderate pravalence,\nannual treatment", "Moderate pravalence,\nsemi-annual treatment", "High pravalence,\nsemi-annual treatment"))
print(p)
ggsave(paste0(outpath, "summaryFigure1.jpeg"), p, dpi=300)

###################################################################################################################


##################################################################################
# Moderate prevalence, annual SB
##################################################################################

meanSB_ICL <- 0.0394
meanWRA_ICL <- 0.0234

sdSB_ICL <- 0.00288
sdWRA_ICL <- 0.00184


meanSB_EMC <- 0.0571 
meanWRA_EMC <- 0.0458

sdSB_EMC <- 0.00308
sdWRA_EMC <- 0.00255


reduction_ICL <- meanSB_ICL - meanWRA_ICL
sd_reduction_ICL <- sdSum(sdSB_ICL, sdWRA_ICL)

ratio_mod_1_ICL <- reduction_ICL / meanSB_ICL
sd_mod_1_ICL <- sdRatio(reduction_ICL, meanSB_ICL, sd_reduction_ICL, sdSB_ICL)


reduction_EMC <- meanSB_EMC - meanWRA_EMC
sd_reduction_EMC <- sdSum(sdSB_EMC, sdWRA_EMC)

ratio_mod_1_EMC <- reduction_EMC / meanSB_EMC
sd_mod_1_EMC <- sdRatio(reduction_EMC, meanSB_EMC, sd_reduction_EMC, sdSB_EMC)



##################################################################################
# Moderate prevalence, semi-annual SB
##################################################################################

meanSB_ICL <- 0.0322
meanWRA_ICL <- 0.0232

sdSB_ICL <- 0.00222
sdWRA_ICL <- 0.00197


meanSB_EMC <- 0.0574 
meanWRA_EMC <- 0.0446

sdSB_EMC <- 0.00314
sdWRA_EMC <- 0.00245



reduction_ICL <- meanSB_ICL - meanWRA_ICL
sd_reduction_ICL <- sdSum(sdSB_ICL, sdWRA_ICL)

ratio_mod_2_ICL <- reduction_ICL / meanSB_ICL
sd_mod_2_ICL <- sdRatio(reduction_ICL, meanSB_ICL, sd_reduction_ICL, sdSB_ICL)


reduction_EMC <- meanSB_EMC - meanWRA_EMC
sd_reduction_EMC <- sdSum(sdSB_EMC, sdWRA_EMC)

ratio_mod_2_EMC <- reduction_EMC / meanSB_EMC
sd_mod_2_EMC <- sdRatio(reduction_EMC, meanSB_EMC, sd_reduction_EMC, sdSB_EMC)



##################################################################################
# High prevalence, semi-annual SB
##################################################################################

meanSB_ICL <- 0.0798
meanWRA_ICL <- 0.0778

sdSB_ICL <- 0.00363
sdWRA_ICL <- 0.00329


meanSB_EMC <- 0.124 
meanWRA_EMC <- 0.110

sdSB_EMC <- 0.00384
sdWRA_EMC <- 0.00382



reduction_ICL <- meanSB_ICL - meanWRA_ICL
sd_reduction_ICL <- sdSum(sdSB_ICL, sdWRA_ICL)

ratio_high_2_ICL <- reduction_ICL / meanSB_ICL
sd_high_2_ICL <- sdRatio(reduction_ICL, meanSB_ICL, sd_reduction_ICL, sdSB_ICL)


reduction_EMC <- meanSB_EMC - meanWRA_EMC
sd_reduction_EMC <- sdSum(sdSB_EMC, sdWRA_EMC)

ratio_high_2_EMC <- reduction_EMC / meanSB_EMC
sd_high_2_EMC <- sdRatio(reduction_EMC, meanSB_EMC, sd_reduction_EMC, sdSB_EMC)


################################################################################

outpath <- "D:\\STH\\ModellingConsortium\\WRAquestion\\Output\\"

#ratios <- c(ratio_mod_1_ICL, ratio_mod_2_ICL, ratio_high_2_ICL, ratio_mod_1_EMC, ratio_mod_2_EMC, ratio_high_2_EMC)*100
ratios <- c(ratio_mod_1_ICL, ratio_mod_1_EMC, ratio_mod_2_ICL, ratio_mod_2_EMC, ratio_high_2_ICL, ratio_high_2_EMC)*100
sds <- c(sd_mod_1_ICL, sd_mod_2_ICL, sd_high_2_ICL, sd_mod_1_EMC, sd_mod_2_EMC, sd_high_2_EMC)*100
scenarios <- seq(1, 12, 2)
unis <- c(rep(c("ICL", "EMC"), 3))

df <- data.frame(ratios=ratios, scenarios=scenarios, sds=sds, unis=unis)

p <- ggplot(data=df, aes(x=scenarios, y=ratios, colour=unis))
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.title=element_text(size=20), axis.text=element_text(size=16), axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank(), legend.text=element_text(size=16), legend.position="bottom")
p <- p + geom_point(size=4)
p <- p + geom_errorbar(aes(ymin=ratios-2*sds, ymax=ratios+2*sds))
p <- p + scale_y_continuous(breaks=seq(0, 110, 20)) + coord_cartesian(ylim=c(0,100)) 
p <- p + xlab("Scenarios") + ylab("Relative reduction in \nM&HI infection prevalence (%)")
p <- p + scale_x_continuous(breaks=c(2, 6, 10), labels=c("Moderate pravalence,\nannual treatment", "Moderate pravalence,\nsemi-annual treatment", "High pravalence,\nsemi-annual treatment"))
print(p)
ggsave(paste0(outpath, "summaryFigure2.jpeg"), p, dpi=300)





