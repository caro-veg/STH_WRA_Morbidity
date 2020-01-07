#############################################################################################
# CODE TO SIMULATE A COHORT OF WOMEN OVER 50/100 (?) YEARS IN AN AREA WHERE STH ARE ENDEMIC
# AN AGE-SPECIFIC FORCE OF INFECTION (FOI) ACTS ON THE COHORT EVERY YEAR
# TREATMENT SCENARIOS:
# PSAC/SAC (2-15 YEARS): 75% COVERAGE, TREATMENT 1x OR 2x PER YEAR DURING SCHOOL-BASED 
# DEWORMING, DEPENDING ON PREVALENCE MEASURED BY KATO-KATZ (KK)
# ADOLESCENT GIRLS (15-18 YEARS): 75% COVERAGE, TREATMENT 1x PER YEAR DURING SCHOOL-BASED
# HPV VACCINE PROGRAMME
# WRA (19-49 YEARS): DURING PREGNANCY WOMEN GET DEWORMED 2x WHEN VISITING MATERNITY HEALTH
# CENTRES, DURING LACTATION WOMEN GET DEWORMED 1x, ASSUME 4 PREGNANCIES PER WOMAN OVER THE 
# REPRODUCTIVE PERIOD AND 80% COVERAGE
# WOMEN 50+ YEARS: NO TREATMENT
#############################################################################################

library(tidyr)

# number of women in cohort
n <- 500

# aggregation parameter 
k <- 0.35

# worm death rate 
sigma <- 0.5


# construct cohort of n women starting at age 0
worms <- rep(0, n)
pregnant <- rep(FALSE, n)
yearAfterPregnancy <- rep(FALSE, n)
treated <- rep(FALSE, n)

cohort <- data.frame(worms=worms, pregnant=pregnant, yearAfterPregnancy=yearAfterPregnancy)


# read in FoI values and calculate mean for every year
path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\L_50-70 and 20-50_KK\\"
file <- "L_20-50KK_pc_treat1x.csv"
data <- read.csv(paste0(path, file), header=TRUE)
sub <- data[, c("iteration", "time", "FoI")]
subWide <- spread(sub, iteration, FoI)
FoI <- rowMeans(subWide[, -1])


# matrix to record worm number for each woman in cohort every half a year
records <- matrix(0, nrow=length(FoI), ncol=n)


# simulate cohort for 50/100(?) years
# in half-year steps from age 0-1.5 years 
# no treatment
for(i in 1:4)
{
	meanWorms <- (FoI[i] - sigma * mean(cohort$worms)) * 0.5

	newWorms <- rnbinom(500, size=k, mu=meanWorms)
	cohort$worms <- cohort$worms + newWorms	
}

# in half-year steps from age 2-14.5 years 
# treatment 1x or 2x per year, 75% coverage
for(i in 5:30)
{
	# worm life-cycle
	meanWorms <- (FoI[i] - sigma * mean(cohort$worms)) * 0.5
	newWorms <- rnbinom(500, size=k, mu=meanWorms)
	cohort$worms <- cohort$worms + newWorms	

	# treatment 
	treated <- rbinom(500, size=1, prob=0.75) 
	wormsTreated <- cohort$worms * treated
	wormsToDie <- rbinom(n, wormsTreated, 0.95)
	cohort$worms <- cohort$worms - wormsToDie
}

# in half-year steps from age 15-18.5 years 
# treatment 1x per year, 75% coverage
for(i in 31:38)
{

}

# in half-year steps from age 19-49.5 years 
# treatment 2x per year if pregnant and 1x per year after pregnancy
# 80% coverage
# assume 4 pregnancies on average per woman
for(i in 39:101)
{
	
}

# in half-year steps from age 50-100 years
# no treatment
for(i in 102:length(FoI))
{

}






