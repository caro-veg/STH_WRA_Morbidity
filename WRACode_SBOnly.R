#############################################################################################
# CODE TO SIMULATE A COHORT OF WOMEN OVER 50/100 (?) YEARS IN AN AREA WHERE STH ARE ENDEMIC
# AN AGE-SPECIFIC FORCE OF INFECTION (FOI) ACTS ON THE COHORT EVERY YEAR
# TREATMENT SCENARIO:
# PSAC/SAC (2-15 YEARS): 75% COVERAGE, TREATMENT 1x OR 2x PER YEAR DURING SCHOOL-BASED 
# DEWORMING, DEPENDING ON PREVALENCE MEASURED BY KATO-KATZ (KK)
# WOMEN 15+ YEARS: NO TREATMENT
#
# OUTPUT: MATRIX OF NUMBERS OF FERTILISED FEMALE WORMS, ROWS: TIMESTEPS, COLS: WOMEN
# MATRICES FOR ITERATIONS ARE APPENDED AND MARKED BY ITERATION NUMBER (FIRST COLUMN)
#############################################################################################

library(tidyr)
library(doParallel)

set.seed(234)

#############################################################################################
# clean up before new run
rm(list=ls())
gc()

closeAllConnections()

#############################################################################################


## command line arguments. 
args <- commandArgs(trailingOnly = TRUE)

path <- args[1]	# path to file with FoIs
file <- args[2]	# file with FoIs
twice <- args[3]	# treatment frequency per year

setwd(path)
twice <- as.logical(twice)


runtime = proc.time()

# number of women in cohort
n <- 500
# aggregation parameter 
k <- 0.35
# worm death rate (per year)
sigma <- 0.5

stub <- substring(file, first=2, last=9)

if(twice==FALSE)
{
	stub <- paste0(stub, "_treat1x")
}else{
	stub <- paste0(stub, "_treat2x")
}

# read in FoI values and calculate mean for every year
data <- read.csv(file, header=TRUE)
sub <- data[, c("iteration", "time", "FoI")]
subWide <- spread(sub, iteration, FoI)
FoI <- rowMeans(subWide[, -1])

resultsList <- list()

deltaT <- 1/12

#################################################################################################################
# REALISATION FUNCTION
#################################################################################################################

doRealisation <- function(n, k, sigma, twice, FoI)
{

# construct cohort of n women starting at age 0
maleWorms <- rep(0, n)
femaleWorms <- rep(0, n)
pregnant <- rep(0, n)
yearAfterPregnancy <- rep(0, n)

cohort <- data.frame(maleWorms=maleWorms, femaleWorms=femaleWorms, pregnant=pregnant, yearAfterPregnancy=yearAfterPregnancy)

deviates <- rgamma(n, shape=k, scale=1/k)

# matrices to record number of male ane female worms for each woman in cohort every half a year
recordMale <- matrix(0, nrow=length(FoI), ncol=n)
recordFemale <- matrix(0, nrow=length(FoI), ncol=n)
recordFertilisedFemales <- matrix(0, nrow=length(FoI), ncol=n)


# simulate cohort for 50/100(?) years
# in half-year steps from age 0 - <2 years 
# no treatment
for(i in 1:24)
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2	
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * deltaT))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * deltaT)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * deltaT))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * deltaT)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms
	
	# store worm number
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms
}

# in half-year steps from age 2 - <15 years 
# treatment 1x or 2x per year, 75% coverage
for(i in 25:179)
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * deltaT))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * deltaT)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * deltaT))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * deltaT)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms

	# store worm number before treatment
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms

	# treatment 1x or 2x per year
	if(twice==TRUE)
	{
		if(i%%6==0)
		{
			treated <- rbinom(500, size=1, prob=0.75) 
			maleWormsTreated <- cohort$maleWorms * treated
			femaleWormsTreated <- cohort$femaleWorms * treated
			maleWormsToDie <- rbinom(n, maleWormsTreated, 0.95)
			femaleWormsToDie <- rbinom(n, femaleWormsTreated, 0.95)
			cohort$maleWorms <- cohort$maleWorms - maleWormsToDie
			cohort$femaleWorms <- cohort$femaleWorms - femaleWormsToDie
		}

	}else
	{
		if(i%%12==0)
		{
			treated <- rbinom(500, size=1, prob=0.75) 
			maleWormsTreated <- cohort$maleWorms * treated
			femaleWormsTreated <- cohort$femaleWorms * treated
			maleWormsToDie <- rbinom(n, maleWormsTreated, 0.95)
			femaleWormsToDie <- rbinom(n, femaleWormsTreated, 0.95)
			cohort$maleWorms <- cohort$maleWorms - maleWormsToDie
			cohort$femaleWorms <- cohort$femaleWorms - femaleWormsToDie
		}
	}
}

# in half-year steps from age 15-100 years 
# no treatment
for(i in 180:length(FoI))
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * deltaT))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * deltaT)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * deltaT))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * deltaT)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms

	# store worm number
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms
}

# determine number of fertilised female worms for each woman at each timestep
for(i in 1:nrow(recordMale))
{
	for(j in 1:ncol(recordMale))
	{
		if(recordMale[i, j] > 0)
		{
			recordFertilisedFemales[i, j] <- recordFemale[i, j]
		}
		else
		{
			recordFertilisedFemales[i, j] <- 0
		}
	}
}

# store results in list
colnames(recordFertilisedFemales) <- c(paste0("W", seq(1, n)))

return(recordFertilisedFemales)

}

cl <- makeCluster(2)
registerDoParallel(cl)

iterations <- 100
resultsList <- foreach(i=1:iterations) %do% doRealisation(n, k, sigma, twice, FoI)

stopImplicitCluster() 

results <- do.call(rbind, resultsList)

reps <- c()
for(i in 1:iterations)
{
	reps <- c(reps, rep(i, 1201))
}

results <- cbind(reps, results)

save(file=paste0(path, "fertilisedFemaleWorms_old", stub, ".RData"), results)

runtime = proc.time() - runtime
print(runtime)


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

#lambda <- 3
#gamma <- 0.02
#k_epg <- 0.35


# get egg counts based on number of female worms
#eggCounts <- apply(results[, -1], c(1,2), getSetOfEggCounts, lambda, gamma, k_epg)

#save(file=paste0(path, "eggCounts_old", stub, ".RData"), results)

#runtime = proc.time() - runtime
#print(runtime)



