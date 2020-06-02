#############################################################################################
# CODE TO SIMULATE A COHORT OF WOMEN OVER 50/100 (?) YEARS IN AN AREA WHERE STH ARE ENDEMIC
# AN AGE-SPECIFIC FORCE OF INFECTION (FOI) ACTS ON THE COHORT EVERY YEAR
# TREATMENT SCENARIO:
# PSAC/SAC (2-15 YEARS): 75% COVERAGE, TREATMENT 1x OR 2x PER YEAR DURING SCHOOL-BASED 
# DEWORMING, DEPENDING ON PREVALENCE MEASURED BY KATO-KATZ (KK)
# ADOLESCENT GIRLS (15-18 YEARS): 75% COVERAGE, TREATMENT 1x PER YEAR DURING SCHOOL-BASED
# HPV VACCINE PROGRAMME
# WRA (19-49 YEARS): DURING PREGNANCY WOMEN GET DEWORMED 2x WHEN VISITING MATERNITY HEALTH
# CENTRES, DURING LACTATION WOMEN GET DEWORMED 1x, ASSUME 4 PREGNANCIES PER WOMAN OVER THE 
# REPRODUCTIVE PERIOD AND 80% COVERAGE
# WOMEN 50+ YEARS: NO TREATMENT
#
# OUTPUT: MATRIX OF NUMBERS OF FERTILISED FEMALE WORMS, ROWS: TIMESTEPS, COLS: WOMEN
# MATRICES FOR ITERATIONS ARE APPENDED AND MARKED BY ITERATION NUMBER (FIRST COLUMN)
#############################################################################################

library(tidyr)
library(doParallel)
library(doRNG)


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


seed <- 234
set.seed(seed)


runtime = proc.time()

# number of women in cohort
n <- 500
# aggregation parameter 
k <- 0.35
# worm death rate (per year)
sigma <- 0.5
# drug efficacy
efficacy <- 0.95


stub <- substring(file, first=2, last=9)


if(twice==FALSE)
{
	stub <- paste0(stub, "_treat1x")
}else{
	stub <- paste0(stub, "_treat2x")
}


data <- read.csv(file, header=TRUE)
sub <- data[, c("iteration", "time", "FoI")]
subWide <- spread(sub, iteration, FoI)
FoI <- rowMeans(subWide[, -1])

deltaT = 1/12

resultsList <- list()


#################################################################################################################
# REALISATION FUNCTION
#################################################################################################################

doRealisation <- function(n, k, sigma, twice, FoI, deltaT)
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
# in monthly steps from age 0 - <2 years 
# no treatment
for(i in 1:24)
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma  * deltaT))
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

# in monthly steps from age 2 - <15 years 
# treatment 1x or 2x per year, 75% coverage
for(i in 25:180)
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
			maleWormsToDie <- rbinom(n, maleWormsTreated, efficacy)
			femaleWormsToDie <- rbinom(n, femaleWormsTreated, efficacy)
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
			maleWormsToDie <- rbinom(n, maleWormsTreated, efficacy)
			femaleWormsToDie <- rbinom(n, femaleWormsTreated, efficacy)
			cohort$maleWorms <- cohort$maleWorms - maleWormsToDie
			cohort$femaleWorms <- cohort$femaleWorms - femaleWormsToDie
		}
	}
}

# in monthly steps from age 15 - < 19 years 
# treatment 1x per year, 75% coverage
for(i in 181:228)
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


	# treat 1x per year
	if(i%%12==0)
	{
		treated <- rbinom(500, size=1, prob=0.75) 
		maleWormsTreated <- cohort$maleWorms * treated
		femaleWormsTreated <- cohort$femaleWorms * treated
		maleWormsToDie <- rbinom(n, maleWormsTreated, efficacy)
		femaleWormsToDie <- rbinom(n, femaleWormsTreated, efficacy)
		cohort$maleWorms <- cohort$maleWorms - maleWormsToDie
		cohort$femaleWorms <- cohort$femaleWorms - femaleWormsToDie
	}
}

# in monthly steps from age 19 - <50 years 
# treatment 2x per year if pregnant and 1x per year after pregnancy
# 80% coverage
# assume 4 pregnancies on average per woman --> 
# i.e. 0.1333 probability per year of getting pregnant or 0.06667 per half year
for(i in 229:600)
{
	# reset pregnancy status and determine who is lactating
	cohort$yearAfterPregnancy[cohort$pregnancy == 12] <- 0
	for(m in 11:1)
	{
		cohort$yearAfterPregnancy[cohort$yearAfterPregnancy == m] <- (m + 1)
	}
	cohort$yearAfterPregnancy[cohort$pregnancy == 9] <- 1

	cohort$pregnant[cohort$pregnant == 9] <- 0
	for(m in 8:1)
	{
		cohort$pregnant[cohort$pregnant == m] <- (m + 1)
	}


	# determine who gets pregnant 
	notPregnant <- which(cohort$pregnant==0 & cohort$yearAfterPregnancy==0)
	pregnant <- rbinom(notPregnant, size=1, 0.01111)
	cohort$pregnant[notPregnant] <- pregnant

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

	# treat 2x per year if pregnant, 1x if lactating, otherwise not
	for(j in 1:n)
	{
		if(cohort$pregnant[j] == 4 | cohort$pregnant[j] == 7 | cohort$yearAfterPregnancy[j] == 6)
		{
			treated <- rbinom(1, size=1, prob=0.8) 
			maleWormsTreated <- cohort$maleWorms[j] * treated
			femaleWormsTreated <- cohort$femaleWorms[j] * treated
			maleWormsToDie <- rbinom(1, maleWormsTreated, efficacy)
			femaleWormsToDie <- rbinom(1, femaleWormsTreated, efficacy)
			cohort$maleWorms[j] <- cohort$maleWorms[j] - maleWormsToDie
			cohort$femaleWorms[j] <- cohort$femaleWorms[j] - femaleWormsToDie
		}
	}
}


# in half-year steps from age 50-100 years
# no treatment
for(i in 601:length(FoI))
{
	# left-over pregnancies from WRA stage get treated
	if(i < 625)
	{
		# treat 2x per year if pregnant, 1x if lactating, otherwise not
		for(j in 1:n)
		{
			if(cohort$pregnant[j] == 4 | cohort$pregnant[j] == 7 | cohort$yearAfterPregnancy[j] == 6)
			{
				treated <- rbinom(1, size=1, prob=0.8) 
				maleWormsTreated <- cohort$maleWorms[j] * treated
				femaleWormsTreated <- cohort$femaleWorms[j] * treated
				maleWormsToDie <- rbinom(1, maleWormsTreated, efficacy)
				femaleWormsToDie <- rbinom(1, femaleWormsTreated, efficacy)
				cohort$maleWorms[j] <- cohort$maleWorms[j] - maleWormsToDie
				cohort$femaleWorms[j] <- cohort$femaleWorms[j] - femaleWormsToDie
			}
		}

		# reset pregnancy status and determine who is lactating
		cohort$yearAfterPregnancy[cohort$pregnancy == 12] <- 0
		for(m in 11:1)
		{
			cohort$yearAfterPregnancy[cohort$yearAfterPregnancy == m] <- (m + 1)
		}
		cohort$yearAfterPregnancy[cohort$pregnancy == 9] <- 1

		cohort$pregnant[cohort$pregnant == 9] <- 0
		for(m in 8:1)
		{
			cohort$pregnant[cohort$pregnant == m] <- (m + 1)
		}
	}


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
#registerDoRNG(seed)

iterations <- 100
resultsList <- foreach(i=1:iterations, .options.RNG=seed) %dorng% doRealisation(n, k, sigma, twice, FoI, deltaT)

stopImplicitCluster() 

results <- do.call(rbind, resultsList)

reps <- c()
for(i in 1:iterations)
{
	reps <- c(reps, rep(i, 1201))
}

results <- cbind(reps, results)

save(file=paste0(path, "fertilisedFemaleWorms_new", stub, ".RData"), results)

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

#save(file=paste0(path, "eggCounts_new", stub, ".RData"), results)

#runtime = proc.time() - runtime
#print(runtime)





