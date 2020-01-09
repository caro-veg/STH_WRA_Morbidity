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

runtime = proc.time()

# number of women in cohort
n <- 500

# aggregation parameter 
k <- 0.35

# worm death rate (per year)
sigma <- 0.5

# read in FoI values and calculate mean for every year
path <- "D:\\STH\\ModellingConsortium\\WRAquestion\\L_50-70 and 20-50_KK_100years\\"
file <- "L_20-50KK_pc_treat1x.csv"
data <- read.csv(paste0(path, file), header=TRUE)
sub <- data[, c("iteration", "time", "FoI")]
subWide <- spread(sub, iteration, FoI)
FoI <- rowMeans(subWide[, -1])

# treatment frequency per year
twice = FALSE

resultsList <- list()

for(k in 1:100){

# construct cohort of n women starting at age 0
maleWorms <- rep(0, n)
femaleWorms <- rep(0, n)
pregnant <- rep(0, n)
yearAfterPregnancy <- rep(0, n)

cohort <- data.frame(maleWorms=maleWorms, femaleWorms=femaleWorms, pregnant=pregnant, yearAfterPregnancy=yearAfterPregnancy)


# matrices to record number of male ane female worms for each woman in cohort every half a year
recordMale <- matrix(0, nrow=length(FoI), ncol=n)
recordFemale <- matrix(0, nrow=length(FoI), ncol=n)
recordFertilisedFemales <- matrix(0, nrow=length(FoI), ncol=n)


# simulate cohort for 50/100(?) years
# in half-year steps from age 0-1.5 years 
# no treatment
for(i in 1:4)
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	deviates <- rgamma(n, shape=k, scale=1/k)
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms
	
	# store worm number
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms
}

# in half-year steps from age 2-14.5 years 
# treatment 1x or 2x per year, 75% coverage
for(i in 5:30)
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	deviates <- rgamma(n, shape=k, scale=1/k)
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms

	# store worm number before treatment
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms


	# treatment 1x or 2x per year
	if(twice==TRUE)
	{
		treated <- rbinom(500, size=1, prob=0.75) 
		maleWormsTreated <- cohort$maleWorms * treated
		femaleWormsTreated <- cohort$femaleWorms * treated
		maleWormsToDie <- rbinom(n, maleWormsTreated, 0.95)
		femaleWormsToDie <- rbinom(n, femaleWormsTreated, 0.95)
		cohort$maleWorms <- cohort$maleWorms - maleWormsToDie
		cohort$femaleWorms <- cohort$femaleWorms - femaleWormsToDie
	}
	else
	{
		if(i%%2==0)
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

# in half-year steps from age 15-18.5 years 
# treatment 1x per year, 75% coverage
for(i in 31:38)
{
	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	deviates <- rgamma(n, shape=k, scale=1/k)
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms

	# store worm number before treatment
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms


	# treat 1x per year
	if(i%%2==0)
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

# in half-year steps from age 19-49.5 years 
# treatment 2x per year if pregnant and 1x per year after pregnancy
# 80% coverage
# assume 4 pregnancies on average per woman --> 
# i.e. 0.1333 probability per year of getting pregnant or 0.06667 per half year
for(i in 39:101)
{
	# reset pregnancy status and determine who is lactating
	cohort$yearAfterPregnancy[cohort$pregnancy == 1] <- 2
	cohort$yearAfterPregnancy[cohort$pregnancy == 2] <- 1
	cohort$pregnant[cohort$pregnant == 2] <- 0
	cohort$pregnant[cohort$pregnant == 1] <- 2

	# determine who gets pregnant 
	notPregnant <- which(cohort$pregnant==0 & cohort$yearAfterPregnancy==0)
	pregnant <- rbinom(notPregnant, size=1, 0.06667)
	cohort$pregnant[notPregnant] <- pregnant

	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	deviates <- rgamma(n, shape=k, scale=1/k)
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$femaleWorms <- femaleWorms

	# store worm number before treatment
	recordMale[i,] <- cohort$maleWorms
	recordFemale[i,] <- cohort$femaleWorms

	# treat 2x per year if pregnant, 1x if lactating, otherwise not
	for(j in 1:n)
	{
		if(cohort$pregnant[j] > 0 | cohort$yearAfterPregnancy[j] == 1)
		{
			treated <- rbinom(1, size=1, prob=0.8) 
			maleWormsTreated <- cohort$maleWorms[j] * treated
			femaleWormsTreated <- cohort$femaleWorms[j] * treated
			maleWormsToDie <- rbinom(1, maleWormsTreated, 0.95)
			femaleWormsToDie <- rbinom(1, femaleWormsTreated, 0.95)
			cohort$maleWorms[j] <- cohort$maleWorms[j] - maleWormsToDie
			cohort$femaleWorms[j] <- cohort$femaleWorms[j] - femaleWormsToDie
		}
	}
}


# in half-year steps from age 50-100 years
# no treatment
for(i in 102:length(FoI))
{
	# left-over pregnancies from WRA stage get treated
	if(i < 105)
	{
		# treat 2x per year if pregnant, 1x if lactating, otherwise not
		for(j in 1:n)
		{
			if(cohort$pregnant[j] > 0 | cohort$yearAfterPregnancy[j] == 1)
			{
				treated <- rbinom(1, size=1, prob=0.8) 
				maleWormsTreated <- cohort$maleWorms[j] * treated
				femaleWormsTreated <- cohort$femaleWorms[j] * treated
				maleWormsToDie <- rbinom(1, maleWormsTreated, 0.95)
				femaleWormsToDie <- rbinom(1, femaleWormsTreated, 0.95)
				cohort$maleWorms[j] <- cohort$maleWorms[j] - maleWormsToDie
				cohort$femaleWorms[j] <- cohort$femaleWorms[j] - femaleWormsToDie
			}
		}

		# reset pregnancy/lactation status 
		cohort$yearAfterPregnancy[cohort$pregnancy == 1] <- 2
		cohort$yearAfterPregnancy[cohort$pregnancy == 2] <- 1
		cohort$pregnant[cohort$pregnant == 2] <- 0
		cohort$pregnant[cohort$pregnant == 1] <- 2
	}


	# calculate current FoI acting on each woman
	meanFoI <- FoI[i] / 2
	deviates <- rgamma(n, shape=k, scale=1/k)
	currentFoIs <- meanFoI * deviates

	# calculate burden of male worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	maleWorms <- rbinom(n, cohort$maleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
	cohort$maleWorms <- maleWorms

	# calculate burden of female worms for each woman
	M <- currentFoIs / sigma * (1 - exp(-sigma * (i-1) * 0.5))
	femaleWorms <- rbinom(n, cohort$femaleWorms, prob=exp(-sigma * (i-1) * 0.5)) + rpois(n, lambda=M)
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
temp <- rep(k, nrow(recordFertilisedFemales))
recordFertilisedFemales <- cbind(temp, recordFertilisedFemales)
colnames(recordFertilisedFemales) <- c("rep", paste0("W", seq(1, n)))

resultsList[[k]] <- recordFertilisedFemales

}

results <- do.call(rbind, resultsList)

write.table(results, file=paste0(path, "fertilisedFemaleWorms.csv"), sep=",", row.names=FALSE, col.names=TRUE)

runtime = proc.time() - runtime
print(runtime)



