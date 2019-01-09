data {
	//Field level
	int Nfield; //Number of fields	
	int numHives[Nfield]; //Number of hives present
	int is2015[Nfield]; //Was field measured in 2015?
	int isGP[Nfield]; //Was field from Grand Prairie?
	int isIrrigated[Nfield]; //Was field irrigated?
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
	vector[Nplot] dist; //Distance from edge
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	int flyVis[Nplot]; //Number of fly visits per plot
	vector[Nplot] totalTime; //Minutes taken for observation/10	
	vector[Nplot] flDens; //Flower density
	//Plant density (stems/m2) - missing and observed
	int Nplot_densObs; //Number of plots where plant density was observed
	int Nplot_densMiss; //Number of plots with plant density missing 
	vector[Nplot_densObs] plDens_obs; //Observed plant density (stems/m2)
	int<lower=1,upper=Nplot> obsPlDens_ind[Nplot_densObs]; //Index for observed plant density
	int<lower=1,upper=Nplot> missPlDens_ind[Nplot_densMiss]; //Index for missing plant density		
	
	//Flower level
	int Nflw; //Number of flowers
	int<lower=1,upper=Nplot>  flowerIndex[Nflw];  //Index for flowers - which plot?
	int pollenCount[Nflw]; //Pollen count

	//Plant level
	int Nplant; //Number of all plants (some measurements missing)	
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?		
	vector[Nplant] plantSize; //Centered mass of vegetative tissue (no seeds) (g)
	vector[Nplant] yield; //Observed yield (g seed per plant)	
	
	//Pod level
	int Npod; //Number of pods
	int seedCount[Npod];  //Seeds per pod
	vector[Npod] seedMass; //Mass of seedCount seeds (g)
	int<lower=1,upper=Nplant> podIndex[Npod]; //Index for pods - which plant?	
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives	
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	vector[Nplot] logFlyVis; //Fly visitation rate
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)	
	vector[Nplant] logFlwCount; //Log flower count
	vector[Npod] logSeedCount; //Log seed count
	
	logHbeeDist=logHbeeDist-mean(logHbeeDist); //Centers distance
	
	for(i in 1:Nfield){
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives		
	}
	
	for(i in 1:Nplot){ //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5); 
		logFlyVis[i] = log((flyVis[i]/totalTime[i])+0.5);
	}
		
	for(i in 1:Nplant){
		logFlwCount[i] = log(flwCount[i]); //Log-transforms flower count per plant
	}
	
	for(i in 1:Npod){
		logSeedCount[i] = log(seedCount[i]); //Log-transforms seed count per pod
	}
	
}

parameters {
	//Claim: logVisRate~plantDens+...
	real slopePlDensVis;

	// Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slope2015PlDens; //Effect of 2015
	real slopeIrrigPlDens; //Effect of irrigation
	real slope2015IrrigPlDens; //Year:irrigation interaction	
	real slopeGPPlDens; //GP effect on plant density	
	real slopeDistPlDens; //Slope of distance into field		
	real<lower=0.01> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field
		
	// hbee Visitation per plot	
	real intVisit; //Intercept 
	real slopeYearVis; //Effect of 2015 
	real slopeGpVis; //Effect of Grand Prairie 
	real slopeYearGpVis; //Effect of year-GP interaction
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real slopeFlDens; //Slope of flower density
	real slopeIrrigVis; //Effect of irrigation
	real<lower=0> sigmaVisField; //SD of field random intercepts
	real<lower=0> visitHbeePhi; //Dispersion parameter	
	vector[Nfield] intVisit_field; //field-level random intercepts		
	real<lower=0> lambdaVisField; //Lambda for skewed random effects	

}

transformed parameters {		
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density	
	vector[Nplot] visitHbeeMu; //Plot-level hbee visits		
	
	// Imputed missing data;
	vector[Nplot] plDens; //Planting density	
	plDens[obsPlDens_ind]=plDens_obs;
	plDens[missPlDens_ind]=plDens_miss;	
	
	for(i in 1:Nplot){ 
		// Plant density per plot
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex[i]] + 
			slope2015PlDens*is2015[plotIndex[i]]+ //Year effect
			slopeIrrigPlDens*isIrrigated[plotIndex[i]]+ //Irrigation effect
			slope2015IrrigPlDens*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]+ //Year:irrigation interaction
			slopeDistPlDens*logHbeeDist[i] + //Distance effect				
			slopeGPPlDens*isGP[plotIndex[i]]; //GP effect
			
		// Hbee Visitation = intercept + random int + time offset + 	
		visitHbeeMu[i] = intVisit + intVisit_field[plotIndex[i]] + logTime[i] + 
			slopeYearVis*is2015[plotIndex[i]] + //Year effect
			slopeGpVis*isGP[plotIndex[i]] + //Grand Prairie effect
			slopeYearGpVis*is2015[plotIndex[i]]*isGP[plotIndex[i]] + //Year:area interaction			
			slopeDistVis*logHbeeDist[i] + //distance from edge 
			slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			slopeFlDens*flDens[i] + //Flower density
			slopeIrrigVis*isIrrigated[plotIndex[i]] + //Irrigation - not sure this should be here, but leaving it for now
			slopePlDensVis*plDens[i]; //Claim
	}
}
	
model {	
	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density	
	hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
			
	//Priors
	//Claim:
	slopePlDensVis ~ normal(0,1);
	
	//Plant density	- informative priors
	intPlDens ~ normal(0,0.1); //Global intercept
	slope2015PlDens ~ normal(0,1); //Year effect
	slopeIrrigPlDens ~ normal(0,1); //Irrigation effect
	slope2015IrrigPlDens ~ normal(0,1); //Year:irrigation interaction
	slopeDistPlDens ~ normal(0,0.05); //Slope of distance into field	
	slopeGPPlDens ~ normal(0,1); // Grand Prairie effect 
	sigmaPlDens ~ gamma(3,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(3,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field
			
	// Visitation - informative priors
	intVisit ~ normal(-1,1); //Intercept	
	slopeDistVis ~ normal(-0.3,0.3); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0.6,0.5); //Slope of hive effect on visits 
	slopeYearVis ~ normal(0.5,1); //2015 effect
	slopeGpVis ~ normal(0,1); //Effect of Grand Prairie	
	slopeYearGpVis ~ normal(0,1); // GP-year interaction
	slopeFlDens ~ normal(0,0.5); //Flower density 
	sigmaVisField ~ gamma(4,4); //Sigma for random field 		
	intVisit_field ~ exp_mod_normal(0,sigmaVisField,lambdaVisField); //Skewed random effects - slightly better than standard normal
	lambdaVisField ~ gamma(4,2); //Lambda for skewed random effects			
	visitHbeePhi ~ gamma(4,10); //Dispersion parameter for NegBin		
}

generated quantities {	
}
