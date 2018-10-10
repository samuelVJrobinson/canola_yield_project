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
	int Nplant_obs; //Number of "complete" (observed) plants
	int Nplant_miss; //Number of missing plants		
	vector[Nplant_obs] plantSize_obs; //Mass of vegetative tissue (no seeds) (g)
	// vector[Nplant_obs] totalSeedMass; //Mass of seeds from plant (g)
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?		
	int<lower=1,upper=Nplant> obsPl_ind[Nplant_obs]; //Index for observed plants
	int<lower=1,upper=Nplant> missPl_ind[Nplant_miss]; //Index for missing plants	
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	
	//Pod level
	int Npod; //Number of pods
	int seedCount[Npod];  //Seeds per pod
	vector[Npod] seedMass; //Mass of seedCount seeds (g)
	int<lower=1,upper=Nplant> podIndex[Npod]; //Index for pods - which plant?	
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distances	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	for(i in 1:Nfield)
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives
	
	for(i in 1:Nplot) //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5); 
}

parameters {
	//Claim
	real slopeFlDensVis;
					
	//hbee Visitation
	real intVisit; //Intercept 
	real slopeYearVis; //Effect of 2015 
	real slopeGpVis; //Effect of Grand Prairie 
	real slopeYearGpVis; //Effect of year-GP interaction
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real slopeFlDens; //Slope of flower density
	real<lower=0> sigmaVisField; //SD of field random intercepts
	real<lower=0> visitHbeePhi; //Dispersion parameter	
	vector[Nfield] intVisit_field; //field-level random intercepts		
}

transformed parameters {		
	//Expected values	
	vector[Nplot] visitHbeeMu; //Plot-level hbee visits		
	
	for(i in 1:Nplot){ 	
		// Hbee Visitation = intercept + random int + time offset + 	
		visitHbeeMu[i] = intVisit + intVisit_field[plotIndex[i]] + logTime[i] + 
			slopeYearVis*is2015[plotIndex[i]] + //Year effect
			slopeGpVis* +isGP[plotIndex[i]] + //Grand Prairie effect
			slopeYearGpVis*is2015[plotIndex[i]]*isGP[plotIndex[i]] + //Year:area interaction			
			slopeDistVis*logHbeeDist[i] + //distance from edge 
			slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			slopeFlDensVis*flDens[i]; //Flower density		
	}	
}
	
model {		
	//Likelihood
	hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)	
		
	//Priors	
	//Claim
	slopeFlDensVis ~ normal(0,1);
	
	//Visitation - informative priors
	intVisit ~ normal(-1,1); //Intercept	
	slopeDistVis ~ normal(-0.3,0.3); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0.6,0.5); //Slope of hive effect on visits 
	slopeYearVis ~ normal(0.5,1); //2015 effect
	slopeGpVis ~ normal(0,1); //Effect of Grand Prairie	
	slopeYearGpVis ~ normal(0,1); // GP-year interaction
	slopeFlDens ~ normal(0,5); //Flower density 
	sigmaVisField ~ gamma(4,2); //Sigma for random field 	
	intVisit_field ~ normal(0,sigmaVisField); //Random field int
	visitHbeePhi ~ gamma(2,2); //Dispersion parameter			
}

generated quantities{
//Plot-level quantities	
	//hbeeVis
	int predHbeeVis[Nplot]; //Generated
	real hbeeVis_resid[Nplot]; //Residual 
	real predHbeeVis_resid[Nplot]; //Residual of generated 			
		
	for(i in 1:Nplot){
		// bee visits (NB version)
		hbeeVis_resid[i]=hbeeVis[i]-(exp(visitHbeeMu[i])); //Residual for actual value
		predHbeeVis[i] = neg_binomial_2_log_rng(visitHbeeMu[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		predHbeeVis_resid[i]=predHbeeVis[i]-(exp(visitHbeeMu[i])); //residual for predicted value					
	}			
}
