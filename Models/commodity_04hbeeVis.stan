//Commodity canola pollination model used for publication
//  Uses plant-level data

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
	vector[Nplant] seedCount; //Average seeds per pod
	vector[Nplant] seedMass; //Average weight per seed (mg)
	vector[Nplant] yield; //Observed yield per plant (g)	
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives	
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	vector[Nplot] logFlyVis; //Fly visitation rate
	vector[Nplant] logitFlwSurv; //(logit) proportion flower survival	
	vector[Nplant] logFlwCount; //Log flower count
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)
	vector[Nplant] calcYield; //Predicted yield per plant (g seed per plant)	
	vector[Nplant] logCalcYield; //Predicted log yield
	
	logHbeeDist=logHbeeDist-mean(logHbeeDist); //Centers distance
	
	for(i in 1:Nfield){
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives		
	}
	
	for(i in 1:Nplot){ //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+1); 
		logFlyVis[i] = log((flyVis[i]/totalTime[i])+1);
	}
		
	for(i in 1:Nplant){
		logFlwCount[i] = log(flwCount[i]); //Log flower count per plant	
		//Necessary for promoting integers to reals. Otherwise does integer division.
		logitFlwSurv[i] = podCount[i]; 
		logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods		
		if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
			logitFlwSurv[i]=0.01;
		else if(logitFlwSurv[i]>=1)
			logitFlwSurv[i]=0.99;	
		//calcYield = log(pod count x seed weight x seed count)
		calcYield[i] = podCount[i]*seedCount[i]*(seedMass[i]/1000);
		logCalcYield[i] = log(calcYield[i]);
	}		
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
	
}

parameters {
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
	
	//Visitation
	vector[Nplot] visitHbeeMu; //Expected hbee visits
	
	for(i in 1:Nplot){ 
		// Honeybee Visitation
		visitHbeeMu[i] = intVisit + //Intercept
		  intVisit_field[plotIndex[i]] + //Field-level random intercept
		  logTime[i] + //Time offset
			slopeYearVis*is2015[plotIndex[i]] + //Year effect
			slopeGpVis*isGP[plotIndex[i]] + //Grand Prairie effect
			slopeYearGpVis*is2015[plotIndex[i]]*isGP[plotIndex[i]] + //Year:area interaction
			slopeDistVis*logHbeeDist[i] + //distance from edge
			slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			slopeFlDens*flDens[i] + //Flower density
			slopeIrrigVis*isIrrigated[plotIndex[i]]; //Irrigation
	}
}
	
model {	
  
  hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
	
	// Priors
	// Visitation - informative priors
	intVisit ~ normal(0,1); //Intercept
	slopeDistVis ~ normal(0,1); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0,1); //Slope of hive effect on visits
	slopeYearVis ~ normal(0,1); //2015 effect
	slopeGpVis ~ normal(0,1); //Effect of Grand Prairie
	slopeYearGpVis ~ normal(0,1); // GP-year interaction
	slopeFlDens ~ normal(0,1); //Flower density
	sigmaVisField ~ gamma(1.1,1); //Sigma for random field
	intVisit_field ~ exp_mod_normal(0,sigmaVisField,lambdaVisField); //Skewed random effects - slightly better than standard normal
	lambdaVisField ~ gamma(1.1,1); //Lambda for skewed random effects
	visitHbeePhi ~ gamma(1.1,1); //Dispersion parameter for NegBin
}

generated quantities {
	// hbeeVis
	int predHbeeVis[Nplot];
	real hbeeVis_resid[Nplot];
	real predHbeeVis_resid[Nplot];
	for(i in 1:Nplot){
		// bee visits (NB version)
		hbeeVis_resid[i]=hbeeVis[i]-exp(visitHbeeMu[i]); //Residual for actual value
		predHbeeVis[i] = neg_binomial_2_log_rng(visitHbeeMu[i],visitHbeePhi); //Predicted value drawn from neg.bin
		predHbeeVis_resid[i]=predHbeeVis[i]-exp(visitHbeeMu[i]); //residual for predicted value
	}
}
