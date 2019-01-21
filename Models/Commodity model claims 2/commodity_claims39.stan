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
	vector[Nplant] logitFlwSurv; //(logit) proportion flower survival	
	vector[Nplant] logFlwCount; //Log flower count
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)		
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
		logFlwCount[i] = log(flwCount[i]); //Log flower count per plant	
		//Necessary for promoting integers to reals. Otherwise does integer division.
		logitFlwSurv[i] = podCount[i]; 
		logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods		
		if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
			logitFlwSurv[i]=0.01;
		else if(logitFlwSurv[i]>=1)
			logitFlwSurv[i]=0.99;				
	}
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
	
	for(i in 1:Npod){
		logSeedCount[i] = log(seedCount[i]); //Log-transforms seed count per pod
	}
	
}

parameters {
	//Claim: flwCount~plantDens+...
	real slopePlDensFlwCount;

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
		
	// Flower count (per plant) 
	real intFlwCount; //Intercept
	real slopePlSizeFlwCount; //Slope of plant size	
	real slopeSurvFlwCount; //Slope of flower survival rate
	real slope2015FlwCount; //Slope of 2015 effect
	real<lower=0> phiFlwCount_field; //SD of field-level random effect	
	vector[Nfield] intFlwCount_field; //Field-level random effect	
	real intPhiFlwCount; //Intercept for sigma	
	real slopePlSizePhiFlwCount; //Effect of plant size on sigma	
	real<lower=0> sigmaPhiFlwCount_field; //Sigma for field level sigma
	vector[Nfield] intPhiFlwCount_field; //Field-level random effect for sigma	
}

transformed parameters {		
}
	
model {	
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplot] flwCountPlot; //Plot-level flowers per plant
	vector[Nplant] flwCountMu; //Expected flowers per plant
	vector[Nplant] phiFlwCount; //Phi for flowers per plant	
	
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
			
		// Flower count per plant (plot level) = intercept + random field int + random plot int
		flwCountPlot[i] = intFlwCount + intFlwCount_field[plotIndex[i]] +
			slope2015FlwCount*is2015[plotIndex[i]] + //2015 effect
			slopePlDensFlwCount*plDens[i]; //Claim
	}
		
	for(i in 1:Nplant){ //For each plant 			
		// Flower count per plant
		flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count 
			slopePlSizeFlwCount*plantSize[i] + //individual plant size
			slopeSurvFlwCount*logitFlwSurv[i]; //Flower survival rate experienced by plant
			
		// Phi (dispersion) for flower count
		phiFlwCount[i] = exp(intPhiFlwCount + intPhiFlwCount_field[plotIndex[plantIndex[i]]] + //intPhiFlwCount_plot[plantIndex[i]] + 
			slopePlSizePhiFlwCount*plantSize[i]); // Term for sigma				
	}		

	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density	
	flwCount ~ neg_binomial_2_log(flwCountMu,phiFlwCount); //Flower count per plant (attempted pods)
			
	//Priors
	slopePlDensFlwCount ~ normal(0,2); //Claim
	
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
	
	//Flower count (per plant)	
	intFlwCount ~ normal(5,1); //Intercept
	slopePlSizeFlwCount ~ normal(1,1); //Slope of plant size
	slopeSurvFlwCount ~ normal(0,1); //Slope of survival rate 
	slope2015FlwCount ~ normal(0,1); //Slope of 2015 effect
	phiFlwCount_field ~ gamma(2,10); //SD of field-level random effect		
	intFlwCount_field ~ normal(0,phiFlwCount_field); //Field-level random effect		
	intPhiFlwCount ~ normal(5,2); //Terms for variance
	slopePlSizePhiFlwCount ~ normal(1,1);
	sigmaPhiFlwCount_field ~ gamma(1,1); //Sigma for field level sigma
	intPhiFlwCount_field ~ normal(0,sigmaPhiFlwCount_field); //Field-level random effect for sigma			
}

generated quantities {
}
