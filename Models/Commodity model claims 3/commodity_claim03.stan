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
	// Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slopeDistPlDens; //Slope of distance into field		
	real<lower=0> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field
	
	// Plant size
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density
	real slopeDistPlSize; //Slope of distance
	real<lower=0> sigmaPlSize_field; //SD of field-level intercept
	real<lower=0> sigmaPlSize_plot; //SD of plot-level intercept
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	vector[Nfield] intPlSize_field; //Random intercept for field - most overlap zero, but not all: probably should keep it
	vector[Nplot] intPlSize_plot; //Random intercept for plot - not great n_eff or Rhat, but looic is worse without it

	// Flower density per plot
	real claim03_slopeNumHiveFlDens; //Claim
	real intFlDens; //Global intercept
	real slopePlSizeFlDens; //Slope of plant size on flower density
	real slopeHbeeDistFlDens; //Slope of distance into field
	real<lower=0> sigmaFlDens; //Sigma for within-field (residual)
	real<lower=0> sigmaFlDens_field; //Sigma for field
	vector[Nfield] intFlDens_field; //Random intercept for field
}

transformed parameters {		
	//Expected values
	
	//Plant density
	vector[Nplot] plDensMu; //Expected plant density
	vector[Nplot] plDens; //Planting density - imputed

	//Plant size
	vector[Nplot] plSizePlotMu; //Plot-level plant size
	vector[Nplant] plSizeMu; //Expected plant size

	//Flower density
	vector[Nplot] flDensMu; //Expected flower density
	
	//Assign imputed data
	plDens[obsPlDens_ind]=plDens_obs; //Observed data
	plDens[missPlDens_ind]=plDens_miss;	//Missing data

	for(i in 1:Nplot){ 
		// Plant density per plot
		plDensMu[i] = intPlDens + //Intercept
		  intPlDens_field[plotIndex[i]] + //Field level random intercept
			slopeDistPlDens*logHbeeDist[i]; //Distance effect
			

		// Plant size (plot-level)
		plSizePlotMu[i] = intPlSize + //Intercept
		  intPlSize_field[plotIndex[i]] + //Field-level random effect
		  intPlSize_plot[i] + //Plot-level random effect
			slopePlDensPlSize*plDens[i] +  //Plant density
			slopeDistPlSize*logHbeeDist[i]; //Distance effect (edge of field has smaller plants)
			// slopeGpPlSize*isGP[plotIndex[i]] + // Location effect
			// slope2015PlSize*is2015[plotIndex[i]] + //Year effect
			// slopeIrrigPlSize*isIrrigated[plotIndex[i]]; //Irrigation effect

		// Flower density
		flDensMu[i] = intFlDens	+ //Intercept
		  intFlDens_field[plotIndex[i]] + //Field-level random effect
			slopePlSizeFlDens*plSizePlotMu[i] + //Plant size effect
			slopeHbeeDistFlDens*logHbeeDist[i] + //Distance effect
			claim03_slopeNumHiveFlDens*logNumHives[plotIndex[i]]; //Claim
	}
	
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]];
	}	
	
}
	
model {	
  
  //Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size
	flDens ~ normal(flDensMu,sigmaFlDens); //Flower density per plot
	
	// Priors
	//Plant density	- informative priors
	intPlDens ~ normal(3.8,5); //Global intercept
	slopeDistPlDens ~ normal(0,5); //Slope of distance into field	
	sigmaPlDens ~ gamma(1,1); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(1,1); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field
	
	// Plant size - informative priors
	intPlSize ~ normal(4.6,5); //Intercept
	slopePlDensPlSize ~ normal(0,5); //Plant density
	slopeDistPlSize ~ normal(0,5); //Distance effect
	sigmaPlSize_field ~ gamma(1,1); //Sigma for random field
	sigmaPlSize_plot ~ gamma(1,1); //Sigma for random plot
	sigmaPlSize ~ gamma(1,1); //Sigma for residual
	intPlSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intPlSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot

	// Flower density per plot
	
	claim03_slopeNumHiveFlDens ~ normal(0,5); //Claim
	intFlDens ~ normal(4.1,5); //Global intercept
	slopePlSizeFlDens ~ normal(0,5); //plant size effect
	slopeHbeeDistFlDens ~ normal(0,5); //distance into field
	sigmaFlDens ~ gamma(1,1); //Sigma for within-field (residual)
	sigmaFlDens_field ~ gamma(1,1); //Sigma for field
	intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field
	
}