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
	vector[Nplant] avgSeedCount; //Average seeds per pod
	vector[Nplant] avgSeedMass; //Average weight per seed (mg)
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
		//calcYield = log(pod count x seed weight x seed count)
		calcYield[i] = podCount[i]*avgSeedCount[i]*(avgSeedMass[i]/1000);
		logCalcYield[i] = log(calcYield[i]);
	}		
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
	
}

parameters {
	// Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slope2015PlDens; //Effect of 2015
	real slopeIrrigPlDens; //Effect of irrigation
	real slope2015IrrigPlDens; //Year:irrigation interaction	
	real slopeGPPlDens; //GP effect on plant density	
	real slopeDistPlDens; //Slope of distance into field		
	real<lower=0> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field

	// Pollen deposition
	// Plot level random effect has bad trace, strongly correlated with lp__
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits	
	real slopeHbeeDistPollen; //Effect of distance into field - p=0.076 (p=0.17 if stocking and stocking:dist interaction included - see below)
	real<lower=0> sigmaPolField; //SD of field random intercepts	
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts	
	// real<lower=0> sigmaPolPlot; //SD of plot random intercepts
	// vector[Nplot] intPollen_plot; //plot-level random intercepts

	// Seed count
  real intSeedCount; //Intercept
  real slopeVisitSeedCount; //Slope of hbee visits - p=0.81
  real slopePolSeedCount; //Slope of pollen deposition - p=0.74
  real slopePlSizeCount; //Slope of plant size - p=0.22
  real slope2015SeedCount; //Year effect - p=0.0003
  real<lower=0> sigmaSeedCount_plot; //SD of plot random effect - not converging well, Rhat 1.1, small n_eff.
  real<lower=0> sigmaSeedCount_field; //SD of field random effect - OK
  vector[Nplot] intSeedCount_plot; //plot-level random intercepts
  vector[Nfield] intSeedCount_field; //field-level random intercepts
  real<lower=0> sigmaSeedCount; //SD of seed count
  real<lower=0> lambdaSeedCount; //Lambda term for exponential process
 
}

transformed parameters {		
	//Expected values
	
	//Plant density
	vector[Nplot] plDensMu; //Expected plant density
	vector[Nplot] plDens; //Planting density - imputed

	//Pollen deposition
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Expected pollen per stigma		
	
  //Average seeds per pod
  vector[Nplot] seedCountMuPlot; //Plot-level seed count
  vector[Nplant] seedCountMu; //Plant-level seed count
	
	//Assign imputed data
	plDens[obsPlDens_ind]=plDens_obs; //Observed data
	plDens[missPlDens_ind]=plDens_miss;	//Missing data

	for(i in 1:Nplot){ 
		// Plant density per plot
		plDensMu[i] = intPlDens + //Intercept
		  intPlDens_field[plotIndex[i]] + //Field level random intercept
			slope2015PlDens*is2015[plotIndex[i]]+ //Year effect
			slopeIrrigPlDens*isIrrigated[plotIndex[i]]+ //Irrigation effect
			slope2015IrrigPlDens*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]+ //Year:irrigation interaction
			slopeDistPlDens*logHbeeDist[i] + //Distance effect
			slopeGPPlDens*isGP[plotIndex[i]]; //Location effect

		// Plot-level pollen deposition
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //Field-level random intercept
		  // intPollen_plot[i] + //Plot-level random intercept
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect						
		// Global intercept is within flower-level term in order to center plot-level variable
	
		// Plot-level seed count
		seedCountMuPlot[i] = intSeedCount + //Intercept
		  intSeedCount_field[plotIndex[i]] + //Field-level random intercept
		  intSeedCount_plot[i] + //Plot-level random intercept
			slopeVisitSeedCount*logHbeeVis[i] + //(log) hbee visits
			slopePolSeedCount*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedCount and intFlwSurv
			slope2015SeedCount*is2015[plotIndex[i]]; //Year effect

	}
		
	for(i in 1:Nflw){ //For each flower stigma
		pollenMu[i] = intPollen + //Intercept
		  pollenPlot[flowerIndex[i]]; //Plot-level pollen 
	}
		
	for(i in 1:Nplant){ //For each plant 	
		// Seed count per pod
		seedCountMu[i] = seedCountMuPlot[plantIndex[i]] + //Plot-level seed count
			slopePlSizeCount*plantSize[i]; //plant size effect
	}	
	
}
	
model {	
	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	avgSeedCount ~ exp_mod_normal(seedCountMu,sigmaSeedCount,lambdaSeedCount); //Average seeds per pod
		
	// Priors
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
	
	// Pollen deposition - informative priors	
	intPollen ~ normal(5.5,1); //Intercept	
	slopeVisitPol ~ normal(0,0.1); //hbee Visitation effect	
	slopeHbeeDistPollen ~ normal(0,0.1); //hbee distance effect	
	sigmaPolField ~ gamma(1.25,3); //Sigma for random field	
	pollenPhi ~ gamma(1.25,3); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	// sigmaPolPlot ~ gamma(1.05,1); //Sigma for random plot - bad Rhat, poor traces   
	// intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int - not a lot of info at plot level
	
  // Average seed count - informative priors
  intSeedCount ~ normal(3.5,1); //Intercept
  slopeVisitSeedCount ~ normal(0,0.1); //Slope of hbee visits
  slopePolSeedCount ~ normal(0,0.5); //Slope of pollen deposition+
  slopePlSizeCount ~ normal(0,0.05); //Slope of plant size
  slope2015SeedCount ~ normal(0,0.5); //Year effect
  // slopeIrrigSeedCount ~ normal(0,0.5); //Irrigation
  sigmaSeedCount_field ~ gamma(2,10); //SD of field random effect
  sigmaSeedCount_plot ~ gamma(2,10); //SD of plot random effect
  intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts
  intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts
  sigmaSeedCount ~ gamma(1,1); //SD of seed count
  lambdaSeedCount ~ gamma(1,1); //Lambda for exp-normal distribution
}

generated quantities {
  // Average seeds per pod
  real predSeedCount[Nplant];
  real seedCount_resid[Nplant];
  real predSeedCount_resid[Nplant];

	for(i in 1:Nplant){
		// Seed count per pod - doesn't work well due to weird generating process (I think)
		seedCount_resid[i] = avgSeedCount[i] - (seedCountMu[i]+(1/lambdaSeedCount));
		predSeedCount[i] = exp_mod_normal_rng(seedCountMu[i],sigmaSeedCount,lambdaSeedCount);
		predSeedCount_resid[i] = predSeedCount[i] - (seedCountMu[i]+(1/lambdaSeedCount));
	}
}
