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
	
	// Flower count (per plant)
	real intFlwCount; //Intercept
	real slopePlSizeFlwCount; //Slope of plant size
	real slopeFlwSurvFlwCount; //Slope of flower survival rate
	// real slope2015FlwCount; //Slope of 2015 effect
	real<lower=0> phiFlwCount_field; //SD of field-level random effect
	vector[Nfield] intFlwCount_field; //Field-level random effect
	// real<lower=0> phiFlwCount_plot; //SD of plot-level random effect
	// vector[Nplot] intFlwCount_plot; //Plot-level random effect
	//Sigma (variance) modeled as a function of plant size
	real intPhiFlwCount; //Intercept for sigma
	real slopePlSizePhiFlwCount; //Effect of plant size on sigma
	real<lower=0> sigmaPhiFlwCount_field; //Sigma for field level sigma
	vector[Nfield] intPhiFlwCount_field; //Field-level random intercept
	// real<lower=0> sigmaPhiFlwCount_plot; //Sigma for plot level sigma - bad traces, high Rhat
	// vector[Nplot] intPhiFlwCount_plot; //Plot-level random intercept - all overlap zero
}

transformed parameters {		
	//Expected values
	
	//Flower count per plant
	vector[Nplot] flwCountPlot; //Plot-level flowers per plant
	vector[Nplant] flwCountMu; //Expected flowers per plant
	vector<lower=0>[Nplant] phiFlwCount; //Phi (variance) for flowers per plant

	for(i in 1:Nplot){ 
		// Flower count per plant (plot level)
    flwCountPlot[i] = intFlwCount + //Intercept
      intFlwCount_field[plotIndex[i]]; //Field-level random intercept
      // intFlwCount_plot[i]; //Plot-level random intercept
    	// slope2015FlwCount*is2015[plotIndex[i]]; //Year effect
	}
	
	for(i in 1:Nplant){ //For each plant 	
    // Flower count per plant (attempted pods)
    flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count
    	slopePlSizeFlwCount*plantSize[i] + //Plant size effect
    	slopeFlwSurvFlwCount*logitFlwSurv[i]; //Flower survival effect
    // Phi (dispersion) for flower count
    phiFlwCount[i] = exp(intPhiFlwCount + //Intercept
      intPhiFlwCount_field[plotIndex[plantIndex[i]]] + //Field-level random intercept
      // intPhiFlwCount_plot[plantIndex[i]] + //Plot-level random intercept
    	slopePlSizePhiFlwCount*plantSize[i]); // Term for sigma
	}	
	
}
	
model {	
  
	//Likelihood		
	flwCount ~ neg_binomial_2_log(flwCountMu,phiFlwCount); //Flower count per plant (attempted pods)
	
	// Priors
	
  //Flower count (per plant)
  intFlwCount ~ normal(2.6,5); //Intercept
  slopePlSizeFlwCount ~ normal(0,5); //Slope of plant size
  slopeFlwSurvFlwCount ~ normal(0,5); //Slope of survival rate
  // slope2015FlwCount ~ normal(0,5); //Slope of 2015 effect
  phiFlwCount_field ~ gamma(1,1); //SD of field-level random effect
  intFlwCount_field ~ normal(0,phiFlwCount_field); //Field-level random intercept
  // phiFlwCount_plot ~ gamma(1,1); //SD of plot-level random effect
  // intFlwCount_plot ~ normal(0,phiFlwCount_plot); //Plot-level random intercept
  //Variance (sigma) terms
  intPhiFlwCount ~ normal(0,5); //Intercept
  slopePlSizePhiFlwCount ~ normal(0,5);
  sigmaPhiFlwCount_field ~ gamma(1,1); //Sigma for field-level sigma
  intPhiFlwCount_field ~ normal(0,sigmaPhiFlwCount_field); //Field-level random intercept
  // sigmaPhiFlwCount_plot ~ gamma(1,1); //Sigma for plot-level sigma - bad traces, high Rhat
  // intPhiFlwCount_plot ~ normal(0, sigmaPhiFlwCount_plot); //Plot-level random intercept - all overlap zero
}

generated quantities {

	// Plant-level
	// flower count per plant (potential pods)
	real predFlwCount[Nplant];
	real flwCount_resid[Nplant];
	real predFlwCount_resid[Nplant];

	for(i in 1:Nplant){
		// flower count per plant
		flwCount_resid[i] = flwCount[i] - exp(flwCountMu[i]); //Residual for actual
		predFlwCount[i] = neg_binomial_2_log_rng(flwCountMu[i],phiFlwCount[i]); //Generates new value from neg. bin.
		predFlwCount_resid[i] = predFlwCount[i] - exp(flwCountMu[i]); //Residual for new value
	}
}
