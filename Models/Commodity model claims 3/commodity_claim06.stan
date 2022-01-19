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
	vector<lower=1.8,upper=5>[Nplot_densMiss] plDens_miss;  
	
	// Plant size
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density
	real slopeDistPlSize; //Slope of distance
	real<lower=0> sigmaPlSize_field; //SD of field-level intercept
	real<lower=0> sigmaPlSize_plot; //SD of plot-level intercept
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	vector[Nfield] intPlSize_field; //Random intercept for field - most overlap zero, but not all: probably should keep it
	vector[Nplot] intPlSize_plot; //Random intercept for plot - not great n_eff or Rhat, but looic is worse without it
  
	// hbee Visitation per plot
	real claim06_slopePlSizeHbee; //Claim
	real slopePlDensHbee; //Other claim
	
	real intVisit; //Intercept
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real slopeFlDens; //Slope of flower density
	real<lower=0> sigmaVisField; //SD of field random intercepts
	real<lower=0> visitHbeePhi; //Dispersion parameter
	vector[Nfield] intVisit_field; //field-level random intercepts
	real<lower=0> lambdaVisField; //Lambda for skewed random effects
}

transformed parameters {		
	//Expected values
	
	//Plant density
	vector[Nplot] plDens; //Planting density - imputed

	//Plant size
	vector[Nplot] plSizePlotMu; //Plot-level plant size
	vector[Nplant] plSizeMu; //Expected plant size

	//Visitation
	vector[Nplot] visitHbeeMu; //Expected hbee visits
	
	//Assign imputed data
	plDens[obsPlDens_ind]=plDens_obs; //Observed data
	plDens[missPlDens_ind]=plDens_miss;	//Missing data
	
	for(i in 1:Nplot){ 
	  
		// Plant size (plot-level)
		plSizePlotMu[i] = intPlSize + //Intercept
		  intPlSize_field[plotIndex[i]] + //Field-level random effect
		  intPlSize_plot[i] + //Plot-level random effect
			slopePlDensPlSize*plDens[i] +  //Plant density
			slopeDistPlSize*logHbeeDist[i]; //Distance effect (edge of field has smaller plants)
			
		// Honeybee Visitation
		visitHbeeMu[i] = intVisit + //Intercept
		  intVisit_field[plotIndex[i]] + //Field-level random intercept
		  logTime[i] + //Time offset
			slopePlDensHbee*plDens[i] + //Other claim
		  claim06_slopePlSizeHbee * plSizePlotMu[i] + //Claim
			slopeDistVis*logHbeeDist[i] + //distance from edge
			slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			slopeFlDens*flDens[i]; //Flower density
	}
	
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]];
	}	
	
}
	
model {	
  
  plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size
  hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
	
	// Priors
	
	// Plant size - informative priors
	intPlSize ~ normal(4.6,5); //Intercept
	slopePlDensPlSize ~ normal(0,5); //Plant density
	slopeDistPlSize ~ normal(0,5); //Distance effect
	sigmaPlSize_field ~ gamma(1,1); //Sigma for random field
	sigmaPlSize_plot ~ gamma(1,1); //Sigma for random plot
	sigmaPlSize ~ gamma(1,1); //Sigma for residual
	intPlSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intPlSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot
	
	// Visitation - informative priors
	claim06_slopePlSizeHbee ~ normal(0,5); //Claim
	slopePlDensHbee ~ normal(0,5); //Other claim
	
	intVisit ~ normal(-1.4,5); //Intercept
	slopeDistVis ~ normal(0,5); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0,5); //Slope of hive effect on visits
	slopeFlDens ~ normal(0,5); //Flower density
	sigmaVisField ~ gamma(1,1); //Sigma for random field
	intVisit_field ~ exp_mod_normal(0,sigmaVisField,lambdaVisField); //Skewed random effects
	lambdaVisField ~ gamma(1,1); //Lambda for skewed random effects
	visitHbeePhi ~ gamma(1,1); //Dispersion parameter for NegBin
}