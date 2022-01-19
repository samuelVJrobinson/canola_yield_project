//Commodity canola pollination model
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
	
	// Pollen deposition
	// Plot level random effect has bad trace, strongly correlated with lp__
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits	
	real slopeHbeeDistPollen; //Effect of distance into field 
	real<lower=0> sigmaPolField; //SD of field random intercepts	
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts	

	// Weight per seed
	real claim32_slopeFlwSurvSeedWeight; //Claim
	
  real intSeedWeight; //Intercept
  real slopeVisitSeedWeight; //Slope of hbee visits - p=0.24
  real slopePolSeedWeight; //Slope of pollen deposition - p=0.50
  real slopeSeedCount; //Slope of seed count - p<0.0001
  real slopePlSizeWeight; //Slope of plant size - p=0.13
  real<lower=0> sigmaSeedWeight; //SD of seed weight
  real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect
  real<lower=0> sigmaSeedWeight_field; //SD of field random effect
  vector[Nplot] intSeedWeight_plot; //plot-level random intercepts
  vector[Nfield] intSeedWeight_field; //field-level random intercepts
  real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
}

transformed parameters {		
  //Expected values
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma		
	vector[Nplot] seedWeightMuPlot; //Plot-level weight per seed
  vector[Nplant] seedWeightMu; //Plant-level weight per seed
	
	// Imputed missing data;
	vector[Nplot] plDens; //Planting density	
	plDens[obsPlDens_ind]=plDens_obs;
	plDens[missPlDens_ind]=plDens_miss;	
	
	for(i in 1:Nplot){ 
	  // Plot-level pollen deposition
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //Field-level random intercept
		 	slopeVisitPol*logHbeeVis[i] + //(log) hbee visits
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect						
			
		// Plot-level seed weight
		seedWeightMuPlot[i] = intSeedWeight + //Intercept
			intSeedWeight_field[plotIndex[i]] + //Field-level random intercept
			intSeedWeight_plot[i] + //Plot-level random intercept
			slopeVisitSeedWeight*logHbeeVis[i] + //(log) hbee visits
			slopePolSeedWeight*pollenPlot[i]; //pollen deposition
	}
		
	for(i in 1:Nflw) //For each flower stigma
			pollenMu[i] = intPollen + //Intercept
		  pollenPlot[flowerIndex[i]]; //Plot-level pollen 
		
	for(i in 1:Nplant){ //For each plant 	
		// Average weight per seed
		seedWeightMu[i] = seedWeightMuPlot[plantIndex[i]] + //Plot-level seed weight
		  claim32_slopeFlwSurvSeedWeight*logitFlwSurv[i] + //Claim
			slopePlSizeWeight*plantSize[i] + //Plant size effect
			slopeSeedCount*seedCount[i]; //Seed count effect (do plants with many seeds/pod have bigger seeds?)
	}	

}
	
model {	
  
	//Likelihood		
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed	(mg)
		
	// Priors
	// Pollen deposition - informative priors	
	intPollen ~ normal(5.6,5); //Intercept	
	slopeVisitPol ~ normal(0,5); //hbee Visitation effect	
	slopeHbeeDistPollen ~ normal(0,5); //hbee distance effect	
	sigmaPolField ~ gamma(1,1); //Sigma for random field	
	pollenPhi ~ gamma(1,1); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	
	// Average weight per seed - informative priors
	
	claim32_slopeFlwSurvSeedWeight ~ normal(0,5); //Claim
	
  intSeedWeight ~ normal(2.0,5); //Intercept
  slopeVisitSeedWeight ~ normal(0,5); //Slope of hbee visits
  slopePolSeedWeight ~ normal(0,5); //Slope of pollen deposition
  slopeSeedCount ~ normal(0,5); //Slope of (log) seed count
  slopePlSizeWeight ~ normal(0,5); //Slope of plant size
  sigmaSeedWeight ~ gamma(1,1); //SD of seed weight
  sigmaSeedWeight_field ~ gamma(1,1); //SD of field random effect
  sigmaSeedWeight_plot ~ gamma(1,1); //SD of plot random effect
  intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts
  intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts
  lambdaSeedWeight ~ gamma(1,1); //Lambda for exp-normal distribution
}
