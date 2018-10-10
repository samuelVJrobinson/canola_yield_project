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
	//Claim: pollen ~ hbee distance
	real slopeHbeeDistPollen;
		
	//Pollen deposition
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits
	real<lower=0> sigmaPolField; //SD of field random intercepts
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts		
}

transformed parameters {			
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma		
	
	for(i in 1:Nplot){ 		
		// Plot-level pollen deposition = random int field + 
		pollenPlot[i] = intPollen_field[plotIndex[i]] + 
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits				
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //Matches to plot-level pollen + global intercept			
}
	
model {		
	//Likelihood
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
		
	//Priors	
	//Claim
	slopeHbeeDistPollen ~ normal(0,1);
		
	// Pollen deposition - informative priors
	intPollen ~ normal(5.5,1); //Intercept	
	slopeVisitPol ~ normal(0,0.1); //hbee Visitation effect
	sigmaPolField ~ gamma(5,10); //Sigma for random field	
	pollenPhi ~ gamma(7,10); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int		
}

generated quantities{
//Plot-level quantities
	//Flower-level
	//pollen deposition
	int predPollenCount[Nflw]; //Generated 
	real pollen_resid[Nflw]; //residual
	real predPollen_resid[Nflw]; //residual of generated				
	
	for(i in 1:Nflw){
		//pollen deposition
		pollen_resid[i]= exp(pollenMu[i]) - pollenCount[i]; //Residual for actual value
		predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulate pollen counts
		predPollen_resid[i] = exp(pollenMu[i]) - predPollenCount[i]; //Residual for predicted		
	}
}
