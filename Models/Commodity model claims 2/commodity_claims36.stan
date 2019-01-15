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
	vector[Nplant] flwSurv;
	vector[Nplant] logitFlwSurv; //Logit flower survival rate
	
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
		
		//Necessary for promoting integers to reals. Otherwise does integer division.
		logitFlwSurv[i] = podCount[i]; 
		logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods		
		if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
			logitFlwSurv[i]=0.01;
		else if(logitFlwSurv[i]>=1)
			logitFlwSurv[i]=0.99;
		logitFlwSurv[i]=logit(logitFlwSurv[i]); //Logit transformed survival
	}	
}

parameters {
	//Claim: avgSeedWeight~plDens+...
	real slopePlDensSeedWeight;

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
	
	// Pollen deposition
	// Plot level random effect has bad trace, strongly correlated with lp__
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits	
	real slopeHbeeDistPollen; //Effect of distance into field - p=0.076 (p=0.17 if stocking and stocking:dist interaction included - see below)
	real<lower=0> sigmaPolField; //SD of field random intercepts	
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts	
	real<lower=0> sigmaPolPlot; //SD of plot random intercepts - Rhat 1.4, small n_eff, strongly correlated with lp__
	vector[Nplot] intPollen_plot; //plot-level random intercepts	
	
	// Weight per seed	
	real intSeedWeight; //Intercept
	real slopeVisitSeedWeight; //Slope of hbee visits - p=0.24
	real slopePolSeedWeight; //Slope of pollen deposition - p=0.50
	real slopeSeedCount; //Slope of seed count - p<0.0001
	real slopePlSizeWeight; //Slope of plant size - p=0.13
	real slopeIrrigSeedWeight; //Irrigation effect - p=0.04	
	// Interactions	
	real<lower=0> sigmaSeedWeight; //SD of seed weight
	real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect - OK
	real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect - not converging well, Rhat 1.08, small n_eff
	real<lower=0> sigmaSeedWeight_field; //SD of field random effect - OK	
	vector[Nplant] intSeedWeight_plant; //plant-level random intercepts		
	vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
}

transformed parameters {		
}
	
model {	
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma			
	vector[Nplot] seedWeightMuPlot; //Plot-level weight per seed
	vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	vector[Npod] seedWeightMu; //Pod-level weight per seed					
	
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
				
		// Plot-level pollen deposition = random int field + random int plot + 
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //intPollen_plot[i] + 
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits			
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect								
			
		// Plot-level seed weight = intercept + random int field + random int plot + 	
		seedWeightMuPlot[i] = intSeedWeight + intSeedWeight_field[plotIndex[i]] + intSeedWeight_plot[i] + 
			slopeVisitSeedWeight*logHbeeVis[i] + //(log) hbee visits 
			slopePolSeedWeight*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedWeight and intFlwSurv
			slopeIrrigSeedWeight*isIrrigated[plotIndex[i]] + //Irrigation effect			
			slopePlDensSeedWeight*plDens[i]; //Claim
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //global intercept + plot-level pollen 
		
	for(i in 1:Nplant){ //For each plant 						
		// Weight per seed = plot-level effect + random int plant + 
		seedWeightPlantMu[i] = seedWeightMuPlot[plantIndex[i]] + intSeedWeight_plant[i] +			
			slopePlSizeWeight*plantSize[i]; //Plant size			
	}	
	
	for(i in 1:Npod){ //For each pod		
		//Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+
			slopeSeedCount*seedCount[i];  //seed count effect						
	}		

	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density	
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed	(mg)	
		
	//Priors
	//Claim
	slopePlDensSeedWeight ~ normal(0,2);
	
	// Plant density	- informative priors
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
	sigmaPolField ~ gamma(5,10); //Sigma for random field	
	pollenPhi ~ gamma(7,10); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int		
	sigmaPolPlot ~ gamma(3,10); //Sigma for random plot 

	// Weight per seed - informative priors
	intSeedWeight ~ normal(1.5,1); //Intercept
	slopeVisitSeedWeight ~ normal(0,0.5); //Slope of hbee visits
	slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	slopeSeedCount ~ normal(0.015,0.015); //Slope of seed count	
	slopePlSizeWeight ~ normal(0,0.5); //Slope of plant size		
	slopeIrrigSeedWeight ~ normal(0,1); //Slope of irrigation
	// Interactions			
	sigmaSeedWeight ~ gamma(3,6); //SD of seed weight
	sigmaSeedWeight_field ~ gamma(4,10); //SD of field random effect	
	sigmaSeedWeight_plot ~ gamma(3,10); //SD of plot random effect
	sigmaSeedWeight_plant ~ gamma(6,10); //SD of plant random effect		
	intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts	
	lambdaSeedWeight ~ gamma(15,10); //Lambda for exp-normal distribution
}

generated quantities {
}
