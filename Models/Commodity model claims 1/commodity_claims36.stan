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
	vector[Nplant] propFlwSurv; //(logit) proportion flower survival
	
	for(i in 1:Nfield)
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives
	
	for(i in 1:Nplot) //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5);
		
	for(i in 1:Nplant){
		//Necessary for promoting integers to reals. Otherwise does integer division.
		propFlwSurv[i] = podCount[i]; 
		propFlwSurv[i] = propFlwSurv[i]/flwCount[i]; //Proportion surviving pods		
		if(propFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
			propFlwSurv[i]=0.01;
		else if(propFlwSurv[i]>=1)
			propFlwSurv[i]=0.99;		
	}
	//Logit transform 
	propFlwSurv=logit(propFlwSurv);		
}

parameters {
	//Claim: seed weight ~ hbee distance
	real slopeHbeeDistSeedWeight;

	//Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slopeDistPlDens; //Slope of distance into field	
	real<lower=0.01> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field
	
	//Plant size
	vector[Nplant_miss] plantSize_miss; //Imputed data for missing values	
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density
	real slopeGpPlSize; //Slope of Grand Prairie effect
	real slopeIrrigPlSize; //Slope of irrigation effect	
	real slope2015PlSize; //Slope of 2015 effect
	real<lower=0> sigmaPlSize_field; //Sigma for field
	real<lower=0> sigmaPlSize_plot; //Sigma for plot - small n_eff, higher
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	vector[Nfield] intPlSize_field; //Random intercept for field
	vector[Nplot] intPlSize_plot; //Random intercept for plot - not converging well	
				
	//Pollen deposition
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits
	real<lower=0> sigmaPolField; //SD of field random intercepts
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts	
		
	//Weight per seed
	real intSeedWeight; //Intercept
	real slopeVisitSeedWeight; //Slope of hbee visits
	real slopePolSeedWeight; //Slope of pollen deposition
	real slopeSeedCount; //Slope of seed count
	real slopePlSizeWeight; //Slope of plant size
	real<lower=0> sigmaSeedWeight; //SD of seed weight
	real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect - OK
	real<lower=0> sigmaSeedWeight_field; //SD of field random effect - OK	
	vector[Nplant] intSeedWeight_plant; //plant-level random intercepts			
	vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
}

transformed parameters {		
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplant] plSizeMu; //Plant size	
	vector[Nplot] plSizePlotMu; //Plot-level plant size		
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma				
	vector[Nplot] seedWeightMuPlot; //Plot-level weight per seed
	vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	vector[Npod] seedWeightMu; //Pod-level weight per seed			
	
	//Imputed missing data;
	vector[Nplant] plantSize; //Plant size
	vector[Nplot] plDens; //Planting density
	plantSize[obsPl_ind]=plantSize_obs;
	plantSize[missPl_ind]=plantSize_miss;	
	plDens[obsPlDens_ind]=plDens_obs;
	plDens[missPlDens_ind]=plDens_miss;
	
	for(i in 1:Nplot){ 
		//Plant density per plot
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex[i]] + 
			slopeDistPlDens*logHbeeDist[i]; //Distance effect				
			
		//Plant size (plot-level) = intercept + random field int + random plot int + distance + planting density effect 
		//Density distance interaction is basically 0, so leaving it out
		plSizePlotMu[i] = intPlSize + intPlSize_field[plotIndex[i]] + intPlSize_plot[i] + 						
			slopePlDensPlSize*plDens[i] +  //Planting density			
			slopeGpPlSize*isGP[plotIndex[i]] + // Grand Prairie
			slope2015PlSize*is2015[plotIndex[i]] + //2015 	
			slopeIrrigPlSize*isIrrigated[plotIndex[i]]; //Irrigation

		// Plot-level pollen deposition = random int field + random int plot + 
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //intPollen_plot[i] + 
			slopeVisitPol*logHbeeVis[i]; //(log) hbee visits							
			
		//Plot-level seed weight = intercept + random int field + random int plot + 	
		seedWeightMuPlot[i] = intSeedWeight + intSeedWeight_field[plotIndex[i]] + 
			slopeVisitSeedWeight*logHbeeVis[i] + //(log) hbee visits 
			slopePolSeedWeight*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedWeight and intFlwSurv
			slopeHbeeDistSeedWeight*logHbeeDist[i]; //Hbee distance
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //Matches to plot-level pollen + global intercept	
		
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]]; 			
			
		// Weight per seed = plot-level effect + random int plant + 
		seedWeightPlantMu[i] = seedWeightMuPlot[plantIndex[i]] + intSeedWeight_plant[i] +			
			slopePlSizeWeight*plantSize[i]; //Plant size 		
	}
	
	for(i in 1:Npod){ //For each pod		
		//Plant-level effect + seedCount 
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]] + slopeSeedCount*seedCount[i]; 
	}
}
	
model {		
	//Likelihood
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed		
		
	//Priors	
	//Claim
	slopeHbeeDistSeedWeight ~ normal(0,1);
	
	//Plant density	- informative priors
	intPlDens ~ normal(0,0.1); //Global intercept
	slopeDistPlDens ~ normal(0,0.05); //Slope of distance into field	
	sigmaPlDens ~ gamma(3,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(3,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field
	
	//Plant size - informative priors
	intPlSize ~ normal(0,0.5); //Intercept
	slopePlDensPlSize ~ normal(0,0.5); //Plant density	
	slopeGpPlSize ~ normal(0,0.5); //Grand Prairie effect
	slopeIrrigPlSize ~ normal(0,0.5); //Irrigation effect
	slope2015PlSize ~ normal(0.3,0.5); //2015 effect	
	sigmaPlSize_field ~ gamma(3,10); //Sigma for random field 
	sigmaPlSize_plot ~ gamma(3.5,10); //Sigma for random plot
	sigmaPlSize ~ gamma(7,10); //Sigma for residual	
	intPlSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intPlSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot	
	
	// Pollen deposition - informative priors
	intPollen ~ normal(5.5,1); //Intercept	
	slopeVisitPol ~ normal(0,0.1); //hbee Visitation effect
	sigmaPolField ~ gamma(5,10); //Sigma for random field	
	pollenPhi ~ gamma(7,10); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int	
		
	//Weight per seed - informative priors
	intSeedWeight ~ normal(1.5,1); //Intercept
	slopeVisitSeedWeight ~ normal(0,0.1); //Slope of hbee visits
	slopePolSeedWeight ~ normal(0,0.2); //Slope of pollen deposition
	slopeSeedCount ~ normal(0.015,0.01); //Slope of seed count
	slopePlSizeWeight ~ normal(0,0.05); //Slope of plant size
	sigmaSeedWeight ~ gamma(3,6); //SD of seed weight
	sigmaSeedWeight_field ~ gamma(4,10); //SD of field random effect	
	sigmaSeedWeight_plot ~ gamma(3,10); //SD of plot random effect
	sigmaSeedWeight_plant ~ gamma(6,10); //SD of plant random effect		
	intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts		
	intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts	
	lambdaSeedWeight ~ gamma(15,10); //Lambda for exp-normal distribution
}

generated quantities{
//Plot-level quantities
	//planting density
	real predPlDens[Nplot]; //Generated
	real plDens_resid[Nplot]; //Residual
	real predPlDens_resid[Nplot]; //Residual of generated
		
	//Flower-level
	//pollen deposition
	int predPollenCount[Nflw]; //Generated 
	real pollen_resid[Nflw]; //residual
	real predPollen_resid[Nflw]; //residual of generated
	
	//Plant-level	
	//plantSize
	real predPlSize[Nplant]; //Generated
	real plSize_resid[Nplant]; //Residual
	real predPlSize_resid[Nplant]; //Residual of generated
		
	//Pod-level
	//weight per seed
	real predSeedMass[Npod]; //Generated
	real seedMass_resid[Npod]; //Residual
	real predSeedMass_resid[Npod]; //Residual of generated
		
	for(i in 1:Nplot){
		// plant density		
		plDens_resid[i] = plDens[i] - plDensMu[i]; //Residual for actual value
		predPlDens[i]= normal_rng(plDensMu[i],sigmaPlDens); //Generated value from normal
		predPlDens_resid[i] = predPlDens[i] - plDensMu[i]; //Residual for predicted value									
	}		
	
	for(i in 1:Nflw){
		//pollen deposition
		pollen_resid[i]= exp(pollenMu[i]) - pollenCount[i]; //Residual for actual value
		predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulate pollen counts
		predPollen_resid[i] = exp(pollenMu[i]) - predPollenCount[i]; //Residual for predicted		
	}
			
	for(i in 1:Nplant){
		//plant size
		plSize_resid[i]= plantSize[i] - plSizeMu[i]; //Residual for actual
		predPlSize[i] = normal_rng(plSizeMu[i],sigmaPlSize); //Generates new value from normal dist.
		predPlSize_resid[i] = predPlSize[i] - plSizeMu[i]; //Residual for new value			
	}
	
	for(i in 1:Npod){ //For each pod		
		// weight per seed - exp-normal works well		
		seedMass_resid[i] = seedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight)); 
		predSeedMass[i] = exp_mod_normal_rng(seedWeightMu[i],sigmaSeedWeight,lambdaSeedWeight); 
		predSeedMass_resid[i] = predSeedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));
	}	
}
