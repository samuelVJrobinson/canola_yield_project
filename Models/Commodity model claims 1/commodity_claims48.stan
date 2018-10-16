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
	//Claim: seedCount ~ flw count per plant
	real slopeFlwCountSeedCount;

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
		
	// Seed count
	real intSeedCount; //Intercept
	real slopeVisitSeedCount; //Slope of hbee visits
	real slopePolSeedCount; //Slope of pollen deposition
	real slopePlSizeCount; //Slope of plant size	
	real<lower=0> sigmaSeedCount_plant; //SD of plant random effect - OK	
	real<lower=0> sigmaSeedCount_field; //SD of field random effect - OK
	vector[Nplant] intSeedCount_plant; //plant-level random intercepts	
	vector[Nfield] intSeedCount_field; //field-level random intercepts 
	real<lower=0> seedCountPhi; //Dispersion parameter	
}

transformed parameters {		
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplant] plSizeMu; //Plant size	
	vector[Nplot] plSizePlotMu; //Plot-level plant size		
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma			
	vector[Nplot] seedCountMuPlot; //Plot-level seed count
	vector[Nplant] seedCountMuPlant; //Plant-level seed count
	vector[Npod] seedCountMu; //Pod-level seed counts				
	
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
		
		//Plot-level seed count	 = intercept + random int field + random int plot + random int plant + 
		seedCountMuPlot[i] = intSeedCount + intSeedCount_field[plotIndex[i]] + 
			slopeVisitSeedCount*logHbeeVis[i] + //(log) hbee visits 
			slopePolSeedCount*pollenPlot[i]; //pollen deposition - large correlation b/w slopePolSeedCount and intFlwSurv
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //Matches to plot-level pollen + global intercept	
		
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]]; 			
			
		// Seed count per pod = plot-level effect + random plant int +
		seedCountMuPlant[i] = seedCountMuPlot[plantIndex[i]] + intSeedCount_plant[i] + 			
			slopePlSizeCount*plantSize[i] + //plant size
			slopeFlwCountSeedCount*log(flwCount[i]); //Flower count per plant effect		
	}
	
	for(i in 1:Npod){ //For each pod
		seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 		
	}
}
	
model {		
	//Likelihood
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate			
	seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod	
		
	//Priors	
	//Claim
	slopeFlwCountSeedCount ~ normal(0,1);
	
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

	//Seed count - informative priors
	intSeedCount ~ normal(3.5,1); //Intercept
	slopeVisitSeedCount ~ normal(0.01,0.05); //Slope of hbee visits
	slopePolSeedCount ~ normal(0,0.5); //Slope of pollen deposition
	slopePlSizeCount ~ normal(0,0.05); //Slope of plant size
	seedCountPhi ~ normal(22,1); //Dispersion parameter
	sigmaSeedCount_field ~ gamma(2,10); //SD of field random effect	
	sigmaSeedCount_plant ~ gamma(2,10); //SD of plant random effect	
	intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts		
	intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts		
}

generated quantities{
//Plot-level quantities
	//Pod-level
	//seeds per pod
	int predSeedCount[Npod]; //Generated
	real seedCount_resid[Npod]; //Residual
	real predSeedCount_resid[Npod]; //Residual of generated	
		
	for(i in 1:Npod){ //For each pod
		//Seed count per pod - doesn't work well due to weird generating process
		seedCount_resid[i] = seedCount[i] - exp(seedCountMu[i]);
		predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); 
		predSeedCount_resid[i] = predSeedCount[i] - exp(seedCountMu[i]);		
	}	
}
