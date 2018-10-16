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
	//Claim
	real XXX;

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
				
	// // Flower density per plot
	// real intFlDens; //Global intercept
	// real slopePlSizeFlDens; //Slope of plant size on flower density
	// real<lower=0.01> sigmaFlDens; //Sigma for within-field (residual)	
	// real<lower=0.01> sigmaFlDens_field; //Sigma for field
	// vector[Nfield] intFlDens_field; //Random intercept for field		

	// //hbee Visitation
	// real intVisit; //Intercept 
	// real slopeYearVis; //Effect of 2015 
	// real slopeGpVis; //Effect of Grand Prairie 
	// real slopeYearGpVis; //Effect of year-GP interaction
	// real slopeDistVis; //Slope of distance
	// real slopeHiveVis; //Slope of hive number
	// real slopeFlDens; //Slope of flower density
	// real<lower=0> sigmaVisField; //SD of field random intercepts
	// real<lower=0> visitHbeePhi; //Dispersion parameter	
	// vector[Nfield] intVisit_field; //field-level random intercepts	
	
	// //Pollen deposition
	// real intPollen; //Intercept
	// real slopeVisitPol; //Slope of hbee visits
	// real<lower=0> sigmaPolField; //SD of field random intercepts
	// real<lower=0> pollenPhi; //Dispersion parameter
	// vector[Nfield] intPollen_field; //field-level random intercepts	
	
	// //Flower count (per plant) 
	// real intFlwCount; //Intercept
	// real slopePlSizeFlwCount; //Slope of plant size
	// real<lower=0.01> sigmaFlwCount_field; //SD of field-level random effect
	// real<lower=0.01> sigmaFlwCount_plot; //SD of plot-level random effect
	// vector[Nfield] intFlwCount_field; //Field-level random effect
	// vector[Nplot] intFlwCount_plot; //Plot-level random effect	
	// real<lower=0.01> flwCountPhi; //Dispersion parameter
	
	// //Flower survival
	// real intFlwSurv; //Intercept
	// real slopeVisitSurv; //Slope of hbee visits
	// real slopePolSurv; //Slope of pollen deposition
	// real slopePlSizeSurv; //Slope of plant size
	// real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts	
	// real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	// vector[Nfield] intFlwSurv_field; //field-level random intercepts
	// vector[Nplot] intFlwSurv_plot; //plot-level random intercepts - both random intercepts are left skewed, and ahave low n_eff; driven by a few plots/fields (probably McKee5)
	// real<lower=0> flwSurvPhi; //Dispersion parameter for beta-binomial
	
	// // Seed count
	// real intSeedCount; //Intercept
	// real slopeVisitSeedCount; //Slope of hbee visits
	// real slopePolSeedCount; //Slope of pollen deposition
	// real slopePlSizeCount; //Slope of plant size	
	// real<lower=0> sigmaSeedCount_plant; //SD of plant random effect - OK	
	// real<lower=0> sigmaSeedCount_field; //SD of field random effect - OK
	// vector[Nplant] intSeedCount_plant; //plant-level random intercepts	
	// vector[Nfield] intSeedCount_field; //field-level random intercepts 
	// real<lower=0> seedCountPhi; //Dispersion parameter
	
	// //Weight per seed
	// real intSeedWeight; //Intercept
	// real slopeVisitSeedWeight; //Slope of hbee visits
	// real slopePolSeedWeight; //Slope of pollen deposition
	// real slopeSeedCount; //Slope of seed count
	// real slopePlSizeWeight; //Slope of plant size
	// real<lower=0> sigmaSeedWeight; //SD of seed weight
	// real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect - OK
	// real<lower=0> sigmaSeedWeight_field; //SD of field random effect - OK	
	// vector[Nplant] intSeedWeight_plant; //plant-level random intercepts			
	// vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	// real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
}

transformed parameters {		
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplant] plSizeMu; //Plant size	
	vector[Nplot] plSizePlotMu; //Plot-level plant size	
	// vector[Nplot] flDensMu; //Predicted flower density	
	// vector[Nplot] visitHbeeMu; //Plot-level hbee visits	
	// vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	// vector[Nflw] pollenMu; //Flower-level pollen per stigma		
	// vector[Nplot] flwCountPlot; //Plot-level flowers per plant
	// vector[Nplant] flwCountMu; //Expected flowers per plant
	// vector[Nplot] flwSurvPlot; //Plot-level flower survival
	// vector[Nplant] flwSurv; //Flower survival rate (logit)
	// vector[Nplot] seedCountMuPlot; //Plot-level seed count
	// vector[Nplant] seedCountMuPlant; //Plant-level seed count
	// vector[Npod] seedCountMu; //Pod-level seed counts	
	// vector[Nplot] seedWeightMuPlot; //Plot-level weight per seed
	// vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	// vector[Npod] seedWeightMu; //Pod-level weight per seed			
	
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

		// // Flower density = intercept + random field int + 
		// flDensMu[i] = intFlDens	+ intFlDens_field[plotIndex[i]] + 
			// slopePlSizeFlDens*plSizePlotMu[i]; //plant size effect 			
	
		// // Hbee Visitation = intercept + random int + time offset + 	
		// visitHbeeMu[i] = intVisit + intVisit_field[plotIndex[i]] + logTime[i] + 
			// slopeYearVis*is2015[plotIndex[i]] + //Year effect
			// slopeGpVis* +isGP[plotIndex[i]] + //Grand Prairie effect
			// slopeYearGpVis*is2015[plotIndex[i]]*isGP[plotIndex[i]] + //Year:area interaction			
			// slopeDistVis*logHbeeDist[i] + //distance from edge 
			// slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			// slopeFlDens*flDens[i]; //Flower density
		
		// // Plot-level pollen deposition = random int field + random int plot + 
		// pollenPlot[i] = intPollen_field[plotIndex[i]] + //intPollen_plot[i] + 
			// slopeVisitPol*logHbeeVis[i]; //(log) hbee visits		
			
		// // Flower count per plant (plot level) = intercept + random field int + random plot int
		// flwCountPlot[i] = intFlwCount + intFlwCount_field[plotIndex[i]] + intFlwCount_plot[i];
		
		// //Plot-level flower survival = intercept + random int field + random int plot + 
		// flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + intFlwSurv_plot[i] + 
			// slopeVisitSurv*logHbeeVis[i] + //hbee visits 
			// slopePolSurv*pollenPlot[i]; //(log) pollen deposition - large correlation b/w slopePolSurv and intFlwSurv
			
		// //Plot-level seed count	 = intercept + random int field + random int plot + random int plant + 
		// seedCountMuPlot[i] = intSeedCount + intSeedCount_field[plotIndex[i]] + 
			// slopeVisitSeedCount*logHbeeVis[i] + //(log) hbee visits 
			// slopePolSeedCount*pollenPlot[i]; //pollen deposition - large correlation b/w slopePolSeedCount and intFlwSurv
			
		// //Plot-level seed weight = intercept + random int field + random int plot + 	
		// seedWeightMuPlot[i] = intSeedWeight + intSeedWeight_field[plotIndex[i]] + 
			// slopeVisitSeedWeight*logHbeeVis[i] + //(log) hbee visits 
			// slopePolSeedWeight*pollenPlot[i]; //pollen deposition - large correlation b/w slopePolSeedWeight and intFlwSurv
	}
		
	// for(i in 1:Nflw) //For each flower stigma
		// pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //Matches to plot-level pollen + global intercept	
		
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]]; 			
		
		// // Flower count per plant
		// flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count 
			// slopePlSizeFlwCount*plantSize[i]; //individual size effect	
	
		// // Flower survival per plant
		// flwSurv[i] = flwSurvPlot[plantIndex[i]] + slopePlSizeSurv*plantSize[i]; //Plot-level plant survival + size effect		
	
		// // Seed count per pod = plot-level effect + random plant int +
		// seedCountMuPlant[i] = seedCountMuPlot[plantIndex[i]] + intSeedCount_plant[i] + 			
			// slopePlSizeCount*plantSize[i]; //plant size
			
		// // Weight per seed = plot-level effect + random int plant + 
		// seedWeightPlantMu[i] = seedWeightMuPlot[plantIndex[i]] + intSeedWeight_plant[i] +			
			// slopePlSizeWeight*plantSize[i]; //Plant size 		
	}
	
	// for(i in 1:Npod){ //For each pod
		// seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
		// //Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
		// seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+slopeSeedCount*seedCount[i]; 
	// }
}
	
model {		
	//Likelihood
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	// flDens ~ normal(flDensMu,sigmaFlDens); //Flower density per plot
	// hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
	// pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	// flwCount ~ neg_binomial_2_log(flwCountMu,flwCountPhi); //Flower count per plant (attempted pods)
	// podCount ~ beta_binomial(flwCount,inv_logit(flwSurv)*flwSurvPhi,(1-inv_logit(flwSurv))*flwSurvPhi); //Flower survival (surviving pods) - beta binomial	
	// seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	// seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed		
		
	//Priors	
	//Claim
	XXX ~ normal(0,1);
	
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
	
	// // Flower density per plot
	// intFlDens ~ normal(0,1); //Global intercept
	// slopePlSizeFlDens ~ normal(2,1); //plant size effect
	// sigmaFlDens ~ gamma(7,2); //Sigma for within-field (residual)	
	// sigmaFlDens_field ~ gamma(4,2); //Sigma for field
	// intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field			
	
	// //Visitation - informative priors
	// intVisit ~ normal(-1,1); //Intercept	
	// slopeDistVis ~ normal(-0.3,0.3); //Slope of distance effect on hbee visits
	// slopeHiveVis ~ normal(0.6,0.5); //Slope of hive effect on visits 
	// slopeYearVis ~ normal(0.5,1); //2015 effect
	// slopeGpVis ~ normal(0,1); //Effect of Grand Prairie	
	// slopeYearGpVis ~ normal(0,1); // GP-year interaction
	// slopeFlDens ~ normal(0,5); //Flower density 
	// sigmaVisField ~ gamma(4,2); //Sigma for random field 	
	// intVisit_field ~ normal(0,sigmaVisField); //Random field int
	// visitHbeePhi ~ gamma(2,2); //Dispersion parameter		
		
	// // Pollen deposition - informative priors
	// intPollen ~ normal(5.5,1); //Intercept	
	// slopeVisitPol ~ normal(0,0.1); //hbee Visitation effect
	// sigmaPolField ~ gamma(5,10); //Sigma for random field	
	// pollenPhi ~ gamma(7,10); //Dispersion parameter
	// intPollen_field ~ normal(0,sigmaPolField); //Random field int	

	// //Flower count (per plant)
	// intFlwCount ~ normal(6,0.5); //Intercept
	// slopePlSizeFlwCount ~ normal(0.9,0.5); //Slope of plant size
	// sigmaFlwCount_field ~ gamma(1,10); //SD of field-level random effect	
	// sigmaFlwCount_plot ~ gamma(1,10); //SD of plot-level random effect
	// intFlwCount_field ~ normal(0,sigmaFlwCount_field); //Field-level random effect	
	// intFlwCount_plot ~ normal(0,sigmaFlwCount_plot); //Plot-level random effect	
	// flwCountPhi ~ normal(35,1); //Dispersion parameter	
	
	// //Flower survival - informative priors
	// intFlwSurv ~ normal(1,1); //Intercept
	// slopeVisitSurv ~ normal(0,0.05); //Slope of hbee visits
	// slopePolSurv ~ normal(0,0.5); //Slope of pollen deposition
	// slopePlSizeSurv ~ normal(0.02,0.05); //Slope of plant size
	// sigmaFlwSurv_field ~ gamma(3,10); //SD of field random effect
	// sigmaFlwSurv_plot ~ gamma(2,10); //SD of plot random effect	
	// intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	// intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts
	// flwSurvPhi ~ normal(20,5); //Dispersion parameter
	
	// //Seed count - informative priors
	// intSeedCount ~ normal(3.5,1); //Intercept
	// slopeVisitSeedCount ~ normal(0.01,0.05); //Slope of hbee visits
	// slopePolSeedCount ~ normal(0,0.5); //Slope of pollen deposition
	// slopePlSizeCount ~ normal(0,0.05); //Slope of plant size
	// seedCountPhi ~ normal(22,1); //Dispersion parameter
	// sigmaSeedCount_field ~ gamma(2,10); //SD of field random effect
	// sigmaSeedCount_plant ~ gamma(2,10); //SD of plant random effect	
	// intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts		
	// intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts	
		
	// //Weight per seed - informative priors
	// intSeedWeight ~ normal(1.5,1); //Intercept
	// slopeVisitSeedWeight ~ normal(0,0.1); //Slope of hbee visits
	// slopePolSeedWeight ~ normal(0,0.2); //Slope of pollen deposition
	// slopeSeedCount ~ normal(0.015,0.01); //Slope of seed count
	// slopePlSizeWeight ~ normal(0,0.05); //Slope of plant size
	// sigmaSeedWeight ~ gamma(3,6); //SD of seed weight
	// sigmaSeedWeight_field ~ gamma(4,10); //SD of field random effect		
	// sigmaSeedWeight_plant ~ gamma(6,10); //SD of plant random effect		
	// intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts		
	// intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts	
	// lambdaSeedWeight ~ gamma(15,10); //Lambda for exp-normal distribution
}

generated quantities{
//Plot-level quantities
	//planting density
	real predPlDens[Nplot]; //Generated
	real plDens_resid[Nplot]; //Residual
	real predPlDens_resid[Nplot]; //Residual of generated
	// //flower density
	// real predFlDens[Nplot]; //Generated
	// real flDens_resid[Nplot]; //Residual
	// real predFlDens_resid[Nplot]; //Residual of generated	
	// //hbeeVis
	// int predHbeeVis[Nplot]; //Generated
	// real hbeeVis_resid[Nplot]; //Residual 
	// real predHbeeVis_resid[Nplot]; //Residual of generated 	
	
	// //Flower-level
	// //pollen deposition
	// int predPollenCount[Nflw]; //Generated 
	// real pollen_resid[Nflw]; //residual
	// real predPollen_resid[Nflw]; //residual of generated
	
	//Plant-level	
	//plantSize
	real predPlSize[Nplant]; //Generated
	real plSize_resid[Nplant]; //Residual
	real predPlSize_resid[Nplant]; //Residual of generated
	// //flower count per plant (potential pods)
	// int predFlwCount[Nplant]; //Generated
	// real flwCount_resid[Nplant]; //Residual
	// real predFlwCount_resid[Nplant]; //Residual of generated	
	// //flower survival (surviving pods)
	// int predPodCount[Nplant]; //Generated
	// real podCount_resid[Nplant]; //Residual
	// real predPodCount_resid[Nplant]; //Residual of generated
	
	// //Pod-level
	// //seeds per pod
	// int predSeedCount[Npod]; //Generated
	// real seedCount_resid[Npod]; //Residual
	// real predSeedCount_resid[Npod]; //Residual of generated	
	// //weight per seed
	// real predSeedMass[Npod]; //Generated
	// real seedMass_resid[Npod]; //Residual
	// real predSeedMass_resid[Npod]; //Residual of generated
		
	for(i in 1:Nplot){
		// plant density		
		plDens_resid[i] = plDens[i] - plDensMu[i]; //Residual for actual value
		predPlDens[i]= normal_rng(plDensMu[i],sigmaPlDens); //Generated value from normal
		predPlDens_resid[i] = predPlDens[i] - plDensMu[i]; //Residual for predicted value							
	
		// // flower density
		// flDens_resid[i] = flDens[i]-flDensMu[i]; //Residual for actual value
		// predFlDens[i] = normal_rng(flDensMu[i],sigmaFlDens); //Generated value from normal
		// predFlDens_resid[i] = predFlDens[i] - flDensMu[i]; //Residual for predicted value		
		
		// // bee visits (NB version)
		// hbeeVis_resid[i]=hbeeVis[i]-(exp(visitHbeeMu[i])); //Residual for actual value
		// predHbeeVis[i] = neg_binomial_2_log_rng(visitHbeeMu[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		// predHbeeVis_resid[i]=predHbeeVis[i]-(exp(visitHbeeMu[i])); //residual for predicted value					
	}		
	
	// for(i in 1:Nflw){
		// //pollen deposition
		// pollen_resid[i]= exp(pollenMu[i]) - pollenCount[i]; //Residual for actual value
		// predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulate pollen counts
		// predPollen_resid[i] = exp(pollenMu[i]) - predPollenCount[i]; //Residual for predicted		
	// }
			
	for(i in 1:Nplant){
		//plant size
		plSize_resid[i]= plantSize[i] - plSizeMu[i]; //Residual for actual
		predPlSize[i] = normal_rng(plSizeMu[i],sigmaPlSize); //Generates new value from normal dist.
		predPlSize_resid[i] = predPlSize[i] - plSizeMu[i]; //Residual for new value	
		// //flower count per plant
		// flwCount_resid[i] = flwCount[i] - exp(flwCountMu[i]); //Residual for actual
		// predFlwCount[i] = neg_binomial_2_log_rng(flwCountMu[i],flwCountPhi); //Generates new value from neg. bin.
		// predFlwCount_resid[i] = predFlwCount[i] - exp(flwCountMu[i]); //Residual for new value
		// // pod count (surviving pods)
		// podCount_resid[i] = podCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for actual		
		// predPodCount[i] = beta_binomial_rng(flwCount[i],inv_logit(flwSurv[i])*flwSurvPhi,(1-inv_logit(flwSurv[i]))*flwSurvPhi); //Generates new value from beta-binomial
		// predPodCount_resid[i] = predPodCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for new value		
	}
	
	// for(i in 1:Npod){ //For each pod
		// //Seed count per pod - doesn't work well due to weird generating process
		// seedCount_resid[i] = seedCount[i] - exp(seedCountMu[i]);
		// predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); 
		// predSeedCount_resid[i] = predSeedCount[i] - exp(seedCountMu[i]);
		// // weight per seed - exp-normal works well		
		// seedMass_resid[i] = seedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight)); 
		// predSeedMass[i] = exp_mod_normal_rng(seedWeightMu[i],sigmaSeedWeight,lambdaSeedWeight); 
		// predSeedMass_resid[i] = predSeedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));
	// }	
}
