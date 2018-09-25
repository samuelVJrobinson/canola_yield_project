data {
	//Field level
	int Nfield; //Number of fields	
	real numHives[Nfield]; //Number of hives present
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
	vector[Nplot] hbee_dist; //Distance from edge
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	real totalTime[Nplot]; //Time taken for observation/10
	//Missing data from plot level
	int Nplot_flsObs; //Number of plots where flower density observed
	int Nplot_flsMiss; //Number of plots where flower density missing 
	vector[Nplot_flsObs] flDens_obs; //Observed flower density (flowers/m2)
	int<lower=1,upper=Nplot> obsFls_ind[Nplot_flsObs]; //Index for observed flower density
	int<lower=1,upper=Nplot> missFls_ind[Nplot_flsMiss]; //Index for missing flower density	
	//Plant density (stems/m2)
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
	int podCount[Nplant_obs]; //Number of pods per plant
	int flwCount[Nplant_obs]; //Number of total flower (pods + missing) per plant
	vector[Nplant_obs] plantSize_obs; //Mass of vegetative tissue (no seeds) (g)
	vector[Nplant_obs] totalSeedMass; //Mass of seeds from plant (g)
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?	
	int<lower=1,upper=Nplot> plantSurvIndex[Nplant_obs]; //Index for "complete" plants 
	int<lower=1,upper=Nplant> obs_ind[Nplant_obs]; //Index for observed plants
	int<lower=1,upper=Nplant> miss_ind[Nplant_miss]; //Index for missing plants	
	
	//Pod level
	int Npod; //Number of pods
	int seedCount[Npod];  //Seeds per pod
	vector[Npod] seedMass; //Mass of seedCount seeds (g)
	int<lower=1> podIndex[Npod]; //Index for pods - which plant?	
}

transformed parameters {
	//Transformations
	logHbeeDist=log(hbee_dist); //Log-transform distances	
	logTime=log(totalTime); //Log-transform time	
}

parameters {
	//Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slopeDistPlDens; //Slope of distance into field	
	real<lower=0.01> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaPlDens_field; //Sigma for field
	vector[Nfield_all] intPlDens_field; //Random intercept for field
		
	//Flower density per plot
	vector<lower=0.01>[Nplot_flsMiss] flDens_miss; //Vector for missing values	
	real intFlDens; //Global intercept
	real slopePlSizeFlDens; //Slope of plant size on flower density
	real<lower=0.01> sigmaFlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaFlDens_field; //Sigma for field
	vector[Nfield_all] intFlDens_field; //Random intercept for field	

	//Plant size
	vector[Nplant_miss] plantSize_miss; //Imputed data for missing values	
	real intSize; //Global intercept
	real slopePolPlSize; //Effect of pollination on plant size
	real<lower=0> sigmaPlSize_field; //Sigma for field
	real<lower=0> sigmaPlSize_plot; //Sigma for plot
	real<lower=0> sigmaPlSize; //Sigma for within-plot	(residual)
	vector[Nfield] intSize_field; //Random intercept for field
	vector[Nplot] intSize_plot; //Random intercept for plot	

	//hbee Visitation
	real intVisit; //Intercept 
	real slopeYear2015Vis; //Effect of 2015 - probably negative
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real<lower=0> sigmaVisField; //SD of field random intercepts
	real<lower=0> visitPhi; //Dispersion parameter	
	vector[Nfield] intVisit_field; //field-level random intercepts
	real<lower=0,upper=1> zeroVisTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. - about 0.2, but doesn't improve traces
	
	//Pollen deposition
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits
	real<lower=0> sigmaPolField; //SD of field random intercepts
	real<lower=0> sigmaPolPlot; //SD of plot random intercepts
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts
	vector[Nplot] intPollen_plot; //plot-level random intercepts
	
	//Flower survival
	real intFlwSurv; //Intercept
	real slopeVisitSurv; //Slope of hbee visits
	real slopePolSurv; //Slope of pollen deposition
	real slopePlSizeSurv; //Slope of plant size
	real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts	
	real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	vector[Nfield] intFlwSurv_field; //field-level random intercepts
	vector[Nplot] intFlwSurv_plot; //plot-level random intercepts
	
	//Seed count
	real intSeedCount; //Intercept
	real slopeVisitSeedCount; //Slope of hbee visits
	real slopePolSeedCount; //Slope of pollen deposition
	real slopePlSizeCount; //Slope of plant size
	real<lower=0> seedCountPhi; //Dispersion parameter
	real<lower=0> sigmaSeedCount_plant; //SD of plant random effect
	real<lower=0> sigmaSeedCount_plot; //SD of plot random effect
	real<lower=0> sigmaSeedCount_field; //SD of field random effect
	vector[Nfield] intSeedCount_field; //field-level random intercepts	
	vector[Nplot] intSeedCount_plot; //plot-level random intercepts	
	vector[Nplant] intSeedCount_plant; //plant-level random intercepts
	
	//Weight per seed
	real intSeedWeight; //Intercept
	real slopeVisitSeedWeight; //Slope of hbee visits
	real slopePolSeedWeight; //Slope of pollen deposition
	real slopeSeedCount; //Slope of seed count
	real slopePlSizeWeight; //Slope of plant size
	real<lower=0> sigmaSeedWeight; //SD of seed weight
	real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect
	real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect
	real<lower=0> sigmaSeedWeight_field; //SD of field random effect	
	vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	vector[Nplant] intSeedWeight_plant; //plant-level random intercepts		
}

transformed parameters {		
	//Expected values
	vector[Nplot] visitMu; //Plot-level hbee visits	
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma		
	vector[Nplot] flwSurvPlot; //Plot-level survival
	vector[Nplant_obs] flwSurv; //Flower survival rate (logit)
	vector[Nplant] seedCountMuPlant; //Plant-level seed count
	vector[Npod] seedCountMu; //Pod-level seed counts	
	vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	vector[Npod] seedWeightMu; //Pod-level weight per seed		
	vector[Nplant] plSizeMu; //Plant size
	
	//Imputed missing data;
	vector[Nplant] plantSize;
	plantSize[obs_ind]=plantSize_obs;
	plantSize[miss_ind]=plantSize_miss;	
	
	for(i in 1:Nplot){ //Expected value for visitation = intercept + random int + distance + numHives + time offset
		visitMu[i] = intVisit + intVisit_field[plotIndex[i]] + slopeDistVis*dist[i] + slopeHiveVis*numHives[plotIndex[i]] + log(totalTime[i]); 
		
		//Expected value for pollen = intercept + random int field + random int plot + hbee visits/10
		pollenPlot[i] = intPollen + intPollen_field[plotIndex[i]] + intPollen_plot[i] + slopeVisitPol*(hbeeVis[i]/totalTime[i]); //Plot level pollen
		
		//Flower survival = intercept + random int field + random int plot + hbee visits + pollen deposition/1000
		flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + intFlwSurv_plot[i] + 
			slopeVisitSurv*(hbeeVis[i]/totalTime[i]) + slopePolSurv*(exp(pollenPlot[i])/1000);
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = pollenPlot[flowerIndex[i]]; //Matches to plot-level pollen		
		
	for(i in 1:Nplant_obs)	
		flwSurv[i] = flwSurvPlot[plantSurvIndex[i]] + slopePlSizeSurv*plantSize[i]; //Plot-level plant survival + size effect
		
	for(i in 1:Nplant){ //For each plant (complete measurements or not)		
		//Seed count per pod = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000 + plant size
		seedCountMuPlant[i] = intSeedCount + intSeedCount_field[plotIndex[plantIndex[i]]] + intSeedCount_plot[plantIndex[i]] + intSeedCount_plant[i] + 
			slopeVisitSeedCount*(visitMu[plantIndex[i]]+log(1-zeroVisHbeeTheta)) + 
			slopePolSeedCount*pollenPlot[plantIndex[i]] + 
			slopePlSizeCount*plantSize[i];
			
		//Weight per seed = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000
		seedWeightPlantMu[i] = intSeedWeight + intSeedWeight_field[plotIndex[plantIndex[i]]] + intSeedWeight_plot[plantIndex[i]] + intSeedWeight_plant[i] +
			slopeVisitSeedWeight*(visitMu[plantIndex[i]]+log(1-zeroVisHbeeTheta)) + slopePolSeedWeight*pollenPlot[plantIndex[i]]+
			slopePlSizeWeight*plantSize[i]; //Plant-level effect		
			
		//Predicted plant size = intercept + random int field + random int plot	
		plSizeMu[i] = intSize + intSize_field[plotIndex[plantIndex[i]]] + intSize_plot[plantIndex[i]] + slopePolPlSize*pollenPlot[plantIndex[i]]; 
	}
	
	for(i in 1:Npod){ //For each pod
		seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+slopeSeedCount*seedCount[i]; //Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
	}
	}
	
model {
	vector[2] bernLL; //pre-calculate LL for zero inflation process
	bernLL[1]=bernoulli_lpmf(0|zeroVisTheta); //LL of no extra zero
	bernLL[2]=bernoulli_lpmf(1|zeroVisTheta); //LL of extra zero
	for(i in 1:Nplot){ //Zero-inflated negbin for visitation frequency
		if(hbeeVis[i]==0)
			target += log_sum_exp(bernLL[2],bernLL[1]+neg_binomial_2_log_lpmf(0|visitMu[i],visitPhi));
		else
			target += bernLL[1]+neg_binomial_2_log_lpmf(hbeeVis[i]|visitMu[i],visitPhi);	
	}		
		
	//Likelihood		
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	podCount ~ binomial_logit(flwCount,flwSurv); //Flower survival
	seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	seedMass ~ lognormal(seedWeightMu,sigmaSeedWeight); //Weight per seed
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size
		
	//Priors
	//Visitation - informative priors
	intVisit ~ normal(-1,1); //Intercept
	sigmaVisField ~ gamma(4,2); //Sigma for random field 
	slopeDistVis ~ normal(0,1); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0,1); //Slope of hive effect on visits 
	visitPhi ~ gamma(4,10); //Dispersion parameter	
	intVisit_field ~ normal(0,sigmaVisField); //Random field int
	// zeroVisTheta ~ beta(2,2); //Zero-inflation parameter
		
	//Pollen deposition - informative priors
	intPollen ~ normal(5.5,1); //Intercept	
	sigmaPolField ~ gamma(2,4); //Sigma for random field
	sigmaPolPlot ~ gamma(4,10); //Sigma for random plot	
	slopeVisitPol ~ normal(0,1); //hbee Visitation effect
	pollenPhi ~ gamma(7,10); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int
	
	//Flower survival - informative priors
	intFlwSurv ~ normal(1.2,1); //Intercept
	slopeVisitSurv ~ normal(0,0.02); //Slope of hbee visits
	slopePolSurv ~ normal(-0.5,1); //Slope of pollen deposition
	slopePlSizeSurv ~ normal(0,3); //Slope of plant size
	sigmaFlwSurv_field ~ gamma(1.5,5); //SD of field random effect
	sigmaFlwSurv_plot ~ gamma(1.5,5); //SD of plot random effect	
	intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts	
	
	//Seed count - informative priors
	intSeedCount ~ normal(3,1); //Intercept
	slopeVisitSeedCount ~ normal(0,0.01); //Slope of hbee visits
	slopePolSeedCount ~ normal(0,1); //Slope of pollen deposition
	slopePlSizeCount ~ normal(0.04,0.05); //Slope of plant size
	seedCountPhi ~ gamma(21,1); //Dispersion parameter
	sigmaSeedCount_field ~ gamma(2,20); //SD of field random effect
	sigmaSeedCount_plot ~ gamma(2,20); //SD of plot random effect
	sigmaSeedCount_plant ~ gamma(2,20); //SD of plant random effect	
	intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts
	
	//Weight per seed - informative priors
	intSeedWeight ~ normal(1,1); //Intercept
	slopeVisitSeedWeight ~ normal(0,0.01); //Slope of hbee visits
	slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	slopeSeedCount ~ normal(0.005,.01); //Slope of seed count
	slopePlSizeWeight ~ normal(0.02,0.02); //Slope of plant size
	sigmaSeedWeight ~ gamma(3,10); //SD of seed weight
	sigmaSeedWeight_field ~ gamma(2,10); //SD of field random effect	
	sigmaSeedWeight_plot ~ gamma(2,10); //SD of plot random effect
	sigmaSeedWeight_plant ~ gamma(2,10); //SD of plant random effect		
	intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts

	//Plant size - informative priors
	intSize ~ normal(0,1); //Intercept
	slopePolPlSize ~ normal(0,1); //Effect of pollen on plant size
	sigmaPlSize_field ~ gamma(2,6); //Sigma for random field 
	sigmaPlSize_plot ~ gamma(2,7); //Sigma for random plot
	sigmaPlSize ~ gamma(3,5); //Sigma for residual	
	intSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot	
}

generated quantities{
//Plot-level quantities
	//plant density
	real predPlDens[Nplot]; //Generated
	real plDens_resid[Nplot]; //Residual
	real predPlDens_resid[Nplot]; //Residual of generated
	//flower density
	real predFlDens[Nplot]; //Generated
	real flDens_resid[Nplot]; //Residual
	real predFlDens_resid[Nplot]; //Residual of generated	
	//hbeeVis
	int predHbeeVis_all[Nplot]; //Generated
	real hbeeVis_resid[Nplot]; //Residual 
	real predHbeeVis_resid[Nplot]; //Residual of generated 	
	
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
	// int predFlwCount[Nplant_obs]; //Generated
	// real flwCount_resid[Nplant_obs]; //Residual
	// real predFlwCount_resid[Nplant_obs]; //Residual of generated	
	// //flower survival (surviving pods)
	// int predPodCount[Nplant_obs]; //Generated
	// real podCount_resid[Nplant_obs]; //Residual
	// real predPodCount_resid[Nplant_obs]; //Residual of generated
	
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
		//plant density		
		plDens_resid[i] = plDens[i] - plDensMu[i]; //Residual for actual value
		predPlDens[i]= normal_rng(plDensMu[i],sigmaPlDens); //Generated value from normal
		predPlDens_resid[i] = predPlDens[i] - plDensMu[i]; //Residual for predicted value							
	
		//flower density
		flDens_resid[i] = flDens[i]-flDensMu[i]; //Residual for actual value
		predFlDens[i] = normal_rng(flDensMu[i],sigmaFlDens); //Generated value from normal
		predFlDens_resid[i] = predFlDens[i] - flDensMu[i]; //Residual for predicted value	
	
		//bee visits
		hbeeVis_resid[i]=hbeeVis_all[i]-(exp(visitMu_hbee[i])*(1-zeroVisHbeeTheta)); //Residual for actual value x offset
		if(bernoulli_rng(zeroVisHbeeTheta)==1) //If zeroVisHbeeTheta generates an extra zero
			predHbeeVis_all[i] = 0; //Predicted value is automatically zero
		else //Otherwise
			predHbeeVis_all[i] = neg_binomial_2_log_rng(visitMu_hbee[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		predHbeeVis_resid[i]=predHbeeVis_all[i]-(exp(visitMu_hbee[i])*(1-zeroVisHbeeTheta)); //Residual for predicted value			
		
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
	}
	// for(i in 1:Nplant_obs){
		// //flower number per plant
		// flwCount_resid[i] = flwCount[i] - exp(flwCountMu[i]); //Residual for actual
		// predFlwCount[i] = neg_binomial_2_log_rng(flwCountMu[i],flwCountPhi); //Generates new value from neg. bin.
		// predFlwCount_resid[i] = predFlwCount[i] - exp(flwCountMu[i]); //Residual for new value
		// //pod count (surviving pods)
		// podCount_resid[i] = podCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for actual
		// predPodCount[i] = binomial_rng(flwCount[i],inv_logit(flwSurv[i])); //Generates new value from binomial
		// predPodCount_resid[i] = predPodCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for new value		
	// }
	
	// for(i in 1:Npod){ //For each pod
		// //Seed count per pod
		// seedCount_resid[i] = seedCount[i] - exp(seedCountMu[i]);
		// predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); 
		// predSeedCount_resid[i] = predSeedCount[i] - exp(seedCountMu[i]);
		// //weight per seed
		// seedMass_resid[i] = seedMass[i] - exp(seedWeightMu[i]+(pow(sigmaSeedWeight,2)/2)); 
		// predSeedMass[i] = lognormal_rng(seedWeightMu[i],sigmaSeedWeight); 
		// predSeedMass_resid[i] = predSeedMass[i] - exp(seedWeightMu[i]+(pow(sigmaSeedWeight,2)/2));
	// }	
}
