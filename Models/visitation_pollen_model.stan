data {
	//Field level
	int Nfield; //Number of fields
	int<lower=1,upper=Nfield> fieldIndex[Nfield]; //Index for fields
	real numHives[Nfield]; //Number of hives present
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
	real dist[Nplot]; //Distance from edge - centered
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	real totalTime[Nplot]; //Time taken for observation
	//int plantDens[Nplot]; //Plants per m2 - some missing values, so might want to do data imputation
	
	//Flower level
	int Nflw; //Number of flowers
	int<lower=1,upper=Nplot>  flowerIndex[Nflw]; 
	int pollenCount[Nflw]; //Number of pollen grains per stigma

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

parameters {
	//hbee Visitation
	real intVisit; //Intercept 
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real<lower=0.2> sigmaVisField; //SD of field random intercepts
	real<lower=0.2> visitPhi; //Dispersion parameter	
	vector[Nfield] intVisit_field; //field-level random intercepts
	// real<lower=0,upper=1> zeroVisTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. - about 0.2, but doesn't improve traces
	
	//Pollen deposition
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits
	real<lower=0.2> sigmaPolField; //SD of field random intercepts
	real<lower=0.2> sigmaPolPlot; //SD of plot random intercepts
	real<lower=0.2> pollenPhi; //Dispersion parameter
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
	
	//Plant size
	vector[Nplant_miss] plantSize_miss; //Imputed data for missing values	
	real intSize; //Global intercept
	real slopePolPlSize; //Effect of pollination on plant size
	real<lower=0> sigmaPlSize_field; //Sigma for field
	real<lower=0> sigmaPlSize_plot; //Sigma for plot
	real<lower=0> sigmaPlSize; //Sigma for within-plot	(residual)
	vector[Nfield] intSize_field; //Random intercept for field
	vector[Nplot] intSize_plot; //Random intercept for plot	
	
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
			slopeVisitSeedCount*(hbeeVis[plantIndex[i]]/totalTime[plantIndex[i]]) + slopePolSeedCount*(exp(pollenPlot[plantIndex[i]])/1000) + 
			slopePlSizeCount*plantSize[i];
			
		//Weight per seed = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000
		seedWeightPlantMu[i] = intSeedWeight + intSeedWeight_field[plotIndex[plantIndex[i]]] + intSeedWeight_plot[plantIndex[i]] + intSeedWeight_plant[i] +
			slopeVisitSeedWeight*(hbeeVis[plantIndex[i]]/totalTime[plantIndex[i]]) + slopePolSeedWeight*(exp(pollenPlot[plantIndex[i]])/1000)+
			slopePlSizeWeight*plantSize[i]; //Plant-level effect		
			
		//Predicted plant size = intercept + random int field + random int plot	
		plSizeMu[i] = intSize + intSize_field[plotIndex[plantIndex[i]]] + intSize_plot[plantIndex[i]] + slopePolPlSize*(exp(pollenPlot[plantIndex[i]])/1000); 
	}
	
	for(i in 1:Npod){ //For each pod
		seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+slopeSeedCount*seedCount[i]; //Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
	}
	}
	
model {
	// vector[2] bernLL; //pre-calculate LL for zero inflation process
	// bernLL[1]=bernoulli_lpmf(0|zeroVisTheta); //LL of no extra zero
	// bernLL[2]=bernoulli_lpmf(1|zeroVisTheta); //LL of extra zero
	// for(i in 1:Nplot){ //Zero-inflated negbin for visitation frequency
		// if(hbeeVis[i]==0)
			// target += log_sum_exp(bernLL[2],bernLL[1]+neg_binomial_2_log_lpmf(0|visitMu[i],visitPhi));
		// else
			// target += bernLL[1]+neg_binomial_2_log_lpmf(hbeeVis[i]|visitMu[i],visitPhi);	
	// }		
		
	//Likelihood
	hbeeVis ~ neg_binomial_2_log(visitMu,visitPhi);	//Visitation rate	
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
int predHbeeVis[Nplot]; //Predicted visits
int predPollenCount[Nflw]; //Predicted pollen counts
int predPodCount[Nplant_obs]; //Simulated surviving pods
int predSeedCount[Npod]; //Simulated seeds per pod
vector[Npod] predSeedMass; //Simulated seed weight

	for(i in 1:Nplot){		
		predHbeeVis[i] = neg_binomial_2_log_rng(visitMu[i],visitPhi);	 //Simulated hbee visits			
	}
	
	for(i in 1:Nflw)
		predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulated pollen counts
		
	for(i in 1:Nplant_obs)
		predPodCount[i] = binomial_rng(flwCount[i],inv_logit(flwSurv[i])); //Simulated surviving pods
	
	for(i in 1:Npod){ //For each pod
		predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); //Seed count per pod
		predSeedMass[i] = lognormal_rng(seedWeightMu[i],sigmaSeedWeight); //Weight per seed
	}	
}
