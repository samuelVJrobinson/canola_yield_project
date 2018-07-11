data {
	//Field level
	int Nfield; //Number of fields
	int<lower=1,upper=Nfield> fieldIndex[Nfield]; //Index for fields
	int numHives[Nfield]; //Number of hives present
	
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
	int Nplant; //Number of plants
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	vector[Nplant] vegMass; //Mass of vegetative tissue (no seeds) (g)
	vector[Nplant] totalSeedMass; //Mass of seeds from plant (g)
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for plants - which plot?
	
	//Pod level
	int Npod; //Number of pods
	int seedCount[Npod];  //Seeds per pod
	vector[Npod] seedMass; //Mass of seedCount seeds (g)
	int<lower=1,upper=Nplant> podIndex[Npod]; //Index for pods - which plant?	
}

parameters {
	//hbee Visitation
	real intVisit; //Intercept 
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real<lower=0.2> sigmaVisField; //SD of field random intercepts
	real<lower=0.2> visitPhi; //Dispersion parameter
	// real<lower=0,upper=1> zeroVisTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. - about 0.2, but doesn't improve traces
	vector[Nfield] intVisit_field; //field-level random intercepts
	
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
	real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts	
	real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	vector[Nfield] intFlwSurv_field; //field-level random intercepts
	vector[Nplot] intFlwSurv_plot; //plot-level random intercepts
	
	//Seed count
	real intSeedCount; //Intercept
	real slopeVisitSeedCount; //Slope of hbee visits
	real slopePolSeedCount; //Slope of pollen deposition
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
	vector[Nplant] flwSurv; //Flower survival rate		
	vector[Nplant] seedCountMuPlant; //Plant-level seed count
	vector[Npod] seedCountMu; //Pod-level seed counts	
	vector[Npod] seedWeightPlantMu; //Plant-level weight per seed
	vector[Npod] seedWeightMu; //Pod-level weight per seed
	
	for(i in 1:Nplot){ //Expected value for visitation = intercept + random int + distance + numHives + time offset
		visitMu[i] = intVisit + intVisit_field[plotIndex[i]] + slopeDistVis*dist[i] + slopeHiveVis*numHives[plotIndex[i]] + log(totalTime[i]); 
		
		//Expected value for pollen = intercept + random int field + random int plot + hbee visits
		pollenPlot[i] = intPollen + intPollen_field[plotIndex[i]] + intPollen_plot[i] + slopeVisitPol*exp(visitMu[i]); //Plot level pollen
		
		//Flower survival = intercept + random int field + random int plot + hbee visits + pollen deposition/1000
		flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + intFlwSurv_plot[i] + 
			slopeVisitSurv*exp(visitMu[i]) + slopePolSurv*(exp(pollenPlot[i])/1000);
	}
		
	for(i in 1:Nflw) //For each flower
		pollenMu[i] = pollenPlot[flowerIndex[i]]; //Matches to plot-level pollen		
		
	for(i in 1:Nplant){ //For each plant
		flwSurv[i] = flwSurvPlot[plantIndex[i]]; //Matches to plot-level plant survival
		//Seed count per pod = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000
		seedCountMuPlant[i] = intSeedCount + intSeedCount_field[plotIndex[plantIndex[i]]] + intSeedCount_plot[plantIndex[i]] + intSeedCount_plant[i] + 
			slopeVisitSeedCount*exp(visitMu[plantIndex[i]]) + slopePolSeedCount*(exp(pollenPlot[plantIndex[i]])/1000);
		//Weight per seed = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000
		seedWeightPlantMu[i] = intSeedWeight + intSeedWeight_field[plotIndex[plantIndex[i]]] + intSeedWeight_plot[plantIndex[i]] + intSeedWeight_plant[i] +
			slopeVisitSeedWeight*exp(visitMu[plantIndex[i]]) + slopePolSeedWeight*(exp(pollenPlot[plantIndex[i]])/1000); //Plant-level effect		
	}
	
	for(i in 1:Npod){ //For each pod
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+slopeSeedCount*seedCount[i]; //Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)
		seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
	}
}

model {		
	// vector[2] bernLL; //pre-calculate LL for zero inflation process
	// bernLL[1]=bernoulli_lpmf(0|zeroVisTheta);
	// bernLL[2]=bernoulli_lpmf(1|zeroVisTheta); 
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
		
	//Priors
	//Visitation
	intVisit ~ normal(0,3); //Intercept
	sigmaVisField ~ gamma(2,1); //Sigma for random field
	slopeDistVis ~ normal(0,5); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0,5); //Slope of hive effect on visits
	visitPhi ~ gamma(2,5); //Dispersion parameter
	// zeroVisTheta ~ beta(2,2); //Zero-inflation parameter
	intVisit_field ~ normal(0,sigmaVisField); //Random field int
		
	//Pollen deposition
	intPollen ~ normal(5,10); //Intercept	
	sigmaPolField ~ gamma(2,1); //Sigma for random field
	sigmaPolPlot ~ gamma(2,1); //Sigma for random plot	
	slopeVisitPol ~ normal(0,5); //hbee Visitation effect
	pollenPhi ~ gamma(3.5,5); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int
	
	//Flower survival
	intFlwSurv ~ normal(0,3); //Intercept
	slopeVisitSurv ~ normal(0,3); //Slope of hbee visits
	slopePolSurv ~ normal(0,3); //Slope of pollen deposition
	sigmaFlwSurv_field ~ gamma(2,1); //SD of field random effect
	sigmaFlwSurv_plot ~ gamma(2,1); //SD of plot random effect	
	intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts	
	
	//Seed count
	intSeedCount ~ normal(0,3); //Intercept
	slopeVisitSeedCount ~ normal(0,3); //Slope of hbee visits
	slopePolSeedCount ~ normal(0,3); //Slope of pollen deposition
	seedCountPhi ~ gamma(2,2); //Dispersion parameter
	sigmaSeedCount_field ~ gamma(2,1); //SD of field random effect
	sigmaSeedCount_plot ~ gamma(2,1); //SD of plot random effect
	sigmaSeedCount_plant ~ gamma(2,1); //SD of plant random effect	
	intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts
	
	//Weight per seed
	intSeedWeight ~ normal(0,3); //Intercept
	slopeVisitSeedWeight ~ normal(0,3); //Slope of hbee visits
	slopePolSeedWeight ~ normal(0,3); //Slope of pollen deposition
	slopeSeedCount ~ normal(0,3); //Slope of seed count
	sigmaSeedWeight ~ gamma(2,1); //SD of seed weight
	sigmaSeedWeight_field ~ gamma(2,1); //SD of field random effect	
	sigmaSeedWeight_plot ~ gamma(2,1); //SD of plot random effect
	sigmaSeedWeight_plant ~ gamma(2,1); //SD of plant random effect		
	intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts	
	
}

// generated quantities{
// int predHbeeVis[Nplot]; //Predicted visits
// int predPollenCount[Nflw]; //Predicted pollen counts

	// for(i in 1:Nplot)
		// predHbeeVis[i] = neg_binomial_2_log_rng(visitMu[i],visitPhi);	
	
	// for(i in 1:Nflw)
		// predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi);
// }
