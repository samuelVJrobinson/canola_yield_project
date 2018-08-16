data {
	//Field level
	int Nfield; //Number of fields	
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 	
	int lbeeStocking[Nplot]; //Half leafcutter stocking? (double-tent treatment)
	vector[Nplot] hbee_dist; //Distance from hbee hives
	int hbeeVis[Nplot]; //Number of hbee visits to plot
	vector[Nplot] lbee_dist; //Distance to lbee shelters
	int lbeeVis[Nplot]; //Number of lbee visits to plot
	int isCent[Nplot]; //Is plot in center of bay?
	int isFBay[Nplot]; //Is plot in female bay?
	vector[Nplot] totalTime; //Minutes spent observing/10 
	//Extra plots/fields from Riley (visitation data)
	int Nfield_extra;	
	int Nplot_extra; 
	int plotIndex_extra[Nplot_extra]; 
	int lbeeStocking_extra[Nplot_extra];
	vector[Nplot_extra] hbee_dist_extra; 
	int hbeeVis_extra[Nplot_extra]; 
	vector[Nplot_extra] lbee_dist_extra;
	int lbeeVis_extra[Nplot_extra];
	int isCent_extra[Nplot_extra];
	int isFBay_extra[Nplot_extra];
	vector[Nplot_extra] totalTime_extra; 
	
	//Flower level
	int Nflw; //Number of flowers (stigmas)
	int flowerIndex[Nflw]; //Index for flowers - which plot?	
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
	int<lower=1,upper=Nplant> podIndex[Npod]; //Index for pods - which plant?	
	
}
transformed data {
	//Amalgamated visitation data
	int Nfield_all = Nfield + Nfield_extra; //Number of fields
	int Nplot_all = Nplot + Nplot_extra; //Number of plots
	int plotIndex_all[Nplot_all]; //Plot index - which field is plot from?
	int lbeeStocking_all[Nplot_all]; //Half leafcutter stocking?
	vector[Nplot_all] hbee_dist_all; //Distance from hbee hives (edge)
	int hbeeVis_all[Nplot_all]; //Visits from hbees
	vector[Nplot_all] lbee_dist_all; //Distance from lbee shelters
	int lbeeVis_all[Nplot_all]; //Visits from lbees
	int isCent_all[Nplot_all]; //Is plot from center of bay?
	int isFBay_all[Nplot_all]; //Is plot from female bay?
	vector[Nplot_all] totalTime_all; //Minutes of observation/10
	vector[Nplot_all] logTime_all; //Log transform of observation time
	vector[Nplot_all] logHbeeDist_all; //Log-distance from hbee hives (edge)
	vector[Nplot_all] logLbeeDist_all; //Log-distance from lbee shelters
		
	plotIndex_all[1:Nplot] = plotIndex;
	plotIndex_all[Nplot+1:Nplot_all] = plotIndex_extra;
	lbeeStocking_all[1:Nplot] = lbeeStocking;
	lbeeStocking_all[Nplot+1:Nplot_all] = lbeeStocking_extra;
	hbee_dist_all[1:Nplot] = hbee_dist/100; //Scales by 100
	hbee_dist_all[Nplot+1:Nplot_all] = hbee_dist_extra/100;
	hbeeVis_all[1:Nplot] = hbeeVis;
	hbeeVis_all[Nplot+1:Nplot_all] = hbeeVis_extra;
	lbee_dist_all[1:Nplot] = lbee_dist;
	lbee_dist_all[Nplot+1:Nplot_all] = lbee_dist_extra;
	lbeeVis_all[1:Nplot] = lbeeVis;
	lbeeVis_all[Nplot+1:Nplot_all] = lbeeVis_extra;
	isCent_all[1:Nplot] = isCent;
	isCent_all[Nplot+1:Nplot_all] = isCent_extra;
	isFBay_all[1:Nplot] = isFBay;
	isFBay_all[Nplot+1:Nplot_all] = isFBay_extra;
	totalTime_all[1:Nplot] = totalTime;
	totalTime_all[Nplot+1:Nplot_all] = totalTime_extra;
	logHbeeDist_all=log(hbee_dist_all); //Log-transform distances
	logLbeeDist_all=log(lbee_dist_all);
	logTime_all=log(totalTime_all); //Log-transform time	
}

parameters {
	//Claim 2: hbee~predPlSize+...
	real slopePlsizeHbee;

	//hbee Visitation
	real intVisitHbee; //Intercept 
	real slopeHbeeDistHbee; //Slope of distance/100	
	real slopeLbeeDistHbee; //Effect of leafcutter distance - competition
	real slopeCentHbee; //Effect of bay position (center)
	real slopeFBayHbee; //Effect of female bay	
	real<lower=0.2> sigmaHbeeVisField; //SD of field random intercepts
	real<lower=0.2> visitHbeePhi; //Dispersion parameter	
	vector[Nfield_all] intVisitHbee_field; //field-level random intercepts
	real<lower=0,upper=1> zeroVisHbeeTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. 
			
	//Plant size
	vector[Nplant_miss] plantSize_miss; //Vector for imputing missing values	
	real intSize; //Global intercept
	real<lower=0> sigmaPlSize_field; //Sigma for field
	real<lower=0> sigmaPlSize_plot; //Sigma for plot
	real<lower=0> sigmaPlSize; //Sigma for within-plot	(residual)
	vector[Nfield_all] intSize_field; //Random intercept for field
	vector[Nplot_all] intSize_plot; //Random intercept for plot		
}

transformed parameters {			
	//Expected values
	//Plot-level
	vector[Nplot_all] visitMu_hbee; //hbee visits - all plot	
	vector[Nplant] plSizeMu; //Plant size

	//Imputed missing data;
	vector[Nplant] plantSize;
	plantSize[obs_ind]=plantSize_obs;
	plantSize[miss_ind]=plantSize_miss;		
	
	for(i in 1:Nplant){			
		//Predicted plant size = intercept + random int field + random int plot	
		plSizeMu[i] = intSize + intSize_field[plotIndex[plantIndex[i]]] + intSize_plot[plantIndex[i]]; 
	}
	
	for(i in 1:Nplot_all){ 		
		//Expected value for hbee visits = intercept + random int + distance + bay position + bay type + time offset
		visitMu_hbee[i] = intVisitHbee + intVisitHbee_field[plotIndex_all[i]] + logTime_all[i] + //Intercepts + offset
			slopeHbeeDistHbee*logHbeeDist_all[i] +  //hbee distance
			slopeLbeeDistHbee*logLbeeDist_all[i] + //lbee distance
			slopeCentHbee*isCent_all[i] + //bay center effect
			slopeFBayHbee*isFBay_all[i] + //F bay effect 				
			slopePlsizeHbee*(intSize + intSize_field[plotIndex_all[i]] + intSize_plot[i]);
	}	

}
	
model {
	vector[2] bernLL_hbee; //pre-calculate LL for zero inflation process
	bernLL_hbee[1]=bernoulli_lpmf(0|zeroVisHbeeTheta); //LL of no extra zero
	bernLL_hbee[2]=bernoulli_lpmf(1|zeroVisHbeeTheta); //LL of extra zero		
	
	//Likelihood
	for(i in 1:Nplot_all){ //Zero-inflated negbin for hbee visitation frequency
		if(hbeeVis_all[i]==0)
			target += log_sum_exp(bernLL_hbee[2],bernLL_hbee[1]+neg_binomial_2_log_lpmf(0|visitMu_hbee[i],visitHbeePhi));
		else
			target += bernLL_hbee[1]+neg_binomial_2_log_lpmf(hbeeVis_all[i]|visitMu_hbee[i],visitHbeePhi);			
		//NB appears to work fine for lbee vis (zero inflation term adds nothing)
	}		
	// lbeeVis_all ~ neg_binomial_2_log(visitMu_lbee,visitLbeePhi); //Lbee visitation rate
	// pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	// podCount ~ binomial_logit(flwCount,flwSurv); //Flower survival
	// seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	// seedMass ~ lognormal(seedWeightMu,sigmaSeedWeight); //Weight per seed
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size
			
	// Priors
	slopePlsizeHbee ~ normal(0,0.5); //Slope of claim
	// Hbee Visitation - informative priors
	intVisitHbee ~ normal(2.5,1); //Intercept	
	slopeHbeeDistHbee ~ normal(-0.1,0.5); //Slope of distance effect on hbee visits
	slopeLbeeDistHbee ~ normal(0.4,0.5); //Effect of leafcutter shelter distance	
	slopeCentHbee ~ normal(-0.3,0.5); //Effect of center of bay
	slopeFBayHbee ~ normal(-0.1,0.5); //Effect of female bay
	sigmaHbeeVisField ~ gamma(1.5,5); //Sigma for random field 
	visitHbeePhi ~ gamma(3.5,5); //Dispersion parameter		
	zeroVisHbeeTheta ~ beta(3,7); //Zero-inflation parameter
	intVisitHbee_field ~ normal(0,sigmaHbeeVisField); //Random field int
		
	//Plant size - informative priors
	intSize ~ normal(1,1); //Intercept
	sigmaPlSize_field ~ gamma(3,10); //Sigma for random field 
	sigmaPlSize_plot ~ gamma(2,10); //Sigma for random plot
	sigmaPlSize ~ gamma(6,10); //Sigma for residual	
	intSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot	
}

generated quantities{	
}
