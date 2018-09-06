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
	//Claim: lbee ~ plSize
	real slopePlsizeLbee;
	
	//lbee Visitation
	real intVisitLbee; //Intercept 
	real slopeHbeeDistLbee; //Slope of honeybee distance (field edge)
	real slopeLbeeDistLbee; //Slope of leafcutter distance (shelter)		
	real slopeCentLbee; //Effect of bay position (center)
	real slopeFBayLbee; //Effect of female bay	
	real slopeStocking; //Effect of half-stocking leafcutter bees
	real slopeCentHbeeDistLbee; //Bay position : honeybee distance interaction term
	real slopeStockingHbeeDistLbee; //Half-stocking:hbee distance interaction	
	real<lower=0.1> sigmaLbeeVisField; //SD of field random intercepts
	real<lower=0.1> visitLbeePhi; //Dispersion parameter	
	vector[Nfield_all] intVisitLbee_field; //field-level random intercepts	
	
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
	vector[Nplot_all] visitMu_lbee; //lbee visits - all plots
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
		//Expected value for lbee visits = intercept + random int + distance to shelter + distance to honeybees (edge) + bay position + bay:hbee dist + stocking:hbee dist + bay type + time offset	
		visitMu_lbee[i] = intVisitLbee + intVisitLbee_field[plotIndex_all[i]] + logTime_all[i] + //intercepts + offset
			slopeLbeeDistLbee*logLbeeDist_all[i] + //lbee distance
			slopeHbeeDistLbee*logHbeeDist_all[i] + //hbee distance
			slopeCentLbee*isCent_all[i] + //bay center
			slopeStocking*lbeeStocking_all[i]+	//half-stocking
			slopeFBayLbee*isFBay_all[i] + //F bay
			slopeCentHbeeDistLbee*isCent_all[i]*logHbeeDist_all[i] + //hbee dist: bay center interaction			
			slopeStockingHbeeDistLbee*lbeeStocking_all[i]*logHbeeDist_all[i] +  //hbee dist: half stocking interaction
			slopePlsizeLbee*(intSize + intSize_field[plotIndex_all[i]] + intSize_plot[i]); //Claim
	}
				
		
}
	
model {
	//Likelihood
	lbeeVis_all ~ neg_binomial_2_log(visitMu_lbee,visitLbeePhi); //Lbee visitation rate	
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	
			
	// Priors		
	//Claim
	slopePlsizeLbee ~ normal(0,1);
	
	// Lbee Visitation - informative priors	
	intVisitLbee ~ normal(4,4); //Intercept	
	slopeHbeeDistLbee ~ normal(-0.15,0.3); //Slope of honeybee distance on lbee visits
	slopeLbeeDistLbee ~ normal(-0.8,0.3); //Slope of shelter distance on lbee visits
	slopeCentLbee ~ normal(-0.6,0.3); //Effect of center of bay
	slopeFBayLbee ~ normal(0,0.5); //Effect of female bay
	slopeStocking ~ normal(0,0.5); //Effect of half-stocking
	slopeCentHbeeDistLbee ~ normal(-0.2,0.3); //Bay center: hbee distance interaction
	slopeStockingHbeeDistLbee ~ normal(0.2,0.3); //Half-stocking: hbee distance interaction		
	sigmaLbeeVisField ~ gamma(2,2); //Sigma for random field 
	visitLbeePhi ~ gamma(2,5); //Dispersion parameter	
	intVisitLbee_field ~ normal(0,sigmaLbeeVisField); //Random field intercepts	
	
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
