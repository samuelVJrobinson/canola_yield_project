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
	vector[Nplot_all] hbee_dist_all; //Distance from hbee hives
	int hbeeVis_all[Nplot_all]; //Visits from hbees
	vector[Nplot_all] lbee_dist_all; //Distance from lbee shelters
	int lbeeVis_all[Nplot_all]; //Visits from lbees
	int isCent_all[Nplot_all]; //Is plot from center of bay?
	int isFBay_all[Nplot_all]; //Is plot from female bay?
	vector[Nplot_all] totalTime_all; //Minutes of observation/10
	
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
}

parameters {
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
	
	//lbee Visitation
	real intVisitLbee; //Intercept - R0 in SSasymp
	real slopeHbeeDistLbee; //Slope of honeybee distance (field edge)
	real slopeLbeeDistLbee; //Slope of leafcutter distance (shelter)	
	// real AsymLbee; //Horizontal asymptote - nonlinear form
	// real lrcLbee; //ln(rate constant)	
	real slopeCentLbee; //Effect of bay position (center)
	real slopeCentHbeeDistLbee; //Bay position : honeybee distance interaction term
	real slopeFBayLbee; //Effect of female bay	
	real<lower=0.1> sigmaLbeeVisField; //SD of field random intercepts
	real<lower=0.1> visitLbeePhi; //Dispersion parameter	
	vector[Nfield_all] intVisitLbee_field; //field-level random intercepts	
	
	// // Pollen deposition
	// real intPol; //Intercept	
	// real slopeHbeePol; //Slope of hbee visits 
	// real slopeLbeePol; //Slope of lbee visits 
	// real pollenPhi; //Dispersion parameter for pollen deposition
	// real<lower=0> sigmaPolField; //Sigma for field-level intercept
	// real<lower=0> sigmaPolPlot; //Sigma for plot-level intercept
	// vector[Nfield] intPol_field; //Field-level random intercept
	// vector[Nplot] intPol_plot; //Plot-level random intercept
	
	// //Flower survival
	// real intFlwSurv; //Intercept	
	// real slopePolSurv; //Slope of pollen deposition
	// real slopePlSizeSurv; //Slope of plant size
	// real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts	
	// real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	// vector[Nfield] intFlwSurv_field; //field-level random intercepts
	// vector[Nplot] intFlwSurv_plot; //plot-level random intercepts
	
	// //Seed count
	// real intSeedCount; //Intercept	
	// real slopePolSeedCount; //Slope of pollen deposition
	// real slopePlSizeCount; //Slope of plant size
	// real<lower=0> seedCountPhi; //Dispersion parameter
	// real<lower=0> sigmaSeedCount_plant; //SD of plant random effect
	// real<lower=0> sigmaSeedCount_plot; //SD of plot random effect
	// real<lower=0> sigmaSeedCount_field; //SD of field random effect
	// vector[Nfield] intSeedCount_field; //field-level random intercepts	
	// vector[Nplot] intSeedCount_plot; //plot-level random intercepts	
	// vector[Nplant] intSeedCount_plant; //plant-level random intercepts
	
	// //Weight per seed
	// real intSeedWeight; //Intercept	
	// real slopePolSeedWeight; //Slope of pollen deposition
	// real slopeSeedCount; //Slope of seed count
	// real slopePlSizeWeight; //Slope of plant size
	// real<lower=0> sigmaSeedWeight; //SD of seed weight
	// real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect
	// real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect
	// real<lower=0> sigmaSeedWeight_field; //SD of field random effect	
	// vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	// vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	// vector[Nplant] intSeedWeight_plant; //plant-level random intercepts		
	
	// //Plant size
	// vector[Nplant_miss] plantSize_miss; //Vector for imputing missing values	
	// real intSize; //Global intercept
	// real<lower=0> sigmaPlSize_field; //Sigma for field
	// real<lower=0> sigmaPlSize_plot; //Sigma for plot
	// real<lower=0> sigmaPlSize; //Sigma for within-plot	(residual)
	// vector[Nfield] intSize_field; //Random intercept for field
	// vector[Nplot] intSize_plot; //Random intercept for plot		
}

transformed parameters {			
	//Expected values
	//Plot-level
	vector[Nplot_all] visitMu_lbee; //lbee visits - all plots
	vector[Nplot_all] visitMu_hbee; //hbee visits - all plot
	// vector[Nplot] pollenMu_plot; //Plot level pollen
	// vector[Nflw] pollenMu; //Expected pollen - flower level
	// vector[Nplot] flwSurvPlot; //Plot-level flower survival
	// vector[Nplant_obs] flwSurv; //Flower survival rate (logit)
	// vector[Nplant] seedCountMuPlant; //Plant-level seed count
	// vector[Npod] seedCountMu; //Pod-level seed counts	
	// vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	// vector[Npod] seedWeightMu; //Pod-level weight per seed	
	// vector[Nplant] plSizeMu; //Plant size

	//Imputed missing data;
	// vector[Nplant] plantSize;
	// plantSize[obs_ind]=plantSize_obs;
	// plantSize[miss_ind]=plantSize_miss;		
	
	for(i in 1:Nplot_all){ 		
		//Expected value for hbee visits = intercept + random int + distance + bay position + bay type + time offset
		visitMu_hbee[i] = intVisitHbee + intVisitHbee_field[plotIndex_all[i]] + 
			slopeHbeeDistHbee*log(hbee_dist_all[i]) + 
			slopeLbeeDistHbee*log(lbee_dist_all[i]) +
			slopeCentHbee*isCent_all[i] + slopeFBayHbee*isFBay_all[i] + log(totalTime_all[i]); 			
		//Expected value for lbee visits = intercept + random int + distance to shelter + distance to honeybees (edge) + bay position + bay:hbee dist + bay type + time offset	
		visitMu_lbee[i] = intVisitLbee + intVisitLbee_field[plotIndex_all[i]] + slopeLbeeDistLbee*log(lbee_dist_all[i]) +
			slopeHbeeDistLbee*log(hbee_dist_all[i]) + slopeCentLbee*isCent_all[i] + slopeCentHbeeDistLbee*isCent_all[i]*log(hbee_dist_all[i]) + 
			slopeFBayLbee*isFBay_all[i] + log(totalTime_all[i]); 
		// visitMu_lbee[i] = //Random intercept for each field is rolled in with AsymLbee parameter
			// ((AsymLbee+intVisitLbee_field[plotIndex_all[i]])+((intVisitLbee)-(AsymLbee+intVisitLbee_field[plotIndex_all[i]]))*
			// exp(-exp(lrcLbee)*sqrt(lbee_dist_all[i]))) + //Nonlinear slope
			// slopeHbeeDistLbee*log(hbee_dist_all[i]) + //Honeybee distance
			// slopeCentLbee*isCent_all[i] + slopeFBayLbee*isFBay_all[i] + log(totalTime_all[i]); 	
	}
	
	// for(i in 1:Nplot){
		// //Expected value for plot-level pollination. Models for visitation aren't very good, so using actual visits for now.	
		// pollenMu_plot[i] = intPol + intPol_field[plotIndex[i]] + intPol_plot[i] + 
			// slopeHbeePol*(hbeeVis_all[i]/totalTime_all[i]) + slopeLbeePol*(lbeeVis_all[i]/totalTime_all[i]);
			
		// //Flower survival = intercept + random int field + random int plot + pollen deposition/1000
		// flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + intFlwSurv_plot[i] + 
			// slopePolSurv*(exp(pollenMu_plot[i])/1000);
	// }
	
	// for(i in 1:Nflw) 
		// pollenMu[i] = pollenMu_plot[flowerIndex[i]]; //Assigns plot level pollen mu to Nflw long vector
		
	// for(i in 1:Nplant_obs) //For each observed plant 		
		// flwSurv[i] = flwSurvPlot[plantIndex[i]] + slopePlSizeSurv*plantSize[i]; //Plot-level plant survival + individual size effect	
	
	// for(i in 1:Nplant){
		// //Seed count per pod = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000 + plant size
		// seedCountMuPlant[i] = intSeedCount + intSeedCount_field[plotIndex[plantIndex[i]]] + intSeedCount_plot[plantIndex[i]] + intSeedCount_plant[i] + 
			// slopePolSeedCount*(exp(pollenMu_plot[plantIndex[i]])/1000) + slopePlSizeCount*plantSize[i];
			
		// //Weight per seed = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition/1000
		// seedWeightPlantMu[i] = intSeedWeight + intSeedWeight_field[plotIndex[plantIndex[i]]] + intSeedWeight_plot[plantIndex[i]] + intSeedWeight_plant[i] +
			// slopePolSeedWeight*(exp(pollenMu_plot[plantIndex[i]])/1000)+ slopePlSizeWeight*plantSize[i];
			
		// //Predicted plant size = intercept + random int field + random int plot	
		// plSizeMu[i] = intSize + intSize_field[plotIndex[plantIndex[i]]] + intSize_plot[plantIndex[i]]; 
	// }
	
	// for(i in 1:Npod){ //For each pod
		// seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
		// seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+slopeSeedCount*seedCount[i]; //Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
	// }		
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
	lbeeVis_all ~ neg_binomial_2_log(visitMu_lbee,visitLbeePhi); //Lbee visitation rate
	// pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	// podCount ~ binomial_logit(flwCount,flwSurv); //Flower survival
	// seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	// seedMass ~ lognormal(seedWeightMu,sigmaSeedWeight); //Weight per seed
	// plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size
			
	// Priors
	// Hbee Visitation - informative priors
	intVisitHbee ~ normal(2.5,1); //Intercept
	sigmaHbeeVisField ~ gamma(1.5,5); //Sigma for random field 
	slopeHbeeDistHbee ~ normal(-0.1,0.5); //Slope of distance effect on hbee visits
	slopeCentHbee ~ normal(-0.3,0.5); //Effect of center of bay
	slopeFBayHbee ~ normal(-0.1,0.5); //Effect of female bay
	slopeLbeeDistHbee ~ normal(0.4,0.5); //Effect of leafcutter shelter distance
	visitHbeePhi ~ gamma(3.5,5); //Dispersion parameter	
	intVisitHbee_field ~ normal(0,sigmaHbeeVisField); //Random field int
	zeroVisHbeeTheta ~ beta(3,7); //Zero-inflation parameter
	
	// Lbee Visitation - informative priors
	intVisitLbee ~ normal(4,4); //Intercept
	// AsymLbee ~ normal(0.5,1); //Horizontal asymptote	
	// lrcLbee ~ normal(-1,1); //ln(rate constant)	
	sigmaLbeeVisField ~ gamma(2,2); //Sigma for random field 
	slopeHbeeDistLbee ~ normal(-0.15,0.3); //Slope of honeybee distance on lbee visits
	slopeLbeeDistLbee ~ normal(-0.8,0.3); //Slope of shelter distance on lbee visits
	slopeCentLbee ~ normal(-0.6,0.3); //Effect of center of bay
	slopeCentHbeeDistLbee ~ normal(-0.2,0.3); //Bay center: hbee distance interaction
	slopeFBayLbee ~ normal(0,0.5); //Effect of female bay
	visitLbeePhi ~ gamma(2,5); //Dispersion parameter	
	intVisitLbee_field ~ normal(0,sigmaLbeeVisField); //Random field int
	zeroVisLbeeTheta ~ beta(1.1,9); //Zero-inflation parameter	
	
	// // Pollen deposition - informative priors
	// intPol ~ normal(2.5,1); //Intercept	
	// sigmaPolField ~ gamma(9,10); //Sigma for random field
	// sigmaPolPlot ~ gamma(4,4); //Sigma for random plot	
	// slopeHbeePol ~ normal(0,0.05); //hbee Visitation effect
	// slopeLbeePol ~ normal(0.015,0.05); //lbee Visitation effect
	// pollenPhi ~ gamma(4,5); //Dispersion parameter
	// intPol_field ~ normal(0,sigmaPolField); //Random field int
	// intPol_plot ~ normal(0,sigmaPolPlot); //Random plot int
	
	// // Flower survival 
	// intFlwSurv ~ normal(1,1); //Intercept	
	// slopePolSurv ~ normal(2,1); //Slope of pollen deposition
	// slopePlSizeSurv ~ normal(0.22,0.5); //Slope of plant size
	// sigmaFlwSurv_field ~ gamma(3,10); //SD of field random effect
	// sigmaFlwSurv_plot ~ gamma(5,10); //SD of plot random effect	
	// intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	// intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts			
	
	// // Seed count 
	// intSeedCount ~ normal(2.7,1); //Intercept	
	// slopePolSeedCount ~ normal(2,1); //Slope of pollen deposition
	// slopePlSizeCount ~ normal(0.1,0.05); //Slope of plant size
	// seedCountPhi ~ gamma(10,4); //Dispersion parameter
	// sigmaSeedCount_field ~ gamma(1.6,10); //SD of field random effect
	// sigmaSeedCount_plot ~ gamma(2,10); //SD of plot random effect
	// sigmaSeedCount_plant ~ gamma(1,10); //SD of plant random effect	
	// intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	// intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	// intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts
	
	// // Weight per seed 
	// intSeedWeight ~ normal(1.3,1); //Intercept	
	// slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	// slopeSeedCount ~ normal(-0.01,0.05); //Slope of seed count
	// slopePlSizeWeight ~ normal(0.05,0.1); //Slope of plant size
	// sigmaSeedWeight ~ gamma(4,10); //SD of seed weight
	// sigmaSeedWeight_field ~ gamma(1.5,10); //SD of field random effect	
	// sigmaSeedWeight_plot ~ gamma(1.2,10); //SD of plot random effect
	// sigmaSeedWeight_plant ~ gamma(1.8,10); //SD of plant random effect		
	// intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	// intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	// intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts

	// //Plant size - informative priors
	// intSize ~ normal(1,1); //Intercept
	// sigmaPlSize_field ~ gamma(3,10); //Sigma for random field 
	// sigmaPlSize_plot ~ gamma(2,10); //Sigma for random plot
	// sigmaPlSize ~ gamma(6,10); //Sigma for residual	
	// intSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	// intSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot	
}

generated quantities{
	int predHbeeVis_all[Nplot_all]; //Predicted hbee visits
	real hbeeVis_resid=0; //Residual of hbeeVis
	real predHbeeVis_resid=0; //Residual of generated hbeeVis
	int predLbeeVis_all[Nplot_all]; //Predicted lbee visits
	real lbeeVis_resid=0; //Residual of lbeeVis
	real predLbeeVis_resid=0; //Residual of generated lbeeVis	
	
	// int predPollenCount[Nflw]; //Predicted pollen counts
	// int predPodCount[Nplant]; //Simulated surviving pods
	// int predSeedCount[Npod]; //Simulated seeds per pod
	// vector[Npod] predSeedMass; //Simulated seed weight
	
	for(i in 1:Nplot_all){
		//Predicted hbee visits
		hbeeVis_resid=hbeeVis_resid+fabs((exp(visitMu_hbee[i])*(1-zeroVisHbeeTheta))-hbeeVis_all[i]); //Residual for actual value
		if(bernoulli_rng(zeroVisHbeeTheta)==1) //If zeroVisHbeeTheta generates an extra zero
			predHbeeVis_all[i] = 0; //Predicted value is automatically zero
		else //Otherwise
			predHbeeVis_all[i] = neg_binomial_2_log_rng(visitMu_hbee[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		predHbeeVis_resid=predHbeeVis_resid+fabs((exp(visitMu_hbee[i])*(1-zeroVisHbeeTheta))-predHbeeVis_all[i]); //Residual for predicted value
		
		//Predicted lbee visits
		lbeeVis_resid=lbeeVis_resid+fabs((exp(visitMu_lbee[i]))-lbeeVis_all[i]); //Residual for actual value
		predLbeeVis_all[i] = neg_binomial_2_log_rng(visitMu_lbee[i],visitLbeePhi); //Predicted value drawn from neg.bin		
		predLbeeVis_resid=predLbeeVis_resid+fabs((exp(visitMu_lbee[i]))-predLbeeVis_all[i]); //Residual for predicted value		
	}
	
	// for(i in 1:Nflw)
		// predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulated pollen counts
		
	// for(i in 1:Nplant)
		// predPodCount[i] = binomial_rng(flwCount[i],inv_logit(flwSurv[i])); //Simulated surviving pods
	
	// for(i in 1:Npod){ //For each pod
		// predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); //Seed count per pod
		// predSeedMass[i] = lognormal_rng(seedWeightMu[i],sigmaSeedWeight); //Weight per seed
	// }	
}
