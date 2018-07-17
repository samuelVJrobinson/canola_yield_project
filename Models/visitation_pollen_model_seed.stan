data {
	//Field level
	int Nfield; //Number of fields	
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 	
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
	vector[Nplot_extra] hbee_dist_extra; 
	int hbeeVis_extra[Nplot_extra]; 
	vector[Nplot_extra] lbee_dist_extra;
	int lbeeVis_extra[Nplot_extra];
	int isCent_extra[Nplot_extra];
	int isFBay_extra[Nplot_extra];
	vector[Nplot_extra] totalTime_extra; 
	
	//Flower level
	int Nflw; //Number of flowers
	int flowerIndex[Nflw]; //Index for flowers - which plot?	
	int pollenCount[Nflw]; //Pollen count
}
transformed data {
	//Amalgamated visitation data
	int Nfield_all = Nfield + Nfield_extra; //Number of fields
	int Nplot_all = Nplot + Nplot_extra; //Number of plots
	int plotIndex_all[Nplot_all]; //Plot index - which field is plot from?
	vector[Nplot_all] hbee_dist_all; //Distance from hbee hives
	int hbeeVis_all[Nplot_all]; //Visits from hbees
	vector[Nplot_all] lbee_dist_all; //Distance from lbee shelters
	int lbeeVis_all[Nplot_all]; //Visits from lbees
	int isCent_all[Nplot_all]; //Is plot from center of bay?
	int isFBay_all[Nplot_all]; //Is plot from female bay?
	vector[Nplot_all] totalTime_all; //Minutes of observation/10
	
	plotIndex_all[1:Nplot] = plotIndex;
	plotIndex_all[Nplot+1:Nplot_all] = plotIndex_extra;
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
	real slopeDistHbee; //Slope of distance/100
	real slopeCentHbee; //Effect of bay position (center)
	real slopeFBayHbee; //Effect of female bay
	real slopeDistLbee; //Effect of leafcutter distance - competition
	real<lower=0.2> sigmaHbeeVisField; //SD of field random intercepts
	real<lower=0.2> visitHbeePhi; //Dispersion parameter	
	vector[Nfield_all] intVisitHbee_field; //field-level random intercepts
	real<lower=0,upper=1> zeroVisHbeeTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. 
	
	//lbee Visitation
	real intVisitLbee; //Intercept - R0 in SSasymp
	real slopeHbeeDistLbee; //Slope of distance
	real AsymLbee; //Horizontal asymptote		
	real lrcLbee; //ln(rate constant)
	// real pwrLbee; //pwr term for Weibull	
	real slopeCentLbee; //Effect of bay position (center)
	real slopeFBayLbee; //Effect of female bay	
	real<lower=0.2> sigmaLbeeVisField; //SD of field random intercepts
	real<lower=0.2> visitLbeePhi; //Dispersion parameter	
	vector[Nfield_all] intVisitLbee_field; //field-level random intercepts for Asym
	real<lower=0,upper=1> zeroVisLbeeTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. 	
	
	//Pollen deposition
	real intPol; //Intercept
	vector<Nfield> intPol_field; //Field-level random intercept
	vector<Nplot> intPol_plot; //Plot-level random intercept
	real slopeHbeePol; //Slope of hbee visits 
	real slopeLbeePol; //Slope of lbee visits 
	real pollenPhi; //Dispersion parameter for pollen deposition
}

transformed parameters {			
	//Expected values
	//Plot-level
	vector[Nplot_all] visitMu_lbee; //lbee visits - all plots
	vector[Nplot_all] visitMu_hbee; //hbee visits - all plot
	vector[Nplot] pollenMu_plot; //Plot level pollen
	vector[Nflw] pollenMu; //Expected pollen - flower level
	
	for(i in 1:Nplot_all){ 
		//Expected value for lbee visits = nonlinear distance to shelter + random int + distance to honeybees (edge) + bay position + bay type + time offset	
		// visitMu_lbee[i] = intVisitLbee + intVisitLbee_field[plotIndex_all[i]] + slopeDistLbee*log(lbee_dist_all[i]) + 
			// slopeCentLbee*isCent_all[i] + slopeFBayLbee*isFBay_all[i] + log(totalTime_all[i]); 
		visitMu_lbee[i] = intVisitLbee_field[plotIndex_all[i]] + //Random intercept for each field
			(AsymLbee+((intVisitLbee)-AsymLbee)*exp(-exp(lrcLbee)*sqrt(lbee_dist_all[i]))) + //Nonlinear slope
			slopeHbeeDistLbee*hbee_dist_all[i] +
			slopeCentLbee*isCent_all[i] + slopeFBayLbee*isFBay_all[i] + log(totalTime_all[i]); 		
		//Expected value for hbee visits = intercept + random int + distance + bay position + bay type + time offset
		visitMu_hbee[i] = intVisitHbee + intVisitHbee_field[plotIndex_all[i]] + slopeDistHbee*log(hbee_dist_all)[i] + slopeDistLbee*log(lbee_dist_all[i]) +
			slopeCentHbee*isCent_all[i] + slopeFBayHbee*isFBay_all[i] + log(totalTime_all[i]); 			
	}
	
	for(i in 1:Nplot){
		//Expected value for plot-level pollination. Models for visitation aren't very good, so using actual visits for now.	
		polMu_plot[i] = intPol + intPol_field[plotIndex[i]] + intPol_plot[i] + 
			slopeHbeePol*(hbeeVis_all[i]/totalTime_all[i]) + slopeLbeePol*(lbeeVis_all[i]/totalTime_all[i]);
	}
	
	for(i in 1:Nflw)
		pollenMu[i] = polMu_plot[flowerIndex[i]]; //Assigns to vector
}
	
model {
	vector[2] bernLL_hbee; //pre-calculate LL for zero inflation process
	vector[2] bernLL_lbee; //pre-calculate LL for zero inflation process	
	bernLL_hbee[1]=bernoulli_lpmf(0|zeroVisHbeeTheta); //LL of no extra zero
	bernLL_hbee[2]=bernoulli_lpmf(1|zeroVisHbeeTheta); //LL of extra zero	
	bernLL_lbee[1]=bernoulli_lpmf(0|zeroVisLbeeTheta); //LL of no extra zero
	bernLL_lbee[2]=bernoulli_lpmf(1|zeroVisLbeeTheta); //LL of extra zero	
	
	//Likelihood
	for(i in 1:Nplot_all){ //Zero-inflated negbin for hbee visitation frequency
		if(hbeeVis_all[i]==0)
			target += log_sum_exp(bernLL_hbee[2],bernLL_hbee[1]+neg_binomial_2_log_lpmf(0|visitMu_hbee[i],visitHbeePhi));
		else
			target += bernLL_hbee[1]+neg_binomial_2_log_lpmf(hbeeVis_all[i]|visitMu_hbee[i],visitHbeePhi);	
		
		if(lbeeVis_all[i]==0) //Zero-inflated negbin for lbee visitation frequency
			target += log_sum_exp(bernLL_lbee[2],bernLL_lbee[1]+neg_binomial_2_log_lpmf(0|visitMu_lbee[i],visitLbeePhi));
		else
			target += bernLL_lbee[1]+neg_binomial_2_log_lpmf(lbeeVis_all[i]|visitMu_lbee[i],visitLbeePhi);			
	}	
	// lbeeVis_all ~ neg_binomial_2_log(visitMu_lbee,visitLbeePhi);
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	// podCount ~ binomial_logit(flwCount,flwSurv); //Flower survival
	// seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	// seedMass ~ lognormal(seedWeightMu,sigmaSeedWeight); //Weight per seed
	// plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size
		
	//Priors
	//Hbee Visitation - informative priors
	intVisitHbee ~ normal(2.5,1); //Intercept
	sigmaHbeeVisField ~ gamma(1.5,5); //Sigma for random field 
	slopeDistHbee ~ normal(-0.1,0.5); //Slope of distance effect on hbee visits
	slopeCentHbee ~ normal(-0.3,0.5); //Effect of center of bay
	slopeFBayHbee ~ normal(-0.1,0.5); //Effect of female bay
	slopeDistLbee ~ normal(0.3,1); //Effect of leafcutter shelter distance
	visitHbeePhi ~ gamma(3.5,5); //Dispersion parameter	
	intVisitHbee_field ~ normal(0,sigmaHbeeVisField); //Random field int
	zeroVisHbeeTheta ~ beta(3,7); //Zero-inflation parameter
	
	//Lbee Visitation - informative priors
	intVisitLbee ~ normal(4,1); //Intercept
	AsymLbee ~ normal(0.5,1); //Horizontal asymptote	
	lrcLbee ~ normal(-1,1); //ln(rate constant)	
	sigmaLbeeVisField ~ gamma(2,2); //Sigma for random field 
	slopeHbeeDistLbee ~ normal(-0.1,0.5); //Slope of honeybee distance on lbee visits
	slopeCentLbee ~ normal(1,1); //Effect of center of bay
	slopeFBayLbee ~ normal(0.2,1); //Effect of female bay
	visitLbeePhi ~ gamma(2.5,5); //Dispersion parameter	
	intVisitLbee_field ~ normal(0,sigmaLbeeVisField); //Random field int
	zeroVisLbeeTheta ~ beta(1.1,9); //Zero-inflation parameter	
		
	// //Pollen deposition - informative priors
	// intPollen ~ normal(5.5,1); //Intercept	
	// sigmaPolField ~ gamma(2,4); //Sigma for random field
	// sigmaPolPlot ~ gamma(4,10); //Sigma for random plot	
	// slopeVisitPol ~ normal(0,1); //hbee Visitation effect
	// pollenPhi ~ gamma(7,10); //Dispersion parameter
	// intPollen_field ~ normal(0,sigmaPolField); //Random field int
	// intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int
	
	// //Flower survival - informative priors
	// intFlwSurv ~ normal(1.2,1); //Intercept
	// slopeVisitSurv ~ normal(0,1); //Slope of hbee visits
	// slopePolSurv ~ normal(-0.5,1); //Slope of pollen deposition
	// slopePlSizeSurv ~ normal(0,3); //Slope of plant size
	// sigmaFlwSurv_field ~ gamma(1.5,5); //SD of field random effect
	// sigmaFlwSurv_plot ~ gamma(1.5,5); //SD of plot random effect	
	// intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	// intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts	
	
	// //Seed count - informative priors
	// intSeedCount ~ normal(3,1); //Intercept
	// slopeVisitSeedCount ~ normal(0,1); //Slope of hbee visits
	// slopePolSeedCount ~ normal(0,1); //Slope of pollen deposition
	// slopePlSizeCount ~ normal(0,2); //Slope of plant size
	// seedCountPhi ~ gamma(21,1); //Dispersion parameter
	// sigmaSeedCount_field ~ gamma(2,20); //SD of field random effect
	// sigmaSeedCount_plot ~ gamma(2,20); //SD of plot random effect
	// sigmaSeedCount_plant ~ gamma(2,20); //SD of plant random effect	
	// intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	// intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	// intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts
	
	// //Weight per seed - informative priors
	// intSeedWeight ~ normal(1,1); //Intercept
	// slopeVisitSeedWeight ~ normal(0,1); //Slope of hbee visits
	// slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	// slopeSeedCount ~ normal(0.005,.1); //Slope of seed count
	// slopePlSizeWeight ~ normal(0,2); //Slope of plant size
	// sigmaSeedWeight ~ gamma(3,10); //SD of seed weight
	// sigmaSeedWeight_field ~ gamma(2,10); //SD of field random effect	
	// sigmaSeedWeight_plot ~ gamma(2,10); //SD of plot random effect
	// sigmaSeedWeight_plant ~ gamma(2,10); //SD of plant random effect		
	// intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	// intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	// intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts

	// //Plant size - informative priors
	// intSize ~ normal(0,1); //Intercept
	// sigmaPlSize_field ~ gamma(2,6); //Sigma for random field 
	// sigmaPlSize_plot ~ gamma(2,7); //Sigma for random plot
	// sigmaPlSize ~ gamma(3,5); //Sigma for residual	
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
	// int predPodCount[Nplant_obs]; //Simulated surviving pods
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
		lbeeVis_resid=lbeeVis_resid+fabs((exp(visitMu_lbee[i])*(1-zeroVisLbeeTheta))-lbeeVis_all[i]); //Residual for actual value
		if(bernoulli_rng(zeroVisLbeeTheta)==1) //If zeroVisLbeeTheta generates an extra zero
			predLbeeVis_all[i] = 0; //Predicted value is automatically zero
		else //Otherwise
			predLbeeVis_all[i] = neg_binomial_2_log_rng(visitMu_lbee[i],visitLbeePhi); //Predicted value drawn from neg.bin		
		predLbeeVis_resid=predLbeeVis_resid+fabs((exp(visitMu_lbee[i])*(1-zeroVisLbeeTheta))-predLbeeVis_all[i]); //Residual for predicted value		
	}
	
	// for(i in 1:Nflw)
		// predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulated pollen counts
		
	// for(i in 1:Nplant_obs)
		// predPodCount[i] = binomial_rng(flwCount[i],inv_logit(flwSurv[i])); //Simulated surviving pods
	
	// for(i in 1:Npod){ //For each pod
		// predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); //Seed count per pod
		// predSeedMass[i] = lognormal_rng(seedWeightMu[i],sigmaSeedWeight); //Weight per seed
	// }	
}
