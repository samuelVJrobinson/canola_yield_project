data { 
	//Field level
	int Nfield; //Number of fields		
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 		
	int lbeeStocking[Nplot]; //Half leafcutter stocking? (double-tent treatment)
	matrix[Nplot,2] lbeeStocking2; //Model matrix for leafcutter stocking (2x tent or 2x tent + extra bees)
	int is2016[Nplot]; //Was field from 2016?	
	vector[Nplot] hbee_dist; //Distance from hbee hives
	int hbeeVis[Nplot]; //Number of hbee visits to plot
	vector[Nplot] lbee_dist; //Distance to lbee shelters
	int lbeeVis[Nplot]; //Number of lbee visits to plot
	int isCent[Nplot]; //Is plot in center of bay?
	int isMBay[Nplot]; //Is plot in male bay?
	vector[Nplot] totalTime; //Minutes spent observing/10 	
	vector[Nplot] polCountPlot; //Plot-level average pollen count	
	//Missing data from plot level
	//Flower density (fls/m2)
	int Nplot_flsObs; //Number of plots where flower density observed
	int Nplot_flsMiss; //Number of plots where flower density missing (2 plots)
	vector[Nplot_flsObs] flDens_obs; //Observed flower density (flowers/m2)
	int<lower=1,upper=Nplot> obsFls_ind[Nplot_flsObs]; //Index for observed flower density
	int<lower=1,upper=Nplot> missFls_ind[Nplot_flsMiss]; //Index for missing flower density	
	//Plant density (stems/m2)
	int Nplot_densObs; //Number of plots where plant density was observed
	int Nplot_densMiss; //Number of plots with plant density missing (about 40% - male bays)
	vector[Nplot_densObs] plDens_obs; //Observed plant density (stems/m2)
	int<lower=1,upper=Nplot> obsPlDens_ind[Nplot_densObs]; //Index for observed plant density
	int<lower=1,upper=Nplot> missPlDens_ind[Nplot_densMiss]; //Index for missing plant density	
	//Index for mapping plot index onto female plot index
	int Nplot_F; //Number of female plots
	int<lower=0,upper=Nplot_F> plotIndex_F[Nplot]; //Index for female plots (which female plot j does plot i belong to? - some zero indices)
	int<lower=1,upper=Nplot> plotIndex_F2[Nplot_F]; //Reverse index (which plot i does female plot j belong to?)
	
	//Extra plots/fields from Riley (visitation data)
	int Nfield_extra;		
	int Nplot_extra; 
	int is2016_extra[Nplot_extra]; //Was field from 2016?
	int lbeeStocking_extra[Nplot_extra]; //Half leafcutter stocking? (double-tent treatment)
	matrix[Nplot_extra,2] lbeeStocking2_extra; //Model matrix for leafcutter stocking (2x tent or 2x tent + extra bees)
	int plotIndex_extra[Nplot_extra]; 	
	vector[Nplot_extra] hbee_dist_extra; 
	int hbeeVis_extra[Nplot_extra]; 
	vector[Nplot_extra] lbee_dist_extra;
	int lbeeVis_extra[Nplot_extra];
	int isCent_extra[Nplot_extra];
	int isMBay_extra[Nplot_extra];
	vector[Nplot_extra] totalTime_extra;
	//Missing flower density data
	int Nplot_flsObs_extra; //Number of plots where flower density observed
	int Nplot_flsMiss_extra; //Number of plots where flower density missing
	vector[Nplot_flsObs_extra] flDens_obs_extra; //Observed flower density (flowers/m2)
	int<lower=1,upper=Nplot_extra> obsFls_ind_extra[Nplot_flsObs_extra]; //Index for observed flower density
	int<lower=1,upper=Nplot_extra> missFls_ind_extra[Nplot_flsMiss_extra]; //Index for missing flower density			
	
	//Flower level (pollen counts)
	int Nflw; //Number of flowers (stigmas)
	int flowerIndex[Nflw]; //Index for flowers - which plot?	
	int pollenCount[Nflw]; //Pollen count
	
	//Plant level
	int Nplant; //Number of all plants 
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	vector[Nplant] plantSize; //Mass of vegetative tissue (no seeds) (g)
	vector[Nplant] yield; //Observed yield (g seed per plant)		
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?		
			
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
	matrix[Nplot_all,2] lbeeStocking2_all; //Half leafcutter stocking?
	int is2016_all[Nplot_all]; //Is field from 2016?
	vector[Nplot_all] hbee_dist_all; //Distance from hbee hives (edge)
	int hbeeVis_all[Nplot_all]; //Visits from hbees
	vector[Nplot_all] lbee_dist_all; //Distance from lbee shelters
	int lbeeVis_all[Nplot_all]; //Visits from lbees
	int isCent_all[Nplot_all]; //Is plot from center of bay?
	int isMBay_all[Nplot_all]; //Is plot from male bay?
	vector[Nplot_all] totalTime_all; //Minutes of observation/10
	vector[Nplot_all] logTime_all; //Log transform of observation time
	vector[Nplot_all] logHbeeDist_all; //Log-distance from hbee hives (edge)
	vector[Nplot_all] logLbeeDist_all; //Log-distance from lbee shelters
	vector[Nplot_all] logHbeeVis_all; //Log-visitation rate for hbees
	vector[Nplot_all] logLbeeVis_all; //Log-visitation rate for lbees
	vector[Nplot] logPolCountPlot; //Log-average pollen count at plot level
	vector[Nplant] logFlwCount; //Log-flower count
	vector[Nplant] logitFlwSurv; //(logit) proportion flower survival	
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)	
	
	//Assign values
	plotIndex_all[1:Nplot] = plotIndex; 
	plotIndex_all[Nplot+1:Nplot_all] = plotIndex_extra;
	lbeeStocking_all[1:Nplot] = lbeeStocking; //Lbee stocking
	lbeeStocking_all[Nplot+1:Nplot_all] = lbeeStocking_extra;
	lbeeStocking2_all[1:Nplot,] = lbeeStocking2; //Other lbee stocking model matrix
	lbeeStocking2_all[Nplot+1:Nplot_all,] = lbeeStocking2_extra; 	
	is2016_all[1:Nplot] = is2016; //Year
	is2016_all[Nplot+1:Nplot_all] = is2016_extra;	
	hbee_dist_all[1:Nplot] = hbee_dist; //Hbee distance
	hbee_dist_all[Nplot+1:Nplot_all] = hbee_dist_extra;
	hbeeVis_all[1:Nplot] = hbeeVis;
	hbeeVis_all[Nplot+1:Nplot_all] = hbeeVis_extra;
	lbee_dist_all[1:Nplot] = lbee_dist;
	lbee_dist_all[Nplot+1:Nplot_all] = lbee_dist_extra;
	lbeeVis_all[1:Nplot] = lbeeVis;
	lbeeVis_all[Nplot+1:Nplot_all] = lbeeVis_extra;
	isCent_all[1:Nplot] = isCent;
	isCent_all[Nplot+1:Nplot_all] = isCent_extra;
	isMBay_all[1:Nplot] = isMBay;
	isMBay_all[Nplot+1:Nplot_all] = isMBay_extra;
	totalTime_all[1:Nplot] = totalTime;
	totalTime_all[Nplot+1:Nplot_all] = totalTime_extra;
	//Transformations
	logHbeeDist_all=log(hbee_dist_all)-mean(log(hbee_dist_all)); //Log-transform and center distances
	logLbeeDist_all=log(lbee_dist_all)-mean(log(lbee_dist_all));
	logTime_all=log(totalTime_all); //Log-transform time	
	for(i in 1:Nplot_all){
		logHbeeVis_all[i] = log(1+(hbeeVis_all[i]/totalTime_all[i])); //Log-transforms observed visitation rates 
		logLbeeVis_all[i] = log(1+(lbeeVis_all[i]/totalTime_all[i]));			
	}
	logPolCountPlot=polCountPlot; //Log-transforms average pollen count per plot - not sure if this is needed
	
	for(i in 1:Nplant){
		logFlwCount[i] = log(flwCount[i]); //Log flower count per plant	
		//Necessary for promoting integers to reals. Otherwise does integer division.
		logitFlwSurv[i] = podCount[i]; 
		logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods		
		if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
			logitFlwSurv[i]=0.01;
		else if(logitFlwSurv[i]>=1)
			logitFlwSurv[i]=0.99;				
	}
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
}

parameters {
	//Claim: hbeeVis ~ plDens
	real slopePlSizeHbeeVis;
	real slopePlDensHbeeVis; //Conditioning set
 
	// Plant density - looks OK
	// Vector for imputing missing plant density values (missing values from my data + all of Riley's data)
	vector[Nplot_densMiss] plDens_miss; //My fields
	vector[Nplot_extra] plDens_miss_extra; //Riley's fields	
	real intPlDens; //Global intercept
	real slopeHbeeDistPlDens; //Slope of distance into field	
	real<lower=1e-10> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=1e-10> sigmaPlDens_field; //Sigma for field
	vector[Nfield_all] intPlDens_field; //Random intercept for field
	
	// Plant size - random effects at plot level are very small, and don't converge well
	// Density:distance interaction is basically 0, so leaving it out
	// Normal distribution had marginal PPchecks; t-dist is much better
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density	
	real slopeDistPlSize; //Slope of hbee distance (edge of field has small plants)			
	// real slope2016PlSize; //Effect of 2016 on plant size - p=0.858
	real<lower=1e-05> sigmaPlSize; //Sigma for within-plot (residual)	
	real<lower=1e-05> sigmaPlSize_field; //Sigma for field		
	vector[Nfield_all] intPlSize_field; //Random intercept for field	
	real nuPlSize; //exp(nu) for t-distribution
	
	// Flower density per plot
	// slopePlDensFlDens correlated with slopeDistFlDens (r=0.61), and hard to separate effects, so leaving out
	vector[Nplot_flsMiss] flDens_miss; //Vectors for imputing missing values
	vector[Nplot_flsMiss_extra] flDens_extra_miss;
	real intFlDens; //Global intercept
	real slopePlSizeFlDens; //Slope of plant size on flower density
	real slopeMBayFlDens; //Effect of male bay
	real slope2016FlDens; //Effect of 2016
	real slopeDistFlDens;  //Effect of distance from edge	
	// real slopePlDensFlDens; //Effect of planting density 
	real<lower=1e-10> sigmaFlDens; //Sigma for within-field (residual)
	real<lower=1e-10> sigmaFlDens_field; //Sigma for field
	vector[Nfield_all] intFlDens_field; //Random intercept for field
	real nuFlDens; //exp(nu) for t-distribution
 
	// hbee Visitation - random effects at field level weren't converging
	real intVisitHbee; //Intercept	
	real slopeFlDensHbee; //Slope of flower density
	real slopeHbeeDistHbee; //Slope of (log) distance
	real slopeLbeeDistHbee; //Slope of leafcutter distance
	real slopeLbeeHbeeDistHbee; //Interaction b/w leafcutter & honeybee distance
	real slopeLbeeVisHbee; //Direct effect of leafcutter visitation	
	real slopeMBayHbee; //Effect of male bay		
	real slopeCentHbee; //Effect of bay position (center)	
	real<lower=1e-10> visitHbeePhi; //Dispersion parameter		
	real<lower=0,upper=1> zeroVisHbeeTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. 
	real<lower=1e-10> sigmaHbeeVis_field; //Sigma for field
	vector[Nfield_all] intVisitHbee_field; //Sigma for field
}

transformed parameters {			
	//Expected values
	//Plot-level
	vector[Nplot_all] plDensMu; //Expected plant density	
	vector[Nplot_all] plSizePlotMu; //Plot-level plant size
	vector[Nplant] plSizeMu; //Expected plant size	
	vector[Nplot_all] flDensMu; //Expected flower density	
	vector[Nplot_all] visitMu_hbee; //hbee visits - all plot
	
	//Imputed missing data;
	vector[Nplot_all] plDens;
	vector[Nplot_all] flDens;	
	//Combine observed with imputed		
	// Plant density
	plDens[obsPlDens_ind]=plDens_obs; //Observed plant density from my fields
	plDens[missPlDens_ind]=plDens_miss[1:Nplot_densMiss]; //Missing data from my fields
	plDens[(Nplot+1):Nplot_all] = plDens_miss_extra; //Riley's fields		
	// Flower density
	flDens[obsFls_ind]=flDens_obs; //Observed flower density
	flDens[missFls_ind]=flDens_miss;
 	for(i in 1:Nplot_flsObs_extra) //For each extra observed plot
		flDens[obsFls_ind_extra[i]+Nplot]=flDens_obs_extra[i];	//Add it to index in flDens
	for(i in 1:Nplot_flsMiss_extra) //For each extra missing plot
		flDens[missFls_ind_extra[i]+Nplot]=flDens_extra_miss[i];	
	
	for(i in 1:Nplot_all){	//Parameters for all fields, all plots
		// Plant density = intercept + random field int + hbee distance effect
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex_all[i]] + 
			slopeHbeeDistPlDens*logHbeeDist_all[i]; //Distance effect				
			
		// Plant size (plot-level) = intercept + random field int + random plot int + distance + planting density effect 
		
		plSizePlotMu[i] = intPlSize + intPlSize_field[plotIndex_all[i]] + //intPlSize_plot[i] + 			
			slopeDistPlSize*logHbeeDist_all[i] + //Distance effect (edge of field has smaller plants)			
			slopePlDensPlSize*plDens[i]; //Planting density effect			
			
		// Flower density = intercept + random field int + plant size effect
		flDensMu[i] = intFlDens	+ intFlDens_field[plotIndex_all[i]] + 
			slopeMBayFlDens*isMBay_all[i] + //Male bay effect
			slopePlSizeFlDens*plSizePlotMu[i] + //Plant size
			slope2016FlDens*is2016_all[i] + //Year effect
			slopeDistFlDens*logHbeeDist_all[i]; //Distance from edge effect
			// slopePlDensFlDens*plDens[i]; //Planting density effect
				
		// Expected value for hbee visits = intercept + random int + distance + bay position + bay type + time offset
		visitMu_hbee[i] = intVisitHbee + logTime_all[i] + intVisitHbee_field[plotIndex_all[i]] + //Intercepts + time offset
			slopeHbeeDistHbee*logHbeeDist_all[i] +  //hbee distance			
			slopeLbeeDistHbee*logLbeeDist_all[i] + //lbee distance
			slopeLbeeHbeeDistHbee*logHbeeDist_all[i]*logLbeeDist_all[i] + //Hbee:lbee distance interaction
			slopeLbeeVisHbee*logLbeeVis_all[i] + //Direct effect of (log) leafcutter visitation						
			slopeCentHbee*isCent_all[i] + //bay center effect
			slopeFlDensHbee*flDens[i] + //Flower density effect
			slopeMBayHbee*isMBay_all[i] + //M bay effect 	
			slopePlDensHbeeVis*plDens[i] + //Conditioning set
			slopePlSizeHbeeVis*plSizePlotMu[i]; //Claim
	}	
	
	for(i in 1:Nplant){	
		// Predicted plant size (taken from plot level measurements above)
		plSizeMu[i] = plSizePlotMu[plantIndex[i]];		
	}	
}
	
model {
	vector[2] bernLL_hbee; //pre-calculate LL for zero inflation process in honey bees
	bernLL_hbee[1]=bernoulli_lpmf(0|zeroVisHbeeTheta); //LL of no extra zero
	bernLL_hbee[2]=bernoulli_lpmf(1|zeroVisHbeeTheta); //LL of extra zero
	
	// Likelihood	
	for(i in 1:Nplot_all){ //Zero-inflated negbin for lbee visitation frequency
		if(hbeeVis_all[i]==0)
			target += log_sum_exp(bernLL_hbee[2],bernLL_hbee[1]+neg_binomial_2_log_lpmf(0|visitMu_hbee[i],visitHbeePhi));
		else
			target += bernLL_hbee[1]+neg_binomial_2_log_lpmf(hbeeVis_all[i]|visitMu_hbee[i],visitHbeePhi);
	}	
	
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density per plot		
	plantSize ~ student_t(exp(nuPlSize),plSizeMu,sigmaPlSize); //Plant size			
	flDens ~ student_t(exp(nuFlDens),flDensMu,sigmaFlDens); 			
				
	// Priors
	//Claim
	slopePlSizeHbeeVis ~ normal(0,5);
	slopePlDensHbeeVis ~ normal(0,5); //Conditioning set
	
	// Planting density
	intPlDens ~ normal(0,0.5); //Intercept
	slopeHbeeDistPlDens ~ normal(0,0.1); //Distance into field		
	sigmaPlDens ~ gamma(2,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(4,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field	
	
	//Plant size - informative priors
	intPlSize ~ normal(0,0.1); //Intercept
	slopePlDensPlSize ~ normal(0,1); //Planting density
	slopeDistPlSize ~ normal(0,0.1); //Distance from edge of field	
	sigmaPlSize ~ gamma(6,10); //Sigma for residual
	sigmaPlSize_field ~ gamma(1.1,10); //Sigma for field	
	intPlSize_field	~ normal(0,sigmaPlSize_field); //Random intercept for field	
	nuPlSize ~ normal(0,3); //nu for student's t
	
	// Flower density	
	intFlDens ~ normal(0,1); //Intercept
	slopePlSizeFlDens ~ normal(0,2); //Plant size effect
	slope2016FlDens ~ normal(0,2); //Effect of 2016
	slopeDistFlDens ~ normal(0,2); //Distance from edge
	slopeMBayFlDens ~ normal(0,5); //Effect of male bay
	// slopePlDensFlDens ~ normal(0,1); //Plant density 
	sigmaFlDens ~ gamma(20,4); //Sigma for plot (residual)
	sigmaFlDens_field ~ gamma(16,4); ; //Sigma for field
	intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field	
	nuFlDens ~ normal(0,1); //nu for student's t
	
	// Hbee Visitation - informative priors
	intVisitHbee ~ normal(2,1); //Intercept	
	slopeHbeeDistHbee ~ normal(0,1); //Slope of distance effect on hbee visits	
	slopeLbeeDistHbee ~ normal(0,1); //Effect of leafcutter shelter distance	
	slopeLbeeHbeeDistHbee ~ normal(0,1); //Hbee-lbee distance interaction	
	slopeLbeeVisHbee ~ normal(0,0.5); //Direct effect of (log) leafcutter visitation
	slopeCentHbee ~ normal(0,1); //Effect of center of bay
	slopeMBayHbee ~ normal(0,1); //Effect of male bay
	slopeFlDensHbee ~ normal(0,0.1); //Flower density effect	
	visitHbeePhi ~ gamma(3.5,5); //Dispersion parameter		
	zeroVisHbeeTheta ~ beta(5,7); // Zero-inflation parameter
	sigmaHbeeVis_field ~ gamma(1,1); //Sigma for field
	intVisitHbee_field ~ normal(0,sigmaHbeeVis_field); //Random intercept for field		
}

generated quantities {
}
