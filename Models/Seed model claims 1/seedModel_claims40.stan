data {
	//Field level
	int Nfield; //Number of fields		
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 		
	int lbeeStocking[Nplot]; //Half leafcutter stocking? (double-tent treatment)
	int is2016[Nplot]; //Was field from 2016?	
	vector[Nplot] hbee_dist; //Distance from hbee hives
	int hbeeVis[Nplot]; //Number of hbee visits to plot
	vector[Nplot] lbee_dist; //Distance to lbee shelters
	int lbeeVis[Nplot]; //Number of lbee visits to plot
	int isCent[Nplot]; //Is plot in center of bay?
	int isFBay[Nplot]; //Is plot in female bay?
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
	
	//Extra plots/fields from Riley (visitation data)
	int Nfield_extra;		
	int Nplot_extra; 
	int is2016_extra[Nplot_extra]; //Was field from 2016?
	int lbeeStocking_extra[Nplot_extra]; //Half leafcutter stocking? (double-tent treatment) - FIX
	int plotIndex_extra[Nplot_extra]; 	
	vector[Nplot_extra] hbee_dist_extra; 
	int hbeeVis_extra[Nplot_extra]; 
	vector[Nplot_extra] lbee_dist_extra;
	int lbeeVis_extra[Nplot_extra];
	int isCent_extra[Nplot_extra];
	int isFBay_extra[Nplot_extra];
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
	int Nplant; //Number of all plants (some measurements missing)
	int Nplant_obs; //Number of "complete" (observed) plants
	int Nplant_miss; //Number of missing plants	
	int podCount[Nplant_obs]; //Number of pods per plant
	int flwCount[Nplant_obs]; //Number of total flower (pods + missing) per plant
	vector[Nplant_obs] plantSize_obs; //Mass of vegetative tissue (no seeds) (g)
	vector[Nplant_obs] totalSeedMass; //Mass of seeds from plant (g)
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?	
	int<lower=1,upper=Nplot> plantSurvIndex[Nplant_obs]; //Index for "complete" plants 
	int<lower=1,upper=Nplant> obsPlant_ind[Nplant_obs]; //Index for observed plants
	int<lower=1,upper=Nplant> missPlant_ind[Nplant_miss]; //Index for missing plants
			
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
	int is2016_all[Nplot_all]; //Is field from 2016?
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
	vector[Nplot_all] logHbeeVis_all; //Log-visitation rate for hbees
	vector[Nplot_all] logLbeeVis_all; //Log-visitation rate for lbees
	vector[Nplot] logPolCountPlot; //Log-average pollen count at plot level
	vector[Nplant_obs] propFlwSurv_obs; //(logit) proportion flower survival	
	
	//Assign values
	plotIndex_all[1:Nplot] = plotIndex; 
	plotIndex_all[Nplot+1:Nplot_all] = plotIndex_extra;
	lbeeStocking_all[1:Nplot] = lbeeStocking; //Lbee stocking
	lbeeStocking_all[Nplot+1:Nplot_all] = lbeeStocking_extra;
	is2016_all[1:Nplot] = is2016; //Year
	is2016_all[Nplot+1:Nplot_all] = is2016_extra;	
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
	//Transformations
	logHbeeDist_all=log(hbee_dist_all); //Log-transform distances
	logLbeeDist_all=log(lbee_dist_all);
	logTime_all=log(totalTime_all); //Log-transform time	
	for(i in 1:Nplot_all){
		logHbeeVis_all[i] = log(1+(hbeeVis_all[i]/totalTime_all[i])); //Log-transforms observed visitation rates 
		logLbeeVis_all[i] = log(1+(lbeeVis_all[i]/totalTime_all[i]));			
	}
	logPolCountPlot=polCountPlot; //Log-transforms average pollen counts
	
	for(i in 1:Nplant_obs){
		//Necessary for promoting integers to reals. Otherwise does integer division.
		propFlwSurv_obs[i] = podCount[i]; 
		propFlwSurv_obs[i] = propFlwSurv_obs[i]/flwCount[i]; //Proportion surviving pods		
		if(propFlwSurv_obs[i]<=0) //Deal with weird 100% and 0% plants
			propFlwSurv_obs[i]=0.01;
		else if(propFlwSurv_obs[i]>=1)
			propFlwSurv_obs[i]=0.99;		
	}
	//Logit transform 
	propFlwSurv_obs=logit(propFlwSurv_obs);			
}

parameters { 
	//Claim: FlwCount ~ hbee visitation
	real slopeHbeeVisFlwCount; 

	//Plant density
	//Vector for imputing missing plant density values (missing values from my data + all of Riley's data)
	vector[Nplot_densMiss] plDens_miss; //My fields
	vector[Nplot_extra] plDens_miss_extra; //Riley's fields	
	real intPlDens; //Global intercept
	real slopeHbeeDistPlDens; //Slope of distance into field	
	real<lower=0> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0> sigmaPlDens_field; //Sigma for field
	vector[Nfield_all] intPlDens_field; //Random intercept for field
	
	//Plant size - random effects at plot/field level weren't converging
	vector[Nplant_miss] plantSize_miss; //Vector for imputing missing values	
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density	
	// real slopeDistPlSize; //Slope of distance (edge of field has small plants)		
	real slope2016PlSize; //Effect of 2016 on plant size
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)	
	
	// Flower count (per plant) - random effects at plot level weren't converging
	real intFlwCount; //Intercept
	real slopePlSizeFlwCount; //Slope of plant size
	real<lower=0> sigmaFlwCount_field; //SD of field-level random effect
	vector[Nfield] intFlwCount_field; //Field-level random effect
	real<lower=0> flwCountPhi; //Dispersion parameter	
}

transformed parameters {			
	//Expected values
	//Plot-level
	vector[Nplot_all] plDensMu; //Expected plant density	
	vector[Nplot_all] plSizePlotMu; //Plot-level plant size
	vector[Nplant] plSizeMu; //Expected plant size			
	vector[Nplot] flwCountPlot; //Plot-level flower production (per plant)
	vector[Nplant] flwCountMu; //Expected flower count for plant		
	
	//Imputed missing data;
	vector[Nplant] plantSize; //Vector for all values
	vector[Nplot_all] plDens;	
	//Combine observed with imputed		
	//Plant density
	plDens[obsPlDens_ind]=plDens_obs; //Observed plant density from my fields
	plDens[missPlDens_ind]=plDens_miss[1:Nplot_densMiss]; //Missing data from my fields
	plDens[(Nplot+1):Nplot_all] = plDens_miss_extra; //Riley's fields	
	//Plant size
	plantSize[obsPlant_ind]=plantSize_obs;  //Observed plant size
	plantSize[missPlant_ind]=plantSize_miss; //Imputed plant size				
	
	for(i in 1:Nplot_all){		
		//Plant density = intercept + random field int + hbee distance effect
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex_all[i]] + 
			slopeHbeeDistPlDens*logHbeeDist_all[i]; //Distance effect				
			
		//Plant size (plot-level) = intercept + random field int + random plot int + distance + planting density effect 
		//Density distance interaction is basically 0, so leaving it out
		plSizePlotMu[i] = intPlSize + //intPlSize_field[plotIndex_all[i]] + //intPlSize_plot[i] + 						
			slopePlDensPlSize*plDens[i] + //Planting density effect			
			slope2016PlSize*is2016_all[i];					
	}	
	
	for(i in 1:Nplot){			
		// Flower count per plant (plot level) = intercept + random field int
		flwCountPlot[i] = intFlwCount + intFlwCount_field[plotIndex[i]] + //+ intFlwCount_plot[i];
			slopeHbeeVisFlwCount*hbeeVis_all[i]; //Hbee visitation
	}	
		
	for(i in 1:Nplant){	
		//Predicted plant size (taken from plot level measurements above)
		plSizeMu[i] = plSizePlotMu[plantIndex[i]];		
		
		//Predicted flower count per plant
		flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count 
			slopePlSizeFlwCount*plantSize[i]; //individual size effect		
	}	
}
	
model {	
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density per plot
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size				
	flwCount ~ neg_binomial_2_log(flwCountMu[obsPlant_ind],flwCountPhi); //Flower count per plant (attempted pods)	
			
	// Priors	
	//Claim
	slopeHbeeVisFlwCount ~ normal(0,1); 
	
	//Planting density
	intPlDens ~ normal(0,0.5); //Intercept
	slopeHbeeDistPlDens ~ normal(0.05,0.1); //Distance into field		
	sigmaPlDens ~ gamma(2,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(4,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field	
	
	//Plant size - informative priors
	intPlSize ~ normal(0,0.2); //Intercept
	slopePlDensPlSize ~ normal(-0.75,0.5); //Planting density
	// slopeDistPlSize ~ normal(0.07,.1); //Distance effect
	slope2016PlSize	~ normal(0,1); //Year effect
	sigmaPlSize ~ gamma(6,10); //Sigma for residual			
			
	//Flower count (per plant)
	intFlwCount ~ normal(6,0.5); //Intercept
	slopePlSizeFlwCount ~ normal(0.5,0.1); //Slope of plant size
	sigmaFlwCount_field ~ gamma(1,5); //SD of field-level random effect	
	intFlwCount_field ~ normal(0,sigmaFlwCount_field); //Field-level random effect	
	flwCountPhi ~ gamma(20,4); //Dispersion parameter	
}

generated quantities{	
	//Plant-level	
	//flower count per plant (potential pods)
	int predFlwCount[Nplant_obs]; //Generated
	real flwCount_resid[Nplant_obs]; //Residual
	real predFlwCount_resid[Nplant_obs]; //Residual of generated		

	for(i in 1:Nplant_obs){
		//flower count per plant
		flwCount_resid[i] = flwCount[i] - exp(flwCountMu[i]); //Residual for actual
		predFlwCount[i] = neg_binomial_2_log_rng(flwCountMu[i],flwCountPhi); //Generates new value from neg. bin.
		predFlwCount_resid[i] = predFlwCount[i] - exp(flwCountMu[i]); //Residual for new value		
	}
}
