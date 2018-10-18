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
	//Claim: seed count ~ distance from lbee shelter
	real slopeMinDistSeedCount; 	
	
	//Plant size - random effects at plot/field level weren't converging
	vector[Nplant_miss] plantSize_miss; //Vector for imputing missing values	
			
	// Pollen deposition
	real intPol; //Intercept	
	real slopeHbeePol; //Slope of hbee visits 
	real slopeLbeePol; //Slope of lbee visits 	
	real<lower=0> pollenPhi; //Dispersion parameter for pollen deposition
	real<lower=0> sigmaPolField; //Sigma for field-level intercept
	real<lower=0> sigmaPolPlot; //Sigma for plot-level intercept
	vector[Nfield] intPol_field; //Field-level random intercept
	vector[Nplot] intPol_plot; //Plot-level random intercept
					
	// Weight per seed
	real intSeedWeight; //Intercept	
	real slopePolSeedWeight; //Slope of pollen deposition
	real slopeSeedCount; //Slope of seed count
	real slopePlSizeSeedWeight; //Slope of plant size
	real<lower=0> sigmaSeedWeight; //SD of seed weight
	real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect
	real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect
	real<lower=0> sigmaSeedWeight_field; //SD of field random effect	
	vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	vector[Nplant] intSeedWeight_plant; //plant-level random intercepts		
	real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
}

transformed parameters {			
	//Expected values
	// Plot-level	
	vector[Nplot] pollenMu_plot; //Plot level pollen
	vector[Nflw] pollenMu; //Expected pollen - flower level		
	vector[Nplant] seedCountMuPlant; //Plant-level seed count	
	vector[Npod] seedCountMu; //Pod-level seed counts		
	
	//Imputed missing data;
	vector[Nplant] plantSize; //Vector for all values		
	//Combine observed with imputed			
	//Plant size
	plantSize[obsPlant_ind]=plantSize_obs;  //Observed plant size
	plantSize[missPlant_ind]=plantSize_miss; //Imputed plant size							
	
	for(i in 1:Nplot){
		// Pollen per plot = intercept + random field int + random plot int + leafcutter effect + honeybee effect + bay center effect + hbee dist effect
		pollenMu_plot[i] = intPol + intPol_field[plotIndex[i]] + intPol_plot[i] + //Intercept + field/plot level random effects 
			slopeLbeePol*logLbeeVis_all[i] +  //Effect of (log) leafcutter visits			
			slopeHbeePol*logHbeeVis_all[i];  //Effect of (log) honeybee visits 											
	}	
				
	for(i in 1:Nflw) 
		pollenMu[i] = pollenMu_plot[flowerIndex[i]]; //Assigns plot level pollen mu to Nflw long vector
		
	for(i in 1:Nplant){			
		// Weight per seed = intercept + random int field + random int plot + random int plant 
		seedWeightPlantMu[i] = intSeedWeight + intSeedWeight_field[plotIndex[plantIndex[i]]] + intSeedWeight_plot[plantIndex[i]] + intSeedWeight_plant[i] +
			slopePolSeedWeight*pollenMu_plot[plantIndex[i]] + //Pollen deposition			
			slopePlSizeSeedWeight*plantSize[i] + //Plant size
			slopeMinDistSeedWeight*logLbeeDist_all[plantIndex[i]]; //Lbee distance effect			
	}
	
	for(i in 1:Npod){ //For each pod
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]] + slopeSeedCount*seedCount[i]; //Plant-level effect + effect of seedCount 
	}		
}
	
model {	
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed				
	// Priors	
	//Claim
	slopeMinDistSeedWeight ~ normal(0,1); 
	
	// Pollen deposition - informative priors
	intPol ~ normal(3,1); //Intercept		
	slopeHbeePol ~ normal(-0.4,0.5); //hbee Visitation effect
	slopeLbeePol ~ normal(0.2,0.5); //lbee Visitation effect
	sigmaPolField ~ gamma(9,10); //Sigma for random field
	sigmaPolPlot ~ gamma(6,10); //Sigma for random plot	
	pollenPhi ~ gamma(4,5); //Dispersion parameter
	intPol_field ~ normal(0,sigmaPolField); //Random field int
	intPol_plot ~ normal(0,sigmaPolPlot); //Random plot int
	
	// Weight per seed 
	intSeedWeight ~ normal(1.4,0.4); //Intercept	
	slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	slopeSeedCount ~ normal(-0.01,0.02); //Slope of seed count
	slopePlSizeSeedWeight ~ normal(0.06,0.05); //Slope of plant size
	sigmaSeedWeight ~ gamma(4,10); //SD of seed weight
	sigmaSeedWeight_field ~ gamma(1,10); //SD of field random effect	
	sigmaSeedWeight_plot ~ gamma(1,10); //SD of plot random effect
	sigmaSeedWeight_plant ~ gamma(1,10); //SD of plant random effect		
	intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts
	lambdaSeedWeight ~ gamma(3,1); //Lambda term for exponential process
}

generated quantities{
	//Pod-level
	//weight per seed
	real predSeedMass[Npod]; //Generated
	real seedMass_resid[Npod]; //Residual
	real predSeedMass_resid[Npod]; //Residual of generated	
		
	for(i in 1:Npod){ //For each pod
		//weight per seed - exp-normal works well		
		seedMass_resid[i] = seedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight)); 
		predSeedMass[i] = exp_mod_normal_rng(seedWeightMu[i],sigmaSeedWeight,lambdaSeedWeight); 
		predSeedMass_resid[i] = predSeedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));	
	}	
}
