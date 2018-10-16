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
	//Claim
	real XXX; 

	//Plant density
	//Vector for imputing missing plant density values (missing values from my data + all of Riley's data)
	vector[Nplot_densMiss] plDens_miss; //My fields
	vector[Nplot_extra] plDens_miss_extra; //Riley's fields	
	real intPlDens; //Global intercept
	real slopeHbeeDistPlDens; //Slope of distance into field	
	real<lower=0.01> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaPlDens_field; //Sigma for field
	vector[Nfield_all] intPlDens_field; //Random intercept for field
	
	//Plant size - random effects at plot/field level weren't converging
	vector[Nplant_miss] plantSize_miss; //Vector for imputing missing values	
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density	
	// real slopeDistPlSize; //Slope of distance (edge of field has small plants)		
	real slope2016PlSize; //Effect of 2016 on plant size
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	
	// Flower density per plot
	vector<lower=0.01>[Nplot_flsMiss] flDens_miss; //Vectors for imputing missing values
	vector<lower=0.01>[Nplot_flsMiss_extra] flDens_extra_miss;
	real intFlDens; //Global intercept
	real slopePlSizeFlDens; //Slope of plant size on flower density
	real<lower=0.01> sigmaFlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaFlDens_field; //Sigma for field
	vector[Nfield_all] intFlDens_field; //Random intercept for field
 
	// hbee Visitation - random effects at field level weren't converging
	real intVisitHbee; //Intercept	
	real slopeFlDensHbee; //Slope of flower density
	real slopeHbeeDistHbee; //Slope of (log) distance
	real slopeLbeeDistHbee; //Slope of leafcutter distance
	// real slopeLbeeHbeeDistHbee; //Interaction b/w leafcutter & honeybee distance
	// real slopeLbeeVisHbee; //Direct effect of leafcutter visitation	
	real slopeFBayHbee; //Effect of female bay		
	real slopeCentHbee; //Effect of bay position (center)	
	real<lower=0.01> visitHbeePhi; //Dispersion parameter		
	real<lower=0,upper=1> zeroVisHbeeTheta; //Zero-inflation parameter - chance that zero is not from neg.bin. 
	
	// lbee Visitation
	real intVisitLbee; //Intercept
	real slopeFlDensLbee; //Slope of flower density
	real slopeLbeeDistLbee; //Slope of leafcutter distance (shelter)		
	real slopeHbeeDistLbee; //Slope of honeybee distance (field edge)
	real slopeFBayLbee; //Effect of female bay	
	real slopeCentLbee; //Effect of bay position (center)	
	real slopeStockingLbee; //Effect of half-stocking leafcutter bees
	// real slopeCentHbeeDistLbee; //Bay position : honeybee distance interaction term
	// real slopeStockingHbeeDistLbee; //Half-stocking:hbee distance interaction			
	real<lower=0.01> sigmaLbeeVisField; //SD of field random intercepts
	real<lower=0.01> visitLbeePhi; //Dispersion parameter	
	vector[Nfield_all] intVisitLbee_field; //field-level random intercepts	
	
	// // Pollen deposition
	// real intPol; //Intercept	
	// real slopeHbeePol; //Slope of hbee visits 
	// real slopeLbeePol; //Slope of lbee visits 
	// // real slopeCentPol; //Bay center effect
	// // real slopeHbeeDistPol; //(log) hbee distance effect
	// real<lower=0.01> pollenPhi; //Dispersion parameter for pollen deposition
	// real<lower=0.01> sigmaPolField; //Sigma for field-level intercept
	// real<lower=0.01> sigmaPolPlot; //Sigma for plot-level intercept
	// vector[Nfield] intPol_field; //Field-level random intercept
	// vector[Nplot] intPol_plot; //Plot-level random intercept
		
	// // Flower count (per plant) - random effects at plot level weren't converging
	// real intFlwCount; //Intercept
	// real slopePlSizeFlwCount; //Slope of plant size
	// real<lower=0.01> sigmaFlwCount_field; //SD of field-level random effect
	// vector[Nfield] intFlwCount_field; //Field-level random effect
	// real<lower=0.01> flwCountPhi; //Dispersion parameter
	
	// // Flower survival
	// real intFlwSurv; //Intercept	
	// real slopePolSurv; //Slope of pollen deposition
	// real slopePlSizeSurv; //Slope of plant size
	// // real slopeEdgeCentSurv; //Slope of edge effect	
	// // real slopeSeedSizeSurv; //Slope of seed size	
	// real<lower=0.01> sigmaFlwSurv_plot; //SD of plot random intercepts	
	// real<lower=0.01> sigmaFlwSurv_field; //SD of field random intercepts
	// vector[Nfield] intFlwSurv_field; //field-level random intercepts
	// vector[Nplot] intFlwSurv_plot; //plot-level random intercepts	
		
	// // Seed count - random effects at plot level weren't converging
	// real intSeedCount; //Intercept	
	// real slopePolSeedCount; //Slope of pollen deposition
	// real slopePlSizeCount; //Slope of plant size
	// // real slopeEdgeCentSeedCount; // Slope of edge effect on seed count 
	// // real slopeHbeeSeedCount; //Slope of honeybee visitation on seed count
	// // real slopeLbeeSeedCount; //Slope of leafcutter visitation on seed count
	// // real slopeHbeeDistSeedCount; //Slope of leafcutter distance on seed count - correlated 	
	// // real slopeSurvSeedCount; //Slope of plant-level survival on seed count
	// real<lower=0.01> seedCountPhi; //Dispersion parameter
	// real<lower=0.01> sigmaSeedCount_plant; //SD of plant random effect
	// real<lower=0.01> sigmaSeedCount_field; //SD of field random effect
	// vector[Nfield] intSeedCount_field; //field-level random intercepts		
	// vector[Nplant] intSeedCount_plant; //plant-level random intercepts
	
	// // Weight per seed
	// real intSeedWeight; //Intercept	
	// real slopePolSeedWeight; //Slope of pollen deposition
	// real slopeSeedCount; //Slope of seed count
	// real slopePlSizeSeedWeight; //Slope of plant size
	// real<lower=0.01> sigmaSeedWeight; //SD of seed weight
	// real<lower=0.01> sigmaSeedWeight_plant; //SD of plant random effect
	// real<lower=0.01> sigmaSeedWeight_plot; //SD of plot random effect
	// real<lower=0.01> sigmaSeedWeight_field; //SD of field random effect	
	// vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	// vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	// vector[Nplant] intSeedWeight_plant; //plant-level random intercepts		
	// real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
}

transformed parameters {			
	//Expected values
	//Plot-level
	vector[Nplot_all] plDensMu; //Expected plant density	
	vector[Nplot_all] plSizePlotMu; //Plot-level plant size
	vector[Nplant] plSizeMu; //Expected plant size	
	vector[Nplot_all] flDensMu; //Expected flower density	
	// vector[Nplot_all] visitMu_hbee; //hbee visits - all plot
	// vector[Nplot_all] visitMu_lbee; //lbee visits - all plots	
	// vector[Nplot] flwCountPlot; //Plot-level flower production (per plant)
	// vector[Nplant] flwCountMu; //Expected flower count for plant	
	// vector[Nplot] pollenMu_plot; //Plot level pollen
	// vector[Nflw] pollenMu; //Expected pollen - flower level
	// vector[Nplot] flwSurvPlot; //Plot-level flower survival
	// vector[Nplant] flwSurv; //Flower survival rate (logit)
	// vector[Nplant] propFlwSurv; //(logit) observed flower survival	
	// vector[Nplant] seedCountMuPlant; //Plant-level seed count	
	// vector[Npod] seedCountMu; //Pod-level seed counts	
	// vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	// vector[Npod] seedWeightMu; //Pod-level weight per seed		
	
	//Imputed missing data;
	vector[Nplant] plantSize; //Vector for all values
	vector[Nplot_all] plDens;
	vector[Nplot_all] flDens;	
	//Combine observed with imputed		
	//Plant density
	plDens[obsPlDens_ind]=plDens_obs; //Observed plant density from my fields
	plDens[missPlDens_ind]=plDens_miss[1:Nplot_densMiss]; //Missing data from my fields
	plDens[(Nplot+1):Nplot_all] = plDens_miss_extra; //Riley's fields	
	//Plant size
	plantSize[obsPlant_ind]=plantSize_obs;  //Observed plant size
	plantSize[missPlant_ind]=plantSize_miss; //Imputed plant size			
	// Flower density
	flDens[obsFls_ind]=flDens_obs; //Observed flower density
	flDens[missFls_ind]=flDens_miss;
 	for(i in 1:Nplot_flsObs_extra) //For each extra observed plot
		flDens[obsFls_ind_extra[i]+Nplot]=flDens_obs_extra[i];	//Add it to index in flDens
	for(i in 1:Nplot_flsMiss_extra) //For each extra missing plot
		flDens[missFls_ind_extra[i]+Nplot]=flDens_extra_miss[i];
	// //(logit) proportion observed flower survival
	// propFlwSurv[obsPlant_ind]=propFlwSurv_obs;
	// propFlwSurv[missPlant_ind]=flwSurv[missPlant_ind]; //Direct imputation from predicted values
	
	for(i in 1:Nplot_all){		
		//Plant density = intercept + random field int + hbee distance effect
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex_all[i]] + 
			slopeHbeeDistPlDens*logHbeeDist_all[i]; //Distance effect				
			
		//Plant size (plot-level) = intercept + random field int + random plot int + distance + planting density effect 
		//Density distance interaction is basically 0, so leaving it out
		plSizePlotMu[i] = intPlSize + //intPlSize_field[plotIndex_all[i]] + //intPlSize_plot[i] + 			
			// slopeDistPlSize*logHbeeDist_all[i] + //Distance effect (edge of field has smaller plants)			
			slopePlDensPlSize*plDens[i] + //Planting density effect			
			slope2016PlSize*is2016_all[i];
			// slopePolPlSize*pollenMu_plot[i] + //Effect of pollen on plant size (plot level) - creates cyclical association			
			
		//Flower density = intercept + random field int + plant size effect
		flDensMu[i] = intFlDens	+ intFlDens_field[plotIndex_all[i]] + 
			slopePlSizeFlDens*plSizePlotMu[i]; 
	
		// // Expected value for lbee visits = intercept + random int + distance to shelter + distance to honeybees (edge) + bay position + bay:hbee dist + stocking:hbee dist + bay type + time offset	
		// visitMu_lbee[i] = intVisitLbee + intVisitLbee_field[plotIndex_all[i]] + logTime_all[i] + //intercepts + time offset
			// slopeLbeeDistLbee*logLbeeDist_all[i] + //lbee distance
			// slopeHbeeDistLbee*logHbeeDist_all[i] + //hbee distance
			// slopeCentLbee*isCent_all[i] + //bay center
			// slopeStockingLbee*lbeeStocking_all[i] +	//half-stocking
			// slopeFBayLbee*isFBay_all[i] + //F bay
			// // slopeCentHbeeDistLbee*isCent_all[i]*logHbeeDist_all[i] + //hbee dist: bay center interaction			
			// // slopeStockingHbeeDistLbee*lbeeStocking_all[i]*logHbeeDist_all[i] +  //hbee dist: half stocking interaction
			// slopeFlDensLbee*flDens[i]; //Flower density effect
			// // slopePlsizeLbee*plSizePlot[i]; //Plant size
			
		// // Expected value for hbee visits = intercept + random int + distance + bay position + bay type + time offset
		// visitMu_hbee[i] = intVisitHbee + logTime_all[i] + //intVisitHbee_field[plotIndex_all[i]] + //Intercepts + time offset
			// // slopeHbeeDistHbee*logHbeeDist_all[i] +  //hbee distance			
			// slopeLbeeDistHbee*logLbeeDist_all[i] + //lbee distance
			// // slopeLbeeHbeeDistHbee*logHbeeDist_all[i]*logLbeeDist_all[i] + //Hbee:lbee distance interaction
			// // slopeLbeeVisHbee*logLbeeVis_all[i] + //Direct effect of (log) leafcutter visitation						
			// slopeCentHbee*isCent_all[i] + //bay center effect
			// slopeFlDensHbee*flDens[i] + //Flower density effect
			// slopeFBayHbee*isFBay_all[i]; //F bay effect 	
	}	
	
	// for(i in 1:Nplot){
		// // Pollen per plot = intercept + random field int + random plot int + leafcutter effect + honeybee effect + bay center effect + hbee dist effect
		// pollenMu_plot[i] = intPol + intPol_field[plotIndex[i]] + intPol_plot[i] + //Intercept + field/plot level random effects 
			// slopeLbeePol*logLbeeVis_all[i] +  //Effect of (log) leafcutter visits			
			// slopeHbeePol*logHbeeVis_all[i] +  //Effect of (log) honeybee visits 
			// // slopeCentPol*isCent_all[i] + //Bay center effect
			// // slopeHbeeDistPol*logHbeeDist_all[i]; //(log) hbee distance effect				
			
		// //Flower survival = intercept + random int field + random int plot +
		// flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + intFlwSurv_plot[i] + 			
			// // slopeEdgeCentSurv*isCent_all[i] + //Bay center
			// slopePolSurv*logPolCountPlot[i]; //Pollen deposition (plot level average)
			
		// // Flower count per plant (plot level) = intercept + random field int
		// flwCountPlot[i] = intFlwCount + intFlwCount_field[plotIndex[i]]; //+ intFlwCount_plot[i];
	// }	
				
	// for(i in 1:Nflw) 
		// pollenMu[i] = pollenMu_plot[flowerIndex[i]]; //Assigns plot level pollen mu to Nflw long vector
		
	for(i in 1:Nplant){	
		//Predicted plant size (taken from plot level measurements above)
		plSizeMu[i] = plSizePlotMu[plantIndex[i]];		
		
		// //Predicted flower count per plant
		// flwCountMu[i] = flwCountPlot[plantIndex[obsPlant_ind[i]]] + //Plot level flower count 
			// slopePlSizeFlwCount*plantSize[obsPlant_ind[i]]; //individual size effect
		
		// flwSurv[i] = flwSurvPlot[plantIndex[obsPlant_ind[i]]] + //Plot-level plant survival
			// slopePlSizeSurv*plantSize[obsPlant_ind[i]] +  //plant size effect 	
			// // slopeSeedCountSurv*seedCountMuPlant[obsPlant_ind[i]] + //Slope of (log) seed count - NEW VARIABLE
			// // slopeSeedSizeSurv*seedWeightPlantMu[obsPlant_ind[i]]; //Slope of seed size - NEW VARIABLE	
	
		// // Seed count per pod = intercept + random int field + random int plot + random int plant + hbee visits + pollen deposition
		// seedCountMuPlant[i] = intSeedCount + intSeedCount_field[plotIndex[plantIndex[i]]] + intSeedCount_plant[i] + //intSeedCount_plot[plantIndex[i]] +
			// slopePlSizeCount*plantSize[i] + //Plant size
			// // slopeSurvSeedCount*flwSurv[i] + //Flower survival 
			// slopePolSeedCount*logPolCountPlot[plantIndex[i]] + //Pollen deposition (plot level average)
			// // slopeEdgeCentSeedCount*isCent_all[plantIndex[i]] + //Bay center effect
			// // slopeHbeeSeedCount*logHbeeVis_all[plantIndex[i]] + //Honeybee visitation
			// // slopeLbeeSeedCount*logLbeeVis_all[plantIndex[i]] + //Leafcutter visitation
			// // slopeHbeeDistSeedCount*logLbeeDist_all[plantIndex[i]]; //(log) Leafcutter distance 			
			
		// // Weight per seed = intercept + random int field + random int plot + random int plant 
		// seedWeightPlantMu[i] = intSeedWeight + intSeedWeight_field[plotIndex[plantIndex[i]]] + intSeedWeight_plot[plantIndex[i]] + intSeedWeight_plant[i] +
			// slopePolSeedWeight*logPolCountPlot[plantIndex[i]] + //Pollen deposition			
			// slopePlSizeSeedWeight*plantSize[i]; //Plant size
	}
	
	// for(i in 1:Npod){ //For each pod
		// seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
		// seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+slopeSeedCount*seedCount[i]; //Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
	// }		
}
	
model {
	// vector[2] bernLL_hbee; //pre-calculate LL for zero inflation process
	// bernLL_hbee[1]=bernoulli_lpmf(0|zeroVisHbeeTheta); //LL of no extra zero
	// bernLL_hbee[2]=bernoulli_lpmf(1|zeroVisHbeeTheta); //LL of extra zero			
	// //Likelihood	
	// for(i in 1:Nplot_all){ //Zero-inflated negbin for hbee visitation frequency
		// if(hbeeVis_all[i]==0)
			// target += log_sum_exp(bernLL_hbee[2],bernLL_hbee[1]+neg_binomial_2_log_lpmf(0|visitMu_hbee[i],visitHbeePhi));
		// else
			// target += bernLL_hbee[1]+neg_binomial_2_log_lpmf(hbeeVis_all[i]|visitMu_hbee[i],visitHbeePhi);					
	// }		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density per plot
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size		
	flDens ~ normal(flDensMu,sigmaFlDens); //Flower density per plot
	// lbeeVis_all ~ neg_binomial_2_log(visitMu_lbee,visitLbeePhi); //Lbee visitation rate
	// pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate
	// flwCount ~ neg_binomial_2_log(flwCountMu[obsPlant_ind],flwCountPhi); //Flower count per plant (attempted pods)
	// podCount ~ binomial_logit(flwCount,flwSurv[obsPlant_ind]); //Flower survival
	// seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	// seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed	
			
	// Priors	
	//Claim
	XXX ~ normal(0,1); 
	
	//Planting density
	intPlDens ~ normal(0,0.5); //Intercept
	slopeHbeeDistPlDens ~ normal(0.05,0.1); //Distance into field	
	// slopeHbeeDistSqPlDens ~ normal(-0.01,0.1); //Distance into field squared	- doesn't appear to add anything
	sigmaPlDens ~ gamma(2,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(4,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field	
	
	//Plant size - informative priors
	intPlSize ~ normal(0,0.2); //Intercept
	slopePlDensPlSize ~ normal(-0.75,0.5); //Planting density
	// slopeDistPlSize ~ normal(0.07,.1); //Distance effect
	slope2016PlSize	~ normal(0,1); //Year effect
	sigmaPlSize ~ gamma(6,10); //Sigma for residual			
	
	// Flower density	
	intFlDens ~ normal(0,1); //Intercept
	slopePlSizeFlDens ~ normal(0,1); //Plant size effect
	sigmaFlDens ~ gamma(20,4); //Sigma for plot (residual)
	sigmaFlDens_field ~ gamma(16,4); ; //Sigma for field
	intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field	
	
	// // Hbee Visitation - informative priors
	// intVisitHbee ~ normal(2,1); //Intercept	
	// // slopeHbeeDistHbee ~ normal(-0.3,0.5); //Slope of distance effect on hbee visits	
	// slopeLbeeDistHbee ~ normal(0.3,0.5); //Effect of leafcutter shelter distance	
	// // slopeLbeeHbeeDistHbee ~ normal(0.1,0.5); //Hbee-lbee distance interaction	
	// // slopeLbeeVisHbee ~ normal(-0.05,0.15); //Direct effect of (log) leafcutter visitation
	// slopeCentHbee ~ normal(-0.4,0.5); //Effect of center of bay
	// slopeFBayHbee ~ normal(0.1,0.5); //Effect of female bay
	// slopeFlDensHbee ~ normal(0.01,0.05); //Flower density effect	
	// visitHbeePhi ~ gamma(3.5,5); //Dispersion parameter		
	// zeroVisHbeeTheta ~ beta(3,7); //Zero-inflation parameter	
	
	// // Lbee Visitation - informative priors
	// intVisitLbee ~ normal(4,1); //Intercept	
	// slopeHbeeDistLbee ~ normal(-0.2,0.2); //Slope of honeybee distance on lbee visits
	// slopeLbeeDistLbee ~ normal(-0.8,0.2); //Slope of shelter distance on lbee visits
	// slopeCentLbee ~ normal(-0.6,0.5); //Effect of center of bay
	// slopeFBayLbee ~ normal(0,0.5); //Effect of female bay
	// slopeStockingLbee ~ normal(0,0.5); //Effect of half-stocking
	// // slopeCentHbeeDistLbee ~ normal(-0.2,0.2); //Bay center: hbee distance interaction
	// // slopeStockingHbeeDistLbee ~ normal(0.25,0.2); //Half-stocking: hbee distance interaction			
	// slopeFlDensLbee ~ normal(0.03,0.05); //Flower density effect
	// sigmaLbeeVisField ~ gamma(2,2); //Sigma for random field 
	// visitLbeePhi ~ gamma(4,10); //Dispersion parameter	
	// intVisitLbee_field ~ normal(0,sigmaLbeeVisField); //Random field intercepts	
	
	// // Pollen deposition - informative priors
	// intPol ~ normal(3,1); //Intercept		
	// slopeHbeePol ~ normal(-0.4,0.5); //hbee Visitation effect
	// slopeLbeePol ~ normal(0.2,0.5); //lbee Visitation effect
	// // slopeCentPol~ normal(-0.4,0.5); //Bay center effect
	// // slopeHbeeDistPol ~ normal(-0.2,0.5); //(log) hbee distance effect	
	// sigmaPolField ~ gamma(9,10); //Sigma for random field
	// sigmaPolPlot ~ gamma(6,10); //Sigma for random plot	
	// pollenPhi ~ gamma(4,5); //Dispersion parameter
	// intPol_field ~ normal(0,sigmaPolField); //Random field int
	// intPol_plot ~ normal(0,sigmaPolPlot); //Random plot int
	
	// //Flower count (per plant)
	// intFlwCount ~ normal(6,0.5); //Intercept
	// slopePlSizeFlwCount ~ normal(0.5,0.1); //Slope of plant size
	// sigmaFlwCount_field ~ gamma(1,5); //SD of field-level random effect	
	// intFlwCount_field ~ normal(0,sigmaFlwCount_field); //Field-level random effect	
	// flwCountPhi ~ gamma(20,4); //Dispersion parameter
		
	// // Flower survival 
	// intFlwSurv ~ normal(-5,1); //Intercept	
	// slopePolSurv ~ normal(2,1); //Slope of pollen deposition
	// slopePlSizeSurv ~ normal(0,0.5); //Slope of plant size
	// // slopeEdgeCentSurv ~ normal(0.4,0.4); //Slope of edge effect
	// // slopeSeedCountSurv ~ normal(2.5,0.5); //Slope of seed count  - I think this is in the wrong direction
	// // slopeSeedSizeSurv ~ normal(0,0.5); //Slope of seed size
	// sigmaFlwSurv_field ~ gamma(1.5,5); //SD of field random effect
	// sigmaFlwSurv_plot ~ gamma(2.5,5); //SD of plot random effect	
	// intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	// intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts			
	
	// // Seed count 
	// intSeedCount ~ normal(2.5,0.5); //Intercept	
	// slopePolSeedCount ~ normal(2,1); //Slope of pollen deposition
	// slopePlSizeCount ~ normal(0.1,0.05); //Slope of plant size
	// // slopeEdgeCentSeedCount ~ normal(-0.2,0.5); // Slope of edge effect on seed count 
	// // slopeHbeeSeedCount ~ normal(0.1,0.5); //Slope of honeybee visitation on seed count
	// // slopeLbeeSeedCount ~ normal(0.02,0.2); //Slope of leafcutter visitation on seed count
	// // slopeHbeeDistSeedCount ~ normal(0,0.5); //Slope of leafcutter distance on seed count	- correlated with other predictors, and not very strong
	// // slopeSurvSeedCount ~ normal(0,1); //Slope of survival
	// seedCountPhi ~ gamma(12,4); //Dispersion parameter
	// sigmaSeedCount_field ~ gamma(1,10); //SD of field random effect
	// // sigmaSeedCount_plot ~ gamma(1,10); //SD of plot random effect
	// sigmaSeedCount_plant ~ gamma(1,10); //SD of plant random effect	
	// intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	// // intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	// intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts
	
	// // Weight per seed 
	// intSeedWeight ~ normal(1.4,0.4); //Intercept	
	// slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	// slopeSeedCount ~ normal(-0.01,0.02); //Slope of seed count
	// slopePlSizeSeedWeight ~ normal(0.06,0.05); //Slope of plant size
	// sigmaSeedWeight ~ gamma(4,10); //SD of seed weight
	// sigmaSeedWeight_field ~ gamma(1,10); //SD of field random effect	
	// sigmaSeedWeight_plot ~ gamma(1,10); //SD of plot random effect
	// sigmaSeedWeight_plant ~ gamma(1,10); //SD of plant random effect		
	// intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	// intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	// intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts
	// lambdaSeedWeight ~ gamma(3,1); //Lambda term for exponential process
}

generated quantities{
	//Plot-level quantities
	//plant density
	real predPlDens[Nplot]; //Generated
	real plDens_resid[Nplot]; //Residual
	real predPlDens_resid[Nplot]; //Residual of generated
	//flower density
	real predFlDens[Nplot_all]; //Generated
	real flDens_resid[Nplot_all]; //Residual
	real predFlDens_resid[Nplot_all]; //Residual of generated	
	// // hbeeVis
	// int predHbeeVis_all[Nplot_all]; //Generated
	// real hbeeVis_resid[Nplot_all]; //Residual 
	// real predHbeeVis_resid[Nplot_all]; //Residual of generated 
	// // lbee visits
	// int predLbeeVis_all[Nplot_all]; //Generated 
	// real lbeeVis_resid[Nplot_all]; //Residual
	// real predLbeeVis_resid[Nplot_all]; //Residual of generated		
	
	// // Flower-level
	// // pollen deposition
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
		
	for(i in 1:Nplot_all){
		// //bee visits
		// hbeeVis_resid[i]=hbeeVis_all[i]-(exp(visitMu_hbee[i])*(1-zeroVisHbeeTheta)); //Residual for actual value x offset
		// if(bernoulli_rng(zeroVisHbeeTheta)==1) //If zeroVisHbeeTheta generates an extra zero
			// predHbeeVis_all[i] = 0; //Predicted value is automatically zero
		// else //Otherwise
			// predHbeeVis_all[i] = neg_binomial_2_log_rng(visitMu_hbee[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		// predHbeeVis_resid[i]=predHbeeVis_all[i]-(exp(visitMu_hbee[i])*(1-zeroVisHbeeTheta)); //Residual for predicted value
		
		// //lbee visits
		// lbeeVis_resid[i]= lbeeVis_all[i]-exp(visitMu_lbee[i]); //Residual for actual value
		// predLbeeVis_all[i] = neg_binomial_2_log_rng(visitMu_lbee[i],visitLbeePhi); //Predicted value drawn from neg.bin		
		// predLbeeVis_resid[i]=predLbeeVis_all[i]-exp(visitMu_lbee[i]); //Residual for predicted value		
		
		//flower density
		flDens_resid[i] = flDens[i]-flDensMu[i]; //Residual for actual value
		predFlDens[i] = normal_rng(flDensMu[i],sigmaFlDens); //Generated value from normal
		predFlDens_resid[i] = predFlDens[i] - flDensMu[i]; //Residual for predicted value	
	}
	
	for(i in 1:Nplot){
		//plant density		
		plDens_resid[i] = plDens[i] - plDensMu[i]; //Residual for actual value
		predPlDens[i]= normal_rng(plDensMu[i],sigmaPlDens); //Generated value from normal
		predPlDens_resid[i] = predPlDens[i] - plDensMu[i]; //Residual for predicted value							
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
		// //weight per seed - exp-normal works well		
		// seedMass_resid[i] = seedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight)); 
		// predSeedMass[i] = exp_mod_normal_rng(seedWeightMu[i],sigmaSeedWeight,lambdaSeedWeight); 
		// predSeedMass_resid[i] = predSeedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));
	// }	
}
