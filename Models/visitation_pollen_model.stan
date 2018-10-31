data {
	//Field level
	int Nfield; //Number of fields	
	int numHives[Nfield]; //Number of hives present
	int is2015[Nfield]; //Was field measured in 2015?
	int isGP[Nfield]; //Was field from Grand Prairie?
	int isIrrigated[Nfield]; //Was field irrigated?
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
	vector[Nplot] dist; //Distance from edge
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	int flyVis[Nplot]; //Number of fly visits per plot
	vector[Nplot] totalTime; //Minutes taken for observation/10	
	vector[Nplot] flDens; //Flower density
	//Plant density (stems/m2) - missing and observed
	int Nplot_densObs; //Number of plots where plant density was observed
	int Nplot_densMiss; //Number of plots with plant density missing 
	vector[Nplot_densObs] plDens_obs; //Observed plant density (stems/m2)
	int<lower=1,upper=Nplot> obsPlDens_ind[Nplot_densObs]; //Index for observed plant density
	int<lower=1,upper=Nplot> missPlDens_ind[Nplot_densMiss]; //Index for missing plant density		
	
	//Flower level
	int Nflw; //Number of flowers
	int<lower=1,upper=Nplot>  flowerIndex[Nflw];  //Index for flowers - which plot?
	int pollenCount[Nflw]; //Pollen count

	//Plant level
	int Nplant; //Number of all plants (some measurements missing)	
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?		
	vector[Nplant] plantSize; //Centered mass of vegetative tissue (no seeds) (g)
	vector[Nplant] yield; //Observed yield (g seed per plant)	
	
	//Pod level
	int Npod; //Number of pods
	int seedCount[Npod];  //Seeds per pod
	vector[Npod] seedMass; //Mass of seedCount seeds (g)
	int<lower=1,upper=Nplant> podIndex[Npod]; //Index for pods - which plant?	
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives
	vector[Nfield] scalNumHives; //Number of hives/10
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	vector[Nplot] logFlyVis; //Fly visitation rate
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)	
	vector[Nplant] logFlwCount; //Log flower count
	
	logHbeeDist=logHbeeDist-mean(logHbeeDist); //Centers distance
	
	for(i in 1:Nfield){
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives
		scalNumHives[i] = numHives[i]; //Scales number of hives by 10
		scalNumHives[i] = scalNumHives[i]/10;
	}
	
	for(i in 1:Nplot){ //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5); 
		logFlyVis[i] = log((flyVis[i]/totalTime[i])+0.5);
	}
		
	for(i in 1:Nplant){
		logFlwCount[i] = log(flwCount[i]); //Log-transforms flower count per plant
	}
}

parameters {
	// Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slope2015PlDens; //Effect of 2015
	real slopeIrrigPlDens; //Effect of irrigation
	real slope2015IrrigPlDens; //Year:irrigation interaction	
	real slopeGPPlDens; //GP effect on plant density	
	real slopeDistPlDens; //Slope of distance into field		
	real<lower=0.01> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field
	
	// Plant size - Density:distance, Year:GP, Irrigation:2015 interactions are 0. Leaving them out.	
	// Plot level random effect has bad trace, strongly correlated with lp__	
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density	
	real slopeDistPlSize; //Slope of distance 
	real slopeGpPlSize; //Slope of Grand Prairie effect	
	real slope2015PlSize; //Slope of 2015 effect	
	real slopeStockingPlSize; //Stocking effect - number of hives
	real slopeIrrigPlSize; //Slope of irrigation
	// Interactions	
	real slopeDistStockingPlSize; //Distance:Stocking interaction
	real slopeDist2015PlSize; //Distance:year interaction
	// Variance terms
	real<lower=0> sigmaPlSize_field; //Sigma for field
	real<lower=0> sigmaPlSize_plot; //Sigma for plot - small n_eff, higher
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	vector[Nfield] intPlSize_field; //Random intercept for field
	vector[Nplot] intPlSize_plot; //Random intercept for plot - not converging well	
				
	// // Flower density per plot
	// real intFlDens; //Global intercept
	// real slopePlSizeFlDens; //Slope of plant size on flower density
	// real slopeHbeeDistFlDens; //Slope of distance into field
	// real<lower=0.01> sigmaFlDens; //Sigma for within-field (residual)	
	// real<lower=0.01> sigmaFlDens_field; //Sigma for field
	// vector[Nfield] intFlDens_field; //Random intercept for field		

	// // hbee Visitation per plot	
	// real intVisit; //Intercept 
	// real slopeYearVis; //Effect of 2015 
	// real slopeGpVis; //Effect of Grand Prairie 
	// real slopeYearGpVis; //Effect of year-GP interaction
	// real slopeDistVis; //Slope of distance
	// real slopeHiveVis; //Slope of hive number
	// real slopeFlDens; //Slope of flower density
	// real slopeIrrigVis; //Effect of irrigation
	// real<lower=0> sigmaVisField; //SD of field random intercepts
	// real<lower=0> visitHbeePhi; //Dispersion parameter	
	// vector[Nfield] intVisit_field; //field-level random intercepts		
	// real<lower=0> lambdaVisField; //Lambda for skewed random effects	
	
	// Pollen deposition
	// Plot level random effect has bad trace, strongly correlated with lp__
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits	
	real slopeHbeeDistPollen; //Effect of distance into field - p=0.076 (p=0.17 if stocking and stocking:dist interaction included)
	real<lower=0> sigmaPolField; //SD of field random intercepts	
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts
	// Unused terms:
	// real<lower=0> sigmaPolPlot; //SD of plot random intercepts - Rhat 1.4, small n_eff, strongly correlated with lp__
	// vector[Nplot] intPollen_plot; //plot-level random intercepts
	// real slopeFlyVisPol; //Slope of fly visits - p=0.577
	// real slopeStockingPollen; //Hive number effect - p=0.613
	// real slopeStockingHbeeDistPollen; //Hive number:hbee distance interaction - p=0.949
	
	// // Flower count (per plant) 
	// real intFlwCount; //Intercept
	// real slopePlSizeFlwCount; //Slope of plant size
	// real<lower=0.01> sigmaFlwCount_field; //SD of field-level random effect
	// real<lower=0.01> sigmaFlwCount_plot; //SD of plot-level random effect
	// vector[Nfield] intFlwCount_field; //Field-level random effect
	// vector[Nplot] intFlwCount_plot; //Plot-level random effect	
	// real<lower=0.01> flwCountPhi; //Dispersion parameter
	
	// // Flower survival
	// // plot-level random intercepts have low n_eff, and basically all intercepts overlap zero; sigma plot correlated with lp__
	// // slopeFlwCountSurv and slopePlSizeSurv are correlated (-0.77), and basically represent the same thing, so using only plant size
	// real intFlwSurv; //Intercept
	// real slopeVisitSurv; //Slope of hbee visits
	// real slopePolSurv; //Slope of pollen deposition
	// real slopePlSizeSurv; //Slope of plant size
	// real slopePlDensSurv; //Slope of plant density
	// // real slopeFlwCountSurv; //Slope of flower count 
	// real slopeIrrigSurv; //Slope of irrigation - p=0.56
	// real slope2015Surv; //Slope of year - p=0.68
	// real slopeIrrig2015Surv; //Irrigation:year interaction - p=0.51
	// real slopePlSizeIrrigSurv; //Plant size:Irrigation interaction - p=0.58
	// // real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts	
	// real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	// vector[Nfield] intFlwSurv_field; //field-level random intercepts
	// // vector[Nplot] intFlwSurv_plot; 
	// real<lower=0> flwSurvPhi; //Dispersion parameter for beta-binomial
	
	// Seed count
	real intSeedCount; //Intercept
	real slopeVisitSeedCount; //Slope of hbee visits
	real slopePolSeedCount; //Slope of pollen deposition
	real slopePlSizeCount; //Slope of plant size	
	real slope2015SeedCount; //Year effect
	real slopeIrrigSeedCount; //Irrigation 
	//Interactions
	real slopeIrrig2015SeedCount; //Irrigation:year
	real slopePlSizeIrrigSeedCount; //Plant size:irrigation
	real slopePlSize2015SeedCount; //Plant size:year
	real slopePlSizeIrrig2015SeedCount; //Plant size:irrigation:2015	
	real<lower=0> sigmaSeedCount_plant; //SD of plant random effect - OK
	// real<lower=0> sigmaSeedCount_plot; //SD of plot random effect - not converging well, Rhat 1.1, small n_eff. 
	real<lower=0> sigmaSeedCount_field; //SD of field random effect - OK
	vector[Nplant] intSeedCount_plant; //plant-level random intercepts
	// vector[Nplot] intSeedCount_plot; //plot-level random intercepts	
	vector[Nfield] intSeedCount_field; //field-level random intercepts 
	real<lower=0> seedCountPhi; //Dispersion parameter
	
	// Weight per seed
	real intSeedWeight; //Intercept
	real slopeVisitSeedWeight; //Slope of hbee visits
	real slopePolSeedWeight; //Slope of pollen deposition
	real slopeSeedCount; //Slope of seed count
	real slopePlSizeWeight; //Slope of plant size
	real slopeIrrigSeedWeight; //Irrigation effect
	real slope2015SeedWeight; //Slope of 2015
	real slope2015IrrigSeedWeight; //Year:irrigation 
	real slopeSeedCountPlSizeSeedWeight; //SeedCount:plant size 
	real slopePlSizeIrrigSeedWeight; // plant size :irrigation 
	real slopeSeedCount2015SeedWeight; // SeedCount:Year
	real slopePlSize2015SeedWeight; // Plant size:Year	
	real slopePlSizeIrrig2015SeedWeight; // Plant size:Irrigated:Year			
	real<lower=0> sigmaSeedWeight; //SD of seed weight
	real<lower=0> sigmaSeedWeight_plant; //SD of plant random effect - OK
	// real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect - not converging well, Rhat 1.08, small n_eff
	real<lower=0> sigmaSeedWeight_field; //SD of field random effect - OK	
	vector[Nplant] intSeedWeight_plant; //plant-level random intercepts		
	// vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
	
	// Total yield (g/plant)	
	real intYield; //Intercept for predicted yield
	real slopeYield; //Proportion of predicted yield that becomes yield		
	vector<lower=0>[2] sigmaYield_field; //SD of field-level intercept/slopes
	vector<lower=0>[2] sigmaYield_plot; //SD of plot-level intercept/slope	
	real<lower=0> sigmaYield; //SD of plant-level yield	
	cholesky_factor_corr[2] L_field; //Cholesky-decomposed correlation matrices
	cholesky_factor_corr[2] L_plot;
	matrix[2,Nfield] zYield_field; //Unit normals for matrix correlation trick
	matrix[2,Nplot] zYield_plot; 	
}

transformed parameters {		
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplant] plSizeMu; //Plant size	
	vector[Nplot] plSizePlotMu; //Plot-level plant size	
	// vector[Nplot] flDensMu; //Predicted flower density	
	// vector[Nplot] visitHbeeMu; //Plot-level hbee visits	
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma		
	// vector[Nplot] flwCountPlot; //Plot-level flowers per plant
	// vector[Nplant] flwCountMu; //Expected flowers per plant
	// vector[Nplot] flwSurvPlot; //Plot-level flower survival
	// vector[Nplant] flwSurv; //Flower survival rate (logit)
	vector[Nplot] seedCountMuPlot; //Plot-level seed count
	vector[Nplant] seedCountMuPlant; //Plant-level seed count
	vector[Npod] seedCountMu; //Pod-level seed counts	
	vector[Nplot] seedWeightMuPlot; //Plot-level weight per seed
	vector[Nplant] seedWeightPlantMu; //Plant-level weight per seed
	vector[Npod] seedWeightMu; //Pod-level weight per seed				
	vector[Nplant] calcYield; //Calculated yield per plant (from pod count, seed weight, seed count)
	vector[Nplant] logYieldMu; //Predicted log(yield) per plant (calculated x coef)
	//Matrix to store random slopes and intercepts for yield
	matrix[2,Nfield] ranEffYield_field; //Field level random effects
	matrix[2,Nplot] ranEffYield_plot; //Plot level random effects
	
	// Imputed missing data;
	vector[Nplot] plDens; //Planting density	
	plDens[obsPlDens_ind]=plDens_obs;
	plDens[missPlDens_ind]=plDens_miss;	
	
	for(i in 1:Nplot){ 
		// Plant density per plot
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex[i]] + 
			slope2015PlDens*is2015[plotIndex[i]]+ //Year effect
			slopeIrrigPlDens*isIrrigated[plotIndex[i]]+ //Irrigation effect
			slope2015IrrigPlDens*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]+ //Year:irrigation interaction
			slopeDistPlDens*logHbeeDist[i] + //Distance effect				
			slopeGPPlDens*isGP[plotIndex[i]]; //GP effect
			
		// Plant size (plot-level) = intercept + random field int + random plot int + distance + planting density effect 		
		plSizePlotMu[i] = intPlSize + intPlSize_field[plotIndex[i]] + intPlSize_plot[i] + 						
			slopePlDensPlSize*plDens[i] +  //Planting density
			slopeDistPlSize*logHbeeDist[i] + //Distance effect (edge of field has smaller plants)			
			slopeGpPlSize*isGP[plotIndex[i]] + // Grand Prairie
			slope2015PlSize*is2015[plotIndex[i]] + //2015
			slopeIrrigPlSize*isIrrigated[plotIndex[i]] + //Irrigation effect						
			slopeStockingPlSize*numHives[plotIndex[i]] + //Stocking
			slopeDistStockingPlSize*numHives[plotIndex[i]]*logHbeeDist[i] + //Stocking:Distance interaction		
			slopeDist2015PlSize*logHbeeDist[i]*is2015[plotIndex[i]]; //Distance:year interaction
			
		// // Flower density = intercept + random field int + 
		// flDensMu[i] = intFlDens	+ intFlDens_field[plotIndex[i]] + 
			// slopePlSizeFlDens*plSizePlotMu[i] + //plant size effect 			
			// slopeHbeeDistFlDens*logHbeeDist[i]; //distance effect
	
		// // Hbee Visitation = intercept + random int + time offset + 	
		// visitHbeeMu[i] = intVisit + intVisit_field[plotIndex[i]] + logTime[i] + 
			// slopeYearVis*is2015[plotIndex[i]] + //Year effect
			// slopeGpVis* +isGP[plotIndex[i]] + //Grand Prairie effect
			// slopeYearGpVis*is2015[plotIndex[i]]*isGP[plotIndex[i]] + //Year:area interaction			
			// slopeDistVis*logHbeeDist[i] + //distance from edge 
			// slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			// slopeFlDens*flDens[i] + //Flower density
			// slopeIrrigVis*isIrrigated[plotIndex[i]]; //Irrigation
		
		// Plot-level pollen deposition = random int field + random int plot + 
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //intPollen_plot[i] + 
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits
			// slopeFlyVisPol*logFlyVis[i] + //(log) fly visits
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect						
		//Switched global intercept to flower level term to "center" plot level measurements
			
		// // Flower count per plant (plot level) = intercept + random field int + random plot int
		// flwCountPlot[i] = intFlwCount + intFlwCount_field[plotIndex[i]] + intFlwCount_plot[i];
		
		// // Plot-level flower survival = intercept + random int field + random int plot + 
		// flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + //intFlwSurv_plot[i] + 
			// slopeVisitSurv*logHbeeVis[i] + //hbee visits 
			// slopePolSurv*pollenPlot[i] + //(log) pollen deposition - large correlation b/w slopePolSurv and intFlwSurv
			// slopePlDensSurv*plDens[i]; //Plant density			
			// // slopeIrrigSurv*isIrrigated[plotIndex[i]] + //Slope of irrigation
			// // slope2015Surv*is2015[plotIndex[i]] + //Slope of year
			// // slopeIrrig2015Surv*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]; //Irrigation:year interaction
			
		// Plot-level seed count = intercept + random int field + random int plot + random int plant + 
		seedCountMuPlot[i] = intSeedCount + intSeedCount_field[plotIndex[i]] + //intSeedCount_plot[i] + 
			slopeVisitSeedCount*logHbeeVis[i] + //(log) hbee visits 
			slopePolSeedCount*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedCount and intFlwSurv
			slope2015SeedCount*is2015[plotIndex[i]] + //Year effect
			slopeIrrigSeedCount*isIrrigated[plotIndex[i]] + //Irrigation
			slopeIrrig2015SeedCount*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]; //Irrigation:year
			
		// Plot-level seed weight = intercept + random int field + random int plot + 	
		seedWeightMuPlot[i] = intSeedWeight + intSeedWeight_field[plotIndex[i]] + //intSeedWeight_plot[i] + 
			slopeVisitSeedWeight*logHbeeVis[i] + //(log) hbee visits 
			slopePolSeedWeight*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedWeight and intFlwSurv
			slopeIrrigSeedWeight*isIrrigated[plotIndex[i]] + //Irrigation effect
			slope2015SeedWeight*is2015[plotIndex[i]] + // Year effect
			slope2015IrrigSeedWeight*is2015[plotIndex[i]]*isIrrigated[plotIndex[i]]; //Year:irrigation effect
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //global intercept + plot-level pollen 
		
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]]; 			
		
		// // Flower count per plant
		// flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count 
			// slopePlSizeFlwCount*plantSize[i]; //individual size effect	
	
		// // Flower survival per plant
		// flwSurv[i] = flwSurvPlot[plantIndex[i]] + slopePlSizeSurv*plantSize[i]; //Plot-level plant survival + size effect		
			// // slopePlSizeIrrigSurv*(plantSize[i]*isIrrigated[plotIndex[plantIndex[i]]]) + //Plant size:Irrigation interaction
			// // slopeFlwCountSurv*logFlwCount[i]; //log-Flower count effect
			
		// Seed count per pod = plot-level effect + random plant int +		
		seedCountMuPlant[i] = seedCountMuPlot[plantIndex[i]] + 
			intSeedCount_plant[i] + 			
			slopePlSizeCount*plantSize[i] + //plant size
			slopePlSizeIrrigSeedCount*plantSize[i]*isIrrigated[plotIndex[plantIndex[i]]] + //Plant size:irrigation
			slopePlSize2015SeedCount*plantSize[i]*is2015[plotIndex[plantIndex[i]]] + //Plant size:year
			slopePlSizeIrrig2015SeedCount*plantSize[i]*isIrrigated[plotIndex[plantIndex[i]]]*is2015[plotIndex[plantIndex[i]]]; //Plant size:irrigation:year
			
		// Weight per seed = plot-level effect + random int plant + 
		seedWeightPlantMu[i] = seedWeightMuPlot[plantIndex[i]] + intSeedWeight_plant[i] +			
			slopePlSizeWeight*plantSize[i] + //Plant size 		
			slopePlSizeIrrigSeedWeight*isIrrigated[plotIndex[plantIndex[i]]]*plantSize[i] + // Plant size :irrigation 
			slopePlSize2015SeedWeight*is2015[plotIndex[plantIndex[i]]]*plantSize[i] + //Plant size: year
			slopePlSizeIrrig2015SeedWeight*isIrrigated[plotIndex[plantIndex[i]]]*is2015[plotIndex[plantIndex[i]]]*plantSize[i]; //Plant size:Irrigated:Year			
	}	
	
	for(i in 1:Npod){ //For each pod
		seedCountMu[i] = seedCountMuPlant[podIndex[i]]; 
		//Plant-level effect + effect of seedCount (do pods with lots of seed also have bigger seeds?)		
		seedWeightMu[i] = seedWeightPlantMu[podIndex[i]]+
			slopeSeedCount*seedCount[i] + //Seed count effect
			slopeSeedCount2015SeedWeight*seedCount[i]*is2015[plotIndex[plantIndex[podIndex[i]]]] + // SeedCount:Year
			slopeSeedCountPlSizeSeedWeight*seedCount[i]*plantSize[podIndex[i]]; //Seed count: plant size		
	}		
	
	//Generate correlated random effects matrices
	ranEffYield_field = diag_pre_multiply(sigmaYield_field,L_field) * zYield_field;
	ranEffYield_plot = diag_pre_multiply(sigmaYield_plot,L_plot) * zYield_plot; 	
	
	for(i in 1:Nplant){
		// Calculated yield per plant (g) = weight per seed (mg)/1000 * seeds per pod * pods per plant 		
		calcYield[i] = (seedWeightPlantMu[i]+slopeSeedCount*exp(seedCountMuPlant[i])+(1/lambdaSeedWeight))/1000 * //weight per seed(g)
			exp(seedCountMuPlant[i]) * //seeds per pod
			podCount[i]; //pods per plant
		// Predicted yield = intercept + random effects + 
		logYieldMu[i] = (intYield+ranEffYield_field[1,plotIndex[plantIndex[i]]]+ranEffYield_plot[1,plantIndex[i]]) + 
			log(calcYield[i])*(slopeYield+ranEffYield_field[2,plotIndex[plantIndex[i]]]+ranEffYield_plot[2,plantIndex[i]]); //Slope + random effects
	}
}
	
model {	
	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	// flDens ~ normal(flDensMu,sigmaFlDens); //Flower density per plot
	// hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	// flwCount ~ neg_binomial_2_log(flwCountMu,flwCountPhi); //Flower count per plant (attempted pods)
	// podCount ~ beta_binomial(flwCount,inv_logit(flwSurv)*flwSurvPhi,(1-inv_logit(flwSurv))*flwSurvPhi); //Flower survival (surviving pods)
	seedCount ~ neg_binomial_2_log(seedCountMu,seedCountPhi); //Seed count per pod
	seedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Weight per seed	(mg)
	logYield ~ normal(logYieldMu,sigmaYield); //Seed yield per plant
		
	//Priors
	//Plant density	- informative priors
	intPlDens ~ normal(0,0.1); //Global intercept
	slope2015PlDens ~ normal(0,1); //Year effect
	slopeIrrigPlDens ~ normal(0,1); //Irrigation effect
	slope2015IrrigPlDens ~ normal(0,1); //Year:irrigation interaction
	slopeDistPlDens ~ normal(0,0.05); //Slope of distance into field	
	slopeGPPlDens ~ normal(0,1); // Grand Prairie effect 
	sigmaPlDens ~ gamma(3,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(3,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field
	
	// Plant size - informative priors
	intPlSize ~ normal(0,0.5); //Intercept
	slopePlDensPlSize ~ normal(0,0.5); //Plant density
	slopeDistPlSize ~ normal(0,0.05); //Distance effect
	slopeGpPlSize ~ normal(0,0.5); //Grand Prairie effect
	slopeIrrigPlSize ~ normal(0,0.5); //Irrigation effect
	slope2015PlSize ~ normal(0.3,0.5); //2015 effect		
	slopeStockingPlSize ~ normal(0,0.1); //Stocking effect
	slopeDistStockingPlSize ~ normal(0,0.1); //Distance:Stocking interaction
	slopeDist2015PlSize ~ normal(0,0.5); //Distance:Year interaction
	sigmaPlSize_field ~ gamma(1,1); //Sigma for random field 
	sigmaPlSize_plot ~ gamma(1,1); //Sigma for random plot
	sigmaPlSize ~ gamma(7,10); //Sigma for residual	
	intPlSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intPlSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot	
	
	// // Flower density per plot
	// intFlDens ~ normal(0,1); //Global intercept
	// slopePlSizeFlDens ~ normal(2,1); //plant size effect
	// slopeHbeeDistFlDens ~ normal(0,1); //distance into field
	// sigmaFlDens ~ gamma(7,2); //Sigma for within-field (residual)	
	// sigmaFlDens_field ~ gamma(4,2); //Sigma for field
	// intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field			
	
	// // Visitation - informative priors
	// intVisit ~ normal(-1,1); //Intercept	
	// slopeDistVis ~ normal(-0.3,0.3); //Slope of distance effect on hbee visits
	// slopeHiveVis ~ normal(0.6,0.5); //Slope of hive effect on visits 
	// slopeYearVis ~ normal(0.5,1); //2015 effect
	// slopeGpVis ~ normal(0,1); //Effect of Grand Prairie	
	// slopeYearGpVis ~ normal(0,1); // GP-year interaction
	// slopeFlDens ~ normal(0,0.5); //Flower density 
	// sigmaVisField ~ gamma(4,4); //Sigma for random field 		
	// intVisit_field ~ exp_mod_normal(0,sigmaVisField,lambdaVisField); //Skewed random effects - slightly better than standard normal
	// lambdaVisField ~ gamma(4,2); //Lambda for skewed random effects			
	// visitHbeePhi ~ gamma(4,10); //Dispersion parameter for NegBin		
		
	// Pollen deposition - informative priors	
	intPollen ~ normal(5.5,1); //Intercept	
	slopeVisitPol ~ normal(0,0.1); //hbee Visitation effect	
	slopeHbeeDistPollen ~ normal(0,0.1); //hbee distance effect	
	sigmaPolField ~ gamma(5,10); //Sigma for random field	
	pollenPhi ~ gamma(7,10); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	// intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int		
	// sigmaPolPlot ~ gamma(3,10); //Sigma for random plot 

	// //Flower count (per plant)
	// intFlwCount ~ normal(6,0.5); //Intercept
	// slopePlSizeFlwCount ~ normal(0.9,0.5); //Slope of plant size
	// sigmaFlwCount_field ~ gamma(1,10); //SD of field-level random effect	
	// sigmaFlwCount_plot ~ gamma(1,10); //SD of plot-level random effect
	// intFlwCount_field ~ normal(0,sigmaFlwCount_field); //Field-level random effect	
	// intFlwCount_plot ~ normal(0,sigmaFlwCount_plot); //Plot-level random effect	
	// flwCountPhi ~ normal(35,1); //Dispersion parameter	
	
	// //Flower survival - informative priors
	// intFlwSurv ~ normal(1,1); //Intercept
	// slopeVisitSurv ~ normal(0,0.05); //Slope of hbee visits
	// slopePolSurv ~ normal(0,0.5); //Slope of pollen deposition
	// slopePlSizeSurv ~ normal(0.02,0.05); //Slope of plant size
	// slopePlDensSurv ~ normal(0,1); //Slope of planting density
	// // slopeFlwCountSurv ~ normal(0,1); //Slope of flower count 
	// slopeIrrigSurv ~ normal(0,0.5); //Slope of irrigation
	// slope2015Surv ~ normal(0,0.5); //Slope of year
	// slopeIrrig2015Surv ~ normal(0,0.5); //Irrigation:year interaction
	// slopePlSizeIrrigSurv ~ normal(0,0.5); //Plant size:Irrigation interaction	
	// sigmaFlwSurv_field ~ gamma(3,10); //SD of field random effect
	// // sigmaFlwSurv_plot ~ gamma(2,10); //SD of plot random effect	
	// intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	// // intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts
	// flwSurvPhi ~ gamma(2,2); //Dispersion parameter
	
	// Seed count - informative priors
	intSeedCount ~ normal(3.5,1); //Intercept
	slopeVisitSeedCount ~ normal(0,0.1); //Slope of hbee visits
	slopePolSeedCount ~ normal(0,0.5); //Slope of pollen deposition
	slopePlSizeCount ~ normal(0,0.05); //Slope of plant size
	slope2015SeedCount ~ normal(0,0.5); //Year effect
	slopeIrrigSeedCount ~ normal(0,0.5); //Irrigation
	//Interactions
	slopeIrrig2015SeedCount ~ normal(0,0.5); //Irrigation:year
	slopePlSizeIrrigSeedCount ~ normal(0,0.1); //Plant size:irrigation
	slopePlSize2015SeedCount ~ normal(0,0.1); //Plant size:year
	slopePlSizeIrrig2015SeedCount ~ normal(0,0.5); //Plant size:irrigation:year	
	seedCountPhi ~ normal(22,1); //Dispersion parameter
	sigmaSeedCount_field ~ gamma(2,10); //SD of field random effect
	// sigmaSeedCount_plot ~ gamma(2,10); //SD of plot random effect
	sigmaSeedCount_plant ~ gamma(2,10); //SD of plant random effect	
	intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	// intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	intSeedCount_plant ~ normal(0,sigmaSeedCount_plant); //plant-level random intercepts	
		
	// Weight per seed - informative priors
	intSeedWeight ~ normal(1.5,1); //Intercept
	slopeVisitSeedWeight ~ normal(0,0.1); //Slope of hbee visits
	slopePolSeedWeight ~ normal(0,0.5); //Slope of pollen deposition
	slopeSeedCount ~ normal(0,0.05); //Slope of seed count
	slopePlSizeWeight ~ normal(0,0.5); //Slope of plant size
	slopeIrrigSeedWeight ~ normal(0,0.5); //Slope of irrigation
	slope2015SeedWeight ~ normal(0,0.5); //Slope of 2015
	//Interactions	
	slope2015IrrigSeedWeight ~ normal(0,0.5); //Year:irrigation 
	slopeSeedCountPlSizeSeedWeight ~ normal(0,0.1); //SeedCount:plant size 
	slopePlSizeIrrigSeedWeight ~ normal(0,0.5); // plant size :irrigation 
	slopeSeedCount2015SeedWeight ~ normal(0,0.1); // SeedCount:Year
	slopePlSize2015SeedWeight ~ normal(0,0.5); // VegMass:Year	
	slopePlSizeIrrig2015SeedWeight ~ normal(0,0.5); // VegMass:Irrigated:Year		
	sigmaSeedWeight ~ gamma(3,6); //SD of seed weight
	sigmaSeedWeight_field ~ gamma(4,10); //SD of field random effect	
	// sigmaSeedWeight_plot ~ gamma(3,10); //SD of plot random effect
	sigmaSeedWeight_plant ~ gamma(6,10); //SD of plant random effect		
	intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	// intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	intSeedWeight_plant ~ normal(0,sigmaSeedWeight_plant); //plant-level random intercepts	
	lambdaSeedWeight ~ gamma(15,10); //Lambda for exp-normal distribution
			
	//Yield per plant
	intYield ~ normal(0,1); //Intercept
	slopeYield ~ normal(1,1); //Slope of calculated yield	
	sigmaYield ~ gamma(2,5); //Sigma for yield			
	//Correlated random effects:
	sigmaYield_field ~ gamma(2,2); //SD of field-level yield intercepts/slopes
	sigmaYield_plot ~ gamma(2,2); //SD of plot-level yield intercepts/slopes
	to_vector(zYield_field) ~ normal(0,1); //Unit normals for correlated random effects
	to_vector(zYield_plot) ~ normal(0,1);
	L_field ~ lkj_corr_cholesky(2); //Standard prior for lkj cholesky matrix - higher values make extreme correlations less likely (see p.394 in McElreath 2016)	
	L_plot ~ lkj_corr_cholesky(2); 
}

generated quantities {
	//Plot-level quantities
	// // planting density
	// real predPlDens[Nplot]; //Generated
	// real plDens_resid[Nplot]; //Residual
	// real predPlDens_resid[Nplot]; //Residual of generated
	// // flower density
	// real predFlDens[Nplot]; //Generated
	// real flDens_resid[Nplot]; //Residual
	// real predFlDens_resid[Nplot]; //Residual of generated	
	// //hbeeVis
	// int predHbeeVis[Nplot]; //Generated
	// real hbeeVis_resid[Nplot]; //Residual 
	// real predHbeeVis_resid[Nplot]; //Residual of generated 	
	
	// //Flower-level
	// //pollen deposition
	// int predPollenCount[Nflw]; //Generated 
	// real pollen_resid[Nflw]; //residual
	// real predPollen_resid[Nflw]; //residual of generated
	
	//Plant-level	
	// //plantSize
	// real predPlSize[Nplant]; //Generated
	// real plSize_resid[Nplant]; //Residual
	// real predPlSize_resid[Nplant]; //Residual of generated
	// // flower count per plant (potential pods)
	// int predFlwCount[Nplant]; //Generated
	// real flwCount_resid[Nplant]; //Residual
	// real predFlwCount_resid[Nplant]; //Residual of generated	
	// // flower survival (surviving pods)
	// int predPodCount[Nplant]; //Generated
	// real podCount_resid[Nplant]; //Residual
	// real predPodCount_resid[Nplant]; //Residual of generated
	// (log) yield per plant
	real predYield[Nplant]; //Generated
	real yield_resid[Nplant]; //Residual
	real predYield_resid[Nplant]; //Residual of generated
	
	//Pod-level
	//seeds per pod
	int predSeedCount[Npod]; //Generated
	real seedCount_resid[Npod]; //Residual
	real predSeedCount_resid[Npod]; //Residual of generated	
	//weight per seed
	real predSeedMass[Npod]; //Generated
	real seedMass_resid[Npod]; //Residual
	real predSeedMass_resid[Npod]; //Residual of generated
		
	// for(i in 1:Nplot){
		// // plant density		
		// plDens_resid[i] = plDens[i] - plDensMu[i]; //Residual for actual value
		// predPlDens[i]= normal_rng(plDensMu[i],sigmaPlDens); //Generated value from normal
		// predPlDens_resid[i] = predPlDens[i] - plDensMu[i]; //Residual for predicted value							
	
		// // flower density
		// flDens_resid[i] = flDens[i]-flDensMu[i]; //Residual for actual value
		// predFlDens[i] = normal_rng(flDensMu[i],sigmaFlDens); //Generated value from normal
		// predFlDens_resid[i] = predFlDens[i] - flDensMu[i]; //Residual for predicted value	
			
		// bee visits (NB version)
		// hbeeVis_resid[i]=hbeeVis[i]-(exp(visitHbeeMu[i])); //Residual for actual value
		// predHbeeVis[i] = neg_binomial_2_log_rng(visitHbeeMu[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		// predHbeeVis_resid[i]=predHbeeVis[i]-(exp(visitHbeeMu[i])); //residual for predicted value					
		// hbeeVis_resid[i]=logHbeeVis[i]-visitHbeeMu[i]; //Residual for actual value
		// predHbeeVis[i] = neg_binomial_2_log_rng(visitHbeeMu[i],visitHbeePhi); //Predicted value drawn from neg.bin		
		// predHbeeVis_resid[i]=log(predHbeeVis[i]+0.5)-visitHbeeMu[i]; //residual for predicted value							
	// }		
	
	// for(i in 1:Nflw){
		// // pollen deposition
		// pollen_resid[i]= exp(pollenMu[i]) - pollenCount[i]; //Residual for actual value
		// predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulate pollen counts
		// predPollen_resid[i] = exp(pollenMu[i]) - predPollenCount[i]; //Residual for predicted		
	// }
			
	for(i in 1:Nplant){
		// //plant size
		// plSize_resid[i]= plantSize[i] - plSizeMu[i]; //Residual for actual
		// predPlSize[i] = normal_rng(plSizeMu[i],sigmaPlSize); //Generates new value from normal dist.
		// predPlSize_resid[i] = predPlSize[i] - plSizeMu[i]; //Residual for new value	
		// // flower count per plant
		// flwCount_resid[i] = flwCount[i] - exp(flwCountMu[i]); //Residual for actual
		// predFlwCount[i] = neg_binomial_2_log_rng(flwCountMu[i],flwCountPhi); //Generates new value from neg. bin.
		// predFlwCount_resid[i] = predFlwCount[i] - exp(flwCountMu[i]); //Residual for new value
		// // pod count (surviving pods)
		// podCount_resid[i] = podCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for actual		
		// predPodCount[i] = beta_binomial_rng(flwCount[i],inv_logit(flwSurv[i])*flwSurvPhi,(1-inv_logit(flwSurv[i]))*flwSurvPhi); //Generates new value from beta-binomial
		// predPodCount_resid[i] = predPodCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for new value				
		// (log) yield per plant
		yield_resid[i]= logYield[i] - logYieldMu[i]; //Residual for actual
		predYield[i] = normal_rng(logYieldMu[i],sigmaYield); //Generates new value from normal dist.
		predYield_resid[i] = predYield[i] - logYieldMu[i]; //Residual for new value	
	}
	
	for(i in 1:Npod){ //For each pod
		//Seed count per pod - doesn't work well due to weird generating process (I think)
		seedCount_resid[i] = seedCount[i] - exp(seedCountMu[i]);
		predSeedCount[i] = neg_binomial_2_log_rng(seedCountMu[i],seedCountPhi); 
		predSeedCount_resid[i] = predSeedCount[i] - exp(seedCountMu[i]);
		// weight per seed - exp-normal works well		
		seedMass_resid[i] = seedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight)); 
		predSeedMass[i] = exp_mod_normal_rng(seedWeightMu[i],sigmaSeedWeight,lambdaSeedWeight); 
		predSeedMass_resid[i] = predSeedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));
	}	
}
