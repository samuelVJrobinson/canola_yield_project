//Commodity canola pollination model used for thesis
//  Uses plant-level data

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
	vector[Nplant] avgSeedCount; //Average seeds per pod
	vector[Nplant] avgSeedMass; //Average weight per seed (mg)
	vector[Nplant] yield; //Observed yield per plant (g)	
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives	
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	vector[Nplot] logFlyVis; //Fly visitation rate
	vector[Nplant] logitFlwSurv; //(logit) proportion flower survival	
	vector[Nplant] logFlwCount; //Log flower count
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)
	vector[Nplant] calcYield; //Predicted yield per plant (g seed per plant)	
	vector[Nplant] logCalcYield; //Predicted log yield
	
	logHbeeDist=logHbeeDist-mean(logHbeeDist); //Centers distance
	
	for(i in 1:Nfield){
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives		
	}
	
	for(i in 1:Nplot){ //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5); 
		logFlyVis[i] = log((flyVis[i]/totalTime[i])+0.5);
	}
		
	for(i in 1:Nplant){
		logFlwCount[i] = log(flwCount[i]); //Log flower count per plant	
		//Necessary for promoting integers to reals. Otherwise does integer division.
		logitFlwSurv[i] = podCount[i]; 
		logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods		
		if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
			logitFlwSurv[i]=0.01;
		else if(logitFlwSurv[i]>=1)
			logitFlwSurv[i]=0.99;	
		//calcYield = log(pod count x seed weight x seed count)
		calcYield[i] = podCount[i]*avgSeedCount[i]*(avgSeedMass[i]/1000);
		logCalcYield[i] = log(calcYield[i]);
	}		
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
	
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
	real<lower=0> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field
	
	// // Plant size 
	// real intPlSize; //Global intercept
	// real slopePlDensPlSize; //Slope of planting density	
	// real slopeDistPlSize; //Slope of distance 
	// real slopeGpPlSize; //Slope of Grand Prairie effect	
	// real slope2015PlSize; //Slope of 2015 effect		
	// real slopeIrrigPlSize; //Slope of irrigation		
	// real<lower=0> sigmaPlSize_field; //SD of field-level intercept
	// real<lower=0> sigmaPlSize_plot; //SD of plot-level intercept
	// real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	// vector[Nfield] intPlSize_field; //Random intercept for field
	// vector[Nplot] intPlSize_plot; //Random intercept for plot - not converging well, but looic is much worse without it		
	// 
	// // Flower density per plot
	// real intFlDens; //Global intercept
	// real slopePlSizeFlDens; //Slope of plant size on flower density
	// real slopeHbeeDistFlDens; //Slope of distance into field
	// real<lower=0> sigmaFlDens; //Sigma for within-field (residual)	
	// real<lower=0> sigmaFlDens_field; //Sigma for field
	// vector[Nfield] intFlDens_field; //Random intercept for field		
	// 
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
	real slopeHbeeDistPollen; //Effect of distance into field - p=0.076 (p=0.17 if stocking and stocking:dist interaction included - see below)
	real<lower=0> sigmaPolField; //SD of field random intercepts	
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts	
	// real<lower=0> sigmaPolPlot; //SD of plot random intercepts
	// vector[Nplot] intPollen_plot; //plot-level random intercepts
	
	// Flower count (per plant)
	real intFlwCount; //Intercept
	real slopePlSizeFlwCount; //Slope of plant size
	real slopeSurvFlwCount; //Slope of flower survival rate
	real slope2015FlwCount; //Slope of 2015 effect
	real<lower=0> phiFlwCount_field; //SD of field-level random effect
	vector[Nfield] intFlwCount_field; //Field-level random effect
	real<lower=0> phiFlwCount_plot; //SD of plot-level random effect
	vector[Nplot] intFlwCount_plot; //Plot-level random effect
	//Sigma (variance) modeled as a function of plant size
	real intPhiFlwCount; //Intercept for sigma
	real slopePlSizePhiFlwCount; //Effect of plant size on sigma
	real<lower=0> sigmaPhiFlwCount_field; //Sigma for field level sigma
	vector[Nfield] intPhiFlwCount_field; //Field-level random intercept
	// real<lower=0> sigmaPhiFlwCount_plot; //Sigma for plot level sigma
	// vector[Nplot] intPhiFlwCount_plot; //Plot-level random intercept
	
  // Flower survival
	// plot-level random intercepts have low n_eff, and basically all intercepts overlap zero; sigma plot correlated with lp__
	// slopeFlwCountSurv and slopePlSizeSurv are correlated (-0.77), and basically represent the same thing, so using only plant size
	real intFlwSurv; //Intercept
	real slopeVisitSurv; //Slope of hbee visits
	real slopePolSurv; //Slope of pollen deposition
	real slopePlSizeSurv; //Slope of plant size
	real slopePlDensSurv; //Slope of plant density
	real slopeIrrigSurv; //Slope of irrigation - p=0.56
	real slope2015Surv; //Slope of year - p=0.68
	real slopePlSizePlDensSurv; //Plant size:plant density
	real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	vector[Nfield] intFlwSurv_field; //field-level random intercepts
	// real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts
	// vector[Nplot] intFlwSurv_plot; //plot-level random intercepts
	// Sigma (variance) modeled as function of plant size
	real intPhiFlwSurv; //Intercept for sigma
	real slopePlSizePhiFlwSurv; //Effect of plant size on phi
	real<lower=0> sigmaPhiFlwSurv_field; //Sigma for field level sigma
	vector[Nfield] intPhiFlwSurv_field; //Field-level random intercept
	// real<lower=0> sigmaPhiFlwSurv_plot; //Sigma for plot level
	// vector[Nplot] intPhiFlwSurv_plot; //Plot-level random intercept
	
	// // Seed count
	// real intSeedCount; //Intercept
	// real slopeVisitSeedCount; //Slope of hbee visits - p=0.81
	// real slopePolSeedCount; //Slope of pollen deposition - p=0.74
	// real slopePlSizeCount; //Slope of plant size - p=0.22	
	// real slope2015SeedCount; //Year effect - p=0.0003
	// // real slopeIrrigSeedCount; //Irrigation - p=0.88	
	// // slopeIrrig2015SeedCount; //Irrigation:2015 - p=0.4736
	// // slopePlSizeIrrigSeedCount; //Plant size:irrigation - p=0.6772
	// // slopePlSize2015SeedCount; //Plant size:2015 - p=0.2196
	// // slopePlSizeIrrig2015SeedCount; //Plant size:irrigation:2015 - p=0.9639	
	// real<lower=0> sigmaSeedCount_plot; //SD of plot random effect - not converging well, Rhat 1.1, small n_eff. 
	// real<lower=0> sigmaSeedCount_field; //SD of field random effect - OK
	// vector[Nplot] intSeedCount_plot; //plot-level random intercepts	
	// vector[Nfield] intSeedCount_field; //field-level random intercepts 
	// real<lower=0> sigmaSeedCount; //SD of seed count
	// real<lower=0> lambdaSeedCount; //Lambda term for exponential process
	// 
	// // Weight per seed
	// // slopeSeedCountPlSizeSeedWeight highly correlated with slopePlSizeWeight (r=-0.8)	
	// // exp-normal distribution works much better than normal (elpd_diff:-71.3,se:23.7)
	// real intSeedWeight; //Intercept
	// real slopeVisitSeedWeight; //Slope of hbee visits - p=0.24
	// real slopePolSeedWeight; //Slope of pollen deposition - p=0.50
	// real slopeSeedCount; //Slope of seed count - p<0.0001
	// real slopePlSizeWeight; //Slope of plant size - p=0.13
	// real slopeIrrigSeedWeight; //Irrigation effect - p=0.04
	// real slope2015SeedWeight; //Slope of 2015 - p=0.28
	// real<lower=0> sigmaSeedWeight; //SD of seed weight
	// real<lower=0> sigmaSeedWeight_plot; //SD of plot random effect
	// real<lower=0> sigmaSeedWeight_field; //SD of field random effect
	// vector[Nplot] intSeedWeight_plot; //plot-level random intercepts	
	// vector[Nfield] intSeedWeight_field; //field-level random intercepts	
	// real<lower=0> lambdaSeedWeight; //Lambda term for exponential process
	// 
	// // Total yield (g/plant)	
	// real intYield; //Intercept for predicted yield
	// real slopeYield; //Proportion of predicted yield that becomes yield		
	// vector<lower=0>[2] sigmaYield_field; //SD of field-level intercept/slopes
	// vector<lower=0>[2] sigmaYield_plot; //SD of plot-level intercept/slope	
	// real<lower=0> sigmaYield; //SD of plant-level yield	
	// cholesky_factor_corr[2] L_field; //Cholesky-decomposed correlation matrices
	// cholesky_factor_corr[2] L_plot;
	// matrix[2,Nfield] zYield_field; //Unit normals for matrix correlation trick
	// matrix[2,Nplot] zYield_plot; 	
	
}

transformed parameters {		
	//Expected values
	
	//Plant density
	vector[Nplot] plDensMu; //Expected plant density
	vector[Nplot] plDens; //Planting density - imputed

	// //Plant size
	// vector[Nplot] plSizePlotMu; //Plot-level plant size		
	// vector[Nplant] plSizeMu; //Expected plant size	
	// 
	// //Flower density
	// vector[Nplot] flDensMu; //Expected flower density	
	// 
	// //Visitation
	// vector[Nplot] visitHbeeMu; //Expected hbee visits	
	
	//Pollen deposition
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Expected pollen per stigma		
	
	//Flower count per plant
	vector[Nplot] flwCountPlot; //Plot-level flowers per plant
	vector[Nplant] flwCountMu; //Expected flowers per plant
	vector<lower=0>[Nplant] phiFlwCount; //Phi (variance) for flowers per plant

	//Flower survival
  vector[Nplot] flwSurvPlot; //Plot-level flower survival
  vector[Nplant] flwSurv; //Flower survival rate (logit)
  vector<lower=0>[Nplant] flwSurvPhi; //Phi for flower survival

//   //Average seeds per pod
//   vector[Nplot] seedCountMuPlot; //Plot-level seed count
//   vector[Nplant] seedCountMu; //Plant-level seed count
//   
//   //Average seed weight
//   vector[Nplot] seedWeightMuPlot; //Plot-level weight per seed
//   vector[Nplant] seedWeightMu; //Plant-level weight per seed
// 	
// 	//Yield
// 	vector[Nplant] logYieldMu; //Predicted log(yield) per plant (calculated x coef)
// 	// Generate correlated random slopes and intercepts matrix
// 	matrix[2,Nfield] ranEffYield_field = diag_pre_multiply(sigmaYield_field,L_field) * zYield_field; //Field level random effects
// 	matrix[2,Nplot] ranEffYield_plot = diag_pre_multiply(sigmaYield_plot,L_plot) * zYield_plot; //Plot level random effects	
	
	//Assign imputed data
	plDens[obsPlDens_ind]=plDens_obs; //Observed data
	plDens[missPlDens_ind]=plDens_miss;	//Missing data

	for(i in 1:Nplot){ 
		// Plant density per plot
		plDensMu[i] = intPlDens + //Intercept
		  intPlDens_field[plotIndex[i]] + //Field level random intercept
			slope2015PlDens*is2015[plotIndex[i]]+ //Year effect
			slopeIrrigPlDens*isIrrigated[plotIndex[i]]+ //Irrigation effect
			slope2015IrrigPlDens*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]+ //Year:irrigation interaction
			slopeDistPlDens*logHbeeDist[i] + //Distance effect
			slopeGPPlDens*isGP[plotIndex[i]]; //Location effect

		// // Plant size (plot-level)
		// plSizePlotMu[i] = intPlSize + //Intercept
		//   intPlSize_field[plotIndex[i]] + //Field-level random effect
		//   intPlSize_plot[i] + //Plot-level random effect
		// 	slopePlDensPlSize*plDens[i] +  //Plant density
		// 	slopeDistPlSize*logHbeeDist[i] + //Distance effect (edge of field has smaller plants)			
		// 	slopeGpPlSize*isGP[plotIndex[i]] + // Location effect
		// 	slope2015PlSize*is2015[plotIndex[i]] + //Year effect
		// 	slopeIrrigPlSize*isIrrigated[plotIndex[i]]; //Irrigation effect								
		// 
		// // Flower density
		// flDensMu[i] = intFlDens	+ //Intercept
		//   intFlDens_field[plotIndex[i]] + //Field-level random effect
		// 	slopePlSizeFlDens*plSizePlotMu[i] + //Plant size effect 			
		// 	slopeHbeeDistFlDens*logHbeeDist[i]; //Distance effect
		// 
		// // Honeybee Visitation
		// visitHbeeMu[i] = intVisit + //Intercept
		//   intVisit_field[plotIndex[i]] + //Field-level random intercept 
		//   logTime[i] + //Time offset
		// 	slopeYearVis*is2015[plotIndex[i]] + //Year effect
		// 	slopeGpVis*isGP[plotIndex[i]] + //Grand Prairie effect
		// 	slopeYearGpVis*is2015[plotIndex[i]]*isGP[plotIndex[i]] + //Year:area interaction			
		// 	slopeDistVis*logHbeeDist[i] + //distance from edge 
		// 	slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
		// 	slopeFlDens*flDens[i] + //Flower density
		// 	slopeIrrigVis*isIrrigated[plotIndex[i]]; //Irrigation
		
		// Plot-level pollen deposition
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //Field-level random intercept
		  // intPollen_plot[i] + //Plot-level random intercept
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect						
		// Global intercept is within flower-level term in order to center plot-level variable
			
    // Flower count per plant (plot level)
    flwCountPlot[i] = intFlwCount + //Intercept
      intFlwCount_field[plotIndex[i]] + //Field-level random intercept
      intFlwCount_plot[i] + //Plot-level random intercept
    	slope2015FlwCount*is2015[plotIndex[i]]; //Year effect
    
    // Plot-level flower survival
    flwSurvPlot[i] = intFlwSurv + //Intercept
      intFlwSurv_field[plotIndex[i]] + //Field-level random intercept
      // intFlwSurv_plot[i] + //Plot-level random intercept
    	slopeVisitSurv*logHbeeVis[i] + //hbee visits
    	slopePolSurv*pollenPlot[i] + //(log) pollen deposition - large correlation b/w slopePolSurv and intFlwSurv
    	slopePlDensSurv*plDens[i] + //Plant density
    	slopeIrrigSurv*isIrrigated[plotIndex[i]] + //Irrigation effect
    	slope2015Surv*is2015[plotIndex[i]]; //Year effect

		// Plot-level seed count
		// seedCountMuPlot[i] = intSeedCount + //Intercept
		//   intSeedCount_field[plotIndex[i]] + //Field-level random intercept
		//   intSeedCount_plot[i] + //Plot-level random intercept
		// 	slopeVisitSeedCount*logHbeeVis[i] + //(log) hbee visits 
		// 	slopePolSeedCount*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedCount and intFlwSurv
		// 	slope2015SeedCount*is2015[plotIndex[i]]; //Year effect		
		// 	
		// // Plot-level seed weight
		// seedWeightMuPlot[i] = intSeedWeight + //Intercept
		// 	intSeedWeight_field[plotIndex[i]] + //Field-level random intercept
		// 	intSeedWeight_plot[i] + //Plot-level random intercept
		// 	slopeVisitSeedWeight*logHbeeVis[i] + //(log) hbee visits 
		// 	slopePolSeedWeight*pollenPlot[i] + //pollen deposition - large correlation b/w slopePolSeedWeight and intFlwSurv
		// 	slopeIrrigSeedWeight*isIrrigated[plotIndex[i]] + //Irrigation effect
		// 	slope2015SeedWeight*is2015[plotIndex[i]]; // Year effect
	}
		
	for(i in 1:Nflw){ //For each flower stigma
		pollenMu[i] = intPollen + //Intercept
		  pollenPlot[flowerIndex[i]]; //Plot-level pollen 
	}
		
	for(i in 1:Nplant){ //For each plant 	
		// //Plant size = plot-level estimate
		// plSizeMu[i] = plSizePlotMu[plantIndex[i]]; 			
		
    // Flower count per plant (attempted pods)
    flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count
    	slopePlSizeFlwCount*plantSize[i] + //Plant size effect
    	slopeSurvFlwCount*logitFlwSurv[i]; //Flower survival effect - do plants with lower flw survival make more flowers?
    // Phi (dispersion) for flower count
    phiFlwCount[i] = exp(intPhiFlwCount + //Intercept
      intPhiFlwCount_field[plotIndex[plantIndex[i]]] + //Field-level random intercept
      // intPhiFlwCount_plot[plantIndex[i]] + //Plot-level random intercept
    	slopePlSizePhiFlwCount*plantSize[i]); // Term for sigma
    	
    // Flower survival per plant
    flwSurv[i] = flwSurvPlot[plantIndex[i]] + //Plot-level plant survival
    	slopePlSizeSurv*plantSize[i] + //Plant size effect
    	slopePlSizePlDensSurv*plantSize[i]*plDens[plantIndex[i]]; //Plant size:plant density
    //Phi (dispersion) for flower survival
    flwSurvPhi[i] = exp(intPhiFlwSurv + //Intercept
      intPhiFlwSurv_field[plotIndex[plantIndex[i]]] + //Field-level random intercept
      // intPhiFlwSurv_plot[plantIndex[i]] + //Plot-level random intercept
      slopePlSizePhiFlwSurv*plantSize[i]);

		// // Seed count per pod = plot-level effect + random plant int +
		// seedCountMu[i] = seedCountMuPlot[plantIndex[i]] + //Plot-level seed count
		// 	slopePlSizeCount*plantSize[i]; //plant size effect
		// 	
		// // Weight per seed = plot-level effect + 
		// seedWeightMu[i] = seedWeightMuPlot[plantIndex[i]] + //Plot-level seed weight
		// 	slopePlSizeWeight*plantSize[i] + //Plant size effect 						
		// 	slopeSeedCount*avgSeedCount[i]; //Seed count effect (do plants with many seeds/pod have bigger seeds?)		
		// 
		// // Predicted yield = intercept + random effects +
		// logYieldMu[i] = (intYield+ranEffYield_field[1,plotIndex[plantIndex[i]]]+ranEffYield_plot[1,plantIndex[i]]) + //Intercept + random effects
		// logCalcYield[i]*(slopeYield+ranEffYield_field[2,plotIndex[plantIndex[i]]]+ranEffYield_plot[2,plantIndex[i]]); //Slope + random effects
	}	
	
}
	
model {	
  
  // //Debugging
  // print(" logYield=",logYield);
  // print(" logYieldMu=",logYieldMu);
  // print(" sigmaYield=",sigmaYield);
  
	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	// plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	// flDens ~ normal(flDensMu,sigmaFlDens); //Flower density per plot
	// hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate	
	flwCount ~ neg_binomial_2_log(flwCountMu,phiFlwCount); //Flower count per plant (attempted pods)
	// podCount ~ beta_binomial(flwCount,inv_logit(flwSurv).*flwSurvPhi,(1-inv_logit(flwSurv)).*flwSurvPhi); //Flower survival (surviving pods)
	// avgSeedCount ~ exp_mod_normal(seedCountMu,sigmaSeedCount,lambdaSeedCount); //Average seeds per pod
	// avgSeedMass ~ exp_mod_normal(seedWeightMu,sigmaSeedWeight,lambdaSeedWeight); //Average weight per seed (mg) - exp-normal version, works much better	
	// logYield ~ normal(logYieldMu,sigmaYield); //Seed yield per plant
		
	// Priors
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
	
	// // Plant size - informative priors
	// intPlSize ~ normal(0,0.5); //Intercept
	// slopePlDensPlSize ~ normal(0,0.5); //Plant density
	// slopeDistPlSize ~ normal(0,0.05); //Distance effect
	// slopeGpPlSize ~ normal(0,0.5); //Grand Prairie effect
	// slopeIrrigPlSize ~ normal(0,0.5); //Irrigation effect
	// slope2015PlSize ~ normal(0.3,0.5); //2015 effect			
	// // slopeStockingPlSize ~ normal(0,0.1); //Stocking effect	
	// // slopePlDensStockingPlSize ~ normal(0,0.1); //Density:Stocking interaction			
	// sigmaPlSize_field ~ gamma(3,10); //Sigma for random field 
	// sigmaPlSize_plot ~ gamma(3,10); //Sigma for random plot
	// sigmaPlSize ~ gamma(7,10); //Sigma for residual	
	// intPlSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	// intPlSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot		
	// 
	// // Flower density per plot
	// intFlDens ~ normal(0,1); //Global intercept
	// slopePlSizeFlDens ~ normal(2,1); //plant size effect
	// slopeHbeeDistFlDens ~ normal(0,1); //distance into field
	// sigmaFlDens ~ gamma(7,2); //Sigma for within-field (residual)	
	// sigmaFlDens_field ~ gamma(4,2); //Sigma for field
	// intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field			
	// 
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
	sigmaPolField ~ gamma(1.25,3); //Sigma for random field	
	pollenPhi ~ gamma(1.25,3); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	// sigmaPolPlot ~ gamma(1.05,1); //Sigma for random plot - bad Rhat, poor traces   
	// intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int - not a lot of info at plot level
	
  //Flower count (per plant)
  intFlwCount ~ normal(5,0.2); //Intercept
  slopePlSizeFlwCount ~ normal(1,0.1); //Slope of plant size
  slopeSurvFlwCount ~ normal(-0.16,0.05); //Slope of survival rate
  slope2015FlwCount ~ normal(-0.1,0.05); //Slope of 2015 effect
  phiFlwCount_field ~ gamma(1.25,3); //SD of field-level random effect
  intFlwCount_field ~ normal(0,phiFlwCount_field); //Field-level random intercept
  phiFlwCount_plot ~ gamma(1.25,3); //SD of plot-level random effect
  intFlwCount_plot ~ normal(0,phiFlwCount_plot); //Plot-level random intercept
  //Variance (sigma) terms
  intPhiFlwCount ~ normal(3.7,0.5); //Intercept
  slopePlSizePhiFlwCount ~ normal(0.7,0.2);
  sigmaPhiFlwCount_field ~ gamma(1.25,3); //Sigma for field-level sigma
  intPhiFlwCount_field ~ normal(0,sigmaPhiFlwCount_field); //Field-level random intercept
  // sigmaPhiFlwCount_plot ~ gamma(1.25,3); //Sigma for plot-level sigma
  // intPhiFlwCount_plot ~ normal(0, sigmaPhiFlwCount_plot); //Plot-level random intercept
  
  //Flower survival - informative priors
  intFlwSurv ~ normal(1,1); //Intercept
  slopeVisitSurv ~ normal(0,0.05); //Slope of hbee visits
  slopePolSurv ~ normal(0,0.5); //Slope of pollen deposition
  slopePlSizeSurv ~ normal(0.02,0.05); //Slope of plant size
  slopePlDensSurv ~ normal(0,1); //Slope of planting density
  slopeIrrigSurv ~ normal(0,0.5); //Slope of irrigation
  slope2015Surv ~ normal(0,0.5); //Slope of year
  slopePlSizePlDensSurv ~ normal(0,1); //PlantSize:plant density
  sigmaFlwSurv_field ~ gamma(1.25,3); //SD of field-level random intercept
  intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
  // sigmaFlwSurv_plot ~ gamma(1.25,3); //SD of plot-level random intercept
  // intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //Plot-level random intercepts
  //Variance (sigma) terms
  intPhiFlwSurv ~ normal(0,5); //Intercept
  slopePlSizePhiFlwSurv ~ normal(0,5); //Effect of plant size on phi
  sigmaPhiFlwSurv_field ~ gamma(1.25,3); //Sigma for field level sigma
  intPhiFlwSurv_field ~ normal(0,sigmaPhiFlwSurv_field); //Field-level random intercept
  // sigmaPhiFlwSurv_plot ~ gamma(1,1); //Sigma for plot level sigma
  // intPhiFlwSurv_plot ~ normal(0,sigmaPhiFlwSurv_plot); //Plot-level random intercept
  
  
  // // Average seed count - informative priors
	// intSeedCount ~ normal(3.5,1); //Intercept
	// slopeVisitSeedCount ~ normal(0,0.1); //Slope of hbee visits
	// slopePolSeedCount ~ normal(0,0.5); //Slope of pollen deposition+
	// slopePlSizeCount ~ normal(0,0.05); //Slope of plant size
	// slope2015SeedCount ~ normal(0,0.5); //Year effect
	// // slopeIrrigSeedCount ~ normal(0,0.5); //Irrigation
	// sigmaSeedCount_field ~ gamma(2,10); //SD of field random effect
	// sigmaSeedCount_plot ~ gamma(2,10); //SD of plot random effect
	// intSeedCount_field ~ normal(0,sigmaSeedCount_field); //field-level random intercepts	
	// intSeedCount_plot ~ normal(0,sigmaSeedCount_plot); //plot-level random intercepts	
	// sigmaSeedCount ~ gamma(1,1); //SD of seed count
	// lambdaSeedCount ~ gamma(1,1); //Lambda for exp-normal distribution
	// 	
	// // Average weight per seed - informative priors
	// intSeedWeight ~ normal(1.5,1); //Intercept
	// slopeVisitSeedWeight ~ normal(0,0.5); //Slope of hbee visits
	// slopePolSeedWeight ~ normal(0,1); //Slope of pollen deposition
	// // slopeSeedCount ~ normal(0.015,0.015); //Slope of seed count
	// slopeSeedCount ~ normal(0,1); //Slope of (log) seed count
	// slopePlSizeWeight ~ normal(0,0.5); //Slope of plant size
	// slopeIrrigSeedWeight ~ normal(0,1); //Slope of irrigation
	// slope2015SeedWeight ~ normal(0,1); //Slope of 2015
	// sigmaSeedWeight ~ gamma(3,6); //SD of seed weight
	// sigmaSeedWeight_field ~ gamma(4,10); //SD of field random effect	
	// sigmaSeedWeight_plot ~ gamma(3,10); //SD of plot random effect
	// intSeedWeight_field ~ normal(0,sigmaSeedWeight_field); //field-level random intercepts	
	// intSeedWeight_plot ~ normal(0,sigmaSeedWeight_plot); //plot-level random intercepts	
	// lambdaSeedWeight ~ gamma(15,10); //Lambda for exp-normal distribution
	// 
	// // Yield per plant
	// intYield ~ normal(0,1); //Intercept
	// slopeYield ~ normal(1,1); //Slope of calculated yield	
	// sigmaYield ~ gamma(2,5); //Sigma for yield			
	// // Correlated random effects:
	// sigmaYield_field[1] ~ gamma(2,10); //SD of field-level yield intercepts/slopes
	// sigmaYield_field[2] ~ gamma(1.5,25);
	// sigmaYield_plot[1] ~ gamma(6,10); //SD of plot-level yield intercepts/slopes
	// sigmaYield_plot[2] ~ gamma(3,10);
	// to_vector(zYield_field) ~ normal(0,1); //Unit normals for correlated random effects
	// to_vector(zYield_plot) ~ normal(0,1);
	// L_field ~ lkj_corr_cholesky(2); //Standard prior for lkj cholesky matrix - higher values make extreme correlations less likely (see p.394 in McElreath 2016)	
	// L_plot ~ lkj_corr_cholesky(2); 
}

// generated quantities {
// 	//Plot-level quantities
// 	// planting density
// 	real predPlDens[Nplot]; //Generated
// 	real plDens_resid[Nplot]; //Residual
// 	real predPlDens_resid[Nplot]; //Residual of generated
// 	// flower density
// 	real predFlDens[Nplot]; 
// 	real flDens_resid[Nplot]; 
// 	real predFlDens_resid[Nplot]; 
// 	// hbeeVis
// 	int predHbeeVis[Nplot]; 
// 	real hbeeVis_resid[Nplot]; 
// 	real predHbeeVis_resid[Nplot]; 
// 
// 	// Flower-level
// 	// pollen deposition
// 	int predPollenCount[Nflw]; 
// 	real pollen_resid[Nflw]; 
// 	real predPollen_resid[Nflw]; 
// 
// 	// Plant-level
// 	// plantSize
// 	real predPlSize[Nplant]; 
// 	real plSize_resid[Nplant]; 
// 	real predPlSize_resid[Nplant]; 
// 	// flower count per plant (potential pods)
// 	real predFlwCount[Nplant]; 
// 	real flwCount_resid[Nplant]; 
// 	real predFlwCount_resid[Nplant];
// 	// flower survival (surviving pods)
// 	int predPodCount[Nplant]; 
// 	real podCount_resid[Nplant]; 
// 	real predPodCount_resid[Nplant]; 
//   // Average seeds per pod
//   real predSeedCount[Nplant]; 
//   real seedCount_resid[Nplant];
//   real predSeedCount_resid[Nplant];
//   // Average seed weight
//   real predSeedWeight[Nplant]; 
//   real seedWeight_resid[Nplant];
//   real predSeedWeight_resid[Nplant];
// 	// (log) yield per plant
// 	real predYield[Nplant]; 
// 	real yield_resid[Nplant];
// 	real predYield_resid[Nplant];
// 	
// 	
// 	for(i in 1:Nplot){
// 		// plant density
// 		plDens_resid[i] = plDens[i] - plDensMu[i]; //Residual for actual value
// 		predPlDens[i]= normal_rng(plDensMu[i],sigmaPlDens); //Generated value from normal
// 		predPlDens_resid[i] = predPlDens[i] - plDensMu[i]; //Residual for predicted value
// 
// 		// flower density
// 		flDens_resid[i] = flDens[i]-flDensMu[i]; //Residual for actual value
// 		predFlDens[i] = normal_rng(flDensMu[i],sigmaFlDens); //Generated value from normal
// 		predFlDens_resid[i] = predFlDens[i] - flDensMu[i]; //Residual for predicted value
// 
// 		// bee visits (NB version)
// 		hbeeVis_resid[i]=hbeeVis[i]-exp(visitHbeeMu[i]); //Residual for actual value
// 		predHbeeVis[i] = neg_binomial_2_log_rng(visitHbeeMu[i],visitHbeePhi); //Predicted value drawn from neg.bin
// 		predHbeeVis_resid[i]=predHbeeVis[i]-exp(visitHbeeMu[i]); //residual for predicted value
// 	}
// 
// 	for(i in 1:Nflw){
// 		// pollen deposition
// 		pollen_resid[i]= pollenCount[i] - exp(pollenMu[i]); //Residual for actual value
// 		predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi); //Simulate pollen counts
// 		predPollen_resid[i] = predPollenCount[i] - exp(pollenMu[i]); //Residual for predicted
// 	}
// 
// 	for(i in 1:Nplant){
// 		//plant size
// 		plSize_resid[i]= plantSize[i] - plSizeMu[i]; //Residual for actual
// 		predPlSize[i] = normal_rng(plSizeMu[i],sigmaPlSize); //Generates new value from normal dist.
// 		predPlSize_resid[i] = predPlSize[i] - plSizeMu[i]; //Residual for new value
// 
// 		// flower count per plant
// 		flwCount_resid[i] = flwCount[i] - exp(flwCountMu[i]); //Residual for actual
// 		predFlwCount[i] = neg_binomial_2_log_rng(flwCountMu[i],phiFlwCount[i]); //Generates new value from neg. bin.
// 		predFlwCount_resid[i] = predFlwCount[i] - exp(flwCountMu[i]); //Residual for new value
// 
// 		// pod count (surviving pods)
// 		podCount_resid[i] = podCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for actual
// 		predPodCount[i] = beta_binomial_rng(flwCount[i],inv_logit(flwSurv[i])*flwSurvPhi[i],(1-inv_logit(flwSurv[i]))*flwSurvPhi[i]); //Generates new value from beta-binomial
// 		predPodCount_resid[i] = predPodCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for new value
// 
// 
// 		// Seed count per pod - doesn't work well due to weird generating process (I think)
// 		seedCount_resid[i] = avgSeedCount[i] - (seedCountMu[i]+(1/lambdaSeedCount));
// 		predSeedCount[i] = exp_mod_normal_rng(seedCountMu[i],sigmaSeedCount,lambdaSeedCount);
// 		predSeedCount_resid[i] = predSeedCount[i] - (seedCountMu[i]+(1/lambdaSeedCount));
// 		
// 		// weight per seed - exp-normal works well
// 		seedWeight_resid[i] = avgSeedMass[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));
// 		predSeedWeight[i] = exp_mod_normal_rng(seedWeightMu[i],sigmaSeedWeight,lambdaSeedWeight);
// 		predSeedWeight_resid[i] = predSeedWeight[i] - (seedWeightMu[i]+(1/lambdaSeedWeight));
// 
// 		// (log) yield per plant
// 		yield_resid[i]= logYield[i] - logYieldMu[i]; //Residual for actual
// 		predYield[i] = normal_rng(logYieldMu[i],sigmaYield); //Generates new value from normal dist.
// 		predYield_resid[i] = predYield[i] - logYieldMu[i]; //Residual for new value
// 	}
// }
