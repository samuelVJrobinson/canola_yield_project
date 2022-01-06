data { 
	//Field level
	int Nfield; //Number of fields		
	int Nfield_extra;	//Extra fields from Riley	
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 		
	int lbeeStocking[Nplot]; //Half leafcutter stocking? (double-tent treatment)
	matrix<lower=0,upper=1>[Nplot,2] lbeeStocking2; //Double tent, Double tent + Bees Treatments
	int is2016[Nplot]; //Was field from 2016?	
	vector[Nplot] hbee_dist; //Distance from hbee hives
	int hbeeVis[Nplot]; //Number of hbee visits to plot
	vector[Nplot] lbee_dist; //Distance to lbee shelters
	int lbeeVis[Nplot]; //Number of lbee visits to plot
	int isCent[Nplot]; //Is plot in center of bay?
	int isMBay[Nplot]; //Is plot in male bay?
	vector[Nplot] totalTime; //Minutes spent observing/10 	
	
	//Flower density (fls/m2)
	int Nplot_flDensObs; //Number of plots where flower density observed
	int Nplot_flDensMiss; //Number of plots where flower density missing (2 plots)
	vector[Nplot_flDensObs] flDens_obs; //Observed flower density (flowers/m2)
	int<lower=1,upper=Nplot> obsflDens_ind[Nplot_flDensObs]; //Index for observed flower density
	int<lower=1,upper=Nplot> missflDens_ind[Nplot_flDensMiss]; //Index for missing flower density	
	
	//Indices for mapping plot index onto female plot index
	int Nplot_F; //Number of female plots
	int<lower=0,upper=Nplot_F> plotIndex_F[Nplot]; //Index for female plots (which female plot j does plot i belong to? - some zero indices)
	int<lower=1,upper=Nplot> plotIndex_F2[Nplot_F]; //Reverse index (which plot i does female plot j belong to?)
	
	//Plant density (stems/m2)
	int Nplot_plDensObs; //Number of plots where plant density was observed
	int Nplot_plDensMiss; //Number of plots with plant density missing
	vector[Nplot_plDensObs] plDens_obs; //Observed plant density (stems/m2)
	int obsPlDens_ind[Nplot_plDensObs]; //Index for observed plant density
	int missPlDens_ind; //Index for missing plant density	(only 1)
	
	//Extra plots/fields from Riley (visitation data)
	int Nplot_extra; 
	int plotIndex_extra[Nplot_extra]; 	//Index for field (which field?)
	int is2016_extra[Nplot_extra]; //Was field from 2016?
	int lbeeStocking_extra[Nplot_extra]; //Half leafcutter stocking? (double-tent treatment)
	matrix<lower=0,upper=1>[Nplot_extra,2] lbeeStocking2_extra; //Double tent, Double tent + Bees Treatments
	vector[Nplot_extra] hbee_dist_extra; //Hbee distance
	int hbeeVis_extra[Nplot_extra]; //Hbee visitation
	vector[Nplot_extra] lbee_dist_extra; //Lbee distance
	int lbeeVis_extra[Nplot_extra]; //Lbee visitation
	int isCent_extra[Nplot_extra]; //Center plot
	int isMBay_extra[Nplot_extra]; //Male bay
	vector[Nplot_extra] totalTime_extra; //Time
	
	//Missing flower density data
	int Nplot_flDensObs_extra; //Number of plots where flower density observed
	int Nplot_flDensMiss_extra; //Number of plots where flower density missing
	vector[Nplot_flDensObs_extra] flDens_obs_extra; //Observed flower density (flowers/m2)
	int<lower=1,upper=Nplot_extra> obsflDens_ind_extra[Nplot_flDensObs_extra]; //Index for observed flower density
	int<lower=1,upper=Nplot_extra> missflDens_ind_extra[Nplot_flDensMiss_extra]; //Index for missing flower density			
	
	//Flower level (pollen counts)
	int Nflw; //Number of flowers (stigmas)
	int flowerIndex[Nflw]; //Index for flowers - which plot?	
	int pollenCount[Nflw]; //Pollen count
	
	//Plant level
	int Nplant; //Number of all plants 
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?		
	vector[Nplant] plantSize; //Mass of vegetative tissue (no seeds) (g)
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	vector[Nplant] logitFlwSurv; //Logit flower survival
	
	//Average seed counts per plant
  int Nplant_seedCountObs; //N observed 
  int Nplant_seedCountMiss; //N missing
  vector[Nplant_seedCountObs] seedCount_obs; //Observed seed weight
  int<lower=1,upper=Nplant> obsSeedCount_ind[Nplant_seedCountObs]; //Index for observed
  int<lower=1,upper=Nplant> missSeedCount_ind[Nplant_seedCountMiss]; //Index for missing
  
  //Avg weight per seed
  int Nplant_seedMassObs; //N observed
  int Nplant_seedMassMiss; //N missing
  vector[Nplant_seedMassObs] seedMass_obs; //Observed seed weight
  int<lower=1,upper=Nplant> obsSeedMass_ind[Nplant_seedMassObs]; //Index for observed
  int<lower=1,upper=Nplant> missSeedMass_ind[Nplant_seedMassMiss]; //Index for missing
	
	vector[Nplant] yield; //Observed yield (g seed per plant)		
	
}

transformed data {
	//Define values
	int Nfield_all = Nfield + Nfield_extra; //Total number of fields
	int Nplot_all = Nplot + Nplot_extra; //Total number of plots
	int plotIndex_all[Nplot_all]; //Plot index - which field is plot from?
	int lbeeStocking_all[Nplot_all]; //Half leafcutter stocking?
	matrix<lower=0,upper=1>[Nplot_all,2] lbeeStocking2_all; //Leafcutter stocking treatments
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
	// vector[Nplant] logFlwCount; //Log-flower count
	// vector[Nplant] logitFlwSurv; //(logit) proportion flower survival	
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)	
	
	//Assign values
	plotIndex_all[1:Nplot] = plotIndex; 
	plotIndex_all[(Nplot+1):Nplot_all] = plotIndex_extra;
	lbeeStocking_all[1:Nplot] = lbeeStocking; //Lbee stocking
	lbeeStocking_all[Nplot+1:Nplot_all] = lbeeStocking_extra;
	lbeeStocking2_all[1:Nplot,] = lbeeStocking2; //Leafcutter stocking treatments
	lbeeStocking2_all[(Nplot+1):Nplot_all,] = lbeeStocking2_extra;
	is2016_all[1:Nplot] = is2016; //Year
	is2016_all[Nplot+1:Nplot_all] = is2016_extra;	
	hbee_dist_all[1:Nplot] = hbee_dist; //Hbee distance
	hbee_dist_all[Nplot+1:Nplot_all] = hbee_dist_extra;
	hbeeVis_all[1:Nplot] = hbeeVis; //Hbee distance
	hbeeVis_all[Nplot+1:Nplot_all] = hbeeVis_extra;
	lbee_dist_all[1:Nplot] = lbee_dist; //Lbee distance
	lbee_dist_all[Nplot+1:Nplot_all] = lbee_dist_extra;
	lbeeVis_all[1:Nplot] = lbeeVis; //Lbee visitation
	lbeeVis_all[Nplot+1:Nplot_all] = lbeeVis_extra;
	isCent_all[1:Nplot] = isCent; //Center plot
	isCent_all[Nplot+1:Nplot_all] = isCent_extra;
	isMBay_all[1:Nplot] = isMBay; //M bay
	isMBay_all[Nplot+1:Nplot_all] = isMBay_extra;
	totalTime_all[1:Nplot] = totalTime; //Observation time
	totalTime_all[Nplot+1:Nplot_all] = totalTime_extra;
	
	//Transformations
	logHbeeDist_all=log(hbee_dist_all)-mean(log(hbee_dist_all)); //Log-transform and center distances
	logLbeeDist_all=log(lbee_dist_all)-mean(log(lbee_dist_all));
	logTime_all=log(totalTime_all); //Log-transform time	
	for(i in 1:Nplot_all){
		logHbeeVis_all[i] = log(1+(hbeeVis_all[i]/totalTime_all[i])); //Log-transforms observed visitation rates 
		logLbeeVis_all[i] = log(1+(lbeeVis_all[i]/totalTime_all[i]));			
	}
	
	// for(i in 1:Nplant){
	// 	logFlwCount[i] = log(flwCount[i]); //Log flower count per plant
	// 	//Necessary for promoting integers to reals. Otherwise does integer division.
	// 	logitFlwSurv[i] = podCount[i];
	// 	logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods
	// 	if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
	// 		logitFlwSurv[i]=0.01;
	// 	else if(logitFlwSurv[i]>=1)
	// 		logitFlwSurv[i]=0.99;
	// }
	// //Logit transform and center surviving flowers
	// logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
}

parameters {
 
	// Flower density per plot
	vector<lower=0,upper=100>[Nplot_flDensMiss] flDens_miss; //Missing from my fields
	vector<lower=0,upper=100>[Nplot_flDensMiss_extra] flDens_miss_extra; //Missing from Riley's fields
	// real intFlDens; //Global intercept
	// real slopeMBayFlDens; //Effect of male bay
	// real slope2016FlDens; //Effect of 2016
	// real slopeDistFlDens;  //Effect of distance from edge
	// real<lower=1e-10> sigmaFlDens; //Sigma for within-field (residual)
	// real<lower=1e-10> sigmaFlDens_field; //Sigma for field
	// vector[Nfield_all] intFlDens_field; //Random intercept for field
	// real nuFlDens; //exp(nu) for t-distribution

	// // Pollen deposition
	// // sigmaPolPlot correlated with lp__ (r=0.85), and most intercepts overlap zero
	// real intPol; //Intercept
	// real slopeHbeePol; //Slope of hbee visits
	// real slopeLbeePol; //Slope of lbee visits
	// real slopeCentPol; //Bay center effect
	// real slopeHbeeDistPol; //(log) hbee distance effect
	// real slopeFlDensPol; //Flower density
	// real<lower=1e-10> sigmaPolField; //Sigma for field-level intercept
	// vector[Nfield] intPol_field; //Field-level random intercept
	// real<lower=1e-10> sigmaPolPlot; //Sigma for plot-level intercept
	// vector[Nplot_F] intPol_plot; //Plot-level random intercept
	// real<lower=1e-05> pollenPhi; //Dispersion parameter
	
	// Flower survival per plant (pod count)
	real intFlwSurv; //Intercept
	real slopePolSurv; //Slope of pollen deposition
	real slopePlSizeSurv; //Slope of plant size
	real slopeEdgeCentSurv; //Slope of edge effect
	real slopeHbeeDistSurv; //Effect of (log) hbee distance
	real slopeLbeeDistSurv; //Effect of (log) lbee distance
	// real slopeFlwCountSurv; //Effect of flower count
	real slopeFlwDensSurv; //Effect of flower density
	// real slopeSeedSizeSurv; //Slope of seed size
	real<lower=1e-10> sigmaFlwSurv_plot; //SD of plot random intercepts
	real<lower=1e-10> sigmaFlwSurv_field; //SD of field random intercepts
	vector[Nfield] intFlwSurv_field; //field-level random intercepts
	vector[Nplot_F] intFlwSurv_plot; //plot-level random intercepts
	real intPhiFlwSurv; //Intercept for sigma - dispersion term for beta binomial
	// real slopePlSizePhiFlwSurv; //Effect of plant size on phi
	// real<lower=1e-05> sigmaPhiFlwSurv_field; //Sigma for field level sigma
	// vector[Nfield] intPhiFlwSurv_field; //Field-level random effect for sigma
}

transformed parameters {
			
	//Expected values
	//Plot-level
	// vector[Nplot_all] flDensMu; //Expected flower density	
	// vector[Nplot_F] pollenMu_plot; //Plot level pollen
	// vector[Nflw] pollenMu; //Expected pollen - flower level
	vector[Nplot_F] flwSurvPlot; //Plot-level flower production
	vector[Nplant] flwSurv; //Flower production (exp)
	vector<lower=0>[Nplant] flwSurvPhi; //Phi for flower survival
	
	//Imputed missing data;
	vector[Nplot_all] flDens; //Flower density
	
	//Combine observed with imputed		
	flDens[obsflDens_ind]=flDens_obs; //Observed flower density
	flDens[missflDens_ind]=flDens_miss;
 	for(i in 1:Nplot_flDensObs_extra) //For each extra observed plot
		flDens[obsflDens_ind_extra[i]+Nplot]=flDens_obs_extra[i];	//Add it to index in flDens
	for(i in 1:Nplot_flDensMiss_extra) //For each extra missing plot
		flDens[missflDens_ind_extra[i]+Nplot]=flDens_miss_extra[i];
	
	// //Plot-level parameters
	// for(i in 1:Nplot_all){	//Parameters for all fields, all plots
	// 	// Flower density = intercept + random field int + plant size effect
	// 	flDensMu[i] = intFlDens	+ intFlDens_field[plotIndex_all[i]] + 
	// 		slopeMBayFlDens*isMBay_all[i] + //Male bay effect
	// 		// slopePlSizeFlDens*plSizePlotMu[i] + //Plant size
	// 		slope2016FlDens*is2016_all[i] + //Year effect
	// 		slopeDistFlDens*logHbeeDist_all[i]; //Distance from edge effect
	// 		// slopePlDensFlDens*plDens[i]; //Planting density effect
	// }	
	
	for(i in 1:Nplot_F){ //Parameters for F plots only
	
	  int plotI = plotIndex_F2[i]; //Matches F plot i to measurements taken at all plots
	
// 		// Pollen per plot = intercept + random field int + random plot int + leafcutter effect + honeybee effect + bay center effect + hbee dist effect
//   	// Moved intPol to Nflw loop to center plot level data
// 	  pollenMu_plot[i] = intPol_field[plotIndex[plotI]] + //Field random intercept
// 	    intPol_plot[i] + //Plot random intercept
//    	  slopeLbeePol*logLbeeVis_all[plotI] +  //Effect of (log) leafcutter visits
//     	slopeHbeePol*logHbeeVis_all[plotI] +  //Effect of (log) honeybee visits
//     	slopeCentPol*isCent_all[plotI] + //Bay center effect
//     	slopeHbeeDistPol*logHbeeDist_all[plotI] + //(log) hbee distance effect
//     	// slopeStockingHbeeDistPol*logHbeeDist_all[plotI]*lbeeStocking_all[plotI] + //hbee dist:lbee stocking
//    	  slopeFlDensPol*flDens[plotI]; //Flower density
			
		// Flower survival = intercept + random int field + random int plot +
		flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[plotI]] + intFlwSurv_plot[i] +
  		slopeEdgeCentSurv*isCent_all[plotI] + //Bay center
  		// slopePolSurv*pollenMu_plot[i] + //Pollen deposition (plot level average)
  		slopeHbeeDistSurv*logHbeeDist_all[plotI] + //Distance from honeybees
  		slopeLbeeDistSurv*logLbeeDist_all[plotI] + //Distance from leafcutters
  		slopeFlwDensSurv*flDens[plotI]; //Effect of flower density
	}
				
	// for(i in 1:Nflw)
	//   pollenMu[i] = intPol + pollenMu_plot[plotIndex_F[flowerIndex[i]]]; //Assigns plot level pollen mu to Nflw long vector
		
	for(i in 1:Nplant){	
	  
	  int plotI = plotIndex_F[plantIndex[i]];
	
		// Predicted pod count (flower survival)
		flwSurv[i] = flwSurvPlot[plotI]; //Plot-level plant survival
		  // slopePlSizeSurv*plantSize[i];  //plant size effect
		//Phi (dispersion) for flower survival
		flwSurvPhi[i] = exp(intPhiFlwSurv);// + 
		  // intPhiFlwSurv_field[plotIndex[plotI]] +
		  // slopePlSizePhiFlwSurv*plantSize[i]);
	}
	
}
	
model {
	// flDens ~ student_t(exp(nuFlDens),flDensMu,sigmaFlDens); 	
	// pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate
	podCount ~ beta_binomial(flwCount,inv_logit(flwSurv).*flwSurvPhi,(1-inv_logit(flwSurv)).*flwSurvPhi); //Flower survival - betabinomial version
			
	// Priors
	// Flower density	
	// intFlDens ~ normal(0,1); //Intercept
	// slope2016FlDens ~ normal(0,1); //Effect of 2016
	// slopeDistFlDens ~ normal(0,1); //Distance from edge
	// slopeMBayFlDens ~ normal(0,1); //Effect of male bay
	// sigmaFlDens ~ gamma(1,1); //Sigma for plot (residual)
	// sigmaFlDens_field ~ gamma(1,1); ; //Sigma for field
	// intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field	
	// nuFlDens ~ normal(0,1); //nu for student's t

	// // Pollen deposition - informative priors
	// intPol ~ normal(0,1); //Intercept
	// slopeHbeePol ~ normal(0,1); //hbee Visitation effect
	// slopeLbeePol ~ normal(0,1); //lbee Visitation effect
	// slopeCentPol~ normal(0,1); //Bay center effect
	// slopeHbeeDistPol ~ normal(0,1); //(log) hbee distance effect
	// // slopeStockingHbeeDistPol ~ normal(0,1); //Stocking:hbee distance interaction
	// slopeFlDensPol ~ normal(0,1); //Flower density
	// sigmaPolField ~ gamma(1,1); //Sigma for random field
	// sigmaPolPlot ~ gamma(1,1); //Sigma for random plot
	// pollenPhi ~ gamma(1,1); //Dispersion parameter
	// intPol_field ~ normal(0,sigmaPolField); //Random field int
	// intPol_plot ~ normal(0,sigmaPolPlot); //Random plot int
			
	// Flower survival (pod count)
	intFlwSurv ~ normal(0,1); //Intercept
	slopePolSurv ~ normal(0,1); //Slope of pollen deposition
	slopePlSizeSurv ~ normal(0,1); //Slope of plant size
	slopeEdgeCentSurv ~ normal(0,1); //Slope of edge effect
	slopeHbeeDistSurv ~ normal(0,1); //Distance from edge
	slopeLbeeDistSurv ~ normal(0,1); //Distance from lbee shelter
	slopeFlwDensSurv ~ normal(0,1); //Flower density effect
	sigmaFlwSurv_field ~ gamma(1,1); //SD of field random effect
	sigmaFlwSurv_plot ~ gamma(1,1); //SD of plot random effect
	intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
	intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //plot-level random intercepts
	intPhiFlwSurv ~ normal(0,1); //Intercept for sigma
	// slopePlSizePhiFlwSurv ~ normal(0,1); //Effect of plant size on phi
	// sigmaPhiFlwSurv_field ~ gamma(1,1); //Sigma for field level sigma
	// intPhiFlwSurv_field ~ normal(0,sigmaPhiFlwSurv_field); //Field-level random effect for sigma
}

generated quantities {
	// flower survival (surviving pods)
	int<lower=0> predPodCount[Nplant]; //Generated
	real podCount_resid[Nplant]; //Residual
	real predPodCount_resid[Nplant]; //Residual of generated

	for(i in 1:Nplant){
		// pod count (surviving pods) - betabinom version
		podCount_resid[i] = podCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for actual
		predPodCount[i] = beta_binomial_rng(flwCount[i],inv_logit(flwSurv[i])*flwSurvPhi[i],(1-inv_logit(flwSurv[i]))*flwSurvPhi[i]); //Generates new value from beta-binomial
		predPodCount_resid[i] = predPodCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for new value
		
		// // pod count (surviving pods) - negbin version
		// podCount_resid[i] = podCount[i] - exp(flwSurv[i]); //Residual for actual
		// predPodCount[i] = neg_binomial_2_log_rng(flwSurv[i],flwSurvPhi[i]); //Generates new value from negbin
		// predPodCount_resid[i] = predPodCount[i] - exp(flwSurv[i]); //Residual for new value
	}	
}
