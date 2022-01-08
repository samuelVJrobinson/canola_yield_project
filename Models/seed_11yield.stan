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
	vector[Nplant] calcYield; //Calculated yield (# pods x seed count x seed mass)
	
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
	vector[Nplant] lCalcYield; //Log calculated yield (copy)
	
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
	
	for(i in 1:Nplant){
	// 	logFlwCount[i] = log(flwCount[i]); //Log flower count per plant
	// 	//Necessary for promoting integers to reals. Otherwise does integer division.
	// 	logitFlwSurv[i] = podCount[i];
	// 	logitFlwSurv[i] = logitFlwSurv[i]/flwCount[i]; //Proportion surviving pods
	// 	if(logitFlwSurv[i]<=0) //Deal with weird 100% and 0% plants
	// 		logitFlwSurv[i]=0.01;
	// 	else if(logitFlwSurv[i]>=1)
	// 		logitFlwSurv[i]=0.99;
	if(calcYield[i]<0){
	  lCalcYield[i] = -999; //NA values represented by -999
	} else {
	  lCalcYield[i] = log(calcYield[i]);
	}
	
	}
	// //Logit transform and center surviving flowers
	// logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
}

parameters {
  
	vector<lower=1,upper=6>[Nplant_seedMassMiss] seedMass_miss; //Missing seed weights
	vector<lower=1,upper=31>[Nplant_seedCountMiss] seedCount_miss; //Missing seeds per pod

	// Total yield (g/plant)
	real intYield; //Intercept for predicted yield
	real slopeYield; //Proportion of predicted yield that becomes yield
	vector<lower=0>[2] sigmaYield_field; //SD of field-level intercept/slopes
	vector<lower=0>[2] sigmaYield_plot; //SD of plot-level intercept/slope
	real<lower=1e-05> sigmaYield; //SD of plant-level yield
	cholesky_factor_corr[2] L_field; //Cholesky-decomposed correlation matrices
	cholesky_factor_corr[2] L_plot;
	matrix[2,Nfield] zYield_field; //Unit normals for matrix correlation trick
	matrix[2,Nplot] zYield_plot;
}

transformed parameters {
			
	//Expected values
	vector[Nplant] logYieldMu; //Predicted log(yield) per plant (calculated x coef)
	vector[Nplant] logCalcYield; //Log-calculated yield
	// Generate correlated random slopes and intercepts matrix
	matrix[2,Nfield] ranEffYield_field = diag_pre_multiply(sigmaYield_field,L_field) * zYield_field; //Field level random effects
	matrix[2,Nplot] ranEffYield_plot = diag_pre_multiply(sigmaYield_plot,L_plot) * zYield_plot; //Plot level random effects
	
	//Imputed missing data;
	vector[Nplant] seedCount; //Seeds per pod
	vector[Nplant] seedMass; //Seed weight
	
	//Combine observed with imputed		
	seedCount[obsSeedCount_ind]	= seedCount_obs; //Observed seeds per pod
	seedCount[missSeedCount_ind]	= seedCount_miss; //Missing seeds per pod
	seedMass[obsSeedMass_ind]	= seedMass_obs; //Observed seed mass
	seedMass[missSeedMass_ind]	= seedMass_miss; //Missing seed mass
	
	//Plot-level parameters
		
	for(i in 1:Nplant){	
	  
	  int plotI = plotIndex_F[plantIndex[i]];
	  
	  //Calculated yield
	  if(lCalcYield[i] == -999){ //If NA
	    logCalcYield[i] = log(seedCount[i]*(seedMass[i]/1000)*podCount[i]);
	  } else {
	    logCalcYield[i] = lCalcYield[i];
	  }
  		
    // Predicted yield
		logYieldMu[i] = intYield + //Intercept
  		ranEffYield_field[1,plotIndex[plantIndex[i]]]+ //Field level random intercepts
  		ranEffYield_plot[1,plantIndex[i]] + //Plot level random intercepts
  		logCalcYield[i] * ( slopeYield + //Slope
  		ranEffYield_field[2,plotIndex[plantIndex[i]]] + //Field-level random slopes
  		ranEffYield_plot[2,plantIndex[i]]); //Plot-level random slopes
	}
	
}
	
model {
	
	logYield ~ normal(logYieldMu,sigmaYield); //Seed yield per plant
			
	// Priors
	// Yield per plant
	intYield ~ normal(0,5); //Intercept
	slopeYield ~ normal(0,5); //Slope of calculated yield
	sigmaYield ~ gamma(1,1); //Sigma for yield
	// Correlated random effects:
	sigmaYield_field[1] ~ gamma(1,1); //SD of field-level yield intercepts/slopes
	sigmaYield_field[2] ~ gamma(1,1);
	sigmaYield_plot[1] ~ gamma(1,1); //SD of plot-level yield intercepts/slopes
	sigmaYield_plot[2] ~ gamma(1,1);
	to_vector(zYield_field) ~ normal(0,5); //Unit normals for correlated random effects
	to_vector(zYield_plot) ~ normal(0,5);
	L_field ~ lkj_corr_cholesky(2); //Standard prior for lkj cholesky matrix - higher values make extreme correlations less likely (see p.394 in McElreath 2016)
	L_plot ~ lkj_corr_cholesky(2);
}

generated quantities {
	// (log) yield per plant
	real predYield[Nplant]; //Generated
	real yield_resid[Nplant]; //Residual
	real predYield_resid[Nplant]; //Residual of generated
	
	real rhoField = multiply_lower_tri_self_transpose(L_field)[2,1];
	real rhoPlot = multiply_lower_tri_self_transpose(L_plot)[2,1];

	for(i in 1:Nplant){
		// (log) yield per plant
		yield_resid[i]= logYield[i] - logYieldMu[i]; //Residual for actual
		predYield[i] = normal_rng(logYieldMu[i],sigmaYield); //Generates new value from normal dist.
		predYield_resid[i] = predYield[i] - logYieldMu[i]; //Residual for new value
	}	
 
}
