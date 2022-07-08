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
	vector<lower=5,upper=51>[Nplot_flDensMiss] flDens_miss; //Missing from my fields
	vector<lower=5,upper=51>[Nplot_flDensMiss_extra] flDens_miss_extra; //Missing from Riley's fields

	// Pollen deposition
	// sigmaPolPlot correlated with lp__ (r=0.85), and most intercepts overlap zero
	real<lower=-10,upper=10> intPollen; //Intercept
	real<lower=-10,upper=10> slopeHbeeVisPollen; //Slope of hbee visits
	real<lower=-10,upper=10> slopeLbeeVisPollen; //Slope of lbee visits
	real<lower=-10,upper=10> slopeCentPollen; //Bay center effect
	real<lower=-10,upper=10> slopeHbeeDistPollen; //(log) hbee distance effect
	real<lower=-10,upper=10> slopeLbeeDistPollen; //(log) lbee distance effect
	real<lower=-10,upper=10> slopeFlDensPollen; //Flower density
	// real slopeStockingHbeeDistPollen; //hbee distance:stocking interaction
	real<lower=1e-10,upper=10> sigmaPollen_field; //Sigma for field-level intercept
	vector<lower=-10,upper=10>[Nfield] intPollen_field; //Field-level random intercept
	real<lower=1e-10,upper=10> sigmaPollen_plot; //Sigma for plot-level intercept
	vector<lower=-10,upper=10>[Nplot_F] intPollen_plot; //Plot-level random intercept
	real<lower=1e-05,upper=10> phiPollen; //Dispersion parameter
	
}

transformed parameters {
			
	//Expected values
	//Plot-level
	vector[Nplot_F] pollenMu_plot; //Plot level pollen
	vector[Nflw] pollenMu; //Expected pollen - flower level

	//Imputed missing data;
	vector[Nplot_all] flDens; //Flower density
	
	//Combine observed with imputed		
	flDens[obsflDens_ind]=flDens_obs; //Observed flower density
	flDens[missflDens_ind]=flDens_miss;
 	for(i in 1:Nplot_flDensObs_extra) //For each extra observed plot
		flDens[obsflDens_ind_extra[i]+Nplot]=flDens_obs_extra[i];	//Add it to index in flDens
	for(i in 1:Nplot_flDensMiss_extra) //For each extra missing plot
		flDens[missflDens_ind_extra[i]+Nplot]=flDens_miss_extra[i];
		
	for(i in 1:Nplot_F){ //Parameters for F plots only
	  // Pollen per plot = intercept + random field int + random plot int + leafcutter effect + honeybee effect + bay center effect + hbee dist effect
  	// Moved intPol to Nflw loop to center plot level data
	  pollenMu_plot[i] = intPollen_field[plotIndex[plotIndex_F2[i]]] + intPollen_plot[i] + //Intercept + field/plot level random effects
   	  slopeLbeeVisPollen*logLbeeVis_all[plotIndex_F2[i]] +  //Effect of (log) leafcutter visits
    	slopeHbeeVisPollen*logHbeeVis_all[plotIndex_F2[i]] +  //Effect of (log) honeybee visits
    	slopeCentPollen*isCent_all[plotIndex_F2[i]] + //Bay center effect
    	slopeHbeeDistPollen*logHbeeDist_all[plotIndex_F2[i]] + //(log) hbee distance effect
    	slopeLbeeDistPollen*logLbeeDist_all[plotIndex_F2[i]] + //(log) lbee distance effect
   	  slopeFlDensPollen*flDens[plotIndex_F2[i]]; //Flower density
	}
				
	for(i in 1:Nflw)
	  pollenMu[i] = intPollen + pollenMu_plot[plotIndex_F[flowerIndex[i]]]; //Assigns plot level pollen mu to Nflw long vector
	
}
	
model {

	pollenCount ~ neg_binomial_2_log(pollenMu,phiPollen); //Pollination rate
			
	// Priors
	
	// Pollen deposition - informative priors
	intPollen ~ normal(3.1,5); //Intercept
	slopeHbeeVisPollen ~ normal(0,5); //hbee Visitation effect
	slopeLbeeVisPollen ~ normal(0,5); //lbee Visitation effect
	slopeCentPollen~ normal(0,5); //Bay center effect
	slopeHbeeDistPollen ~ normal(0,5); //(log) hbee distance effect
	slopeLbeeDistPollen ~ normal(0,5); //(log) lbee distance effect
	// slopeStockingHbeeDistPollen ~ normal(0,5); //Stocking:hbee distance interaction
	slopeFlDensPollen ~ normal(0,5); //Flower density
	sigmaPollen_field ~ gamma(1,1); //Sigma for random field
	sigmaPollen_plot ~ gamma(1,1); //Sigma for random plot
	phiPollen ~ gamma(1,1); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPollen_field); //Random field int
	intPollen_plot ~ normal(0,sigmaPollen_plot); //Random plot int
}

generated quantities {
	// Flower-level
	// pollen deposition
	int predPollenCount[Nflw]; //Generated
	real pollen_resid[Nflw]; //residual
	real predPollen_resid[Nflw]; //residual of generated

	for(i in 1:Nflw){
  	//pollen deposition
  	pollen_resid[i]= exp(pollenMu[i]) - pollenCount[i]; //Residual for actual value
  	predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],phiPollen); //Simulate pollen counts
  	predPollen_resid[i] = exp(pollenMu[i]) - predPollenCount[i]; //Residual for predicted
	}
}
