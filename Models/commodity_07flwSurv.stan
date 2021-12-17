//Commodity canola pollination model used for publication
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
	// // Plant density
	// vector[Nplot_densMiss] plDens_miss;
	// real intPlDens; //Global intercept
	// real slope2015PlDens; //Effect of 2015
	// real slopeIrrigPlDens; //Effect of irrigation
	// real slope2015IrrigPlDens; //Year:irrigation interaction
	// real slopeGPPlDens; //GP effect on plant density
	// real slopeDistPlDens; //Slope of distance into field
	// real<lower=0> sigmaPlDens; //Sigma for within-field (residual)
	// real<lower=0> sigmaPlDens_field; //Sigma for field
	// vector[Nfield] intPlDens_field; //Random intercept for field

	// Pollen deposition
	// Plot level random effect has bad trace, strongly correlated with lp__
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits
	real slopeHbeeDistPollen; //Effect of distance into field - p=0.076 (p=0.17 if stocking and stocking:dist interaction included - see below)
	real<lower=0> sigmaPolField; //SD of field random intercepts
	vector[Nfield] intPollen_field; //field-level random intercepts
	real<lower=0> pollenPhi; //Dispersion parameter
	// // real<lower=0> sigmaPolPlot; //SD of plot random intercepts
	// // vector[Nplot] intPollen_plot; //plot-level random intercepts
	
	
  // Flower survival
	// slopeFlwCountSurv and slopePlSizeSurv are correlated (-0.77), and basically represent the same thing, so using only plant size
	real intFlwSurv; //Intercept
	real slopeVisitSurv; //Slope of hbee visits
	real slopePlSizeSurv; //Slope of plant size
	real slopeIrrigSurv; //Slope of irrigation - p=0.56
	real slope2015Surv; //Slope of year - p=0.68
	real slopePolSurv; //Slope of pollen deposition - requires other models
	// real slopePlDensSurv; //Slope of plant density - requires other models
	real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts - bad traces, high Rhat
	vector[Nfield] intFlwSurv_field; //field-level random intercepts - all overlap zero
	// real<lower=0> sigmaFlwSurv_plot; //SD of plot random intercepts - bad traces, high Rhat
	// vector[Nplot] intFlwSurv_plot; //plot-level random intercepts - all overlap zero
	// Sigma (variance) modeled as function of plant size
	real intPhiFlwSurv; //Intercept for sigma
	// real slopePlSizePhiFlwSurv; //Effect of plant size on phi
	// real<lower=0> sigmaPhiFlwSurv_field; //Sigma for field level sigma - bad traces, high Rhat
	// vector[Nfield] intPhiFlwSurv_field; //Field-level random intercept - all overlap zero
	// real<lower=0> sigmaPhiFlwSurv_plot; //Sigma for plot level - bad traces, high Rhat
	// vector[Nplot] intPhiFlwSurv_plot; //Plot-level random intercept - all overlap zero
	vector<lower=0, upper=1>[Nplant] flwSurvTheta;
}

transformed parameters {		
	//Expected values
	
	// //Plant density
	// vector[Nplot] plDensMu; //Expected plant density
	// vector[Nplot] plDens; //Planting density - imputed

	//Pollen deposition
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Expected pollen per stigma

	//Flower survival
  vector[Nplot] flwSurvPlot; //Plot-level flower survival
  vector[Nplant] flwSurv; //Flower survival rate (logit)
  vector<lower=0>[Nplant] flwSurvPhi; //Phi for flower survival
  vector<lower=0>[Nplant] flwSurvAlpha;
  vector<lower=0>[Nplant] flwSurvBeta;

  // //Assign imputed data
  // plDens[obsPlDens_ind]=plDens_obs; //Observed data
  // plDens[missPlDens_ind]=plDens_miss;	//Missing data

	for(i in 1:Nplot){ 
		// // Plant density per plot
		// plDensMu[i] = intPlDens + //Intercept
		//   intPlDens_field[plotIndex[i]] + //Field level random intercept
		// 	slope2015PlDens*is2015[plotIndex[i]]+ //Year effect
		// 	slopeIrrigPlDens*isIrrigated[plotIndex[i]]+ //Irrigation effect
		// 	slope2015IrrigPlDens*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]+ //Year:irrigation interaction
		// 	slopeDistPlDens*logHbeeDist[i] + //Distance effect
		// 	slopeGPPlDens*isGP[plotIndex[i]]; //Location effect

		// Plot-level pollen deposition
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //Field-level random intercept
		  // intPollen_plot[i] + //Plot-level random intercept
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect
		// Global intercept is within flower-level term in order to center plot-level variable
		
    // Plot-level flower survival
    flwSurvPlot[i] = intFlwSurv + //Intercept
      intFlwSurv_field[plotIndex[i]] + //Field-level random intercept
      // intFlwSurv_plot[i] + //Plot-level random intercept
    	slopeVisitSurv*logHbeeVis[i] + //hbee visits
    	// slopePolSurv*pollenPlot[i] + //(log) pollen deposition - large correlation b/w slopePolSurv and intFlwSurv
    	// slopePlDensSurv*plDens[i] + //Plant density
    	slopeIrrigSurv*isIrrigated[plotIndex[i]] + //Irrigation effect
    	slope2015Surv*is2015[plotIndex[i]]; //Year effect
	}
		
	for(i in 1:Nflw){ //For each flower stigma
		pollenMu[i] = intPollen + //Intercept
		  pollenPlot[flowerIndex[i]]; //Plot-level pollen
	}
		
	for(i in 1:Nplant){ //For each plant
    // Flower survival per plant
    flwSurv[i] = flwSurvPlot[plantIndex[i]] + //Plot-level plant survival
    	slopePlSizeSurv*plantSize[i]; //Plant size effect
    //Phi (dispersion) for flower survival
    flwSurvPhi[i] = exp(intPhiFlwSurv); //Intercept
      // intPhiFlwSurv_field[plotIndex[plantIndex[i]]] + //Field-level random intercept
      // intPhiFlwSurv_plot[plantIndex[i]] + //Plot-level random intercept
      // slopePlSizePhiFlwSurv*plantSize[i]; //Plant size
    flwSurvAlpha[i] = inv_logit(flwSurv[i])*flwSurvPhi[i]; //Outside of loop, used dot-multiply for vector multiplication
    flwSurvBeta[i] = (1-inv_logit(flwSurv[i]))*flwSurvPhi[i];
    }
}
	
model {	
  
  // print("plDensMu: ",plDensMu);
  // print("sigmaPlDens: ",sigmaPlDens);
  // print("pollenMu: ",pollenMu);
  // print("pollenPhi:",pollenPhi);
  // print("flwSurvTheta:",flwSurvTheta);
  // print("flwSurvAlpha:",flwSurvAlpha);
  // print("flwSurvBeta:",flwSurvBeta);
  
	//Likelihood		
	// plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate
	podCount ~ binomial(flwCount,flwSurvTheta); //Flower survival (surviving pods)
	flwSurvTheta ~ beta(flwSurvAlpha,flwSurvBeta); //For beta-binomial distribution

	// Priors
	// //Plant density	- informative priors
	// intPlDens ~ normal(0,1); //Global intercept
	// slope2015PlDens ~ normal(0,1); //Year effect
	// slopeIrrigPlDens ~ normal(0,1); //Irrigation effect
	// slope2015IrrigPlDens ~ normal(0,1); //Year:irrigation interaction
	// slopeDistPlDens ~ normal(0,1); //Slope of distance into field
	// slopeGPPlDens ~ normal(0,1); // Grand Prairie effect
	// sigmaPlDens ~ gamma(1,1); //Sigma for within-field (residual)
	// sigmaPlDens_field ~ gamma(1,1); //Sigma for field
	// intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field

	// // Pollen deposition - informative priors
	intPollen ~ normal(0,1); //Intercept
	slopeVisitPol ~ normal(0,1); //hbee Visitation effect
	slopeHbeeDistPollen ~ normal(0,1); //hbee distance effect
	pollenPhi ~ gamma(1,1); //Dispersion parameter
	sigmaPolField ~ gamma(1,1); //Sigma for random field
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	// sigmaPolPlot ~ gamma(1,1); //Sigma for random plot - bad Rhat, poor traces
	// // intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int - not a lot of info at plot level
	
  //Flower survival - informative priors
  intFlwSurv ~ normal(0,1); //Intercept
  slopeVisitSurv ~ normal(0,1); //Slope of hbee visits
  slopePolSurv ~ normal(0,1); //Slope of pollen deposition
  slopePlSizeSurv ~ normal(0,1); //Slope of plant size
  // slopePlDensSurv ~ normal(0,1); //Slope of planting density
  slopeIrrigSurv ~ normal(0,1); //Slope of irrigation
  slope2015Surv ~ normal(0,1); //Slope of year
  sigmaFlwSurv_field ~ gamma(1,1); //SD of field-level random intercept
  intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts
  // sigmaFlwSurv_plot ~ gamma(1,1); //SD of plot-level random intercept
  // intFlwSurv_plot ~ normal(0,sigmaFlwSurv_plot); //Plot-level random intercepts
  //Variance (sigma) terms
  intPhiFlwSurv ~ normal(0,5); //Intercept
  // slopePlSizePhiFlwSurv ~ normal(0,0.5); //Effect of plant size on phi
  // sigmaPhiFlwSurv_field ~ gamma(1.25,3); //Sigma for field level sigma
  // intPhiFlwSurv_field ~ normal(0,sigmaPhiFlwSurv_field); //Field-level random intercept
  // sigmaPhiFlwSurv_plot ~ gamma(1,1); //Sigma for plot level sigma
  // intPhiFlwSurv_plot ~ normal(0,sigmaPhiFlwSurv_plot); //Plot-level random intercept
}

generated quantities {
	// flower survival (surviving pods)
	int predPodCount[Nplant];
	real podCount_resid[Nplant];
	real predPodCount_resid[Nplant];

	for(i in 1:Nplant){
		// pod count (surviving pods)
		podCount_resid[i] = podCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for actual
		predPodCount[i] = beta_binomial_rng(flwCount[i],flwSurvAlpha[i],flwSurvBeta[i]); //Generates new value from beta-binomial
		predPodCount_resid[i] = predPodCount[i] - (flwCount[i]*inv_logit(flwSurv[i])); //Residual for new value
	}
}