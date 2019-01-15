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
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	vector[Nplot] logFlyVis; //Fly visitation rate
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)	
	vector[Nplant] logFlwCount; //Log flower count
	vector[Nplant] logitFlwSurv; //(logit) proportion flower survival	
	vector[Npod] logSeedCount; //Log seed count
	
	
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
	}
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
	
	for(i in 1:Npod){
		logSeedCount[i] = log(seedCount[i]); //Log-transforms seed count per pod
	}	
}

parameters {
	//Claim:flwCount~avgSurv+...
	real slopeSurvFlwCount;

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
			
	// Pollen deposition
	// Plot level random effect has bad trace, strongly correlated with lp__
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits	
	real slopeHbeeDistPollen; //Effect of distance into field - p=0.076 (p=0.17 if stocking and stocking:dist interaction included - see below)
	real<lower=0> sigmaPolField; //SD of field random intercepts	
	real<lower=0> pollenPhi; //Dispersion parameter
	vector[Nfield] intPollen_field; //field-level random intercepts	
	real<lower=0> sigmaPolPlot; //SD of plot random intercepts - Rhat 1.4, small n_eff, strongly correlated with lp__
	vector[Nplot] intPollen_plot; //plot-level random intercepts	
	
	// Flower count (per plant) 
	real intFlwCount; //Intercept
	real slopePlSizeFlwCount; //Slope of plant size	
	real<lower=0> phiFlwCount_field; //SD of field-level random effect	
	vector[Nfield] intFlwCount_field; //Field-level random effect	
	real intPhiFlwCount; //Intercept for sigma	
	real slopePlSizePhiFlwCount; //Effect of plant size on sigma	
	real<lower=0> sigmaPhiFlwCount_field; //Sigma for field level sigma
	vector[Nfield] intPhiFlwCount_field; //Field-level random effect for sigma	
	
	// Flower survival (pod production)
	// plot-level random intercepts have low n_eff, and basically all intercepts overlap zero; sigma plot correlated with lp__
	// slopeFlwCountSurv and slopePlSizeSurv are correlated (-0.77), and basically represent the same thing, so using only plant size
	real intFlwSurv; //Intercept
	real slopeVisitSurv; //Slope of hbee visits
	real slopePolSurv; //Slope of pollen deposition
	real slopePlSizeSurv; //Slope of plant size
	real slopePlDensSurv; //Slope of plant density		
	//Interactions	
	real slopePlSizePlDensSurv; //Plant size:plant density
	real<lower=0> sigmaFlwSurv_field; //SD of field random intercepts
	vector[Nfield] intFlwSurv_field; //field-level random intercepts
	real intPhiFlwSurv; //Intercept for sigma
	real slopePlSizePhiFlwSurv; //Effect of plant size on phi
	real<lower=0> sigmaPhiFlwSurv_field; //Sigma for field level sigma
	vector[Nfield] intPhiFlwSurv_field; //Field-level random effect for sigma
}

transformed parameters {		
}
	
model {	
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density	
	vector[Nplot] pollenPlot; //Plot-level pollen per stigma
	vector[Nflw] pollenMu; //Flower-level pollen per stigma
	vector[Nplot] flwCountPlot; //Plot-level flowers per plant
	vector[Nplant] flwCountMu; //Expected flowers per plant
	vector[Nplant] phiFlwCount; //Phi for flowers per plant	
	vector[Nplot] flwSurvPlot; //Plot-level flower survival
	vector[Nplant] flwSurv; //Flower survival rate (logit)
	vector[Nplant] flwSurvPhi; //Phi for flower survival	
	
	
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
		
		// Plot-level pollen deposition = random int field + random int plot + 
		pollenPlot[i] = intPollen_field[plotIndex[i]] + //intPollen_plot[i] + 
			slopeVisitPol*logHbeeVis[i] + //(log) hbee visits			
			slopeHbeeDistPollen*logHbeeDist[i]; //Distance effect								
		
		// Flower count per plant (plot level) = intercept + random field int + random plot int
		flwCountPlot[i] = intFlwCount + intFlwCount_field[plotIndex[i]];// + intFlwCount_plot[i];	
			
		// Plot-level flower survival = intercept + random int field + random int plot + 
		flwSurvPlot[i] = intFlwSurv + intFlwSurv_field[plotIndex[i]] + //intFlwSurv_plot[i] + 
			slopeVisitSurv*logHbeeVis[i] + //hbee visits 
			slopePolSurv*pollenPlot[i] + //(log) pollen deposition - large correlation b/w slopePolSurv and intFlwSurv
			slopePlDensSurv*plDens[i]; //Plant density									
	}
		
	for(i in 1:Nflw) //For each flower stigma
		pollenMu[i] = intPollen + pollenPlot[flowerIndex[i]]; //global intercept + plot-level pollen 
		
	for(i in 1:Nplant){ //For each plant 

	// Flower count per plant
		flwCountMu[i] = flwCountPlot[plantIndex[i]] + //Plot level flower count 
			slopePlSizeFlwCount*plantSize[i] + //individual plant size
			slopeSurvFlwCount*logitFlwSurv[i]; //Claim
			
		// Phi (dispersion) for flower count
		phiFlwCount[i] = exp(intPhiFlwCount + intPhiFlwCount_field[plotIndex[plantIndex[i]]] + //intPhiFlwCount_plot[plantIndex[i]] + 
			slopePlSizePhiFlwCount*plantSize[i]); 				
			
		// Flower survival per plant
		flwSurv[i] = flwSurvPlot[plantIndex[i]] + slopePlSizeSurv*plantSize[i] + //Plot-level plant survival + plant size effect		
			slopePlSizePlDensSurv*plantSize[i]*plDens[plantIndex[i]]; //Plant size:plant density 			
			
		//Phi (dispersion) for flower survival	
		flwSurvPhi[i] = exp(intPhiFlwSurv + intPhiFlwSurv_field[plotIndex[plantIndex[i]]] + slopePlSizePhiFlwSurv*plantSize[i]);		
	}	

	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density		
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi); //Pollination rate		
	flwCount ~ neg_binomial_2_log(flwCountMu,phiFlwCount); //Flower count per plant (attempted pods)
	podCount ~ beta_binomial(flwCount,inv_logit(flwSurv).*flwSurvPhi,(1-inv_logit(flwSurv)).*flwSurvPhi); //Flower survival (surviving pods)	
	
	//Priors
	//Claim
	slopeSurvFlwCount ~ normal(0,2);
	
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
			
	// Pollen deposition - informative priors	
	intPollen ~ normal(5.5,1); //Intercept	
	slopeVisitPol ~ normal(0,0.1); //hbee Visitation effect	
	slopeHbeeDistPollen ~ normal(0,0.1); //hbee distance effect	
	sigmaPolField ~ gamma(5,10); //Sigma for random field	
	pollenPhi ~ gamma(7,10); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int		
	sigmaPolPlot ~ gamma(3,10); //Sigma for random plot

	//Flower count (per plant)	
	intFlwCount ~ normal(5,1); //Intercept
	slopePlSizeFlwCount ~ normal(1,1); //Slope of plant size	
	phiFlwCount_field ~ gamma(2,10); //SD of field-level random effect		
	intFlwCount_field ~ normal(0,phiFlwCount_field); //Field-level random effect		
	intPhiFlwCount ~ normal(5,2); //Terms for variance
	slopePlSizePhiFlwCount ~ normal(1,1);
	sigmaPhiFlwCount_field ~ gamma(1,1); //Sigma for field level sigma
	intPhiFlwCount_field ~ normal(0,sigmaPhiFlwCount_field); //Field-level random effect for sigma			
	
	//Flower survival - informative priors
	intFlwSurv ~ normal(1,1); //Intercept
	slopeVisitSurv ~ normal(0,0.05); //Slope of hbee visits
	slopePolSurv ~ normal(0,0.5); //Slope of pollen deposition
	slopePlSizeSurv ~ normal(0.02,0.05); //Slope of plant size
	slopePlDensSurv ~ normal(0,1); //Slope of planting density		
	slopePlSizePlDensSurv ~ normal(0,1); //PlantSize:plant density
	sigmaFlwSurv_field ~ gamma(3,10); //SD of field random effect	
	intFlwSurv_field ~ normal(0,sigmaFlwSurv_field); //field-level random intercepts		
	intPhiFlwSurv ~ normal(0,5); //Intercept for sigma
	slopePlSizePhiFlwSurv ~ normal(0,5); //Effect of plant size on phi
	sigmaPhiFlwSurv_field ~ gamma(1,1); //Sigma for field level sigma
	intPhiFlwSurv_field ~ normal(0,sigmaPhiFlwSurv_field); //Field-level random effect for sigma
}

generated quantities {
}
