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
	vector[Nplant] seedCount; //Average seeds per pod
	vector[Nplant] seedMass; //Average weight per seed (mg)
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
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+1); 
		logFlyVis[i] = log((flyVis[i]/totalTime[i])+1);
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
		calcYield[i] = podCount[i]*seedCount[i]*(seedMass[i]/1000);
		logCalcYield[i] = log(calcYield[i]);
	}		
	//Logit transform and center surviving flowers
	logitFlwSurv=logit(logitFlwSurv)-mean(logit(logitFlwSurv));
	
}

parameters {
	
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
	//Yield
	vector[Nplant] logYieldMu; //Predicted log(yield) per plant (calculated x coef)
	// Generate correlated random slopes and intercepts matrix
	matrix[2,Nfield] ranEffYield_field = diag_pre_multiply(sigmaYield_field,L_field) * zYield_field; //Field level random effects
	matrix[2,Nplot] ranEffYield_plot = diag_pre_multiply(sigmaYield_plot,L_plot) * zYield_plot; //Plot level random effects
	
	for(i in 1:Nplant){ //For each plant 	
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
  
	//Likelihood		
	logYield ~ normal(logYieldMu,sigmaYield); //Seed yield per plant
		
	// Priors
	// Yield per plant
	intYield ~ normal(0,1); //Intercept
	slopeYield ~ normal(0,1); //Slope of calculated yield
	sigmaYield ~ gamma(1,1); //Sigma for yield
	// Correlated random effects:
	sigmaYield_field[1] ~ gamma(1,1); //SD of field-level yield intercepts/slopes
	sigmaYield_field[2] ~ gamma(1,1);
	sigmaYield_plot[1] ~ gamma(1,1); //SD of plot-level yield intercepts/slopes
	sigmaYield_plot[2] ~ gamma(1,1);
	to_vector(zYield_field) ~ normal(0,1); //Unit normals for correlated random effects
	to_vector(zYield_plot) ~ normal(0,1);
	L_field ~ lkj_corr_cholesky(2); //Standard prior for lkj cholesky matrix - higher values make extreme correlations less likely (see p.394 in McElreath 2016)
	L_plot ~ lkj_corr_cholesky(2);
}

generated quantities {
  
  real rhoField = multiply_lower_tri_self_transpose(L_field)[1,2]; //Correlation for field intercepts/slopes 
  real rhoPlot = multiply_lower_tri_self_transpose(L_plot)[1,2]; //Correlation for plot intercepts/slopes
  
	// (log) yield per plant
	real predYield[Nplant];
	real yield_resid[Nplant];
	real predYield_resid[Nplant];


	for(i in 1:Nplant){
		// (log) yield per plant
		yield_resid[i]= logYield[i] - logYieldMu[i]; //Residual for actual
		predYield[i] = normal_rng(logYieldMu[i],sigmaYield); //Generates new value from normal dist.
		predYield_resid[i] = predYield[i] - logYieldMu[i]; //Residual for new value
	}
}
