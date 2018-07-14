data {
	//Field level
	int Nfield; //Number of fields
	int<lower=1,upper=Nfield> fieldIndex[Nfield]; //Index for fields
		
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
		
	//Plant level
	int Nplant; //Number of all plants (some measurements missing)
	int Nplant_obs; //Number of "complete" plants
	int Nplant_miss; //Number of missing plants	
	vector[Nplant_obs] plantSize_obs; //Mass of vegetative tissue (no seeds) (g)	
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?	
	int<lower=1,upper=Nplot> plantSurvIndex[Nplant_obs]; //Index for "complete" plants 	
	int<lower=1,upper=Nplant> obs_ind[Nplant_obs]; //Index for observed plants
	int<lower=1,upper=Nplant> miss_ind[Nplant_miss]; //Index for missing plants
}

parameters {
	vector[Nplant_miss] plantSize_miss; //Imputed data for missing values	
	real intSize; //Global intercept
	real<lower=0> sigmaField; //Sigma for field
	real<lower=0> sigmaPlot; //Sigma for plot
	real<lower=0> sigmaSize; //Sigma for within-plot	(residual)
	vector[Nfield] intSize_field; //Random intercept for field
	vector[Nplot] intSize_plot; //Random intercept for plot	
}

transformed parameters {	
}
	
model {		
	//Expected plant size values	
	vector[Nplant] sizeMu; //Predicted plant-level size
	
	//Imputed missing data;
	vector[Nplant] plantSize;
	plantSize[obs_ind]=plantSize_obs;
	plantSize[miss_ind]=plantSize_miss;
			
	for(i in 1:Nplant)
		sizeMu[i] = intSize + intSize_field[plotIndex[plantIndex[i]]] + intSize_plot[plantIndex[i]]; 

	//Likelihood
	plantSize ~ normal(sizeMu,sigmaSize);	//Plant size
		
	//Priors
	//Plant size
	intSize ~ normal(0,1); //Intercept
	sigmaField ~ gamma(2,6); //Sigma for random field 
	sigmaPlot ~ gamma(2,7); //Sigma for random plot
	sigmaSize ~ gamma(3,5); //Sigma for residual	
	intSize_field ~ normal(0,sigmaField); //Random field int
	intSize_plot ~ normal(0,sigmaPlot); //Random int plot		
}

generated quantities{
}
