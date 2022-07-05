parameters {
  
 	// Flower density per plot
	vector<lower=5,upper=51>[Nplot_flDensMiss] flDens_miss; //Missing from my fields

	// Pollen deposition
	real<lower=-10,upper=10> intPollen; //Intercept
	real<lower=-10,upper=10> slopeHbeeVisPollen; //Slope of hbee visits
	real<lower=-10,upper=10> slopeLbeeVisPollen; //Slope of lbee visits
	real<lower=-10,upper=10> slopeCentPollen; //Bay center effect
	real<lower=-10,upper=10> slopeHbeeDistPollen; //(log) hbee distance effect
	real<lower=-10,upper=10> slopeFlDensPollen; //Flower density
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
	vector[Nplot] flDens; //Flower density
	
	//Combine observed with imputed		
	flDens[obsflDens_ind]=flDens_obs; //Observed flower density
	flDens[missflDens_ind]=flDens_miss;
	
	for(i in 1:Nplot_F){ //Parameters for F plots only
	
	  int plotI = plotIndex_F2[i]; //Matches F plot i to measurements taken at all plots
	
		// Pollen per plot
	  pollenMu_plot[i] = intPollen_field[plotIndex[plotI]] + //Field random intercept
	    intPollen_plot[i] + //Plot random intercept
   	  slopeLbeeVisPollen*logLbeeVis_all[plotI] +  //Effect of (log) leafcutter visits
    	slopeHbeeVisPollen*logHbeeVis_all[plotI] +  //Effect of (log) honeybee visits
    	slopeCentPollen*isCent_all[plotI] + //Bay center effect
    	slopeHbeeDistPollen*logHbeeDist_all[plotI] + //(log) hbee distance effect
   	  slopeFlDensPollen*flDens[plotI]; //Flower density
   	  
   	CLAIM*pollenMu_plot[i];
	}
				
	//Assigns plot level pollen mu to Nflw long vector			
	for(i in 1:Nflw)
	  pollenMu[i] = intPollen + pollenMu_plot[plotIndex_F[flowerIndex[i]]]; 
	
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
	slopeFlDensPollen ~ normal(0,5); //Flower density
	sigmaPollen_field ~ gamma(1,1); //Sigma for random field
	sigmaPollen_plot ~ gamma(1,1); //Sigma for random plot
	phiPollen ~ gamma(1,1); //Dispersion parameter
	intPollen_field ~ normal(0,sigmaPollen_field); //Random field int
	intPollen_plot ~ normal(0,sigmaPollen_plot); //Random plot int
}
