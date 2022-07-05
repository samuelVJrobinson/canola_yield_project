parameters {
  
	// Plant size - random effects at plot level are very small, and don't converge well
	// Density:distance interaction is basically 0, so leaving it out
	// Normal distribution had marginal PPchecks; t-dist is much better
	real<lower=-10,upper=10> intPlSize; //Global intercept
	real<lower=-10,upper=10> slopePlDensPlSize; //Slope of planting density
	real<lower=-10,upper=10> slopeHbeeDistPlSize; //Slope of hbee distance (edge of field has small plants)
	real<lower=1e-05,upper=10> sigmaPlSize; //Sigma for within-plot (residual)
	real<lower=1e-05,upper=10> sigmaPlSize_field; //Sigma for field
	vector<lower=-10,upper=10>[Nfield] intPlSize_field; //Random intercept for field
	real<lower=-10,upper=10> nuPlSize; //exp(nu) for t-distribution

}

transformed parameters {
			
	//Expected values
	//Plot-level
	vector[Nplot_F] plSizePlotMu; //Plot-level plant size
	vector[Nplant] plSizeMu; //Expected plant size
	
	
	//Plot-level parameters
	for(i in 1:Nplot_F){	//Parameters for all fields, all plots
	
	  //Matches F plot i to measurements taken at all plots
	  // "plotI" used to index measurements from ALL plots, "i" used otherwise
	  int plotI = plotIndex_F2[i]; 
	  
	  // Plant size (plot-level) 
		plSizePlotMu[i] = intPlSize + //Intercept
		  intPlSize_field[plotIndex_all[plotI]] +  //Field level intercept
			slopeHbeeDistPlSize*logHbeeDist_all[plotI] + //Distance effect (edge of field has smaller plants)
			slopePlDensPlSize*plDens[i]; //Planting density effect
			
		CLAIM*plSizePlotMu[i]; 	
	}	
	
	for(i in 1:Nplant){
	  int plotI = plotIndex_F[plantIndex[i]]; //Matches plant to F plot
		// Predicted plant size (taken from plot level measurements above)
		plSizeMu[i] = plSizePlotMu[plotI];
	}
}
	
model {

	plantSize ~ student_t(exp(nuPlSize),plSizeMu,sigmaPlSize); //Plant size - T dist

	// Priors
	
	//Plant size - informative priors
	intPlSize ~ normal(3.2,5); //Intercept
	slopePlDensPlSize ~ normal(0,5); //Planting density
	slopeHbeeDistPlSize ~ normal(0,5); //Distance from edge of field
	sigmaPlSize ~ gamma(1,1); //Sigma for residual
	sigmaPlSize_field ~ gamma(1,1); //Sigma for field
	intPlSize_field	~ normal(0,sigmaPlSize_field); //Random intercept for field
	nuPlSize ~ normal(0,5); //nu for student's t
	
}