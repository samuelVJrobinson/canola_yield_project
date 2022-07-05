parameters {
  
  real<lower=2.4,upper=4.5> plDens_miss; //Imputed plant density for my fields (only 1)

}

transformed parameters {
	
	//Imputed missing data;
	vector[Nplot_F] plDens; //Plant density
	
	//Combine observed with imputed		
	plDens[plotIndex_F[obsPlDens_ind]]=plDens_obs; //Observed plant density from my fields
	plDens[plotIndex_F[missPlDens_ind]]=plDens_miss; //Missing plant density (only 1) 
	
	//Plot-level parameters
	for(i in 1:Nplot_F){	//Parameters for all fields, all plots
	
	  //Matches F plot i to measurements taken at all plots
	  // "plotI" used to index measurements from ALL plots, "i" used otherwise
	  int plotI = plotIndex_F2[i]; 
				
		CLAIM*plDens[i]; //Claim
	}	
}
	
model {
	
}