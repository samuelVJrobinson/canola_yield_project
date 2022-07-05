parameters {
 
	// Flower density per plot
	vector<lower=5,upper=51>[Nplot_flDensMiss] flDens_miss; //Missing from my fields
}

transformed parameters {
			
	//Imputed missing data;
	vector[Nplot] flDens; //Flower density
	
	//Combine observed with imputed		
	flDens[obsflDens_ind]=flDens_obs; //Observed flower density
	flDens[missflDens_ind]=flDens_miss;
		
	for(i in 1:Nplot_F){ //Parameters for F plots only
	
	  //Matches F plot i to measurements taken at all plots
	  // "plotI" used to index measurements from ALL plots, "i" used otherwise
	  int plotI = plotIndex_F2[i]; 

   	  CLAIM*flDens[plotI];
   	  
	}
	
}
	
model {

}
