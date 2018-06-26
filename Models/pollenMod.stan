data {
	int Npollen; //Number of pollen measurements	
	int<lower=0> PollenCount[Npollen]; //Pollen per stigma	
}

parameters {
	//Pollen distribution model - negative binomial
	real<lower=0> pollenMu;
	real<lower=0> pollenPhi;
}

model {
	//Priors
	pollenMu ~ normal(290,40); // Should be about 293
	pollenPhi ~ gamma(2,3); // Should be about 0.61
	
	//Likelihood		
	PollenCount ~ neg_binomial_2(pollenMu,pollenPhi); //Pollen counts			
}
