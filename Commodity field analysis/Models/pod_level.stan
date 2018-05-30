data {
	int Npod; //number of pods measured
	int FertCount[Npod]; //number of fertilized ovules (all NAs)
	int<lower=0> SeedCount[Npod]; //seeds per pod	
}

transformed data {
	int<lower=0> Novules; //Number of ovules (constant)
	Novules=120;
}

parameters {
	real int_Pol;		
	real slope_Fert;	
	//int FertCount[Npod];	
	
}

transformed parameters{
	//Transformed parameters
	real<lower=0> p1[Npod];
	real<lower=0> p2[Npod];
	
	for(pod in 1:Npod){
		p1[Npod]=exp(int_Pol)*Novules;
		p2[Npod]=exp(slope_Fert*FertCount[Npod])*FertCount[Npod];
	}

}

model {

	//Priors
	int_Pol ~ normal(-2,5);
	slope_Fert ~ normal(0,5);
	
	//Likelihood
	FertCount ~ poisson(p1);
	SeedCount ~ poisson(p2); 
	
	/* slope_Fert ~ normal(0.1,1);
	
	for(pod in 1:Npod){
		FertCount[pod] ~ poisson_log(p1*Novules[pod]);
		SeedCount[pod] ~ poisson_log(p2*FertCount[pod]);
	} */
}
