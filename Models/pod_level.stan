data {
	int Npod; //number of pods measured	
	int Nunique; //Length of unique seed counts
	int<lower=1> uniquePodCount[Nunique]; //Unique seed counts
	int<lower=1> uniqueMatches[Npod]; //Matching index for seed counts
	int<lower=1> SeedCount[Npod]; //Seeds per pod
	//int<lower=0> PodWeight[Npod]; //Weight per pod
	
}

transformed data {
	int<lower=0> Novules = 100; //Number of ovules (constant)		
}

parameters {
	//Seed count model	
	real<lower=0,upper=1> intPol; //Proportion of ovules -> fertilized ovules
	//(inv-logit) Proportion of fertilized ovules -> seeds
	real intFert; //Intercept
	real slopeFert; //Slope (effect of fert ov count)
	
	
	
	// //Pod weight model
	// real int_Podmass; //Intercept
	// real slope_PodCountPodmass; //Slope of mass-count relationship
	// real<lower=0> sd_Podmass; //SD		
}

transformed parameters {	
}

model {
	matrix[Nunique,Novules] lpSeedsLookup = rep_matrix(-1000,Nunique,Novules); //Lookup matrix for seed log probability
	vector[Npod] lpSeeds; //Log-prob for seed count model		
	
	//Idea: calculate lp once, then look up number for each observed count
	for(seed in 1:Nunique){ //For each observed seed number
		for(ov in 1:Novules){ //For each potential fertilized ovule number
			if(seed<=ov){ //If seed number <= potential fert. ovule number
			//Log-probability of fertilized ovule production & seed production
				lpSeedsLookup[seed,ov]= binomial_lpmf(ov|Novules,inv_logit(intPol)) + //p(Fert ov)
					binomial_lpmf(seed|ov,inv_logit(intFert+slopeFert*ov));  //p(Seed)					
			}
		}	
	}	
	for(i in 1:Npod)
		lpSeeds[i] = sum(lpSeedsLookup[uniqueMatches[i],]); //Look up and assign log-prob to lpSeeds 
	
	//Priors
	intPol ~ normal(0,10); //Should be about -1.83 +/- 0.05 (JAGS)
	intFert ~ normal(0,10); 
	slopeFert ~ normal(0,5); //Should be about 0.009 +- 0.001 (JAGS)
	
	//Likelihood
	target += log_sum_exp(lpSeeds); //Increment log-prob for seed count model
}

generated quantities {
	// int<lower=1,upper=Novules> simFertCount[Npod]; //Simulated fertilized seeds
	// int<lower=1,upper=Novules> simPodCount[Npod]; //Simulated seeds per pod	
	
	// for(i in 1:Npod){
		// simFertCount[i]=categorical_logit_rng(lpSeeds[i,]');	
		// simPodCount[i]=binomial_rng(simFertCount[i],inv_logit(intFert+slopeFert*simFertCount[i]));
	// }
}
