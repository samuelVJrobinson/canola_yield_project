data {
	int Npod; //number of pods measured	
	int Nunique; //Length of unique seed counts
	int Npollen; //Number of pollen measurements
	int<lower=1> uniquePodCount[Nunique]; //Unique seed counts
	vector[Nunique] NuniquePodCount; //Number of matching pod counts (from contingency table)
	int<lower=1> uniqueMatches[Npod]; //Matching index for seed counts
	int<lower=1> SeedCount[Npod]; //Seeds per pod
	int<lower=0> PollenCount[Npollen]; //Pollen per stigma
	int Nplants; //Number of plants
	int Pods[Nplants]; //Number of observed pods
	int PodsMissing[Nplants]; //Number of missing pods	
}

transformed data {
	//Limits for marginalization procedure
	int<lower=0> maxPolNum = 1800; //Max number of pollen
	int<lower=0> maxOvNum = 52; //Max number of ovules
	int totalFlowers[Nplants]; //Total flowers produced by a plant
	
	for(i in 1:Nplants)
		 totalFlowers[i]= Pods[i]+PodsMissing[i]; 	
}

parameters {
	//Pollen distribution model - negative binomial
	real<lower=0> pollenMu;
	real<lower=0> pollenPhi;

	//Seed count model	
	real<lower=0,upper=1> propPol; //Proportion of pollen reaching ovules
	//(inv-logit) Proportion of pods that become seeds
	real intPodSurv; //Intercept
	real slopePodSurv; //Slope (effect of fert ov count)	
}

model {
	vector[maxPolNum] lpPollen; //Log-prob for pollen counts
	vector[maxOvNum] lpOv; //Log prob for ovule counts		
	matrix[maxOvNum,maxPolNum] lpFertSeed; //Log prob for pollen success | ovule success
	vector[maxOvNum] lpPodAbort; //Log prob for pod abortion | fert seed number
	matrix[maxOvNum,maxPolNum] lpMarg[maxOvNum]; //Array of matrices to store marginal lp
	vector[maxOvNum] lpSeeds; //Lp for each seed count		
	real probZero; //Prob of observing zero seeds per pod
	
	//Calculate log prob for marginal vectors
	for(i in 1:maxPolNum) //Pollination
		lpPollen[i] = neg_binomial_2_lpmf(i|pollenMu,pollenPhi);
	for(i in 1:maxOvNum){ //Ovule production
		lpOv[i] = poisson_lpmf(i|32);
		lpPodAbort[i] = bernoulli_logit_lpmf(1|intPodSurv+slopePodSurv*i); //Pod abortion
		lpMarg[i] = rep_matrix(0,maxOvNum,maxPolNum); //Intial values for lpMarg
		for(j in i:maxPolNum){	//Upper triangular matrix with pollen success log-probs | pollen	 
			lpFertSeed[i,j] = binomial_lpmf(i|j,propPol);
		}		
	}	
	
	//Marginal probability matrix
	for(seed in 1:maxOvNum){ //For each seed count		
		for(pol in seed:maxPolNum){ //For each pollen count >= seed count			
			//Use pre-calculated lps for seed+1:maxOvNum 
			lpMarg[seed,(seed+1):maxOvNum,pol] = rep_vector(lpPollen[pol],maxOvNum-seed) + //lp from pollen process
				lpOv[(seed+1):maxOvNum]+ //lp from ovule process				
				rep_vector(lpFertSeed[seed,pol],maxOvNum-seed)+ //lp from pollination success
				rep_vector(lpPodAbort[seed],maxOvNum-seed); //lp from pod success process
			//Use mixture of lps for pol==seed
			lpMarg[seed,seed,pol] = lpPollen[pol]+ //lp from pollen process
				lpOv[seed]+ //lp from ovule process				
				(pol==seed ? lpFertSeed[seed,pol] : binomial_lccdf(seed-1|pol,propPol))+ //lp from pollination success				
				lpPodAbort[seed]; //lp from pod success process					
		}		
	}
	
	for(i in 1:maxOvNum){ //Sums lp for each unique seed count		
		lpSeeds[i] = log_sum_exp(lpMarg[i,i:maxOvNum,i:maxPolNum]);
	}
	//print("lpSeeds: ",lpSeeds);
	
	//Calculate p(seed!=1:maxOvNum), basically p(seeds==0)
	probZero = 1-sum(exp(lpSeeds));
	print("ProbZero: ",probZero);
		
	for(i in 1:Nunique){ //Multiply log-lik for seed count by number of data in each category 
		target += lpSeeds[uniquePodCount[i]]*NuniquePodCount[i]; //Scales by observed number of seed counts, increments LL
	}	
	// print("Target lpSeeds LL:",sum(lpSeeds));
	
	//Priors
	propPol ~ beta(1.1,3); //Should be low
	intPodSurv ~ normal(0,6); 
	slopePodSurv ~ normal(0,6); //Should be about 0.009 +- 0.001 (JAGS)
	pollenMu ~ normal(290,40); // Should be about 293
	pollenPhi ~ gamma(2,3); // Should be about 0.61
	
	//Likelihood	
	// print("Target lpSeeds LL:",sum(lpSeeds));
	PollenCount ~ neg_binomial_2(pollenMu,pollenPhi); //Pollen counts	
	
	//PodsMissing ~ binomial(totalFlowers,probZero); //Missing pod model
}
