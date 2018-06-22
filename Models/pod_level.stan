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
	int maxPolNum = 1500; //Max number of pollen ~ 98% of nbinom density
	int maxOvNum = 52; //Max number of ovules
	int minOvNum = 12; //Min number of ovules
	int totalFlowers[Nplants]; //Total flowers produced by a plant
	real<lower=-10,upper=0.1> slopePolSurv=0; 
	
	for(i in 1:Nplants) //Calculate total number of flowers (missing+observed pods)
		 totalFlowers[i]= Pods[i]+PodsMissing[i]; 	
}

parameters {
	//Pollen distribution model - negative binomial
	real<lower=0,upper=400> pollenMu;
	real<lower=0,upper=2> pollenPhi;

	//Seed count model	
	//(inv-logit)Proportion of pollen reaching ovules
	// real<lower=0,upper=1> propPol; 
	real<lower=-10,upper=2> intPolSurv; 
	// real<lower=-10,upper=0.1> slopePolSurv; 
	//(inv-logit) Proportion of pods that survive
	real<lower=-10,upper=10> intPodSurv; //Intercept
	real<lower=-0.1,upper=10> slopePodSurv; //Slope (effect of fert ov count)	
}

model {
	vector[maxPolNum] lpPollen; //Log-prob for pollen counts
	vector[maxOvNum-minOvNum+1] lpOv; //Log prob for ovule counts		
	matrix[maxOvNum,maxPolNum] lpFertSeed; //Log prob for pollen success | ovule success
	vector[maxOvNum] lpPodAbort; //Log prob for pod abortion | fert seed number
	matrix[maxOvNum-minOvNum+1,maxPolNum] lpMarg[maxOvNum]; //Array of matrices to store marginal lp
	vector[maxOvNum] lpSeeds; //Lp for each seed count		
	real probZero; //Prob of observing zero seeds per pod
	
	//Calculate log prob for marginal vectors
	for(i in 1:maxPolNum) //Pollination
		lpPollen[i] = neg_binomial_2_lpmf(i|pollenMu,pollenPhi);
	for(i in 1:(maxOvNum-minOvNum+1))
		lpOv[i] = poisson_lpmf((i+minOvNum-1)|32); //Ovule production
	for(i in 1:maxOvNum){ 		
		lpPodAbort[i] = bernoulli_logit_lpmf(1|intPodSurv+slopePodSurv*i); //Pod abortion ~ #seeds
		lpMarg[i] = rep_matrix(0,maxOvNum-minOvNum+1,maxPolNum); //Intial values for lpMarg
		for(j in i:maxPolNum){	//Upper triangular matrix with pollen success log-probs | pollen	 
			lpFertSeed[i,j] = binomial_logit_lpmf(i|j,intPolSurv+slopePolSurv*j); //#Successful ovules~#Pollen
		}		
	}	
	
	//Marginal probability matrix
	for(seed in 1:maxOvNum){ //For each seed count		
		int start=(seed>minOvNum ? seed-minOvNum+1 : 1); //Start at position 1 (minimum ovule count), unless seed > minOvNum			
		int end=maxOvNum-minOvNum+1; //End position
		for(pol in seed:maxPolNum){ //For each pollen count >= seed count			
			//Use pre-calculated lps for seed+1:maxOvNum 
			//print("start:",start," end:",end);
			if(start<end){ //If index not at end of matrix (ie last possible seed size)
				lpMarg[seed,(start+1):end,pol] = rep_vector(lpPollen[pol],end-start) + //lp from pollen process
					lpOv[(start+1):end]+ //lp from ovule process				
					rep_vector(lpFertSeed[seed,pol],end-start)+ //lp from pollination success - not sure if indexing is right??
					rep_vector(lpPodAbort[seed],end-start); //lp from pod success process
			}
			//Use mixture of lps for pol==seed
			lpMarg[seed,start,pol] = lpPollen[pol]+ //lp from pollen process
				lpOv[start]+ //lp from ovule process				
				(pol==seed ? lpFertSeed[seed,pol] : binomial_lccdf(seed-1|pol,inv_logit(intPolSurv+slopePolSurv*pol)))+ //lp from pollination success				
				//(pol==seed ? lpFertSeed[seed,pol] : log_sum_exp(lpFertSeed[,pol]))+ //lp from pollination success				
				lpPodAbort[seed]; //lp from pod success process				
		}		
		lpSeeds[seed] = log_sum_exp(lpMarg[seed,start:end,seed:maxPolNum]); //Sums lp for each unique seed count		
	}	
	// print("lpSeeds: ",lpSeeds);	
	
	//Calculate p(seed!=1:maxOvNum), basically p(seeds==0)
	probZero = 1-sum(exp(lpSeeds));
	// print("ProbZero: ",probZero, " intPolSurv:",intPolSurv," slopePolSurv:",slopePolSurv," intPodSurv:",intPodSurv," slopePodSurv:",slopePodSurv);
		
	for(i in 1:Nunique){ //Multiply log-lik for seed count by number of data in each category 
		target += lpSeeds[uniquePodCount[i]]*NuniquePodCount[i]; //Scales by observed number of seed counts, increments LL
	}		
	
	//Priors
	// propPol ~ beta(1.1,3); //Should be low
	intPolSurv ~ normal(0,3);
	// slopePolSurv ~ normal(0,3);
	intPodSurv ~ normal(0,3); 
	slopePodSurv ~ normal(0,3); //Should be about 0.009 +- 0.001 (JAGS)
	pollenMu ~ normal(290,40); // Should be about 293
	pollenPhi ~ gamma(2,3); // Should be about 0.61
	
	//Likelihood	
	// print("Target lpSeeds LL:",sum(lpSeeds));
	PollenCount ~ neg_binomial_2(pollenMu,pollenPhi); //Pollen counts	
	
	PodsMissing ~ binomial(totalFlowers,probZero); //Missing pod model
}
