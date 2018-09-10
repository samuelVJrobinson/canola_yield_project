/* This chunk of code was written to experiment with creating a PDF for canola seed counts that would predict both pod success and seed count. Ultimately I think this is interesting but it takes way too long to run (>2hrs). Perhaps if the math was simplified or optimized somehow it could be better. R does an OK job at finding the ML estimate, but perhaps this is too complicated for a Bayesian framework. */

functions {

	//Returns log probability for seeds 1:maxOvNum, given:
	// minOvNum: minimum # of ovules to consider (marginalize across)
	// maxOvNum: maximum # of ovules to consider
	// maxPolNum: maximum # of pollen to consider
	// polMu: mu for neg. binomial polliantion process
	// polPhi: phi for neg. binomial pollination process
	// ovLambda: lambda for poisson ovule production process
	// intPolSurv: intercept for pollen survival (pollen -> fertilized seed)
	// slopePolSurv: slope for pollen survival
	// intPodSurv: intercept for pod survival
	// slopePodSuv: slope for pod survival
	
	vector lpSeeds(int minOvNum, int maxOvNum, int maxPolNum, real polMu, real polPhi, real ovLambda, real intPolSurv, real slopePolSurv, real intPodSurv, real slopePodSurv) {	
		vector[maxPolNum] lpPollen; //Log-prob for pollen counts
		vector[maxOvNum-minOvNum+1] lpOv; //Log prob for ovule counts		
		matrix[maxOvNum,maxPolNum] lpFertSeed; //Log prob for pollen success | ovule success
		vector[maxOvNum] lpPodAbort; //Log prob for pod abortion | fert seed number
		matrix[maxOvNum-minOvNum+1,maxPolNum] lpMarg[maxOvNum]; //Array of matrices to store marginal lp
		vector[maxOvNum] lp; //Lp for each seed count				
		
		//Pre-calculate log prob for marginal vectors
		for(i in 1:maxPolNum) //Pollination
			lpPollen[i] = neg_binomial_2_lpmf(i|polMu,polPhi);
		for(i in 1:(maxOvNum-minOvNum+1))
			lpOv[i] = poisson_lpmf((i+minOvNum-1)|ovLambda); //Ovule production
		for(i in 1:maxOvNum){ 		
			lpPodAbort[i] = bernoulli_logit_lpmf(1|intPodSurv+slopePodSurv*i*.1); //Pod abortion ~ #seeds 
			lpMarg[i] = rep_matrix(0,maxOvNum-minOvNum+1,maxPolNum); //Intial values for lpMarg
			for(j in i:maxPolNum){	//Upper triangular matrix with pollen success log-probs | pollen	 
				lpFertSeed[i,j] = binomial_logit_lpmf(i|j,intPolSurv+slopePolSurv*j*.001); //#Successful ovules~#Pollen (/1000)
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
					(pol==seed ? lpFertSeed[seed,pol] : binomial_lccdf(seed-1|pol,inv_logit(intPolSurv+slopePolSurv*pol/1000)))+ //lp from pollination success				
					//(pol==seed ? lpFertSeed[seed,pol] : log_sum_exp(lpFertSeed[,pol]))+ //lp from pollination success				
					lpPodAbort[seed]; //lp from pod success process				
			}		
			lp[seed] = log_sum_exp(lpMarg[seed,start:end,seed:maxPolNum]); //Sums lp for each unique seed count		
		}	
		return lp;	
	}
	
}

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
	real pollenMu = 293.75; //Parameters for pollen arrival
	real pollenPhi = 0.61;	 
	real ovLambda = 32; //Lambda for ovule production process
	
	int totalFlowers[Nplants]; //Total flowers produced by a plant		
	for(i in 1:Nplants) //Calculate total number of flowers (missing+observed pods)
		 totalFlowers[i]= Pods[i]+PodsMissing[i]; 	
		 
	
}

parameters {
	//Pollen distribution model - negative binomial
	// real<lower=0> pollenMu;
	// real<lower=0> pollenPhi;

	//Seed count model	
	//(inv-logit)Proportion of pollen reaching ovules	
	real intPolSurv; 
	real slopePolSurv; 
		
	//(inv-logit) Proportion of pods that survive
	real intPodSurv; //Intercept
	real slopePodSurv; //Slope (effect of fert ov count)	
}

model {
	vector[maxOvNum] logProbSeeds; //Log probability for seed counts from 1:maxOvNum	
	real probZero; //Prob of observing zero seeds per pod
	
	logProbSeeds = lpSeeds(minOvNum,maxOvNum,maxPolNum,pollenMu,pollenPhi,ovLambda,intPolSurv,slopePolSurv,intPodSurv,slopePodSurv);
	
	//Calculate p(seed!=1:maxOvNum), basically p(seeds==0)
	probZero = 1-sum(exp(logProbSeeds));	
		
	for(i in 1:Nunique){ //Multiply log-lik for seed count by number of data in each category 
		target += logProbSeeds[uniquePodCount[i]]*NuniquePodCount[i]; //Scales by observed number of seed counts, increments LL
	}		
	
	//Priors
	//Pollen survival	
	intPolSurv ~ normal(1,3);
	slopePolSurv ~ normal(-3.5,3);
	//Pod survival
	intPodSurv ~ normal(-3,3); 
	slopePodSurv ~ normal(3.5,3);
	//Pollen distribution
	// pollenMu ~ normal(293.75,10); // Should be about 293
	// pollenPhi ~ gamma(8,12); // Should be about 0.61
	
	//Likelihood	
	// print("Target lpSeeds LL:",sum(lpSeeds));
	// PollenCount ~ neg_binomial_2(pollenMu,pollenPhi); //Pollen counts	
	
	PodsMissing ~ binomial(totalFlowers,probZero); //Missing pod model
}
