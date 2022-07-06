
// The input data is a vector 'y' of length 'N'.
data {
  int N;
  int Nbeta;
  matrix[N,Nbeta] X;
  int y[N]; 
}

parameters {
  vector[Nbeta] Beta;
  real<lower=0,upper=1> theta;
  real<lower=0> phi;
}

transformed parameters {
  vector[N] mu = X * Beta;
}

model {
  vector[2] bernLL; //pre-calculate LL for zero inflation process
	bernLL[1]=bernoulli_lpmf(0|theta); //LL of no extra zero
	bernLL[2]=bernoulli_lpmf(1|theta); //LL of extra zero
	
	//Zero-inflated negbin
	for(i in 1:N){ 
		if(y[i]==0) //If y is zero
			target += log_sum_exp(bernLL[2],bernLL[1]+neg_binomial_2_log_lpmf(0|mu[i],phi));
		else //If y is not zero
			target += bernLL[1]+neg_binomial_2_log_lpmf(y[i]|mu[i],phi);
	}
	
	//Priors
	Beta ~ normal(0,5);
	theta ~ beta(2,2);
	phi ~ gamma(1,1);
}

generated quantities {
	
	int generated[N]; //Generated
	real resid[N]; //Residual 
	real generated_resid[N]; //Residual of generated 
	for(i in 1:N){
		//hbee visits - ZI neg bin
		resid[i]=y[i]-(exp(mu[i])*(1-theta)); //Residual for actual value x offset
		if(bernoulli_rng(theta)==1) //If theta generates an extra zero
			generated[i] = 0; //Predicted value is automatically zero
		else //Otherwise
			generated[i] = neg_binomial_2_log_rng(mu[i],phi); //Predicted value drawn from neg.bin		
		generated_resid[i]=generated[i]-(exp(mu[i])*(1-theta)); //Residual for predicted value
	}
}


