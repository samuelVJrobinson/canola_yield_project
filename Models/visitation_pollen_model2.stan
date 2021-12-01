
data {
	//Field level
	int Nfield; //Number of fields	
	int numHives[Nfield]; //Number of hives present

	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
	vector[Nplot] dist; //Distance from edge
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	int flyVis[Nplot]; //Number of fly visits per plot
	vector[Nplot] totalTime; //Minutes taken for observation/10	
	vector[Nplot] flDens; //Flower density
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives	
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	
	logHbeeDist=logHbeeDist-mean(logHbeeDist); //Centers distance
	
	for(i in 1:Nfield){
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives		
	}
	
	for(i in 1:Nplot){ //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5); 
	}
}

parameters {

	// hbee Visitation per plot	
	real intVisit; //Intercept 
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real slopeFlDens; //Slope of flower density
	real<lower=0> sigmaVisField; //SD of field random intercepts
	real<lower=0> visitHbeePhi; //Dispersion parameter	
	vector[Nfield] intVisit_field; //field-level random intercepts		
	real<lower=0> lambdaVisField; //Lambda for skewed random effects	
	
}

transformed parameters {		
	//Expected values
	vector[Nplot] visitHbeeMu; //Plot-level hbee visits	
	
	for(i in 1:Nplot){ 
		
		// Hbee Visitation = intercept + random int + time offset + 	
		visitHbeeMu[i] = intVisit + intVisit_field[plotIndex[i]] + //Intercept + random effect
		  logTime[i] + //Offset
			slopeDistVis*logHbeeDist[i] + //distance from edge 
			slopeHiveVis*logNumHives[plotIndex[i]] + //(log) Number of hives
			slopeFlDens*flDens[i]; //Flower density
	}
}
	
model {	
	//Likelihood		
	hbeeVis ~ neg_binomial_2_log(visitHbeeMu,visitHbeePhi); //Honeybee visitation (no ZI-process)
		
	// Priors
	
	// Visitation - informative priors
	intVisit ~ normal(-1,1); //Intercept	
	slopeDistVis ~ normal(-0.3,0.3); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0.6,0.5); //Slope of hive effect on visits 
	slopeFlDens ~ normal(0,0.5); //Flower density 
	sigmaVisField ~ gamma(4,4); //Sigma for random field 		
	intVisit_field ~ exp_mod_normal(0,sigmaVisField,lambdaVisField); //Skewed random effects - slightly better than standard normal
	lambdaVisField ~ gamma(4,2); //Lambda for skewed random effects			
	visitHbeePhi ~ gamma(4,10); //Dispersion parameter for NegBin		
		
}
 
