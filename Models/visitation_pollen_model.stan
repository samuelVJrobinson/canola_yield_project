data {
	//Field level
	int Nfield; //Number of fields
	int<lower=1,upper=Nfield> fieldIndex[Nfield]; //Index for fields
	int numHives[Nfield]; //Number of hives present
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field 
	real dist[Nplot]; //Distance from edge - centered
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	real totalTime[Nplot]; //Time taken for observation
	
	//Flower level
	int Nflw; //Number of flowers
	int<lower=1,upper=Nplot>  flowerIndex[Nflw]; 
	int pollenCount[Nflw]; //Number of pollen grains per stigma	
}

parameters {
	//hbee Visitation
	real intVisit; //Intercept 
	real slopeDistVis; //Slope of distance
	real slopeHiveVis; //Slope of hive number
	real<lower=0.2> sigmaVisField; //SD of field random effect
	real<lower=0.2> visitPhi; //Dispersion parameter
	
	//Pollen deposition
	real intPollen; //Intercept
	real slopeVisitPol; //Slope of hbee visits
	real<lower=0.2> sigmaPolField; //SD of field random effect
	real<lower=0.2> sigmaPolPlot; //SD of plot random effect
	real<lower=0.2> pollenPhi; //Dispersion parameter
	
	//Field-level random intercepts
	vector[Nfield] intVisit_field; //hbee visitation
	vector[Nfield] intPollen_field; //pollen deposition
	// Plot-level random intercepts
	vector[Nplot] intPollen_plot; //pollen deposition
}

transformed parameters {	
	//Expected values
	vector[Nplot] visitMu;
	vector[Nflw] pollenMu;
	
	for(i in 1:Nplot) //Expected value for visitation = intercept + random int + distance + numHives + time offset
		visitMu[i] = intVisit + intVisit_field[plotIndex[i]] + slopeDistVis*dist[i] + slopeHiveVis*numHives[plotIndex[i]] + log(totalTime[i]); 
		
	for(i in 1:Nflw) //Expected value for pollen = intercept + random int field + random int plot + hbee visits
		pollenMu[i] = intPollen + intPollen_field[plotIndex[flowerIndex[i]]] + intPollen_plot[flowerIndex[i]] + slopeVisitPol*visitMu[flowerIndex[i]];
}

model {	
	//Likelihood
	hbeeVis ~ neg_binomial_2_log(visitMu,visitPhi);	
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi);
		
	//Priors
	//Visitation
	intVisit ~ normal(0,5); //Intercept
	sigmaVisField ~ gamma(2,1); //Sigma for random field
	intVisit_field ~ normal(0,sigmaVisField); //Random field int
	slopeDistVis ~ normal(0,5); //Slope of distance effect on hbee visits
	slopeHiveVis ~ normal(0,5); //Slope of hive effect on visits
	visitPhi ~ gamma(2,5); //Dispersion parameter
		
	//Pollen deposition
	intPollen ~ normal(5,10); //Intercept	
	sigmaPolField ~ gamma(2,1); //Sigma for random field
	sigmaPolPlot ~ gamma(2,1); //Sigma for random plot
	intPollen_field ~ normal(0,sigmaPolField); //Random field int
	intPollen_plot ~ normal(0,sigmaPolPlot); //Random plot int
	slopeVisitPol ~ normal(0,5); //hbee Visitation effect
	pollenPhi ~ gamma(3.5,5); //Dispersion parameter
}

generated quantities{
int predHbeeVis[Nplot]; //Predicted visits
int predPollenCount[Nflw]; //Predicted pollen counts

	for(i in 1:Nplot)
		predHbeeVis[i] = neg_binomial_2_log_rng(visitMu[i],visitPhi);	
	
	for(i in 1:Nflw)
		predPollenCount[i] = neg_binomial_2_log_rng(pollenMu[i],pollenPhi);
}
