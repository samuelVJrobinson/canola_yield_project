data {
	//Field level
	int Nfield; //Number of fields
	int<lower=1,upper=Nfield> fieldIndex[Nfield]; //Index for fields
	int numHives[Nfield]; //Number of hives present
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots	
	real dist[Nplot]; //Distance from edge
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	real totalTime[Nplot]; //Time taken for observation
	
	//Flower level
	int Nflw; //Number of flowers
	int<lower=1,upper=Nplot>  flowerIndex[Nflw]; 
	int pollenCount[Nflw]; //Number of pollen grains per stigma	
}

parameters {
	real intVisit; //Intercept for hbee visitation
	real slopeDistVis; //Effect of distance on hbee visitation
	real slopeHiveVis; //Effect of hive number on hbee visitation
	
	real intPollen; //Intercept for pollen counts
	real slopeVisitPol; //Effect of hbee visit on pollen count
}

transformed parameters{
	//Field-level intercepts
	vector[Nfield] intVisit_field; //hbee visitation
	vector[Nfield] intPollen_field; //pollen deposition
	
	//Plot-level intercepts
	vector[Nplot] intPollen_plot; //pollen deposition
	
	//Predicted values
	vector[Nplot] visitMu;
	vector[Nflw] pollenMu;
	
	for(i in 1:Nfield){
		visitMu[i] = intVisit + slopeDistVis*dist + slopeHiveVis*numHives + exp(totalTime); 
	}
	
	
	pollenMu = intPollen + slopeVisitPol*hbeeVis;	
	
}

model {
	
	hbeeVis ~ neg_binomial_2_log(visitMu,visitPhi);
	
	
	pollenCount ~ neg_binomial_2_log(pollenMu,pollenPhi);
	
	
	//Priors
	visitPhi ~ dgamma(2,1);
	slopeDistVis ~ dnorm(0,5);
	slopeHiveVis ~ dnorm(0,5);
}
