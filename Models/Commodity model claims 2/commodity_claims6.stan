data {
	//Field level
	int Nfield; //Number of fields	
	int numHives[Nfield]; //Number of hives present
	int is2015[Nfield]; //Was field measured in 2015?
	int isGP[Nfield]; //Was field from Grand Prairie?
	int isIrrigated[Nfield]; //Was field irrigated?
	
	//Plot level	
	int Nplot; //Number of plots	
	int<lower=1,upper=Nfield> plotIndex[Nplot]; //Index for plots - which field? 
	vector[Nplot] dist; //Distance from edge
	int hbeeVis[Nplot]; //Number of honeybee visits per plot	
	int flyVis[Nplot]; //Number of fly visits per plot
	vector[Nplot] totalTime; //Minutes taken for observation/10	
	vector[Nplot] flDens; //Flower density
	//Plant density (stems/m2) - missing and observed
	int Nplot_densObs; //Number of plots where plant density was observed
	int Nplot_densMiss; //Number of plots with plant density missing 
	vector[Nplot_densObs] plDens_obs; //Observed plant density (stems/m2)
	int<lower=1,upper=Nplot> obsPlDens_ind[Nplot_densObs]; //Index for observed plant density
	int<lower=1,upper=Nplot> missPlDens_ind[Nplot_densMiss]; //Index for missing plant density		
	
	//Flower level
	int Nflw; //Number of flowers
	int<lower=1,upper=Nplot>  flowerIndex[Nflw];  //Index for flowers - which plot?
	int pollenCount[Nflw]; //Pollen count

	//Plant level
	int Nplant; //Number of all plants (some measurements missing)	
	int podCount[Nplant]; //Number of pods per plant
	int flwCount[Nplant]; //Number of total flower (pods + missing) per plant
	int<lower=1,upper=Nplot> plantIndex[Nplant]; //Index for all plants - which plot?		
	vector[Nplant] plantSize; //Centered mass of vegetative tissue (no seeds) (g)
	vector[Nplant] yield; //Observed yield (g seed per plant)	
	
	//Pod level
	int Npod; //Number of pods
	int seedCount[Npod];  //Seeds per pod
	vector[Npod] seedMass; //Mass of seedCount seeds (g)
	int<lower=1,upper=Nplant> podIndex[Npod]; //Index for pods - which plant?	
}

transformed data {
	//Transformations
	vector[Nfield] logNumHives; //Log number of hives	
	vector[Nplot] logHbeeDist=log(dist); //Log-transform distance	
	vector[Nplot] logTime=log(totalTime); //Log-transform time
	vector[Nplot] logHbeeVis;  //Hbee visitation rate
	vector[Nplot] logFlyVis; //Fly visitation rate
	vector[Nplant] logYield = log(yield); //Log yield (g seed per plant)	
	vector[Nplant] logFlwCount; //Log flower count
	vector[Npod] logSeedCount; //Log seed count
	vector[Nplot] totalFlw_obs = rep_vector(0,Nplot); //Total flowers per plot (observed)
	vector[Nplot] numPlants_obs = rep_vector(0,Nplot); //Number of plants observed 	
	vector[Nplot] avgFlws_obs = rep_vector(0,Nplot); //Average observed flowers
	// int numEmptyPlots = 0; //Number of plots with observed flower counts			
	
	logHbeeDist=logHbeeDist-mean(logHbeeDist); //Centers distance
	
	for(i in 1:Nfield){
		logNumHives[i]=log(numHives[i]+1); //Log transform number of hives		
	}
	
	for(i in 1:Nplot){ //Log transform of honeybee visitation rate (per 10 mins)
		logHbeeVis[i] = log((hbeeVis[i]/totalTime[i])+0.5); 
		logFlyVis[i] = log((flyVis[i]/totalTime[i])+0.5);
	}
		
	for(i in 1:Nplant){
		logFlwCount[i] = log(flwCount[i]); //Log-transforms flower count per plant
	}
	
	for(i in 1:Npod){
		logSeedCount[i] = log(seedCount[i]); //Log-transforms seed count per pod
	}
	
	//Calculates average number of flowers per plant (plot level)
	for(i in 1:Nplant){
		totalFlw_obs[plantIndex[i]] = totalFlw_obs[plantIndex[i]]+flwCount[i]; //Adds flowers from plant
		numPlants_obs[plantIndex[i]] = numPlants_obs[plantIndex[i]]+1; //Adds plant to plot
	}
	
	for(i in 1:Nplot){ //Count number of plots with no data
		if(numPlants_obs[i]>0){ //If plants observed
			// numEmptyPlots = numEmptyPlots+1; //	Increment number of empty plots
		// else
			avgFlws_obs[i] = totalFlw_obs[i]/numPlants_obs[i]; //Average number of flowers
		}
	}
}

parameters {
	//Claim: flDens~flwCount
	real slopeFlwCountFlDens; 

	// Plant density	
	vector[Nplot_densMiss] plDens_miss; 
	real intPlDens; //Global intercept
	real slope2015PlDens; //Effect of 2015
	real slopeIrrigPlDens; //Effect of irrigation
	real slope2015IrrigPlDens; //Year:irrigation interaction	
	real slopeGPPlDens; //GP effect on plant density	
	real slopeDistPlDens; //Slope of distance into field		
	real<lower=0.01> sigmaPlDens; //Sigma for within-field (residual)
	real<lower=0.01> sigmaPlDens_field; //Sigma for field
	vector[Nfield] intPlDens_field; //Random intercept for field
	
	// Plant size 	
	real intPlSize; //Global intercept
	real slopePlDensPlSize; //Slope of planting density	
	real slopeDistPlSize; //Slope of distance 
	real slopeGpPlSize; //Slope of Grand Prairie effect	
	real slope2015PlSize; //Slope of 2015 effect		
	real slopeIrrigPlSize; //Slope of irrigation		
	real<lower=0> sigmaPlSize_field; //SD of field-level intercept
	real<lower=0> sigmaPlSize_plot; //SD of plot-level intercept
	real<lower=0> sigmaPlSize; //Sigma for within-plot (residual)
	vector[Nfield] intPlSize_field; //Random intercept for field
	vector[Nplot] intPlSize_plot; //Random intercept for plot - not converging well, but looic is much worse without it		
	
	// Flower density per plot
	real intFlDens; //Global intercept
	real slopePlSizeFlDens; //Slope of plant size on flower density
	real slopeHbeeDistFlDens; //Slope of distance into field
	real<lower=0.01> sigmaFlDens; //Sigma for within-field (residual)	
	real<lower=0.01> sigmaFlDens_field; //Sigma for field
	vector[Nfield] intFlDens_field; //Random intercept for field			
}

transformed parameters {		
	//Expected values
	vector[Nplot] plDensMu; //Predicted plant density
	vector[Nplant] plSizeMu; //Plant size	
	vector[Nplot] plSizePlotMu; //Plot-level plant size		
	vector[Nplot] flDensMu; //Predicted flower density	
	
	// Imputed missing data;
	vector[Nplot] plDens; //Planting density
	vector[Nplot] avgFlws; //Average flowers per plant (some imputed)		
	//Creates empty plot index for flowers per plant (plot level)
	// vector[numEmptyPlots] emptyPlotIndex; //Index for plots with no flower counts per plant		
	// int j = 1; //Starting number for plot index
	// for(i in 1:Nplot){
		// if(numPlants_obs[i]==0){
			// emptyPlotIndex[j] = i;
			// j = j+1;
		// }				
	// }	
	plDens[obsPlDens_ind]=plDens_obs;
	plDens[missPlDens_ind]=plDens_miss;	
	for(i in 1:Nplot){
		if(avgFlws_obs[i]>0)
			avgFlws[i] = avgFlws_obs[i];		
	}	 
	
	for(i in 1:Nplot){ 
		// Plant density per plot
		plDensMu[i] = intPlDens + intPlDens_field[plotIndex[i]] + 
			slope2015PlDens*is2015[plotIndex[i]]+ //Year effect
			slopeIrrigPlDens*isIrrigated[plotIndex[i]]+ //Irrigation effect
			slope2015IrrigPlDens*isIrrigated[plotIndex[i]]*is2015[plotIndex[i]]+ //Year:irrigation interaction
			slopeDistPlDens*logHbeeDist[i] + //Distance effect				
			slopeGPPlDens*isGP[plotIndex[i]]; //GP effect
			
		// Plant size (plot-level) = intercept + random field int + random plot int + distance + planting density effect 		
		plSizePlotMu[i] = intPlSize + intPlSize_field[plotIndex[i]] + intPlSize_plot[i] + 								
			slopePlDensPlSize*plDens[i] +  //Planting density
			slopeDistPlSize*logHbeeDist[i] + //Distance effect (edge of field has smaller plants)			
			slopeGpPlSize*isGP[plotIndex[i]] + // Grand Prairie
			slope2015PlSize*is2015[plotIndex[i]] + //2015
			slopeIrrigPlSize*isIrrigated[plotIndex[i]]; //Irrigation effect												
	
		// Flower density = intercept + random field int + 
		flDensMu[i] = intFlDens	+ intFlDens_field[plotIndex[i]] + 
			slopePlSizeFlDens*plSizePlotMu[i] + //plant size effect 			
			slopeHbeeDistFlDens*logHbeeDist[i] + //distance effect		
			slopeFlwCountFlDens*avgFlws[i]; //Claim
	}
		
	for(i in 1:Nplant){ //For each plant 	
		//Plant size = plot-level estimate
		plSizeMu[i] = plSizePlotMu[plantIndex[i]]; 			
	}	
}
	
model {	
	//Likelihood		
	plDens ~ normal(plDensMu,sigmaPlDens); //Plant density
	plantSize ~ normal(plSizeMu,sigmaPlSize); //Plant size	
	flDens ~ normal(flDensMu,sigmaFlDens); //Flower density per plot	
	
	//Claim
	slopeFlwCountFlDens ~ normal(0,1); 
		
	//Priors
	//Plant density	- informative priors
	intPlDens ~ normal(0,0.1); //Global intercept
	slope2015PlDens ~ normal(0,1); //Year effect
	slopeIrrigPlDens ~ normal(0,1); //Irrigation effect
	slope2015IrrigPlDens ~ normal(0,1); //Year:irrigation interaction
	slopeDistPlDens ~ normal(0,0.05); //Slope of distance into field	
	slopeGPPlDens ~ normal(0,1); // Grand Prairie effect 
	sigmaPlDens ~ gamma(3,10); //Sigma for within-field (residual)
	sigmaPlDens_field ~ gamma(3,10); //Sigma for field
	intPlDens_field ~ normal(0,sigmaPlDens_field); //Random intercept for field
	
	// Plant size - informative priors
	intPlSize ~ normal(0,0.5); //Intercept
	slopePlDensPlSize ~ normal(0,0.5); //Plant density
	slopeDistPlSize ~ normal(0,0.05); //Distance effect
	slopeGpPlSize ~ normal(0,0.5); //Grand Prairie effect
	slopeIrrigPlSize ~ normal(0,0.5); //Irrigation effect
	slope2015PlSize ~ normal(0.3,0.5); //2015 effect				
	sigmaPlSize_field ~ gamma(3,10); //Sigma for random field 
	sigmaPlSize_plot ~ gamma(3,10); //Sigma for random plot
	sigmaPlSize ~ gamma(7,10); //Sigma for residual	
	intPlSize_field ~ normal(0,sigmaPlSize_field); //Random field int
	intPlSize_plot ~ normal(0,sigmaPlSize_plot); //Random int plot		
	
	// Flower density per plot
	intFlDens ~ normal(0,1); //Global intercept
	slopePlSizeFlDens ~ normal(2,1); //plant size effect
	slopeHbeeDistFlDens ~ normal(0,1); //distance into field
	sigmaFlDens ~ gamma(7,2); //Sigma for within-field (residual)	
	sigmaFlDens_field ~ gamma(4,2); //Sigma for field
	intFlDens_field ~ normal(0,sigmaFlDens_field); //Random intercept for field				
}

generated quantities {
}
