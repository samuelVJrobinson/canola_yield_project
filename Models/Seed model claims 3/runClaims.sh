#!/bin/bash

# chmod u+x runClaims.sh #Run this to make shell script executable

# for i in {28..31}; do Rscript ./testScript.R $i; done #Test
for i in {15}; do Rscript ./runClaimsTemplate.R $i; done #Actual
