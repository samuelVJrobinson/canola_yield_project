#Run claims

for(N in 1:2){
  txt <- readLines('./Seed model claims 3/runClaimsTemplate.R')
  txt <- gsub('MODNUMBER',N,txt)
  writeLines(txt,'./Seed model claims 3/runClaimsCurrent.R')
  tempEnv <- new.env()  
  source('./Seed model claims 3/runClaimsCurrent.R',local=tempEnv,echo=FALSE,chdir = TRUE)
}
