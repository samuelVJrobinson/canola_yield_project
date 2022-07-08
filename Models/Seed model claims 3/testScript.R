#Template for running claim sets

i <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(i)
i <- i[!is.na(i)]
if(length(i)>1) stop('Too many arguments')

print(paste('This is a test: ',i))
