library(BaseSpaceR)
ACCESS_TOKEN <- 'e09b20a266a44225918dc3bd58df3f26'
PROJECT_ID <- '17750733'  ## Get proj ID from url of the project
experiment_name <- "gayther-jae_jung-perin-rong_lu_150721"

aAuth<- AppAuth(access_token = ACCESS_TOKEN)
selProj <- Projects(aAuth, id = PROJECT_ID, simplify = TRUE) 
sampl <- listSamples(selProj, limit= 1000)
sampl <- sampl[unlist(lapply(1:length(sampl), function(i) sampl[i]$ExperimentName))==experiment_name]
inSample <- Samples(aAuth, id = Id(sampl), simplify = TRUE)
for(s in inSample){ 
    f <- listFiles(s, Extensions = ".gz")
    print(Name(f))
    getFiles(aAuth, id= Id(f), destDir = '072215/', verbose = TRUE)
}
