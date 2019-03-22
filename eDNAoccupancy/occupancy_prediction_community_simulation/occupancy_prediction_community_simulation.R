library(eDNAoccupancy)
library(plyr)
library(reshape2)
library(tidyverse)
require(VGAM)
setwd("~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/workdir/")
rm(list=ls())
load(file="~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/Data/grenoble_experiment/output_final_dataset/final_data.Rdata")

OTU <- "HISEQ:267:CAJCDANXX:2:1101:3550:2152_CONS_SUB_SUB"

# function to fit eDNAoccupancy model on a single OTU ####
fitEDNAmodel <- function(OTU=myOTU,modeltype="null", niter = niter, burnin = burnin){
  this_otu <- data.frame(final_data$final_otus[,OTU])
  rownames(this_otu) <-rownames(final_data$final_otus)
  names(this_otu) <- "nuReeds"
  
  #separate the sample columns
  this_otu$dna<-sapply(row.names(this_otu),function(x)strsplit(x,"_")[[1]][1])
  this_otu$Rep<-sapply(row.names(this_otu),function(x)strsplit(x,"_")[[1]][2])
  this_otu$Lake<-sapply(row.names(this_otu),function(x)substr(x,1,2))
  this_otu$Depth<-sapply(row.names(this_otu),function(x)strsplit(x,"\\.")[[1]][2])
  this_otu$Depth<-as.numeric(sapply(this_otu$Depth,function(x)strsplit(x,"_")[[1]][1]))
  
  #add PA data
  this_otu$PA <- ifelse(this_otu$nuReeds>0,1,0)
  #table(this_otu$PA)
  
  #get data in the format for eDNA occupancy detection function
  #recode depth levels - adding level2
  this_otu <- ddply(this_otu,.(Lake),function(x){
    x$level2 <- as.integer(factor(x$Depth))
    return(x)
  })
  this_otu$Lake <- as.factor(this_otu$Lake)
  
  #   #order data my Lake and sample
  # this_otu <- arrange(this_otu,level2,Lake)
  # #str(this_otu)
  
  # format data for occModel input
  this_otu_Cast <- dcast(this_otu[,c("Lake","level2","Rep","PA")],
                         Lake+level2~Rep,value.var="PA",fun=sum)
  formatted4Model <- occData(this_otu_Cast, 
                             siteColName = "Lake", sampleColName = "level2")
  # str(formatted4Model$y)
  # str(formatted4Model$K)
  
  # get proxy data, get site data ####
  final_data$proxies %>%
    select(dna, lat, lon, date_final, prod_min, prod_max, LoI, Pb, DDT, S_Fe, erosion) ->
    ourSiteSampleData
  
  #subset data frame for when we have OTU data ####
  ourSiteSampleDataO <- subset(ourSiteSampleData,dna %in% this_otu$dna)
  #centre/scale the variables
  ourSiteSampleData <- scaleData(ourSiteSampleDataO)
  ourSiteSampleData$original_date_final <- ourSiteSampleDataO$date_final
  
  #formatting year groups ####
  #round to nearest 5
  ourSiteSampleData$Year <-  5 * round(ourSiteSampleData$original_date_final/5)
  #dateSummary <- ddply(ourSiteSampleData,.(Year),summarise,
  #                     nuLakes=length(unique(Lake)))

  #before 1900-1950, group by decade
  ourSiteSampleData$Decade <-  10 * round(ourSiteSampleData$original_date_final/10)

  #pre 1900, 50 years?
  #ourSiteSampleData$Quarter <-  25 * round(ourSiteSampleData$original_date_final/25)
  ourSiteSampleData$Half <-  50 * round(ourSiteSampleData$original_date_final/50)

  #pre 1700, pool all
  ourSiteSampleData$Century <-  100 * round(ourSiteSampleData$original_date_final/100)

  #combine groups
  ourSiteSampleData$yearGroups <- NA
  ourSiteSampleData$yearGroups[ourSiteSampleData$original_date_final<1700] <- 1700
  ourSiteSampleData$yearGroups[ourSiteSampleData$original_date_final>1700] <- ourSiteSampleData$Half[ourSiteSampleData$original_date_final>1700]
  ourSiteSampleData$yearGroups[ourSiteSampleData$original_date_final>1900] <- ourSiteSampleData$Decade[ourSiteSampleData$original_date_final>1900]
  ourSiteSampleData$yearGroups[ourSiteSampleData$original_date_final>1950] <- ourSiteSampleData$Year[ourSiteSampleData$original_date_final>1950]
  #str(ourSiteSampleData)
  
  #add level 2
  ourSiteSampleData$level2 <- this_otu$level2[match(ourSiteSampleData$dna,this_otu$dna)]
  ourSiteSampleData$Lake <- sapply(as.character(ourSiteSampleData$dna),
                                   function(x) substr(x,1,2))
  ourSiteSampleData <- arrange(ourSiteSampleData,level2,Lake)
  ourSiteSampleData$level2 <- as.integer(ourSiteSampleData$level2)
  ourSiteSampleData$Lake <- as.factor(ourSiteSampleData$Lake)

    #for testing effect of ecological covariates on community metrics e,.g. richness
    fitModel <- occModel(formulaSite = ~ 1, #lake occurence probabiliy (lake area?)
                         formulaSiteAndSample = ~ factor(yearGroups)-1, #lake-time occurrence probability
                         formulaReplicate = ~ original_date_final, #detection probability at each lake-time
                         detectionMats=formatted4Model,
                         siteColName="Lake",
                         sampleColName="level2",
                         siteAndSampleData = ourSiteSampleData,
                         niter=niter)
    
  #predicted (95%CI for each parameter)
  output <- posteriorSummary(fitModel,burnin=burnin,mcError=T,outputSummary=T)
  
  #combine into a data frame
  Params <- data.frame(output[[1]])
  SEs <- data.frame(output[[2]])
  names(SEs) <- sapply(names(SEs),function(x)paste("MCSE",x,sep="_"))
  output <- cbind(Params,SEs)
  output$Param <- row.names(output)
  
  #add OTU to model output
  output$OTU <- OTU
  
  #return it
  return(output)
}

# many OTU models ####
OTUlist <- c("HISEQ:267:CAJCDANXX:2:1101:3550:2152_CONS_SUB_SUB",
             "HISEQ:267:CAJCDANXX:2:1101:10865:35945_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:3573:3618_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:19514:2820_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:9155:4054_CONS_SUB_SUB",
             "HISEQ:267:CAJCDANXX:2:1101:4872:2562_CONS_SUB_SUB",
             "HISEQ:267:CAJCDANXX:2:1101:18124:5794_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:8581:3757_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:3946:5901_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:20906:26553_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1101:16520:5482_CONS_SUB_SUB",
             "HISEQ:267:CAJCDANXX:2:1101:20004:9611_CONS_SUB_SUB_CMP",
             "HISEQ:267:CAJCDANXX:2:1102:14534:13167_CONS_SUB_SUB_CMP")

allFits_occupancy0 <- ldply(OTUlist[1:3],
                      function(x)fitEDNAmodel(OTU=x,
                                              modeltype="occupancy",
                                              niter = 100, burnin = 50))

# community sampling ####
# subset to occupancy period estimates
allFits_occupancy <- allFits_occupancy0[grepl("alpha.original",
                                   allFits_occupancy0$Param),]

#pull out year data
allFits_occupancy$Year <- gsub("alpha.factor.yearGroups.","",allFits_occupancy$Param)

#get sds
allFits_occupancy$SD1 <- (allFits_occupancy$X97.5.- allFits_occupancy$Mean)/1.96
allFits_occupancy$SD2 <- (allFits_occupancy$Mean - allFits_occupancy$X2.5.)/1.96
allFits_occupancy$SD <- ifelse(allFits_occupancy$SD1>allFits_occupancy$SD2,allFits_occupancy$SD1,allFits_occupancy$SD2)

#create a 1000 communities - the presence/absence of each species is drawn by its probability of occurence  
myMatrix <- matrix(data=NA,nrow=nrow(allFits_occupancy),ncol=1000) 
dim(myMatrix)

#random assign the P/A for each OTU based on its occurrence probabilty
for(j in 1:ncol(myMatrix)){
  for(i in 1:nrow(myMatrix)){
    #myp <- runif(1,min=allFits_occupancy$X2.5.[i],max=allFits_occupancy$X97.5.[i])#unif distribution
    myp <- rnorm(1,mean=allFits_occupancy$Mean[i],sd=allFits_occupancy$SD[i])#normal distribution
    #back-transform
    myp <- probitlink(myp, inverse=T)
    myMatrix[i,j] <- rbinom(1,1,myp) # this changes probability into PA
  }
}

#community random matrix - each column is a simulation
temp <- data.frame(myMatrix)

#add on year/OTU information
temp <- cbind(allFits_occupancy[,c("OTU","Year")],temp)

#get species richness for each community simulation
out <- ddply(temp,.(Year),function(x){
  numcolwise(sum)(x)})
#sum means add up the number of species present

#get mean and 95% CI across simulated communities
meanRichness<-apply(out[,2:1001],1,median)
lowerRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.025))
upperRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.975))

#combine all
finalDF <- data.frame(Year=as.numeric(as.character(sort(unique(allFits_occupancy$Year)))),
                      meanRichness,lowerRichness,upperRichness)

#plotting
ggplot(finalDF)+
  geom_line(aes(x=Year,y=meanRichness))+
  geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness),alpha=0.5)

################################################################################
# not used ####
# # function to get model summaries ####
# get_model_summaries <- function(fitModel){
#   
#   #predicted (95%CI for each parameter)
#   output <- posteriorSummary(fitModel,burnin=500,mcError = T, outputSummary=T)
#   
#   #combine into a data frame
#   Params <- data.frame(output[[1]])
#   SEs <- data.frame(output[[2]])
#   names(SEs) <- sapply(names(SEs),function(x)paste("MCSE",x,sep="_"))
#   output <- cbind(Params,SEs)
#   output$Param <- row.names(output)
#   
#   #add OTU to model output
#   output$OTU <- OTU
#   return(output)
# }

# covariate and null models ####
# #fit the occupancy model ##
# if(modeltype=="null"){
#   #null model
#   fitModel <- occModel(formulaSite = ~1,
#                        formulaSiteAndSample = ~1,
#                        formulaReplicate = ~1,
#                        detectionMats=formatted4Model,
#                        siteColName="Lake",
#                        sampleColName="level2",
#                        niter = niter)
# } else if (modeltype=="covariates"){
#   #with ecological covariates
#   #testing effect of ecological covariates on OTU occurence
#   fitModel <- occModel(formulaSite = ~ 1,#lake occurence probabiliy (lake area?)
#                        formulaSiteAndSample = ~ Pb+erosion+S_Fe,#lake-time occurrence probability
#                        formulaReplicate = ~ original_date_final,#detection probability at each lake-time
#                        detectionMats=formatted4Model,
#                        siteColName="Lake",
#                        sampleColName="level2",
#                        siteAndSampleData = ourSiteSampleData,
#                        niter=niter)
# }else if (modeltype =="occupancy"){
#   #with time-varying covariates   


