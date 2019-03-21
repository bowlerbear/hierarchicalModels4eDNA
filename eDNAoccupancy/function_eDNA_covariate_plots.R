# setwd("C:/Users/db40fysa/Desktop")
# rm(list=ls())
# load("final_data.RData")
library(eDNAoccupancy)
library(plyr)
library(reshape2)
library(tidyverse)
setwd("~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/workdir/")
rm(list=ls())
load(file="~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/Data/grenoble_experiment/output_final_dataset/final_data.Rdata")

##############################################################################################

myOTU = "HISEQ:267:CAJCDANXX:2:1101:13278:2098_CONS_SUB_SUB"

fitEDNAmodel <- function(OTU=myOTU,modeltype="null"){

# example OTU
most_common_otu <- data.frame(final_data$final_otus[,OTU])
rownames(most_common_otu) <-rownames(final_data$final_otus)
names(most_common_otu) <- "nuReeds"

#separate the sample columns
most_common_otu$dna<-sapply(row.names(most_common_otu),function(x)strsplit(x,"_")[[1]][1])
most_common_otu$Rep<-sapply(row.names(most_common_otu),function(x)strsplit(x,"_")[[1]][2])
most_common_otu$Lake<-sapply(row.names(most_common_otu),function(x)strsplit(x,"\\.")[[1]][1])
most_common_otu$Depth<-sapply(row.names(most_common_otu),function(x)strsplit(x,"\\.")[[1]][2])
most_common_otu$Depth<-as.numeric(sapply(most_common_otu$Depth,function(x)strsplit(x,"_")[[1]][1]))
#head(most_common_otu)

############################################################################################

#add PA data
most_common_otu$PA <- ifelse(most_common_otu$nuReeds>0,1,0)
#table(most_common_otu$PA)

#############################################################################################

#examine distribution of depths across lakes
# library(plyr)
# depthSummary <- ddply(most_common_otu,.(Depth),summarise,
#                       nuLakes=length(unique(Lake)),
#                       nuReps = length(unique(Rep)))
# 
# library(reshape2)
# otuMatrix <- acast(most_common_otu,Lake~original_date_final,
#                    fun=sum,
#                    fill=NA_real_,
#                    value.var="nuReeds")

#########################################################################################
#format without date groups

#how many reps do we have per date
#library(plyr)
#dateSummary <- ddply(most_common_otu,.(original_date_final),summarise,
#                      nuLakes=length(unique(Lake)))
#dateSummary


###############################################################################

#get data in the format for eDNA occupancy detection function

#recode depth levels - adding level2
most_common_otu <- ddply(most_common_otu,.(Lake),function(x){
  x$level2 <- as.integer(factor(x$Depth))
  return(x)
})
most_common_otu$Lake <- as.factor(most_common_otu$Lake)
#order data my Lake and sample
most_common_otu <- arrange(most_common_otu,level2,Lake)
#str(most_common_otu)

######################################################################################
#format on our own
most_common_otu_F <- most_common_otu[,c("Lake","level2","Rep","PA")]
most_common_otu_Cast <- dcast(most_common_otu_F,Lake+level2~Rep,value.var="PA",fun=sum)

#run occ Data data
d <- most_common_otu_Cast
siteColName="Lake"
sampleColName="level2"

M = length(unique(d[, siteColName]))#number of lakes
J = max(d[, sampleColName])#number of samples

y = matrix(nrow = M, ncol = J)
K = matrix(nrow = M, ncol = J)
siteID = unique(d[, siteColName])
for (i in 1:length(siteID)) {
  ind = d[, siteColName] == siteID[i]
  sampleNum = d[ind, sampleColName]
  ymat = d[ind, -c(1:2)]
  y[i, sampleNum] = rowSums(ymat, na.rm = TRUE)
  K[i, sampleNum] = rowSums(!is.na(ymat))
}
indMiss = is.na(K)
K[indMiss] = 0
y[K == 0] = NA
dimnames(y)[[1]] = as.character(siteID)
dimnames(K)[[1]] = as.character(siteID)
formatted4Model <- list(y = y, K = K)

str(formatted4Model$y)
str(formatted4Model$K)

####################################################################################
#get proxy data
# get site data
final_data$proxies %>%
  select(dna, lat, lon, date_final, prod_min, prod_max, Pb, DDT, S_Fe, erosion) ->
  ourSiteSampleData

#subset data frame for when we have OTU data
ourSiteSampleDataO <- subset(ourSiteSampleData,dna %in% most_common_otu$dna)
#centre/scale the variables
ourSiteSampleData <- scaleData(ourSiteSampleDataO)
ourSiteSampleData$original_date_final <- ourSiteSampleDataO$date_final

#add level 2
ourSiteSampleData$level2 <- most_common_otu$level2[match(ourSiteSampleData$dna,most_common_otu$dna)]
ourSiteSampleData$Lake <- sapply(as.character(ourSiteSampleData$dna),function(x)strsplit(x,"\\.")[[1]][1])
ourSiteSampleData <- arrange(ourSiteSampleData,level2,Lake)
ourSiteSampleData$level2 <- as.integer(ourSiteSampleData$level2)
ourSiteSampleData$Lake <- as.factor(ourSiteSampleData$Lake)

#str(ourSiteSampleData)

###########################################################################################

#formatting data groups

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

############################################################################################
#fit the occupancy models

##############################################################################################

if(modeltype=="null"){
#null model
fitModel <- occModel(formulaSite = ~1,
                     formulaSiteAndSample = ~1,
                     formulaReplicate = ~1,
                     detectionMats=formatted4Model,
                     siteColName="Lake",
                     sampleColName="level2")

##############################################################################################
}else if (modeltype=="covariates"){
#with ecological covariates
#testing effect of ecological covariates on OTU occurence
fitModel <- occModel(formulaSite = ~ 1,#lake occurence probabiliy (lake area?)
                     formulaSiteAndSample = ~ Pb+erosion+S_Fe,#lake-time occurrence probability
                     formulaReplicate = ~ original_date_final,#detection probability at each lake-time
                     detectionMats=formatted4Model,
                     siteColName="Lake",
                     sampleColName="level2",
                     siteAndSampleData = ourSiteSampleData,
                     niter=2000)

###############################################################################################

}else if (modeltype =="time"){
#with time-varying covariates    
#for testing effect of ecological covariates on community metrics e,.g. richness
  fitModel <- occModel(formulaSite = ~ 1,#lake occurence probabiliy (lake area?)
                       formulaSiteAndSample = ~ factor(yearGroups)-1,#lake-time occurrence probability
                       formulaReplicate = ~ original_date_final,#detection probability at each lake-time
                       detectionMats=formatted4Model,
                       siteColName="Lake",
                       sampleColName="level2",
                       siteAndSampleData = ourSiteSampleData,
                       niter=2000)
}

##########################################################################################

#predicted (95%CI for each parameter)
output <- posteriorSummary(fitModel,burnin=500,mcError=T,outputSummary=T)

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

#posteriors of each level
posteriorSummaryOfSiteOccupancy(fitModel)#probability of a lake being occupied
#posteriorSummaryOfSampleOccupancy(fitModel)#probabilty at a given time given lake is occupied
#posteriorSummaryOfDetection(fitModel)#probability of detection given present at a given time and lake

#########################################################################################

#model checking
plotTrace(fitModel,paramName="beta..Intercept.")
plotTrace(fitModel,paramName="alpha..Intercept.")
plotTrace(fitModel,paramName="delta..Intercept.")
plotTrace(fitModel,paramName="alpha.Pb")
plotTrace(fitModel,paramName="alpha.erosion")
plotTrace(fitModel,paramName="alpha.S_Fe")
plotTrace(fitModel,paramName="delta.original_date_final")
# plotACF(fitModel,paramName="beta..Intercept.")
# plotACF(fitModel,paramName="alpha..Intercept.")
# 
# ########################################################################################
# #model selection
# posteriorPredictiveLoss(fitModel)
# 
# #model predictive ability
# posteriorSummaryOfAUC(fitModel)
######################################################################################

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


library(plyr)
allFits <- ldply(OTUlist,function(x)fitEDNAmodel(OTU=x,modeltype="covariates"))

# Plot segment plots ####
allFits %>%
  filter(Param == "alpha.S_Fe") %>%
  ggplot(aes(reorder(OTU, X50.), y = X50.)) + # reorder: reorders factor with values from a variable
  geom_crossbar(aes(ymin = X2.5., ymax = X97.5.)) +
  coord_flip() +
  geom_hline(aes(yintercept=0)) +
  theme(legend.position = "none")
###############################################################################

#community sampling

#create a 1000 communities - the presence/absence of each species is drawn by its probability of occurence  
myMatrix <- matrix(data=NA,nrow=nrow(modelSummary_Ad),ncol=1000) 


#random assign the P/A for each OTU based on its occurrence probabilty
for(j in 1:ncol(myMatrix)){
  for(i in 1:nrow(myMatrix)){
    myp <- dnorm(1,modelSummary_Ad$mean[i],modelSummary_Ad$sd[i]) 
    myp[myp>1]<-0.9999
    myp[myp<0]<-0.0001
    myMatrix[i,j] <- rbinom(1,1,myp)
  }
}
randomMatrix<-cbind(modelSummary_Ad[,c("State","Year","Species")],myMatrix)
save(randomMatrix,file="randomMatrix.RData")

#get species richness for each community
out <- ddply(randomMatrix,.(State,Year),function(x){
  numcolwise(sum)(x)})

#get mean and 95% CI across simulated communities
out$meanRichness<-apply(out[,3:1002],1,median)
out$lowerRichness<-apply(out[,3:1002],1,function(x)quantile(x,0.025))
out$upperRichness<-apply(out[,3:1002],1,function(x)quantile(x,0.975))

#####################################################################################