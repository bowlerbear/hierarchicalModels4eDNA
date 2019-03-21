# setwd("C:/Users/db40fysa/Desktop")
# rm(list=ls())
# load("final_data.RData")
library(eDNAoccupancy)
library(plyr)
library(reshape2)

# Miki setup
setwd("~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/workdir/")
rm(list=ls())
load(file="~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/Data/grenoble_experiment/output_final_dataset/final_data.Rdata")


##############################################################################################

# example OTU
most_common_otu <- data.frame(final_data$final_otus[,"HISEQ:267:CAJCDANXX:2:1101:8019:2385_CONS_SUB_SUB"])
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
table(most_common_otu$PA)

#############################################################################################

#examine distribution of depths across lakes
# library(plyr)
# depthSummary <- ddply(most_common_otu,.(Depth),summarise,
#                       nuLakes=length(unique(Lake)),
#                       nuReps = length(unique(Rep)))
# 
# library(reshape2)
# otuMatrix <- acast(most_common_otu,Lake~date_final,
#                    fun=sum,
#                    fill=NA_real_,
#                    value.var="nuReeds")

##############################################################################################
# link dates to depths
most_common_otu$Horiz <- sapply(rownames(most_common_otu),
                                function(x) strsplit(x,"_")[[1]][1])

library(magrittr)
library(dplyr)

most_common_otu <-
most_common_otu %>%
  left_join(final_data$proxies, by = c("Horiz" = "dna")) %>%
  select(nuReeds,Rep,Lake,Depth, date_final)
head(most_common_otu)

# #########################################################################################
# #format without date groups
# 
# #how many reps do we have per date
# library(plyr)
# dateSummary <- ddply(most_common_otu,.(date_final),summarise,
#                      nuLakes=length(unique(Lake)))
# dateSummary

#round to nearest 5
most_common_otu$Year <-  5 * round(most_common_otu$date_final/5)
#dateSummary <- ddply(most_common_otu,.(Year),summarise,
#                     nuLakes=length(unique(Lake)))

#before 1900-1950, group by decade
most_common_otu$Decade <-  10 * round(most_common_otu$date_final/10)

#pre 1900, 50 years?
#most_common_otu$Quarter <-  25 * round(most_common_otu$date_final/25)
most_common_otu$Half <-  50 * round(most_common_otu$date_final/50)

#pre 1700, pool all
most_common_otu$Century <-  100 * round(most_common_otu$date_final/100)

#combine groups
most_common_otu$yearGroups <- NA
most_common_otu$yearGroups[most_common_otu$date_final<1700] <- 1700
most_common_otu$yearGroups[most_common_otu$date_final>1700] <- most_common_otu$Half[most_common_otu$date_final>1700]
most_common_otu$yearGroups[most_common_otu$date_final>1900] <- most_common_otu$Decade[most_common_otu$date_final>1900]
most_common_otu$yearGroups[most_common_otu$date_final>1950] <- most_common_otu$Year[most_common_otu$date_final>1950]

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
str(most_common_otu)

###############################################################################
#format using inbuilt function
# gee <- occData(most_common_otu,siteColName="Lake",sampleColName = "level2")
# gee <- my_occData(most_common_otu,siteColName="Lake",sampleColName = "level2")
# 
# #error checking
# d <- most_common_otu_Cast
# siteColName="Lake"
# sampleColName="level2"
# 
# siteID = unique(d[, siteColName])
# for (i in 1:length(siteID)) {
#   ind = d[, siteColName] == siteID[i]
#   sampleNum = d[ind, sampleColName]
#   if (!identical(sort(sampleNum), 1:length(sampleNum))) {
#     errMsg = paste("Sample numbers of site", siteID[i], 
#                    "do not contain consecutive positive integers in detection file")
#     stop(errMsg)
#   }
# }
# 
# ind = d[, siteColName] == siteID[1]
# sampleNum = d[ind, sampleColName]
# identical(sort(sampleNum), 1:length(sampleNum))
#           
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

#add level 2
ourSiteSampleData$level2 <- most_common_otu$level2[match(ourSiteSampleData$dna,most_common_otu$dna)]
ourSiteSampleData$Lake <- sapply(as.character(ourSiteSampleData$dna),function(x)strsplit(x,"\\.")[[1]][1])
ourSiteSampleData <- arrange(ourSiteSampleData,level2,Lake)
ourSiteSampleData$level2 <- as.integer(ourSiteSampleData$level2)
ourSiteSampleData$Lake <- as.factor(ourSiteSampleData$Lake)
# add yearGroup

#subset data frame for when we have OTU data
ourSiteSampleData <- subset(ourSiteSampleData,dna %in% most_common_otu$dna)
str(ourSiteSampleData)

# center and scale the covariates
ourSiteSampleData_cs <- scaleData(ourSiteSampleData)

############################################################################################
#fit the occupancy model

#null model
fitModel1 <- occModel(formulaSite = ~1,
                     formulaSiteAndSample = ~1,
                     formulaReplicate = ~1,
                     detectionMats=formatted4Model,
                     siteColName="Lake",
                     sampleColName="level2")

#with ecological covariates
#testing effect of ecological covariates on OTU occurence
fitModel2 <- occModel(formulaSite = ~ lat + lon,#lake occurence probabiliy (lake area?)
                     formulaSiteAndSample = ~ prod_max + Pb + DDT + S_Fe + erosion, #lake-time occurrence probability
                     formulaReplicate = ~ date_final,#detection probability at each lake-time
                     detectionMats=formatted4Model,
                     siteColName="Lake",
                     sampleColName="level2",
                     siteAndSampleData = ourSiteSampleData)

psi2 <- posteriorSummaryOfSiteOccupancy(fitModel2, burnin = 1000)
theta <- posteriorSummaryOfSampleOccupancy(fitModel2, burnin = 1000)
p <- posteriorSummaryOfDetection(fitModel2, burnin=1000)
# cbind(psi = psi2$median, theta = theta$median[,1], p = p$median[,1])

#predicted (95%CI for each parameter)
posteriorSummary(fitModel2,burnin=100,mcError=T)

#model checking
plotTrace(fitModel2,paramName=c("beta..Intercept.","beta.lat","beta.lon"))
          
plotTrace(fitModel2,paramName=c("alpha..Intercept.","alpha.prod_max",
                                "alpha.Pb","alpha.DDT"))
          
plotTrace(fitModel2,paramName=c("alpha.S_Fe","alpha.erosion",
                                "delta..Intercept.","delta.date_final"))

plotACF(fitModel2,paramName="beta..Intercept.")

########################################
# testing only time effects
fitModel2 <- occModel(formulaSite = ~ 1,#lake occurence probabiliy (lake area?)
                      formulaSiteAndSample = ~ , #lake-time occurrence probability
                      formulaReplicate = ~ date_final,#detection probability at each lake-time
                      detectionMats=formatted4Model,
                      siteColName="Lake",
                      sampleColName="level2",
                      siteAndSampleData = ourSiteSampleData)

psi2 <- posteriorSummaryOfSiteOccupancy(fitModel2, burnin = 1000)
theta <- posteriorSummaryOfSampleOccupancy(fitModel2, burnin = 1000)
p <- posteriorSummaryOfDetection(fitModel2, burnin=1000)
# cbind(psi = psi2$median, theta = theta$median[,1], p = p$median[,1])

#predicted (95%CI for each parameter)
posteriorSummary(fitModel2,burnin=100,mcError=T)

#model checking
plotTrace(fitModel2,paramName=c("beta..Intercept.","beta.lat","beta.lon"))

plotTrace(fitModel2,paramName=c("alpha..Intercept.","alpha.prod_max",
                                "alpha.Pb","alpha.DDT"))

##############################################################################################
##############################################################################################
siteAndSampleData=ourSiteSampleData
siteColName="Lake"
sampleColName="level2"

if (!is.null(siteAndSampleData)) {
  siteID = unique(siteAndSampleData[, siteColName])
  for (i in 1:length(siteID)) {
    
    ind = siteAndSampleData[, siteColName] == siteID[i]
    sampleNum = siteAndSampleData[ind, sampleColName]
    if (any(sort(sampleNum) != (1:length(sampleNum)))) {
      errMsg = paste("Sample numbers of site", siteID[i], 
                     "do not contain consecutive positive integers in data frame")
      stop(errMsg)
    }
  }
}

gee1 <- sort(sampleNum)
gee2 <- 1:length(sampleNum)
str(gee1)
str(gee2)

###############################################################################################
#with time-varying covariates    
#for testing effect of ecological covariates on community metrics e,.g. richness
fitModel <- occModel(formulaSite = ~1,#lake occurence probabiliy
                     formulaSiteAndSample = ~ factor(YearGroups),#lake-time occurrence probability
                     formulaReplicate = ~ date_final,#detection probability at each lake-time
                     detectionMats=formatted4Model,
                     siteColName="Lake",
                     sampleColName="level2")

##########################################################################################

#predicted (95%CI for each parameter)
posteriorSummary(fitModel,burnin=100,mcError=T)

#posteriors of each level
posteriorSummaryOfSiteOccupancy(fitModel)#probability of a lake being occupied
posteriorSummaryOfSampleOccupancy(fitModel)#probabilty at a given time given lake is occupied
posteriorSummaryOfDetection(fitModel)#probability of detection given present at a given time and lake

#########################################################################################
#model checking
plotTrace(fitModel,paramName="beta..Intercept.")
plotTrace(fitModel,paramName="alpha..Intercept.")
plotACF(fitModel,paramName="beta..Intercept.")
plotACF(fitModel,paramName="alpha..Intercept.")

########################################################################################
#model selection
posteriorPredictiveLoss(fitModel)

#model predictive ability
posteriorSummaryOfAUC(fitModel)
######################################################################################

my_occData<-
function (d, siteColName = "site", sampleColName = "sample") 
{
  if (any(apply(apply(d, 2, is.na), 2, sum) == nrow(d))) {
    warning("Your data file appears to have a column of all NAs. You may want to double check your data. This may have indicated that your file editor (e.g., Excel introduced extra columns in your data.")
  }
  if (!any(grepl(siteColName, names(d), fixed = TRUE))) {
    stop(paste("Column name of sites in detection file does not match", 
               siteColName))
  }
  if (!any(grepl(sampleColName, names(d), fixed = TRUE))) {
    stop(paste("Column name of samples in detection file does not match", 
               sampleColName))
  }
  #siteID = unique(d[, siteColName])
  #for (i in 1:length(siteID)) {
  #  ind = d[, siteColName] == siteID[i]
  #  sampleNum = d[ind, sampleColName]
  #  if (!identical(sort(sampleNum), 1:length(sampleNum))) {
  #    errMsg = paste("Sample numbers of site", siteID[i], 
  #                   "do not contain consecutive positive integers in detection file")
  #    stop(errMsg)
  #  }
  #}
  
  M = length(unique(d[, siteColName]))
  J = max(d[, sampleColName])
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
  list(y = y, K = K)
}

####################################################################################################

#otu surveys
otuSurveys <- unique(most_common_otu[,c("dna","Lake","Depth")])

#proxy surveys
final_data$proxies %>%
  select(dna, lat, lon, date_final, prod_min, prod_max, Pb, DDT, S_Fe, erosion) ->
  ourSiteSampleData

#merge data frames fully
allData <- merge(otuSurveys,ourSiteSampleData,by=c("dna"),all=T)

nrow(subset(allData,is.na(Depth)))
nrow(subset(allData,is.na(Pb)))

####################################################################################################
