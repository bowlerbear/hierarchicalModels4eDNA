library(eDNAoccupancy)
library(plyr)
library(reshape2)
library(tidyverse)
require(VGAM)
setwd("~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/workdir/")
rm(list=ls())
load(file="~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/Data/grenoble_experiment/output_final_dataset/final_data.Rdata")
# load(file="../final_data.Rdata")

# OTUs
final_data$final_assign %>%
  filter(final_data$final_assign$domain == "Eukaryota") %>%
  select(otu) %>%
  flatten() %>%
  unlist() ->
  OTUlist

set1 <- OTUlist[1:2000]
set2 <- OTUlist[2001:4000]
set3 <- OTUlist[4001:6000]
set4 <- OTUlist[6001:length(OTUlist)]

# get environmental covariates ####
final_data$proxies %>%
  select(dna, lat, lon, date_final, prod_min, prod_max, LoI, Pb, DDT, S_Fe, erosion) ->
  ourSiteSampleData

ourSiteSampleData <- subset(ourSiteSampleData, dna %in% rownames(final_data$aggr_replsum_otus))

# # put replicates into site sample data
data.frame(final_data$final_assign, t(final_data$aggr_replsum_otus)) %>%
  gather(AR1.010:WM1.490, key="horizon", value="reps") %>%
  filter(domain == "Eukaryota", !is.na(phylum)) %>%
  group_by(horizon) %>%
  summarize(reps_sum = sum(reps)) ->
  replicates

ourSiteSampleData %>%
  left_join(replicates, by = c("dna" = "horizon")) ->
  ourSiteSampleData

# scaled reps_sum
ourSiteSampleData$reps_sum_scaled <- scale(ourSiteSampleData$reps_sum)

# add lake, depth and upper layers, date polynomials
ourSiteSampleData$Lake <- sapply(as.character(ourSiteSampleData$dna),
                                 function(x) substr(x,1,2))
ourSiteSampleData$Lake <- as.factor(ourSiteSampleData$Lake)
ourSiteSampleData$Depth <- as.numeric(sapply(as.character(ourSiteSampleData$dna),
                                             function(x) substr(x,5,7)))
ourSiteSampleData$upper <- ifelse(ourSiteSampleData$Depth <= 20, "upper", "lower")
ourSiteSampleData$upper <- factor(ourSiteSampleData$upper, levels = c("lower", "upper"))
ourSiteSampleData$date_final_1 <- as.numeric(scale(ourSiteSampleData$date_final))
ourSiteSampleData$date_final_2 <- as.numeric(scale(ourSiteSampleData$date_final^2))
ourSiteSampleData$date_final_3 <- as.numeric(scale(ourSiteSampleData$date_final^3))
ourSiteSampleData$date_final_4 <- as.numeric(scale(ourSiteSampleData$date_final^4))
ourSiteSampleData$date_final_5 <- as.numeric(scale(ourSiteSampleData$date_final^5))
str(ourSiteSampleData)

OTU <- "HISEQ:267:CAJCDANXX:2:1101:4969:25485_CONS_SUB_SUB"

# function to fit eDNAoccupancy model on a single OTU ####
fitEDNAmodel <- function(OTU=myOTU, niter = niter, burnin = burnin){
  this_otu <- data.frame(final_data$final_otus[,OTU])
  rownames(this_otu) <-rownames(final_data$final_otus)
  names(this_otu) <- "nuReeds"
  
  #separate the sample columns
  this_otu$dna<-sapply(row.names(this_otu),function(x)strsplit(x,"_")[[1]][1])
  this_otu$Rep<-sapply(row.names(this_otu),function(x)strsplit(x,"_")[[1]][2])
  this_otu$Lake<-sapply(row.names(this_otu),function(x)substr(x,1,2))
  this_otu$Depth<-sapply(row.names(this_otu),function(x)strsplit(x,"\\.")[[1]][2])
  this_otu$Depth<-as.numeric(sapply(this_otu$Depth,function(x)strsplit(x,"_")[[1]][1]))
  
  #get data in the format for eDNA occupancy detection function
  #recode depth levels - adding level2
  this_otu <- ddply(this_otu,.(Lake),function(x){
    x$level2 <- as.integer(factor(x$Depth))
    return(x)
  })
  this_otu$Lake <- as.factor(this_otu$Lake)
  
  #add PA data ####
  this_otu$PA <- ifelse(this_otu$nuReeds>0,1,0)
  #table(this_otu$PA)
  
  #add level 2 to env. covariates
  ourSiteSampleData$level2 <- this_otu$level2[match(ourSiteSampleData$dna,this_otu$dna)]
  ourSiteSampleData$level2 <- as.integer(ourSiteSampleData$level2)
  # ourSiteSampleData <- arrange(ourSiteSampleData,level2,Lake)
  
  # format data for occModel input ####
  this_otu_Cast <- dcast(this_otu[,c("Lake","level2","Rep","PA")],
                         Lake+level2~Rep,value.var="PA",fun=sum)
  formatted4Model <- occData(this_otu_Cast, 
                             siteColName = "Lake", sampleColName = "level2")
  
  # fit model ####
    #for testing effect of ecological covariates on community metrics e,.g. richness
  fitModel <- occModel(formulaSite = ~ 1,
                       formulaSiteAndSample = ~ date_final_1 + date_final_2 + date_final_3,
                       formulaReplicate = ~ reps_sum_scaled + upper,
                       detectionMats=formatted4Model,
                       siteColName="Lake",
                       sampleColName="level2",
                       siteAndSampleData = ourSiteSampleData,
                       niter=1000)
    
  niter = 1000
    # plot the traces
    post_sum <- posteriorSummary(fitModel, outputSummary=T, mcError = T, 
                                 burnin = niter/2)
    post_sample <- posteriorSummaryOfSampleOccupancy(fitModel, mcError = T, 
                                                     burnin = niter/2)
    post_detect <- posteriorSummaryOfDetection(fitModel, mcError = T, 
                                                     burnin = niter/2)
    my_params <- rownames(data.frame(post_sum[[1]]))
    
    pdf(file=paste(OTU,"_traceplot.pdf"))
    for (i in seq(from = 1, to = length(my_params),by = 4)){
      if (i < 5) {
        plotTrace(fitModel, c(my_params[i],my_params[i+1],
                              my_params[i+2],my_params[i+3]))
        } else {
          plotTrace(fitModel, c(my_params[i],my_params[i+1],
                                my_params[i+2], my_params[i+3]))
        }
      }
    dev.off()
    
    to_return <- list(fitModel, post_sum, post_sample, post_detect)
    return(to_return)
}

# many OTU models ####
allFits_occupancy <- llply(set4,
                           function(x) tryCatch(fitEDNAmodel(OTU=x, niter = 20000),
                                                error = function(t)
                                                  {cat("\n theta overparametrized\n\n")}),
                           .inform = T)
names(allFits_occupancy) <- set4
save(file="set4_models.Rdata", allFits_occupancy)

# this part of the script was run for all OTUs on malloy, start: 190422


model_ok <- !(unlist(llply(allFits_occupancy,
                                 function(x)
                                   all(x == "\n theta overparametrized\n\n"))))
sum(model_ok)

plot(1:59,
allFits_good$`HISEQ:267:CAJCDANXX:2:1101:9028:4557_CONS_SUB_SUB_CMP`[[3]]$mean["AR",])

plot(1:40,
final_data$aggr_replsum_otus[grep("AR", rownames(final_data$aggr_replsum_otus)),
                             "HISEQ:267:CAJCDANXX:2:1101:9028:4557_CONS_SUB_SUB_CMP"]
)


# 
# # models that could fit
allFits_good <- allFits_occupancy[model_ok]
# save(file="allFits_good.Rdata", allFits_good)
# load("allFits_good.Rdata")

# # community sampling ####
# # subset to occupancy period estimates
allFits_occupancy <- allFits_occupancy[grepl("alpha.original",
                                   allFits_occupancy$Param),]
# 
# #pull out year data
# allFits_occupancy$Year <- gsub("alpha.factor.yearGroups.","",allFits_occupancy$Param)
# 
# #get sds
# allFits_occupancy$SD1 <- (allFits_occupancy$X97.5.- allFits_occupancy$Mean)/1.96
# allFits_occupancy$SD2 <- (allFits_occupancy$Mean - allFits_occupancy$X2.5.)/1.96
# allFits_occupancy$SD <- ifelse(allFits_occupancy$SD1>allFits_occupancy$SD2,allFits_occupancy$SD1,allFits_occupancy$SD2)
# 
# #create a 1000 communities - the presence/absence of each species is drawn by its probability of occurence  
# myMatrix <- matrix(data=NA,nrow=nrow(allFits_occupancy),ncol=1000) 
# dim(myMatrix)
# 
# #random assign the P/A for each OTU based on its occurrence probabilty
# for(j in 1:ncol(myMatrix)){
#   for(i in 1:nrow(myMatrix)){
#     #myp <- runif(1,min=allFits_occupancy$X2.5.[i],max=allFits_occupancy$X97.5.[i])#unif distribution
#     myp <- rnorm(1,mean=allFits_occupancy$Mean[i],sd=allFits_occupancy$SD[i])#normal distribution
#     #back-transform
#     myp <- probitlink(myp, inverse=T)
#     myMatrix[i,j] <- rbinom(1,1,myp) # this changes probability into PA
#   }
# }
# 
# #community random matrix - each column is a simulation
# temp <- data.frame(myMatrix)
# 
# #add on year/OTU information
# temp <- cbind(allFits_occupancy[,c("OTU","Year")],temp)
# 
# #get species richness for each community simulation
# out <- ddply(temp,.(Year),function(x){
#   numcolwise(sum)(x)})
# #sum means add up the number of species present
# 
# #get mean and 95% CI across simulated communities
# meanRichness<-apply(out[,2:1001],1,median)
# lowerRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.025))
# upperRichness<-apply(out[,2:1001],1,function(x)quantile(x,0.975))
# 
# #combine all
# finalDF <- data.frame(Year=as.numeric(as.character(sort(unique(allFits_occupancy$Year)))),
#                       meanRichness,lowerRichness,upperRichness)
# 
# #plotting
# ggplot(finalDF)+
#   geom_line(aes(x=Year,y=meanRichness))+
#   geom_ribbon(aes(x=Year,ymin=lowerRichness,ymax=upperRichness),alpha=0.5)

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

#   #predicted (95%CI for each parameter)
#   output <- posteriorSummary(fitModel,burnin=burnin,mcError=T,outputSummary=T)
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
#   
#   #return it
#   return(output)
# }

# #round to nearest 5
# ourSiteSampleData$Year <-  5 * round(ourSiteSampleData$original_date_final/5)
# #dateSummary <- ddply(ourSiteSampleData,.(Year),summarise,
# #                     nuLakes=length(unique(Lake)))
# 
# #before 1900-1950, group by decade
# ourSiteSampleData$Decade <-  10 * round(ourSiteSampleData$original_date_final/10)
# 
# #pre 1900, 50 years?
# #ourSiteSampleData$Quarter <-  25 * round(ourSiteSampleData$original_date_final/25)
# ourSiteSampleData$Half <-  50 * round(ourSiteSampleData$original_date_final/50)
# 
# #pre 1700, pool all
# ourSiteSampleData$Century <-  100 * round(ourSiteSampleData$original_date_final/100)

# #formatting time periods ####
# ourSiteSampleData$yearGroups <- NA
# ourSiteSampleData$yearGroups[ourSiteSampleData$date_final<1950] <- 1
# ourSiteSampleData$yearGroups[ourSiteSampleData$date_final>=1950 &
#                                ourSiteSampleData$date_final<1990] <- 2
# ourSiteSampleData$yearGroups[ourSiteSampleData$date_final>=1990] <- 3
# ourSiteSampleData$yearGroups <- factor(ourSiteSampleData$yearGroups,
#                                        levels = c("1","2","3"))

#subset data frame for when we have OTU data ####
# ourSiteSampleData <- subset(ourSiteSampleData,dna %in% this_otu$dna)
# #centre/scale the date
# ourSiteSampleData$original_date_final <- ourSiteSampleData$date_final
# ourSiteSampleData$date_final <- scale(ourSiteSampleData$date_final)

# fitModel <- occModel(formulaSite = ~ 1, #lake occurence probabiliy (lake area?)
#                      formulaSiteAndSample = ~ yearGroups -1, #lake-time occurrence probability
#                      formulaReplicate = ~ original_date_final + upper, #detection probability at each lake-time
#                      detectionMats=formatted4Model,
#                      siteColName="Lake",
#                      sampleColName="level2",
#                      siteAndSampleData = ourSiteSampleData,
#                      niter=niter)


# OTUlist <- c("HISEQ:267:CAJCDANXX:2:1101:4969:25485_CONS_SUB_SUB",
# #"HISEQ:267:CAJCDANXX:2:1111:2773:5692_CONS_SUB_SUB",
# "HISEQ:267:CAJCDANXX:2:1105:10506:33069_CONS_SUB_SUB_CMP",
# #"HISEQ:267:CAJCDANXX:2:1101:16250:29354_CONS_SUB_SUB_CMP",
# "HISEQ:267:CAJCDANXX:2:1101:13278:2098_CONS_SUB_SUB")


# # fit models to many OTUs
# allFits_occupancy <- llply(OTUlist,
#                            function(x) fitEDNAmodel(OTU=x, niter=100),
#                            .inform = T)
# names(allFits_occupancy) <- unlist(OTUlist)


