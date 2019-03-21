library(eDNAoccupancy)
library(tidyverse)
library(plyr) 

setwd("~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/workdir/")
rm(list=ls())
load(file="~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/Data/grenoble_experiment/output_final_dataset/final_data.Rdata")

#separate the sample columns 
final_data$Rep<-sapply(row.names(final_data),function(x)strsplit(x,"_")[[1]][2]) 
final_data$Lake<-sapply(row.names(final_data),function(x)strsplit(x,"\\.")[[1]][1]) 
final_data$Depth<-sapply(row.names(final_data),function(x)strsplit(x,"\\.")[[1]][2]) 
final_data$Depth<-as.numeric(sapply(final_data$Depth,function(x)strsplit(x,"_")[[1]][1])) 

# example OTU 
most_common_otu <- data.frame(final_data$final_otus[,"HISEQ:267:CAJCDANXX:2:1101:13278:2098_CONS_SUB_SUB"]) 
rownames(most_common_otu) <-rownames(final_data$final_otus) 
names(most_common_otu) <- "nuReeds" 

#examine distribution of depths across lakes 

depthSummary <- ddply(most_common_otu,.(Depth),summarise, 
                      nuLakes=length(unique(Lake)), 
                      nuReps = length(unique(Rep))) 

#format without date groups 
library(reshape2) 
otuMatrix <- acast(most_common_otu,Lake~Depth,fun=function(x)sum(x), 
                   value.var="nuReeds") 


# link dates to depths 
most_common_otu$Horiz <- sapply(rownames(most_common_otu), 
                                function(x) strsplit(x,"_")[[1]][1]) 

library(magrittr) 
library(dplyr) 
most_common_otu <- 
  most_common_otu %>% 
  left_join(final_data$proxies, by = c("Horiz" = "dna")) %>% 
  select(nuReeds,Rep,Lake,Depth, date_final) 


#format without date groups 
library(reshape2) 
otuMatrix <- acast(most_common_otu,Lake~date_final, 
                   fun=sum, 
                   fill=NA_real_, 
                   value.var="nuReeds") 

#how many reps do we have per date 
library(plyr) 
dateSummary <- ddply(most_common_otu,.(date_final),summarise, 
                     nuLakes=length(unique(Lake))) 
dateSummary 

#per year 
most_common_otu$Year <- round(most_common_otu$date_final) 
dateSummary <- ddply(most_common_otu,.(Year),summarise, 
                     nuLakes=length(unique(Lake))) 
dateSummary 

#per 2-year period 
most_common_otu$Year <-  2 * round(most_common_otu$Year/2) 
dateSummary <- ddply(most_common_otu,.(Year),summarise, 
                     nuLakes=length(unique(Lake))) 

library(reshape2) 
otuMatrix <- acast(most_common_otu,Lake~Year, 
                   fun=sum, 
                   fill=NA_real_, 
                   value.var="nuReeds") 

#round to nearest 5 
most_common_otu$Year <-  5 * round(most_common_otu$date_final/5) 
dateSummary <- ddply(most_common_otu,.(Year),summarise, 
                     nuLakes=length(unique(Lake))) 


#before 1900-1950, group by decade 
most_common_otu$Decade <-  10 * round(most_common_otu$date_final/10) 

#pre 1900, 50 years? 
most_common_otu$Quarter <-  25 * round(most_common_otu$date_final/25) 
most_common_otu$Half <-  50 * round(most_common_otu$date_final/50) 

#pre 1700, pool all 
most_common_otu$Century <-  100 * round(most_common_otu$date_final/100) 


#combine groups 
most_common_otu$yearGroups <- NA 
most_common_otu$yearGroups[most_common_otu$date_final<1700] <- 1700 
most_common_otu$yearGroups[most_common_otu$date_final>1700] <- most_common_otu$Half[most_common_otu$date_final>1700] 
most_common_otu$yearGroups[most_common_otu$date_final>1900] <- most_common_otu$Decade[most_common_otu$date_final>1900] 
most_common_otu$yearGroups[most_common_otu$date_final>1950] <- most_common_otu$Year[most_common_otu$date_final>1950] 


#get M x J matrix 

#this is y in our analysis 
#change into Reed pres/abs 

most_common_otu$PA <- ifelse(most_common_otu$nuReeds>0,1,0) 

#aggregate data into year groups 
most_common_otu_Ag <- ddply(most_common_otu,c(yearGroups,Lake), 
                            summarise=) 

otuMatrix <- acast(most_common_otu,Lake~yearGroups, 
                   fun=sum, 
                   fill=NA_real_, 
                   value.var="PA") 


#get data in the format for eDNA occupancy detection function 
head(most_common_otu) 

#check the number of years into each year groups -Miki to check 
most_common_otu_F <- most_common_otu[,c("Lake","date","Rep","PA")] 
most_common_otu_Cast <- dcast(most_common_otu_F,Lake+yearGroups~Rep,value.var="PA",fun=sum) 


#recode depth levels 
most_common_otu <- ddply(most_common_otu,.(Lake),function(x){ 
  x$level2 <- as.numeric(factor(x$date_final)) 
  return(x) 
}) 
most_common_otu_F <- most_common_otu[,c("Lake","level2","Rep","PA")] 
most_common_otu_Cast <- dcast(most_common_otu_F,Lake+level2~Rep,value.var="PA",fun=sum) 

dataMatrix <- occData(most_common_otu_Cast,'Lake','level2') 

d <- most_common_otu_Cast 
siteID = unique(d[, "Lake"]) 
ind = d[, "Lake"] == siteID[1] 
sampleNum = d[ind, "level2"] 
sort(sampleNum), 1:length(sampleNum)) 

# example OTU
most_common_otu <- data.frame(most_common = final_data$final_otus[,"HISEQ:267:CAJCDANXX:2:1101:13278:2098_CONS_SUB_SUB"])
rownames(most_common_otu) <-rownames(final_data$final_otus)

# link dates to depths
most_common_otu$Horiz <- sapply(rownames(most_common_otu), 
                                function(x) strsplit(x,"_")[[1]][1])

most_common_otu %>%
  left_join(final_data$proxies, by = c("Horiz" = "dna")) %>%
  select(most_common, date_final)

### example with gobies
str(gobyDetectionData)
gobyDetection <- occData(gobyDetectionData, siteColName = "site",
                         sampleColName = "sample")

# get site data
final_data$proxies %>%
  select(dna, lat, lon, date_final, prod_min, prod_max, Pb, DDT, S_Fe, erosion) ->
  ourSiteSampleData

