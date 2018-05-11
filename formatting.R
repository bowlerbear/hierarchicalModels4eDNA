#formatting script

setwd("C:/Users/diana.bowler/OneDrive - NINA/eDNA")

#read in all OTUS
load("all_otus.RData")
dim(all_otus)
#1854 44003

#view a sample of the data
all_otus[1:10,1:10]
#columns are OTUS

#get info on assignement
load("all_assign.RData")
dim(all_assign)

#view a sample of the data
all_assign[1:10,1:8]

#look at possible categories
unique(all_assign$domain)# possible values: Bacteria  Archaea   Eukaryota <NA>
unique(all_assign$taxon)

#Miki says restrict the analysis to Eukaryotes
all_otus <- all_otus[,colnames(all_otus)%in%all_assign$otu[all_assign$domain=="Eukaryota"]]

#now pull out metadata on lake, replicate, depth etc...
metaDF<-data.frame(rows=row.names(all_otus))
metaDF$Lake <- sapply(metaDF$rows,function(x) strsplit(as.character(x),"\\.")[[1]][1])
metaDF$Lake <- sapply(metaDF$Lake,function(x) substr(x,1,2))
metaDF$Replicate <- sapply(metaDF$rows,function(x) strsplit(as.character(x),"_")[[1]][2])
metaDF$dna <- sapply(metaDF$rows,function(x) strsplit(as.character(x),"_")[[1]][1])
metaDF$Depth <- sapply(metaDF$dna,function(x) strsplit(as.character(x),"\\.")[[1]][2])

#reformat the data
otuSummary<-data.frame(metaDF,all_otus)
library(reshape2)
otuSummary<-melt(otuSummary,id=c("rows","Lake","Replicate","Depth","dna"))
names(otuSummary)[which(names(otuSummary)=="variable")]<-"OTU"
names(otuSummary)[which(names(otuSummary)=="value")]<-"Count"

#use OTUs seen in all lakes
library(plyr)
unique(metaDF$Lake)
#AR1" "BL1" "CR1" "FH3" "OR2" "PL1" "SL1" "ST8" "TF1" "WM1"
length(unique(metaDF$Lake))#10
otuSummaryLakes<-ddply(otuSummary,.(OTU),summarise,nuLakes=length(unique(Lake[Count>0])))
nrow(subset(otuSummaryLakes,nuLakes==10))#137
otuSummary<-subset(otuSummary,OTU%in%otuSummaryLakes$OTU[otuSummaryLakes$nuLakes==10])

#read in proxy data
load("proxies.Rdata")
dim(proxies)
proxies[1:10,1:10]

#From Mikie: proxies of interest are:
#The dates of the horizons are the date_final. The erosion proxy is the "Al", the metal is "Pb", the DDT is "ddt"
proxies<-proxies[,c("dna","code","depth","lake","date_final","Al","Pb","ddt")]
head(proxies)

unique(proxies$code)
#AR BL CR FH OR PL SL ST TF WM

#merge proxies and occurence data
otuSummary<-merge(otuSummary,proxies,by="dna")

#get occurence matrix
occMatrix<-acast(otuSummary,Lake+OTU+depth+date_final+Al+Pb+ddt~Replicate,value.var="Count")
occMatrix[1:10,1:6]

#turn into presence/absence
occMatrix[occMatrix>0]<-1

#set site covariates
varMatrix<-dcast(otuSummary,Lake+OTU+depth+date_final+Al+Pb+ddt~Replicate,value.var="Count")[,1:7]
#Just use the years from the date_final info
varMatrix$Year<-as.numeric(sapply(varMatrix$date_final,function(x)strsplit(as.character(x),"\\.")[[1]][1]))

#look at relationships among the covariates
library(GGally)
ggpairs(varMatrix[,5:8])

# Bundle and summarize data set for the BUGS model
bugs.data <- list(y = as.matrix(occMatrix),
                  M = nrow(occMatrix), 
                  J = ncol(occMatrix),
                  Year = varMatrix$Year,
                  Al = varMatrix$Al - median(varMatrix$Al), 
                  Pb = varMatrix$Pb - median(varMatrix$Pb),
                  ddt = varMatrix$ddt - median(varMatrix$ddt),
                  OTU = as.numeric(factor(varMatrix$OTU)),
                  Lake = as.numeric(factor(varMatrix$Lake)),
                  nLake = length(unique(varMatrix$Lake)),
                  nOTU = length(unique(varMatrix$OTU)))