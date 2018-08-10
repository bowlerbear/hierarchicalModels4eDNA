#formatting script

# setwd("C:/Users/diana.bowler/OneDrive - NINA/eDNA")

#read in all OTUS
load("/Users/noemikrisztinanagy/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses/grenoble_experiment/data_prepare/all_otus.Rdata")
dim(all_otus)
#1854 44003

#view a sample of the data
all_otus[1:10,1:10]
#columns are OTUS

#get info on assignement
load("/Users/noemikrisztinanagy/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses/grenoble_experiment/data_prepare/all_assign.RData")
dim(all_assign)

#view a sample of the data
all_assign[1:10,1:8]

#look at possible categories
unique(all_assign$domain)# possible values: Bacteria  Archaea   Eukaryota <NA>
unique(all_assign$taxon)

#Miki says restrict the analysis to Eukaryotes
all_otus <- all_otus[,colnames(all_otus)%in%all_assign$otu[all_assign$domain=="Eukaryota"]]

# count OTU occurences in replicates, for the total data, and in each lake
## count OTU occurences in all replicates
replicate_data <- data.frame(otu = colnames(all_otus),
                             replicate_count =
                               specnumber(t(all_otus))
)

## count OTU occurences in replicates of each lake
my_lakes <-
  levels(as.factor(sapply(rownames(all_otus),
                          function(x)strsplit(x,"\\.")[[1]][1])))
n <- ncol(replicate_data)
for (i in 1:length(my_lakes)) {
  n <- n+1
  my_set <- all_otus[grep(my_lakes[i], rownames(all_otus)),]
  my_replicates <- specnumber(t(my_set))
  replicate_data <- cbind(replicate_data,
                          my_name = my_replicates)
  colnames(replicate_data)[n] <-
    paste(my_lakes[i], "repli_count", sep="_")
}
rm(list = c("i", "my_set", "my_lakes", "my_replicates", "n"))

# create a vector with names of OTUs frequent in stechlin
# and filter all_otus with this vector
replicate_data %>%
  filter(ST8_repli_count >= 12) %>%
  select(otu) ->
  otus_frequent_ST_names
otus_frequent_ST <- all_otus[,names(all_otus) %in% 
                               otus_frequent_ST_names$otu]
rm(otus_frequent_ST_names)

#now pull out metadata on lake, replicate, depth etc...
metaDF<-data.frame(rows=row.names(all_otus))
metaDF$Lake <- sapply(metaDF$rows,function(x) strsplit(as.character(x),"\\.")[[1]][1])
metaDF$Lake <- sapply(metaDF$Lake,function(x) substr(x,1,2))
metaDF$Replicate <- sapply(metaDF$rows,function(x) strsplit(as.character(x),"_")[[1]][2])
metaDF$dna <- sapply(metaDF$rows,function(x) strsplit(as.character(x),"_")[[1]][1])
metaDF$Depth <- sapply(metaDF$dna,function(x) strsplit(as.character(x),"\\.")[[1]][2])

# keep only one lake
otus_frequent_ST <- data.frame(metaDF, otus_frequent_ST)
otus_frequent_ST <- filter(otus_frequent_ST, Lake == "ST")
## lecsekkol
summary(apply(otus_frequent_ST[,6:298], 2, sum))
summary(specnumber(t(otus_frequent_ST[,6:298])))

#reformat the data
library(reshape2)
otuSummary_ST<-melt(otus_frequent_ST,id=c("rows","Lake","Replicate","Depth","dna"))
names(otuSummary_ST)[which(names(otuSummary_ST)=="variable")]<-"OTU"
names(otuSummary_ST)[which(names(otuSummary_ST)=="value")]<-"Count"

#read in proxy data
load("/Users/noemikrisztinanagy/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses/grenoble_experiment/multiproxy/proxies.Rdata")
dim(proxies)
proxies[1:10,1:10]

#From Mikie: proxies of interest are:
#The dates of the horizons are the date_final. The erosion proxy is the "Al", the metal is "Pb", the DDT is "ddt"
proxies<-proxies[,c("dna","code","depth","lake","date_final", "date_trusted","Al","Pb","ddt")]
head(proxies)

unique(proxies$code)
#AR BL CR FH OR PL SL ST TF WM

#merge proxies and occurence data
otuSummary_ST<-merge(otuSummary_ST,proxies,by="dna")

#get occurence matrix
occMatrix_ST<-acast(otuSummary_ST,Lake+OTU+depth+date_final+date_trusted+Al+Pb+ddt~Replicate,
                    value.var="Count")
View(occMatrix_ST[1:10,1:6])

#turn into presence/absence
occMatrix_ST[occMatrix_ST>0]<-1

#set site covariates
varMatrix_ST<-dcast(otuSummary_ST,Lake+OTU+date_trusted+depth+date_final+Al+Pb+ddt~Replicate,
                 value.var="Count")[,1:8]
#Just use the years from the date_final info
varMatrix_ST$Year<-as.numeric(sapply(varMatrix_ST$date_final,
                                     function(x)strsplit(as.character(x),
                                                         "\\.")[[1]][1]))

#look at relationships among the covariates
library(GGally)
ggpairs(varMatrix_ST[,4:8])

# Bundle and summarize data set for the BUGS model
bugs.data.ST <- list(y = as.matrix(occMatrix_ST),
                  M = nrow(occMatrix_ST), 
                  J = ncol(occMatrix_ST),
                  Year = varMatrix_ST$Year,
                  Al = varMatrix_ST$Al - median(varMatrix_ST$Al), 
                  Pb = varMatrix_ST$Pb - median(varMatrix_ST$Pb),
                  ddt = varMatrix_ST$ddt - median(varMatrix_ST$ddt),
                  OTU = as.numeric(factor(varMatrix_ST$OTU)),
                  nOTU = length(unique(varMatrix_ST$OTU)))

save(file="bugs.data.ST.Rdata", bugs.data.ST)
