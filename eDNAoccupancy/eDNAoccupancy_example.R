library(eDNAoccupancy)

setwd("~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/workdir/")
rm(list=ls())
load(file="~/Google Drive/Lake diversity/Lake_diversity_Prepared/Analyses Lake Div/Data/grenoble_experiment/output_final_dataset/final_data.Rdata")

# example OTU
most_common_otu <- data.frame(most_common = final_data$final_otus[,"HISEQ:267:CAJCDANXX:2:1101:13278:2098_CONS_SUB_SUB"])
rownames(most_common_otu) <-rownames(final_data$final_otus)
