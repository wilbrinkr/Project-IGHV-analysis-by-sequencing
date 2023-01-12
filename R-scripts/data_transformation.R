# General information ----
# Title: Data transformation
# Author: Rick Wilbrink
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: r.wilbrink01@umcg.nl
# Collaboration: please ask permission from the author before using this script

# importing libraries ----
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(gridExtra))

# reading data ----
setwd("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones")
db <- read_rearrangement(file.path("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones/allclones_final_combined.tsv"))
colnames(db) # show the column names in the database

#assign values to cell subsets, use function grepl(test, yes, no)
db$subset = ifelse(grepl("totb",db$sample), "Total B cells" , 0)
db$subset = ifelse(grepl("CD27m",db$sample), "CD27negativeCD38lowCD21low" , db$subset)
db$subset = ifelse(grepl("CD27p",db$sample), "CD27positiveCD38lowCD21low" , db$subset)
db$subset = ifelse(grepl("plas",db$sample), "plasmablasts" , db$subset)

db$group = ifelse(grepl("AS",db$sample), "AS" , 0)
db$group = ifelse(grepl("HD",db$sample), "HD" , db$group)
db$group = ifelse(grepl("pSS",db$sample), "pSS" , db$group)

# subset data
#Match <- c('AS6B',' AS13B', 'AS26', 'HD1','HD2','HD3','pSS80','pSS87','pSS115')
#Match <- c('AS27',' AS28', 'AS32', 'HD4','HD6','HD7','pSS131','pSS136','pSS145')
#db <- db[grepl(paste(Match, collapse='|'), db$sample),]

#db <- db[!grepl('CD27m', db$sample),]

rm(list=setdiff(ls(), "db"))

# saving the database
write.table(db, file = paste0("allclones_final_transformed.tsv"), row.names=FALSE, sep="\t")
