# General information ----
# Title: mutational frequency
# Author: Rick Wilbrink
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: r.wilbrink01@umcg.nl
# Collaboration: please ask permission from the author before using this script

# importing libraries
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))
library(ggpubr)
library(gridExtra)

setwd(setwd("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones"))

# read data 
ACS <- read_rearrangement(file.path("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones/ACS.tsv"))
colnames(ACS) # show the column names in the database

# subset data
#Match <- c('AS6B',' AS13B', 'AS26', 'HD1','HD2','HD3','pSS80','pSS87','pSS115')
Match <- c('AS27',' AS28', 'AS32', 'HD4','HD6','HD7','pSS131','pSS136','pSS145')
#ACS <- ACS[grepl(paste(Match, collapse='|'), ACS$sample),]

ACS <- ACS[!grepl('CD27m', ACS$sample),]

#assign values to cell subsets, use function grepl(test, yes, no)
ACS$subset = ifelse(grepl("CD27m",ACS$sample), "CD27negativeCD38lowCD21low" , 0)
ACS$subset = ifelse(grepl("CD27p",ACS$sample), "CD27positiveCD38lowCD21low" , ACS$subset)
ACS$subset = ifelse(grepl("plas",ACS$sample), "plasmablasts" , ACS$subset)

ACS$group = ifelse(grepl("AS",ACS$sample), "AS" , 0)
ACS$group = ifelse(grepl("HD",ACS$sample), "HD" , ACS$group)
ACS$group = ifelse(grepl("pSS",ACS$sample), "pSS" , ACS$group)

# Calculate selection scores using the output from expectedMutations
baseline <- calcBaseline(ACS, testStatistic="focused", 
                         regionDefinition=IMGT_V, nproc=8)

rm(list=setdiff(ls(), "baseline"))

grouped2 <- groupBaseline(baseline, groupBy=c("subset", "group"))

#subset_colors <- c('CD27negativeCD38lowCD21low' = "#f3e26c", 'CD27positiveCD38lowCD21low' = "#7adc98", "plasmablasts" = '#dd95ea')
#subset_colors <- c('CD27negativeCD38lowCD21low' = 'firebrick', 'CD27positiveCD38lowCD21low' = "seagreen", "plasmablasts" = 'steelblue')
group_colors <- c('AS' = 'firebrick', 'HD' = "seagreen", "pSS" = 'steelblue')

plotBaselineDensity(grouped2, idColumn = "group", groupColumn = "subset",
                    colorValues=group_colors, sigmaLimits=c(-1.5, 1), size = 1.2)

plotBaselineDensity(grouped2, idColumn = "subset", groupColumn = "group",
                    colorValues=subset_colors, sigmaLimits=c(-1, 1), size = 1.2, facetBy = 'subset')

View(baseline)

test <- testBaseline(grouped2, groupBy = 'subset')
View(test)

##
##


testBaseline(grouped2, groupBy = 'group')

View(grouped2)
