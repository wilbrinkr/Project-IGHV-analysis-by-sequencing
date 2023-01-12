# General information ----
# Title: mutational frequency
# Author: Rick Wilbrink
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: r.wilbrink01@umcg.nl
# Collaboration: please ask permission from the author before using this script

# importing libraries ----
suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(scoper))
library(ggpubr)
library(gridExtra)

# reading data ----
setwd("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones")
db <- read_rearrangement(file.path("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones/allclones_final_transformed.tsv"))
colnames(db) # show the column names in the database

db <- db[!grepl('totalB', db$sample),]

# Collapse clonal groups into single sequences
clones <- collapseClones(db, cloneColumn="clone_id", 
                         sequenceColumn="sequence_alignment", 
                         germlineColumn="germline_alignment_d_mask", 
                         regionDefinition=IMGT_V, 
                         method="thresholdedFreq", minimumFrequency=0.6,
                         includeAmbiguous=FALSE, breakTiesStochastic=FALSE, 
                         nproc=4)

View(db)

# Count observed mutations and append mu_count columns to the output
observed <- observedMutations(clones, 
                              sequenceColumn="clonal_sequence",
                              germlineColumn="clonal_germline",
                              regionDefinition=IMGT_V, nproc=4)

# Count expected mutations and append mu_exptected columns to the output
expected <- expectedMutations(observed, 
                              sequenceColumn="clonal_sequence",
                              germlineColumn="clonal_germline",
                              targetingModel=HH_S5F,
                              regionDefinition=IMGT_V, nproc=4)

# subset data
#Match <- c('AS6B',' AS13B', 'AS26', 'HD1','HD2','HD3','pSS80','pSS87','pSS115')
#Match <- c('AS27',' AS28', 'AS32', 'HD4','HD6','HD7','pSS131','pSS136','pSS145')
#db <- db[grepl(paste(Match, collapse='|'), db$sample),]

CD27m <- db1[!grepl('CD27m', db$sample),]

names <- c('TotalB', "CD27m", "CD27p", 'plas')

for (name in names){

db <- read_rearrangement(file.path("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones/expected.tsv"))

subset <- db[grepl(name, db$sample),]

# Calculate selection scores using the output from expectedMutations
baseline <- calcBaseline(subset, testStatistic="focused", 
                         regionDefinition=IMGT_V, nproc=8)

rm(list=setdiff(ls(), "baseline"))

grouped <- groupBaseline(baseline, groupBy=c("subset", "group"))

saveRDS(grouped, file = paste0(name,".RDS")) 
}


#subset_colors <- c('CD27negativeCD38lowCD21low' = "#f3e26c", 'CD27positiveCD38lowCD21low' = "#7adc98", "plasmablasts" = '#dd95ea')
#subset_colors <- c('CD27negativeCD38lowCD21low' = 'firebrick', 'CD27positiveCD38lowCD21low' = "seagreen", "plasmablasts" = 'steelblue')
group_colors <- c('AS' = 'firebrick', 'HD' = "seagreen", "pSS" = 'steelblue')

plotBaselineDensity(grouped2, idColumn = "group", groupColumn = "subset",
                    colorValues=group_colors, sigmaLimits=c(-1.5, 1), size = 1.2)

plotBaselineDensity(grouped2, idColumn = "subset", groupColumn = "group",
                    colorValues=subset_colors, sigmaLimits=c(-1, 1), size = 1.2, facetBy = 'subset')

test <- testBaseline(grouped2, groupBy = 'subset')
test
