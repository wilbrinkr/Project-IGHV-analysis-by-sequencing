# General information ----
# Title: construction of clonal lineage trees
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
suppressPackageStartupMessages(library(dowser))
library(ggpubr)
library(gridExtra)
library(dplyr)

setwd(setwd("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones"))

# reading data ----
clones <- read_rearrangement(file.path("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones/allclones_final_combined.tsv"))
colnames(clones) # show the column names in the database

#assign values to cell subsets, use function grepl(test, yes, no)
clones$subset = ifelse(grepl("totalB",clones$sample), "Total B cells" , 0)
clones$subset = ifelse(grepl("CD27m",clones$sample), "CD27negativeCD38lowCD21low" , clones$subset)
clones$subset = ifelse(grepl("CD27p",clones$sample), "CD27positiveCD38lowCD21low" , clones$subset)
clones$subset = ifelse(grepl("plas",clones$sample), "plasmablasts" , clones$subset)

clones$group = ifelse(grepl("AS",clones$sample), "AS" , 0)
clones$group = ifelse(grepl("HD",clones$sample), "HD" , clones$group)
clones$group = ifelse(grepl("pSS",clones$sample), "pSS" , clones$group)

#subset sample containing shared clones
sample1 <- clones[grepl('AS35', clones$sample),]

# Building lineage trees

# subset
toMatch <- c('1220', '1228', '3160')
example <- sample1[sample1$clone_id %in% toMatch,]
example$subject_id = 'AS35'

# Process example data keeping samples from different times
# distinct, adding duplicate_count among collapsed sequences,
# and show the sample_id within each clone in the tibble.
newClones = formatClones(data = example, traits=c("sample", 'group', 'subset', "v_call"),minseq = 2, columns=c("subject_id"))

trees <- getTrees(newClones)

# simple tree plotting with ggtree R package with isotypes at tips
plots <- plotTrees(trees, tips="subset",tipsize=2)

# plot tree of largest clone
plots[[1]]








# Build trees using dnapars.
# exec here is set to dnapars position in the Docker image.
igphylm.trees = getTrees(newClones, build="igphyml", 
                  exec="/usr/local/share/igphyml/src/igphyml", nproc=1)

igphylm.plots <- plotTrees(trees, tips="subset",tipsize=2)
plots[[1]]

clones$parameters[[1]]$omega_mle










dup <- sample1[duplicated(sample1$clone_id),]
dupset <- dup[,c('clone_id', 'group', 'subset', 'sample')]




duphuh <- dupset %>% 
  group_by(clone_id, subset) %>%
  filter(n()>4)

View(duphuh)




# Load required packages
library(alakazam)
library(dowser)

# load example AIRR tsv data
data(ExampleAirr)

# subset data for this example
ExampleAirr = ExampleAirr[ExampleAirr$clone_id %in% c("3170", "3184"),]






# subset data
#Match <- c('AS6B',' AS13B', 'AS26', 'HD1','HD2','HD3','pSS80','pSS87','pSS115')
Match <- c('AS27',' AS28', 'AS32', 'HD4','HD6','HD7','pSS131','pSS136','pSS145')
#clones <- clones[grepl(paste(Match, collapse='|'), clones$sample),]

clones <- clones[!grepl('CD27m', clones$sample),]

#assign values to cell subsets, use function grepl(test, yes, no)
clones$subset = ifelse(grepl("CD27m",clones$sample), "CD27negativeCD38lowCD21low" , 0)
clones$subset = ifelse(grepl("CD27p",clones$sample), "CD27positiveCD38lowCD21low" , clones$subset)
clones$subset = ifelse(grepl("plas",clones$sample), "plasmablasts" , clones$subset)

clones$group = ifelse(grepl("AS",clones$sample), "AS" , 0)
clones$group = ifelse(grepl("HD",clones$sample), "HD" , clones$group)
clones$group = ifelse(grepl("pSS",clones$sample), "pSS" , clones$group)

# Calculate selection scores using the output from expectedMutations
baseline <- calcBaseline(clones, testStatistic="focused", 
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

testBaseline(grouped2, groupBy = 'group')

View(grouped2)
