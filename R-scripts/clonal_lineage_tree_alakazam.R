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
suppressPackageStartupMessages(library(dowser))
library(dplyr)
library(readr)

# read data 
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

# biggest clone
View(table(clones$clone_id))

View(clones)

#subset sample containing shared clones
sample1 <- clones[grepl('AS35', clones$sample),]

# Building lineage trees

View(sample1)

table(sample1$clone_id)

x <- as.data.frame(table(sample$clone_id)); x[order(x$Freq, decreasing = TRUE), ]
x

View(table(sample1$clone_id))



library(dplyr)
sample1 %>% 
  group_by(clone_id) %>% 
  summarise_all(sum)

# subset
toMatch <- c('1220', '1228', '3160')
ToMatch1 <- '1220'
example <- sample1[sample1$clone_id %in% ToMatch1,]
example$subject_id = 'AS35'

# trees alakazam
# Load required packages
library(alakazam)
library(igraph)
library(dplyr)

# This example data set does not have ragged ends
# Preprocess clone without ragged end masking (default)
clone <- makeChangeoClone(example, text_fields=c("subject_id", "subset"))

# Show combined annotations
clone@data[, c("subject_id", "subset")]

View(clone@data)

View(clone)

phylip_exec <- 'C:/Users/rickw/Desktop/phylip-3.698/exe/dnapars.exe'
phylip_exec <- 
  
graph <- buildPhylipLineage(clone, phylip_exec = 'C:/Users/rickw/Desktop/phylip-3.698/exe/dnapars.exe', rm_temp=TRUE)

# Modify graph and plot attributes
V(graph)$color <- "steelblue"
V(graph)$color[V(graph)$name == "Germline"] <- "black"
V(graph)$color[grepl("Inferred", V(graph)$name)] <- "white"
V(graph)$label <- V(graph)$subset
E(graph)$label <- ""

# Remove large default margins
par(mar=c(0, 0, 0, 0) + 0.1)
# Plot graph
plot(graph, layout=layout_as_tree, edge.arrow.mode=0, vertex.frame.color="black",
     vertex.label.color="black", vertex.size=40)
# Add legend
legend("topleft", c("Germline", "Inferred", "Sample"), 
       fill=c("black", "white", "steelblue"), cex=0.75)
