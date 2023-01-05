# General information ----
# Title: V genotyping and determening clonal thresholds
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

# creating for loop
files <- list.files(path="D:/Onedrive/PhD/R/Immcantation/results/databases/combined_combined", pattern="*.tsv", full.names=FALSE, recursive=FALSE)

for (file in files){

# set working directory
setwd("D:/Onedrive/PhD/R/Immcantation")

# Load data  
db <- read_rearrangement(file.path("D:/Onedrive/PhD/R/Immcantation/results/databases/combined_combined",paste0(file)))

# Clonal distance ----

# get the distance to nearest neighbors
db <- distToNearest(db, model = "ham", normalize = "len", vCallColumn = "v_call", nproc = 4)

# determine the threshold
threshold <- findThreshold(db$dist_nearest, method = "gmm")
thr <- round(threshold@threshold, 2)
thrg_plot <- plot(threshold, binwidth=0.01, title="Distance to nearest neighbors (Hamming) - gmm") # plot the distribution

# set and create directory
setwd("D:/Onedrive/PhD/R/Immcantation/results/clonal_thresholds/plots")

tiff(paste0("thresholds_gmm", file,".tiff"), units="in", width=8, height=8, res=300, compression = 'lzw')
print(thrg_plot)
dev.off()

}