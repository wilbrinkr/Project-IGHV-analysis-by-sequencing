# General information ----
# Title: clonal lineage tree dowser
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
suppressPackageStartupMessages(library(shazam))
library(dplyr)
library(readr)
library(data.table)
library(ggvenn)

# read data 
db <- read_tsv("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/clones/allclones_final_combined.tsv")
db1 <- read_tsv("D:/Onedrive/PhD/Projecten/project - IGHVseq/databases/final_combined_files/final_combined.tsv")
colnames(db) # show the column names in the database

#assign values to cell subsets, use function grepl(test, yes, no)
db$subset = ifelse(grepl("totalB",db$sample), "Total B cells" , 0)
db$subset = ifelse(grepl("CD27m",db$sample), "CD27negativeCD38lowCD21low" , db$subset)
db$subset = ifelse(grepl("CD27p",db$sample), "CD27positiveCD38lowCD21low" , db$subset)
db$subset = ifelse(grepl("plas",db$sample), "plasmablasts" , db$subset)

db$group = ifelse(grepl("AS",db$sample), "AS" , 0)
db$group = ifelse(grepl("HD",db$sample), "HD" , db$group)
db$group = ifelse(grepl("pSS",db$sample), "pSS" , db$group)

db <- db[!grepl('Total B cells', db$subset),]
#

# Add sample_id
db$sample_id = ifelse(grepl("AS6B-",db$sample), "AS6B" , 0)
db$sample_id = ifelse(grepl("AS13B-",db$sample), "AS13B" , db$sample_id)
db$sample_id = ifelse(grepl("AS26-",db$sample), "AS26" , db$sample_id)
db$sample_id = ifelse(grepl("AS27-",db$sample), "AS27" , db$sample_id)
db$sample_id = ifelse(grepl("AS28-",db$sample), "AS28" , db$sample_id)
db$sample_id = ifelse(grepl("AS32-",db$sample), "AS32" , db$sample_id)
db$sample_id = ifelse(grepl("AS35-",db$sample), "AS35" , db$sample_id)
db$sample_id = ifelse(grepl("AS36-",db$sample), "AS36" , db$sample_id)
db$sample_id = ifelse(grepl("AS38-",db$sample), "AS38" , db$sample_id)
db$sample_id = ifelse(grepl("AS39-",db$sample), "AS39" , db$sample_id)

db$sample_id = ifelse(grepl("HD1-",db$sample), "HD1" , db$sample_id)
db$sample_id = ifelse(grepl("HD2-",db$sample), "HD2" , db$sample_id)
db$sample_id = ifelse(grepl("HD3-",db$sample), "HD3" , db$sample_id)
db$sample_id = ifelse(grepl("HD4-",db$sample), "HD4" , db$sample_id)
db$sample_id = ifelse(grepl("HD6-",db$sample), "HD6" , db$sample_id)
db$sample_id = ifelse(grepl("HD7-",db$sample), "HD7" , db$sample_id)
db$sample_id = ifelse(grepl("HD8-",db$sample), "HD8" , db$sample_id)
db$sample_id = ifelse(grepl("HD9-",db$sample), "HD9" , db$sample_id)
db$sample_id = ifelse(grepl("HD10-",db$sample), "HD10" , db$sample_id)
db$sample_id = ifelse(grepl("HD11-",db$sample), "HD11" , db$sample_id)

db$sample_id = ifelse(grepl("pSS80-",db$sample), "pSS80" , db$sample_id)
db$sample_id = ifelse(grepl("pSS87-",db$sample), "pSS87" , db$sample_id)
db$sample_id = ifelse(grepl("pSS115-",db$sample), "pSS115" , db$sample_id)
db$sample_id = ifelse(grepl("pSS131-",db$sample), "pSS131" , db$sample_id)
db$sample_id = ifelse(grepl("pSS136-",db$sample), "pSS136" , db$sample_id)
db$sample_id = ifelse(grepl("pSS145-",db$sample), "pSS145" , db$sample_id)
db$sample_id = ifelse(grepl("pSS149-",db$sample), "pSS149" , db$sample_id)
db$sample_id = ifelse(grepl("pSS160-",db$sample), "pSS160" , db$sample_id)
db$sample_id = ifelse(grepl("pSS193-",db$sample), "pSS193" , db$sample_id)

rm(list=setdiff(ls(), "db"))

write.table(db, file = paste0("allclones_final_transformed.tsv"), row.names=FALSE, sep="\t")

# ---- 

setwd("D:/Onedrive/PhD/Projecten/project - IGHVseq/github-code/Plots/VennDiagram")

unique_samples <- unique(db$sample_id)

unique_samples

for (samples in unique_samples){

## filter
db.filtered <- db %>%
  subset(db$sample_id == samples) 

db.CD27m <- db.filtered %>%
  filter(subset == 'CD27negativeCD38lowCD21low') 

db.CD27p <- db.filtered %>%
  filter(subset == 'CD27positiveCD38lowCD21low') 

db.plas <- db.filtered %>%
  filter(subset == 'plasmablasts') 

df <- data.frame(matrix(ncol = 3, nrow = 30000))
colnames(df) <- c('CD27m', 'CD27p', 'plas')

df$CD27m <- ifelse(rownames(df) %in% as.character(db.CD27m$clone_id), paste0(rownames(df)), df$CD27m)
df$CD27p <- ifelse(rownames(df) %in% as.character(db.CD27p$clone_id), paste0(rownames(df)), df$CD27p)
df$plas <- ifelse(rownames(df) %in% as.character(db.plas$clone_id), paste0(rownames(df)), df$plas)

venn.input <- list('CD27m' = df$CD27m,
                   'CD27p' = df$CD27p,
                   'plas' = df$plas)

venn <- ggvenn(
            venn.input, 
            fill_color = c("#ffe86f",
                           "#64e5cb",
                           "#f784ff"),
            stroke_size = 0.5, set_name_size = 4)

tiff(paste0(samples, "_VennDiagram.tiff"), units="in", width=5, height=5, res=600, compression = 'lzw')
print(venn)
dev.off()

}

# ----

setwd("D:/Onedrive/PhD/Projecten/project - IGHVseq/github-code/data")