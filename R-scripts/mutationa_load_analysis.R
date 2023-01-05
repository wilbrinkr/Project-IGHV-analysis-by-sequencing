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

# Load data  
db <- read_rearrangement(file.path("D:","Shared folder","share","databases","final","AS13B-CD27m-U_S161.tsv"))
colnames(db) # show the column names in the database

#assign values to cell subsets, use function grepl(test, yes, no) 
db$group = ifelse(grepl("CD27m",db$sample), "CD27negativeCD38lowCD21low" , 0)
db$group = ifelse(grepl("CD27p",db$sample), "CD27positiveCD38lowCD21low" , db$group)
db$group = ifelse(grepl("plas",db$sample), "plasmablasts" , db$group)

# Calculate R and S mutation counts
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=FALSE, 
                            nproc=1)

# Show new mutation count columns
db_obs %>% 
  select(sequence_id, starts_with("mu_count_")) %>%
  head(n=4)

# Calculate R and S mutation frequencies
db_obs <- observedMutations(db_obs, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            nproc=1)

# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)

# Calculate combined R and S mutation frequencies
db_obs <- observedMutations(db, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment_d_mask",
                            regionDefinition=NULL,
                            frequency=TRUE, 
                            combine=TRUE,
                            nproc=1)

# Show new mutation frequency columns
db_obs %>% 
  select(sequence_id, starts_with("mu_freq_")) %>%
  head(n=4)



cols = c('blue', 'green', 'purple')
# creating a ggplot 
g1 <- ggplot(db_obs, aes(x=group, y=mu_freq, fill=group)) +
  theme_bw() + ggtitle("Total mutations") +
  xlab("Subset") + ylab("Mutation frequency") +
  scale_fill_manual(name="groups", values=cols) +
  geom_boxplot()
plot(g1)
