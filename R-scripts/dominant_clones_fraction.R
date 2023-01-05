#### fraction dominant clones per group ####

#### loading packages ####
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(rstatix)

#### import and load data #### 
df <- read.delim("D:/PhD/Projecten/VH-gen project/Dataset/9-5-2022 Combi 3 runs/vjcdr3-clones-mut-Wilbrink-IGH_HUMAN.csv")

## remove clones with a freq of 1, thus singletons that are not clones
#df[df$freq != 1, ] 

#assign values to cell subsets, use function grepl(test, yes, no) 
df$group = ifelse(grepl("totalB",df$Sample), "total B cells" , 0)
df$group = ifelse(grepl("CD27m",df$Sample), "CD27negativeCD38lowCD21low" , df$group)
df$group = ifelse(grepl("CD27p",df$Sample), "CD27positiveCD38lowCD21low" , df$group)
df$group = ifelse(grepl("plas",df$Sample), "plasmablasts" , df$group)

#### filter the subets you want to compare ####
group.name = 'plas'

# Clones were considered dominant when the number of clonally related sequences was
fraction = 0.05  

#### AS patients ####
df2 <- filter(df, grepl('AS', Sample))
df3 <- filter(df2, grepl(group.name , Sample))
df3 <- transform(df3, sample_id = match(Sample, unique(Sample)))

#### filter for every unique sample and apply function ####
n1 <- filter(df3, sample_id == 1)
n1$freq.fraction <- (n1$freq/sum(n1$freq))*100

n2 <- filter(df3, sample_id == 2)
n2$freq.fraction <- (n2$freq/sum(n2$freq))*100

n3 <- filter(df3, sample_id == 3)
n3$freq.fraction <- (n3$freq/sum(n3$freq))*100

n4 <- filter(df3, sample_id == 4)
n4$freq.fraction <- (n4$freq/sum(n4$freq))*100

n5 <- filter(df3, sample_id == 5)
n5$freq.fraction <- (n5$freq/sum(n5$freq))*100

n6 <- filter(df3, sample_id == 6)
n6$freq.fraction <- (n6$freq/sum(n6$freq))*100

n7 <- filter(df3, sample_id == 7)
n7$freq.fraction <- (n7$freq/sum(n7$freq))*100

n8 <- filter(df3, sample_id == 8)
n8$freq.fraction <- (n8$freq/sum(n8$freq))*100

n9 <- filter(df3, sample_id == 9)
n9$freq.fraction <- (n9$freq/sum(n9$freq))*100

n10 <- filter(df3, sample_id == 10)
n10$freq.fraction <- (n10$freq/sum(n1$freq))*100

## bind rows 
df4 <- bind_rows(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)

## what clones are dominant 

df4$clone.type <- ifelse(df4$freq.fraction >= fraction, "dominant", "not dominant")

## what percentage of clones are dominant per sample

clone.type.sample_id <- table(df4$sample_id, df4$clone.type)
df5 <- round(prop.table(clone.type.sample_id,1)*100,digits=2)

##### create data frames ####
df6 <- apply(as.matrix.noquote(df5),2,as.numeric)
df.AS <- as.data.frame(df6)

## create variable group
df.AS$group <- rep("AS", 10)







#### pSS patients ####
df2 <- filter(df, grepl('pSS', Sample))
df3 <- filter(df2, grepl(group.name , Sample))
df3 <- transform(df3, sample_id = match(Sample, unique(Sample)))

#### filter for every unique sample and apply function ####
n1 <- filter(df3, sample_id == 1)
n1$freq.fraction <- (n1$freq/sum(n1$freq))*100

n2 <- filter(df3, sample_id == 2)
n2$freq.fraction <- (n2$freq/sum(n2$freq))*100

n3 <- filter(df3, sample_id == 3)
n3$freq.fraction <- (n3$freq/sum(n3$freq))*100

n4 <- filter(df3, sample_id == 4)
n4$freq.fraction <- (n4$freq/sum(n4$freq))*100

n5 <- filter(df3, sample_id == 5)
n5$freq.fraction <- (n5$freq/sum(n5$freq))*100

n6 <- filter(df3, sample_id == 6)
n6$freq.fraction <- (n6$freq/sum(n6$freq))*100

n7 <- filter(df3, sample_id == 7)
n7$freq.fraction <- (n7$freq/sum(n7$freq))*100

n8 <- filter(df3, sample_id == 8)
n8$freq.fraction <- (n8$freq/sum(n8$freq))*100

n9 <- filter(df3, sample_id == 9)
n9$freq.fraction <- (n9$freq/sum(n9$freq))*100

n10 <- filter(df3, sample_id == 10)
n10$freq.fraction <- (n10$freq/sum(n1$freq))*100

## bind rows 
df4 <- bind_rows(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)

## what clones are dominant 

df4$clone.type <- ifelse(df4$freq.fraction >= fraction, "dominant", "not dominant")

## what percentage of clones are dominant per sample

clone.type.sample_id <- table(df4$sample_id, df4$clone.type)
df5 <- round(prop.table(clone.type.sample_id,1)*100,digits=2)

##### create data frames ####
df6 <- apply(as.matrix.noquote(df5),2,as.numeric)
df.pSS <- as.data.frame(df6)

## create variable group
df.pSS$group <- rep("pSS", 10)



#### AS patients ####
df2 <- filter(df, grepl('HD', Sample))
df3 <- filter(df2, grepl(group.name , Sample))
df3 <- transform(df3, sample_id = match(Sample, unique(Sample)))

#### filter for every unique sample and apply function ####
n1 <- filter(df3, sample_id == 1)
n1$freq.fraction <- (n1$freq/sum(n1$freq))*100

n2 <- filter(df3, sample_id == 2)
n2$freq.fraction <- (n2$freq/sum(n2$freq))*100

n3 <- filter(df3, sample_id == 3)
n3$freq.fraction <- (n3$freq/sum(n3$freq))*100

n4 <- filter(df3, sample_id == 4)
n4$freq.fraction <- (n4$freq/sum(n4$freq))*100

n5 <- filter(df3, sample_id == 5)
n5$freq.fraction <- (n5$freq/sum(n5$freq))*100

n6 <- filter(df3, sample_id == 6)
n6$freq.fraction <- (n6$freq/sum(n6$freq))*100

n7 <- filter(df3, sample_id == 7)
n7$freq.fraction <- (n7$freq/sum(n7$freq))*100

n8 <- filter(df3, sample_id == 8)
n8$freq.fraction <- (n8$freq/sum(n8$freq))*100

n9 <- filter(df3, sample_id == 9)
n9$freq.fraction <- (n9$freq/sum(n9$freq))*100

n10 <- filter(df3, sample_id == 10)
n10$freq.fraction <- (n10$freq/sum(n1$freq))*100

## bind rows 
df4 <- bind_rows(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10)

## what clones are dominant 

df4$clone.type <- ifelse(df4$freq.fraction >= fraction, "dominant", "not dominant")

## what percentage of clones are dominant per sample

clone.type.sample_id <- table(df4$sample_id, df4$clone.type)
df5 <- round(prop.table(clone.type.sample_id,1)*100,digits=2)

##### create data frames ####
df6 <- apply(as.matrix.noquote(df5),2,as.numeric)
df.HD <- as.data.frame(df6)

## create variable group
df.HD$group <- rep("HD", 10)

df.merged <- rbind(df.AS, df.HD, df.pSS)



## Mann-Whitney U test, comparing the three different groups
wth <- df.merged %>%
  pairwise_wilcox_test(dominant ~ group, p.adjust.method = "holm")
wth <- wth %>% add_xy_position(x = "group")
wth


## ggplot
facet.plot <- ggplot(df.merged, aes(group, dominant, color = group)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter() + 
  xlab("Groups") + ylab("% of dominant clones") + 
  ggtitle("Plasmablasts") + 
  theme(
    # axis elements  
    axis.text = element_text(colour = NULL, size = 15, margin = margin(t = 0, r = 0, b = 10, l = 0)),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(colour = NULL, size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(colour = NULL, size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    # strip elements
    strip.text = element_text(face = NULL, size = rel(1.2)),
    strip.background = element_rect(fill = "white", size = 1),
    # legend
    legend.position = "none"
  )
facet.plot



#Visualize data with ggplot 
k <- ggplot(df.merged, aes(group, dominant, color = group)) +
  
  geom_boxplot(outlier.shape = NA
               
  ) +
  
  geom_jitter(
    width = 0.15, # points spread out over 15% of available width
    height = 0, # do not move position on the y-axis
    alpha = 0.5,
    size = 4
  ) +
  
  ggtitle("Plasmablasts") +
  
  xlab(NULL) + ylab("% Dominant clones") +
  
  theme( # theme of the figure 
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(colour = "black", size = 17, margin = margin(t = 0, r = 0, b = 10, l = 0)),
    axis.title.x = element_text(colour = "black", size = 20),
    axis.title.y = element_text(colour = "black", size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    
    legend.position = "none", # legend theme
    
    panel.border = element_blank(), # panel theme
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", size = 1.1)
  )

k

##statistics
wth <- df.merged %>%
  pairwise_wilcox_test(dominant ~ group, p.adjust.method = "holm")
wth <- wth %>% add_xy_position(x = "group")
wth


#### create ggplot with facet_wrap ####
facet.plot <- ggplot(df.merged, aes(group, dominant, color = group)) +
  facet_wrap(~condition, nrow = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  
  theme_half_open(12) +
  xlab("Groups") + ylab("% of clones") +
  theme(
    # axis elements  
    axis.text = element_text(colour = NULL, size = 15, margin = margin(t = 0, r = 0, b = 10, l = 0)),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(colour = NULL, size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(colour = NULL, size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    # strip elements
    strip.text = element_text(face = NULL, size = rel(1.2)),
    strip.background = element_rect(fill = "white", size = 1),
    # legend
    legend.position = "none"
  )

facet.plot
