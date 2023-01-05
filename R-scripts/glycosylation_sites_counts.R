#### Glycosylation sites per cell subset in a group ####

#### loading packages ####
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cowplot)

#### import and load data #### 
df <- read.delim("D:/PhD/Projecten/VH-gen project/Dataset/9-5-2022 Combi 3 runs/vjcdr3-clones-mut-Wilbrink-IGH_HUMAN.csv")


View(df)

## remove clones with a freq of 1, thus singletons that are not clones
#df[df$freq != 1, ] 

#df <- df %>% filter(df$mut.count_x.mode >= 20)

## Assign mutation counts
df$nr_sites.count <- ifelse(df$nr_sites.mean == 0, 1, 0)
df$nr_sites.count <- ifelse(df$nr_sites.mean > 0 & df$nr_sites.mean < 1, 2, df$nr_sites.count)
df$nr_sites.count <- ifelse(df$nr_sites.mean >= 1 & df$nr_sites.mean < 2, 3, df$nr_sites.count)
df$nr_sites.count <- ifelse(df$nr_sites.mean >= 2, 4, df$nr_sites.count)


## Factor
df$nr_sites.count <- factor(df$nr_sites.count, c(1, 2, 3, 4),labels = c("0","<1", "1-2", "2<"))

## Assign groups/cell subsets
#first check unique <- unique(df$Sample)

#assign values to cell subsets, use function grepl(test, yes, no) 
df$group = ifelse(grepl("totalB",df$Sample), "total B cells" , 0)
df$group = ifelse(grepl("CD27m",df$Sample), "CD27negativeCD38lowCD21low" , df$group)
df$group = ifelse(grepl("CD27p",df$Sample), "CD27positiveCD38lowCD21low" , df$group)
df$group = ifelse(grepl("plas",df$Sample), "plasmablasts" , df$group)

#### filter per specific patient group ####
df2 <- df %>% 
  filter(grepl('AS', Sample))

#### transforming data, give every patient/sample an unique value #### 

## total B cells ## 
df.totalB <- df2 %>% 
  filter(grepl('total B cells', group))
df.totalB <- transform(df.totalB, sample_id = match(Sample, unique(Sample)))

## CD27+CD38lowCD21low ##
df.CD27m <- df2 %>% 
  filter(grepl('CD27negative', group))
df.CD27m <- transform(df.CD27m, sample_id = match(Sample, unique(Sample)))

## CD27+CD38lowCD21low ## 
df.CD27p <- df2 %>% 
  filter(grepl('CD27positive', group))
df.CD27p <- transform(df.CD27p, sample_id = match(Sample, unique(Sample)))

## plasmablasts ##
df.p <- df2 %>% 
  filter(grepl('plasmablasts', group))
df.p <- transform(df.p, sample_id = match(Sample, unique(Sample)))

#### create table with percentages of mutations within unique samples per cell subset ####
nr_sites.count.sample_id <- table(df.totalB$sample_id, df.totalB$nr_sites.count)
t.totB <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

nr_sites.count.sample_id <- table(df.CD27m$sample_id, df.CD27m$nr_sites.count)
t.CD27m <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

nr_sites.count.sample_id <- table(df.CD27p$sample_id, df.CD27p$nr_sites.count)
t.CD27p <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

nr_sites.count.sample_id <- table(df.p$sample_id, df.p$nr_sites.count)
t.plas <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

##### create data frames ####
t.totB.1 <- apply(as.matrix.noquote(t.totB),2,as.numeric)
df.totB <- as.data.frame(t.totB.1)

t.CD27m.1 <- apply(as.matrix.noquote(t.CD27m),2,as.numeric)
df.CD27m <- as.data.frame(t.CD27m.1)

t.CD27p.1 <- apply(as.matrix.noquote(t.CD27p),2,as.numeric)
df.CD27p <- as.data.frame(t.CD27p.1)

t.plas.1 <- apply(as.matrix.noquote(t.plas),2,as.numeric)
df.plas <- as.data.frame(t.plas.1)

## create columns that can be used for plotting ##
df.totB <- data.frame("totB" = c(df.totB[,"0"], df.totB[,"<1"], df.totB[,"1-2"], df.totB[,"2<"]))
df.totB$group <- c(rep("0", 10), rep("<1", 10), rep("1-2", 10), rep("2<", 10))
df.totB$group <- fct_relevel(df.totB$group,"0", "<1", "1-2", "2<")

df.CD27m <- data.frame("CD27m" = c(df.CD27m[,"0"], df.CD27m[,"<1"], df.CD27m[,"1-2"], df.CD27m[,"2<"]))
df.CD27m$group <- c(rep("0", 10), rep("<1", 10), rep("1-2", 10), rep("2<", 10))
df.CD27m$group <- fct_relevel(df.CD27m$group,"0", "<1", "1-2", "2<")

df.CD27p <- data.frame("CD27p" = c(df.CD27p[,"0"], df.CD27p[,"<1"], df.CD27p[,"1-2"], df.CD27p[,"2<"]))
df.CD27p$group <- c(rep("0", 10), rep("<1", 10), rep("1-2", 10), rep("2<", 10))
df.CD27p$group <- fct_relevel(df.CD27p$group,"0", "<1", "1-2", "2<")

df.plas <- data.frame("plas" = c(df.plas[,"0"], df.plas[,"<1"], df.plas[,"1-2"], df.plas[,"2<"]))
df.plas$group <- c(rep("0", 10), rep("<1", 10), rep("1-2", 10), rep("2<", 10))
df.plas$group <- fct_relevel(df.plas$group,"0", "<1", "1-2", "2<")

## create dataframe for facet_wrap plot
colnames(df.totB) <- c("percentage", "group")
colnames(df.CD27m) <- c("percentage", "group")
colnames(df.CD27p) <- c("percentage", "group")
colnames(df.plas) <- c("percentage", "group")
df.merged <- rbind(df.totB, df.CD27m, df.CD27p, df.plas)
df.merged$condition <- c(rep("Total B cells", 40), rep("CD27-CD38lowCD21low", 40), rep("CD27+CD38lowCD21low", 40), rep("Plasmablasts", 40))
df.merged$condition <- fct_relevel(df.merged$condition,"Total B cells", "CD27-CD38lowCD21low", "CD27+CD38lowCD21low", "Plasmablasts")

#### create ggplot with facet_wrap ####
facet.plot <- ggplot(df.merged, aes(group, percentage, color = condition)) +
  facet_wrap(~condition, nrow = 1) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  
  theme_half_open(12) +
  xlab("Number of glycosylation sites") + ylab("% of clones") +
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

#### save as tiff file ####
tiff("plot HD.tiff", units="in", width=10.4, height=6.98, res=600)

facet.plot

dev.off()






















# create data table for exporting ----
# matrix 

#### import and load data #### 
df <- read.delim("D:/PhD/Projecten/VH-gen project/Dataset/9-5-2022 Combi 3 runs/vjcdr3-clones-mut-Wilbrink-IGH_HUMAN.csv")

## remove clones with a freq of 1, thus singletons that are not clones
#df[df$freq != 1, ] 

#df <- df %>% filter(df$mut.count_x.mode >= 20)

## Assign mutation counts
df$nr_sites.count <- ifelse(df$nr_sites.mean == 0, 1, 0)
df$nr_sites.count <- ifelse(df$nr_sites.mean > 0 & df$nr_sites.mean < 1, 2, df$nr_sites.count)
df$nr_sites.count <- ifelse(df$nr_sites.mean >= 1 & df$nr_sites.mean < 2, 3, df$nr_sites.count)
df$nr_sites.count <- ifelse(df$nr_sites.mean >= 2, 4, df$nr_sites.count)


## Factor
df$nr_sites.count <- factor(df$nr_sites.count, c(1, 2, 3, 4),labels = c("0","<1", "1-2", "2<"))

## Assign groups/cell subsets
#first check unique <- unique(df$Sample)

#assign values to cell subsets, use function grepl(test, yes, no) 
df$group = ifelse(grepl("totalB",df$Sample), "total B cells" , 0)
df$group = ifelse(grepl("CD27m",df$Sample), "CD27negativeCD38lowCD21low" , df$group)
df$group = ifelse(grepl("CD27p",df$Sample), "CD27positiveCD38lowCD21low" , df$group)
df$group = ifelse(grepl("plas",df$Sample), "plasmablasts" , df$group)

#### filter per specific patient group ####
df2 <- df %>% 
  filter(grepl('pSS', Sample))

df2 <- df2 %>% filter(!grepl('pSS208', Sample))

#### transforming data, give every patient/sample an unique value #### 

## total B cells ## 
df.totalB <- df2 %>% 
  filter(grepl('total B cells', group))
df.totalB <- transform(df.totalB, sample_id = match(Sample, unique(Sample)))

## CD27+CD38lowCD21low ##
df.CD27m <- df2 %>% 
  filter(grepl('CD27negative', group))
df.CD27m <- transform(df.CD27m, sample_id = match(Sample, unique(Sample)))

## CD27+CD38lowCD21low ## 
df.CD27p <- df2 %>% 
  filter(grepl('CD27positive', group))
df.CD27p <- transform(df.CD27p, sample_id = match(Sample, unique(Sample)))

## plasmablasts ##
df.p <- df2 %>% 
  filter(grepl('plasmablasts', group))
df.p <- transform(df.p, sample_id = match(Sample, unique(Sample)))

#### create table with percentages of mutations within unique samples per cell subset ####
nr_sites.count.sample_id <- table(df.totalB$sample_id, df.totalB$nr_sites.count)
t.totB <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

nr_sites.count.sample_id <- table(df.CD27m$sample_id, df.CD27m$nr_sites.count)
t.CD27m <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

nr_sites.count.sample_id <- table(df.CD27p$sample_id, df.CD27p$nr_sites.count)
t.CD27p <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

nr_sites.count.sample_id <- table(df.p$sample_id, df.p$nr_sites.count)
t.plas <- round(prop.table(nr_sites.count.sample_id,1)*100,digits=2)

# create matrix per subset
totB <- as.matrix(t.totB)
colnames(totB)<-paste(colnames(totB),"totB",sep="_")
totB

CD27m <- as.matrix(t.CD27m)
colnames(CD27m)<-paste(colnames(CD27m),"CD27m",sep="_")
CD27m

CD27p <- as.matrix(t.CD27p)
colnames(CD27p)<-paste(colnames(CD27p),"CD27p",sep="_")
CD27p

plas <- as.matrix(t.plas)
colnames(plas)<-paste(colnames(plas),"plas",sep="_")
plas

m <- cbind(totB, CD27m, CD27p, plas)

df.all <- as.data.frame(m , row.names = TRUE)
rownames(df.all) <-rnames


# AS
df.AS <- df.all
rnames <- c('AS13B', 'AS26', 'AS27', 'AS28', 'AS32' ,'AS35', 'AS36', 'AS38', 'AS39', 'AS6B')

# HD
df.HD <- df.all
rnames <- c('HD1', 'HD10', 'HD11' , 'HD2', 'HD3', 'HD4', 'HD6', 'HD7', 'HD8', 'HD9')

# pSS
df.pSS <- df.all
rnames <- c('pSS115', 'pSS131', 'pSS136', 'pSS145', 'pSS149', 'pSS160', 'pSS193', 'pSS80', 'pSS87')

# bind together
df.final <- rbind(df.AS, df.HD, df.pSS)
df.final$sample <- rownames(df.final)

View(df.final)
str(df.final)

library("writexl")
write_xlsx(df.final,"D:/OneDrive/PhD/Projecten/VH-gen project/Excel/Glyc_sites_allsamples.xlsx")






























## create columns that can be used for plotting ##
df.totB <- data.frame("totB" = c(df.totB[,"0-2"], df.totB[,"2-10"], df.totB[,"10<"]))
df.totB$group <- c(rep("0-2", 10), rep("2-10", 10), rep("10<", 10))
df.totB$group <- fct_relevel(df.totB$group,"0-2", "2-10", "10<")

df.CD27m <- data.frame("CD27m" = c(df.CD27m[,"0-2"], df.CD27m[,"2-10"], df.CD27m[,"10<"]))
df.CD27m$group <- c(rep("0-2", 10), rep("2-10", 10), rep("10<", 10))
df.CD27m$group <- fct_relevel(df.CD27m$group,"0-2", "2-10", "10<")

df.CD27p <- data.frame("CD27p" = c(df.CD27p[,"0-2"], df.CD27p[,"2-10"], df.CD27p[,"10<"]))
df.CD27p$group <- c(rep("0-2", 10), rep("2-10", 10), rep("10<", 10))
df.CD27p$group <- fct_relevel(df.CD27p$group,"0-2", "2-10", "10<")

df.plas <- data.frame("plas" = c(df.plas[,"0-2"], df.plas[,"2-10"], df.plas[,"10<"]))
df.plas$group <- c(rep("0-2", 10), rep("2-10", 10), rep("10<", 10))
df.plas$group <- fct_relevel(df.plas$group,"0-2", "2-10", "10<")


#### create single ggplots  ####
plot.totB <- ggplot(df.totB, aes(group, totB, group = group)) + 
  geom_jitter(color = "blue", size = 2) +
  ggtitle("Total B cells") +
  xlab(NULL) + ylab("% of clones") +
  scale_y_continuous(breaks = seq(0,100, by = 20))
  

plot.CD27m <- ggplot(df.CD27m, aes(group, CD27m, group = group)) + 
  geom_jitter(color = "gray", size = 2) +
  ggtitle("CD27m") +
  xlab(NULL) + ylab("% of clones") +
  scale_y_continuous(breaks = seq(0,100, by = 20))

plot.CD27p <- ggplot(df.CD27p, aes(group, CD27p, group = group)) + 
  geom_jitter(color = "green", size = 2) +
  ggtitle("CD27p") +
  xlab(NULL) + ylab("% of clones") +
  scale_y_continuous(breaks = seq(0,100, by = 20))

plot.plas <- ggplot(df.plas, aes(group, plas, group = group)) + 
  geom_jitter(color = "purple", size = 2) +
  ggtitle("plasmablasts") +
  xlab(NULL) + ylab("% of clones") +
  scale_y_continuous(breaks = seq(0,100, by = 20))

merged.plot <- cowplot::plot_grid(plot.totB,
                   plot.CD27m + theme(axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      axis.title.y = element_blank() ) ,
                   plot.CD27p + theme(axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    axis.title.y = element_blank() ),
                   plot.plas + theme(axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank(),
                                  axis.title.y = element_blank() ), 
                   nrow = 1
                   )




































merged.plot.2 <- align_plots(plot.totB, 
                           plot.CD27m + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ), 
            plot.CD27p + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
            plot.plas + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank() ),
          nrow = 1)

plot_grid(plot.totB, plot.CD27m, plot.CD27p, plot.plas, nrow = 1, align = "v")


merged.plot

merged.plot.2










## create nice looking fig

plot.totB2 <- ggplot(df.plas, aes(group, plas, group = group)) +
  
  geom_boxplot(
    lwd = 1.2,
    color = "lightblue",
    outlier.shape = F
    
  ) +
  
  geom_jitter(
    width = 0.15, # points spread out over 15% of available width
    height = 0, # do not move position on the y-axis
    alpha = 0.5,
    size = 4,
    color = "blue"
  ) +
  
  ggtitle("Total B cells") +
  xlab(NULL) + ylab("% of clones") +
  
  theme( # theme of the figure 
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(colour = "black", size = 15),
    axis.title.x = element_text(colour = "black", size = 40),
    axis.title.y = element_text(colour = "black", size = 20),
    legend.position = "none", # legend theme
    panel.border = element_blank(), # panel theme
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", size = 1.1)
  )













## create dataframe for facet_wrap plot
colnames(df.totB) <- c("percentage", "group")
colnames(df.CD27m) <- c("percentage", "group")
colnames(df.CD27p) <- c("percentage", "group")
colnames(df.plas) <- c("percentage", "group")
df.merged <- rbind(df.totB, df.CD27m, df.CD27p, df.plas)
df.merged$condition <- c(rep("Total B cells", 30), rep("CD27-CD38lowCD21low", 30), rep("CD27+CD38lowCD21low", 30), rep("Plasmablasts", 30))
df.merged$condition <- fct_relevel(df.merged$condition,"Total B cells", "CD27-CD38lowCD21low", "CD27+CD38lowCD21low", "Plasmablasts")

## ggplot with facet_wrap
facet.plot <- ggplot(df.merged, aes(group, percentage, color = condition)) +
  facet_wrap(~condition, nrow = 1) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_half_open(12) +
  xlab("Mutations counts") + ylab("Percentage of clones")
  
facet.plot

#### CODE FOR MERGER ####

## create columns that can be used for plotting ##
df.totB <- data.frame("total B cells" = c(df.totB[,"0-2"], df.totB[,"2-10"], df.totB[,"10<"]))
df.totB$id <- c(1:30)

df.CD27m <- data.frame("CD27m" = c(df.CD27m[,"0-2"], df.CD27m[,"2-10"], df.CD27m[,"10<"]))
df.CD27m$id <- c(1:30)

df.CD27p <- data.frame("CD27p" = c(df.CD27p[,"0-2"], df.CD27p[,"2-10"], df.CD27p[,"10<"]))
df.CD27p$id <- c(1:30)

df.plas <- data.frame("plas" = c(df.plas[,"0-2"], df.plas[,"2-10"], df.plas[,"10<"]))
df.plas$id <- c(1:30)

#put all data frames into list
df_list <- list(df.totB, df.CD27m, df.CD27p, df.plas)      

#merge all data frames together
df.merged <- df_list %>% reduce(full_join, by = "id")

#add metadata 
df.merged$group <- c(rep("0-2", 10), rep("2-10", 10), rep("10<", 10))
df.merged$sample <- c(rep(1:10, 3))