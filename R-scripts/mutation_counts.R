#### creating a mutational load script for VH-gen analysis ####

####library####
library(readxl)
library(ggplot2)
library(qqplotr)
library(dplyr)
library(DescTools)
library(tidyverse)
library(ggsignif)
library(rstatix)
library(ggpubr)


# import data 
library(readr)
df <- read_delim("C:/Users/rickw/Downloads/RUN43_vjcdr3-clones-mut-Wilbrink-IGH_HUMAN.csv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)

# add junction length
df$junction_l <- nchar(df$cdr3nuc.min_mode)

# load dplyr package
library(dplyr)

# filter rows that contain the string 'AS' in the Sample column
df_AS <- df %>% filter(grepl('AS', Sample))

# assign values to groups 
df_AS$group = ifelse(grepl("totalB",df_AS$Sample),1,0)
df_AS$group = ifelse(grepl("CD27m",df_AS$Sample),2,df_AS$group)
df_AS$group = ifelse(grepl("CD27p",df_AS$Sample),3,df_AS$group)
df_AS$group = ifelse(grepl("plas",df_AS$Sample),4,df_AS$group)

# attach levels and labels
df_AS$subset <- factor(df_AS$group, 
                   levels = c(1,2,3,4),
                   labels = c("TotalB", "CD27-CD21low", "CD27+CD21low", "plasma"))


####statistical analysis####

# descriptive statistics
ds <- df_AS %>% select(subset, mut.frac_x.mean) %>% group_by(subset) %>% 
  summarise(n = n(), 
            mean = mean(mut.frac_x.mean, na.rm = TRUE), 
            sd = sd(mut.frac_x.mean, na.rm = TRUE),
            stderr = sd/sqrt(n),
            LCL = mean - qt(1 - (0.05 / 2), n - 1) * stderr,
            UCL = mean + qt(1 - (0.05 / 2), n - 1) * stderr,
            median = median(mut.frac_x.mean, na.rm = TRUE),
            min = min(mut.frac_x.mean, na.rm = TRUE), 
            max = max(mut.frac_x.mean, na.rm = TRUE),
            IQR = IQR(mut.frac_x.mean, na.rm = TRUE),
            LCLmed = MedianCI(mut.frac_x.mean, na.rm=TRUE)[2],
            UCLmed = MedianCI(mut.frac_x.mean, na.rm=TRUE)[3])
print(ds)

# Perform QQ plots by group
qqplot <- ggplot(data = df_AS, mapping = aes(sample = mut.frac_x.mean, color = subset, fill = subset)) +
  stat_qq_band(alpha=0.5, conf=0.95, qtype=1, bandType = "boot") +
  stat_qq_line(identity=TRUE) +
  stat_qq_point(col="black") +
  facet_wrap(~ subset, scales = "free") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + theme_bw()

print(qqplot)

# change fill and outline color manually 
ggplot(df_AS, aes(x = mut.frac_x.mean)) +
  geom_histogram(aes(color = subset, fill = "white"), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#de2d26", "#2c7fb8")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#de2d26", "#2c7fb8"))

# perform the kurskal-Wallis test
kt <- kruskal.test(mut.frac_x.mean ~ subset, data = df_AS)

print(kt)

# wilcox.test
pairwise.wilcox.test(df_AS$mut.frac_x.mean, df_AS$subset, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)

# my comparisons for statistical testing 
my_comparisons <- list(c("total B", "CD27-CD21low"), c("total B", "CD27+CD21low"), c("total B", "plasma"),
                       c("CD27-CD21low", "CD27+CD21low"), c("CD27-CD21low", "plasma"), c("CD27+CD21low", "plasma")
                       )


# ggplot subsets
p <- ggplot(df_AS, aes(subset, mut.frac_x.mean, fill = subset)) +
  
geom_boxplot() +
  
  
  ggtitle("Mutational status B cell populations of AS patients") +
  
  xlab(NULL) + ylab("Mutational fraction") +
  
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
   

  theme( # theme of the figure 
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 18),
    axis.title.y = element_text(colour = "black", size = 20),
    
    legend.position = "none", # legend theme
    
    panel.border = element_blank(), # panel theme
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", size = 1.1)
  )

print(p)

