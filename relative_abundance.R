########################################################################################################################
##              https://joey711.github.io/phyloseq/plot_richness-examples.html                                        ##
########################################################################################################################

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("readxl")
library("rlang")
library("dplyr")
library("tibble")
library("reshape2")

#sample_data
#data("GlobalPatterns")

#Import data in txt format. 
otu_mat <- read_excel("~/R/Lab_experiment_2017/Relative Abundance/methanogen_otu.xlsx")
row.names(otu_mat) <- otu_mat$otu
otu_mat <- otu_mat %>% select (-otu)
head(otu_mat)

tax_mat <- read_excel("~/R/Lab_experiment_2017/Relative Abundance/methano_tax.xlsx")
tax_mat
row.names(tax_mat) <- tax_mat$otu
tax_mat <- tax_mat %>% select (-otu) 
View(tax_mat)

samples_df<- read_excel("~/R/Lab_experiment_2017/Diversity/Alpha_diversity/samples_df.xlsx")
samples_df
row.names(samples_df) <- samples_df$Sample
samples_df


#OTU's are the absolute abundances of taxa. Columns are counts, Row 1 is OTU1, OTU2, OTU3 .....ect
otu_mat<- as.matrix(otu_mat)
class(otu_mat)


#Taxmat is the imported taxanomic heiarchy save in .txt format
tax_mat<- as.matrix(tax_mat)
class(tax_mat)

#Taxmat is the imported taxanomic heiarchy save in .txt format
samples_df <- as.data.frame(samples_df)
head(sam_data)
class(sam_data)

#combine otumat and taxmat into a phyloseq object (Must be as data matrix)
library("phyloseq")
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAM = sample_data(samples_df)
head(OTU)
head(TAX)
head(SAM)

#plots abundance vs sample
physeq = phyloseq(OTU, TAX, SAM)
physeq


########################################################################
###            Data pruning                                          ###
########################################################################

# only look into high-abundance/high-prevelance OTUs over 10
physeq.f <- filter_taxa(physeq,function(x) sum(x >= 10) > (0.01*length(x)), 
                        prune = TRUE)
#oringinal
physeq
#post pruning
physeq.f

# Covert to relative abundance
physeq.f.ra <- transform_sample_counts(physeq.f, function(x) x*100/sum(x))
physeq.f.ra

############################################################################
###                 Plotting                                             ###
############################################################################
a <- plot_taxa_bar(physeq.f.ra)+ 
  geom_bar(aes(color=Species, fill=genus), stat="identity", position="stack", colour="black")
  theme_bw()
a  

ra <- a + theme(axis.text.x = element_text(angle = 90),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))
ra

methanogen_rlab_pathway <- ra + labs(title = NULL) + xlab("Samples") + ylab("Relative Abundance (%)")+ 
        guides(fill=guide_legend(title="Methanogen (Genus Level)"))
methanogen_rlab_pathway

############################################ Test for differences in functional groups #####################################

#import raw data (processed in excel seperatly, took the means of each functional group according to top, middle and roots)
library(readxl)
methanogen_fucntional_group_absolute <- read_excel("~/R/Lab_experiment_2017/Relative Abundance/test_function.xlsx")
View(methanogen_fucntional_group_absolute)

#normalise dataset using dylyr
library(dplyr)
my_data <- methanogen_fucntional_group_absolute %>% mutate_at(c("Abundance"), ~(scale(.) %>% as.vector))
my_data

# load packages
library("dplyr")
library("ggpubr")
library("DESeq2")

#Density plot and Q-Q plot can be used to check normality visually.
#Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.

ggdensity(my_data$Abundance, 
          main = "Density plot of functional group",
          xlab = "Values")


# Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution. 
#A 45-degree reference line is also plotted.'

ggqqplot(my_data$Abundance)

# Normality test

#Visual inspection, described in the previous section, is usually unreliable. It's possible to use a significance test comparing the sample distribution to a normal one in order to ascertain whether data show or not a serious deviation from normality.
#There are several methods for normality test such as Kolmogorov-Smirnov (K-S) normality test and Shapiro-Wilk's test.

#The null hypothesis of these tests is that "sample distribution is normal". If the test is significant, the distribution is non-normal.

#Shapiro-Wilk's method is widely recommended for normality test and it provides better power than K-S. It is based on the correlation between the data and the corresponding normal scores.
#Note that, normality test is sensitive to sample size. Small samples most often pass normality tests. Therefore, it's important to combine visual inspection and significance test in order to take the right decision.

shapiro.test(my_data$Abundance)

#this example is not normally distributed therefore we use a non parametric test

#Compute Kurskal wallis test

#Check strings to make sure the treatment column is a factor
str(my_data)

#convert to factore if necessary (from chr to factor)
#one approach it to index with the $ sign and the as.factor function
my_data$`Functional group` <- as.factor(my_data$`Functional group`)

#check strings again to ensure it worked
str(my_data)

# Show the group levels
levels(my_data$`Functional group`)

#If the levels are not automatically in the correct order, re-order them as follow:
  
  #methanogen_fucntional_group_absolute$`Functional group` <- ordered(methanogen_fucntional_group_absolute$`Functional group`,
                           levels = c("Acetoclastic" ,"Hydrogenotrophic","Multiple pathway", "Methylotrophic")

  
#Compute summary statistics  

  library(dplyr)
  group_by(my_data, `Functional group`) %>%
    summarise(
      count = n(),
      mean = mean(Abundance, na.rm = TRUE),
      sd = sd(Abundance, na.rm = TRUE),
      median = median(Abundance, na.rm = TRUE),
      IQR = IQR(Abundance, na.rm = TRUE)
    )

  # Box plots
  # ++++++++++++++++++++
  # Plot abundance by functional group and color by group
  library("ggpubr")
  ggboxplot(my_data, x = "Functional group", y = "Abundance", 
            color = "Functional group", palette = "grey",
            ylab = "Abundance", xlab = "Functional Group")
  
  # Mean plots
  # ++++++++++++++++++++
  # Plot weight by group
  # Add error bars: mean_se
  # (other values include: mean_sd, mean_ci, median_iqr, ....)

  ggline(methanogen_fucntional_group_absolute, x = "Functional group", y = "Abundance", 
         add = c("mean_se", "jitter"),
         ylab = "Abundance", xlab = "Functional group")  
  
  #The test can be performed using the function kruskal.test() as follow:
    
    kruskal.test(Abundance ~ methanogen_fucntional_group_absolute$`Functional group`, data = methanogen_fucntional_group_absolute)
 
    pairwise.wilcox.test(methanogen_fucntional_group_absolute$Abundance, methanogen_fucntional_group_absolute$`Functional group`,
                         p.adjust.method = "BH")
    
    
