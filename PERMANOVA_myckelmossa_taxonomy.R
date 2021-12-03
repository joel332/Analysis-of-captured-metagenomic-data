###################################################################################
# Script for Session 6 of the course: Applied multivariate Statistics with R  	  #
# 		  created by Eduard Szöcs,	 modified by JDW, 	2020				  #
# 						        PERMANOVA					  		  	  		  #
###################################################################################

###set working directory
setwd("~/R/Skogaryd/Taxanomic/PERMANOVA")


### ------------ Load data and package -----------------------------------------
### Load data and package
require(vegan)
library(ggplot2)
library(dplyr)

myckel_sp <-read.csv2("~/R/Skogaryd/Taxanomic/PERMANOVA/myckelmossa_spp_2017_2018.csv", sep = ";",
                       header = TRUE)
head(myckel_sp)
str(myckel_sp)
class(myckel_sp)

myckel_env <- read.table("~/R/Skogaryd/Taxanomic/PERMANOVA/myckelmossa_env.csv", sep = ";", header = TRUE)
head(myckel_env)
str(myckel_env)

myckel_proportions <-read.csv2("~/R/Skogaryd/Taxanomic/PERMANOVA/proportions.csv", sep = ";",
                      header = TRUE)
head(myckel_proportions)
str(myckel_proportions)
class(myckel_proportions)

### ------------------ Proportions of Methanogens to Methanotrophs ---------------------------###

result_proportions <- myckel_proportions %>%
  group_by(Functional_group) %>%
  get_summary_stats(Abundance, type = "full")
result_proportions

result_proportions <- myckel_proportions %>%
  group_by(Year) %>%
  get_summary_stats(Abundance, type = "full")
result_proportions




norm_proportions <- scale(myckel_proportions$Abundance)

#Save the details and add new column to data matrix then re-import
results <- capture.output(print(norm_proportions))
writeLines(results, con = file("normalized_proportions.txt"))

#add the scaled values to the original data matrix and re-import with new colum
myckel_proportions <-read.csv2("~/R/Skogaryd/Taxanomic/PERMANOVA/proportions.csv", sep = ";",
                               header = TRUE)
head(myckel_proportions)

#set columns as factors
myckel_proportions$Year <- as.factor(myckel_proportions$Year)
myckel_proportions$Functional_group <-  as.factor(myckel_proportions$Functional_group)

### Boxplot Results
prop<-ggplot(myckel_proportions, aes(x=Functional_group, y=Normalized, fill=Year)) +
  ylab(expression(paste("Normalized abundances")))+
  scale_color_manual(values = c("steelblue4", "orangered4")) +
  xlab("Year") +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8), width = 0.5)+
  theme_bw()
prop

#save plot
ggsave("AD_year.tiff", units="in", width=10, height=8, dpi=1000, compression = 'lzw')

### ------------ Distance matrix ----------------------------------------------
### Compute distance matrix using Bray-Curtis on double-root transformed abundances

#use this to ensure transformed data lies between 0 - 10. This makes sure highly abundant species dont over influence the results
range(myckel_sp)
range(myckel_sp ^ 0.5)
range(myckel_sp ^ 0.25)
range(myckel_sp ^ 0.22)
range(myckel_sp ^0.2193)

#I have chosen to use 0.25 convertion

# We use a double root transformation to reduce the range of the data

dist_myckel <- vegdist(myckel_sp ^ 0.21, method = "bray")
dist_myckel
# otherwise highly abundant genes would dominate the distance measures

### ------------ NMDS ----------------------------------------------------------
### Run NMDS
nmds <- metaMDS(dist_myckel, distance = "bray", k=2, trymax = 1000)


#Export "metaMDS" values to ggplot
scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs), Management = myckel_env$Year)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Management, data = scrs, FUN = mean)
segs <- merge(scrs, setNames(cent, c('Management','oNMDS1','oNMDS2')),
              by = 'Management', sort = FALSE)

# Plot NMDS using ggplot
p <- ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = factor(Management), alpha = 0.7)) +
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed() +
  scale_color_manual(values = c("steelblue4", "orangered4")) +
  theme_bw()
p


ggsave("nmds_year.tiff", units="in", width=10, height=8, dpi=1000, compression = 'lzw')

#############################################################
# Treatment "year" appears to group seperatly, however there is some overlap between the years
# 2017 appears to have the largest speard when compare to 2018 indicating 
# more variation in gene abundance in 2017.


### ------------ PERMANOVA -----------------------------------------------------

#PERMANOVA Year
pmv_year <- adonis(
  myckel_sp ^ 0.21 ~ Year, data = myckel_env,
  permutations = 999,
  method = "bray")
pmv_year


#PERMANOVA Ecotype
pmv_ecotype <- adonis(
  myckel_sp ^ 0.21 ~ Ecotype, data = myckel_env,
  permutations = 999,
  method = "bray")
pmv_ecotype

#pairwise permtutation Post hoc test with FDR adjstment
Wilks_pairwise_ecotype <- pairwise.perm.manova(dist(myckel_sp,"euclidean"),myckel_env$Ecotype,nperm=999, test = "Wilks")
Wilks_pairwise_ecotype

#Save the details
results <- capture.output(print(pmv_year), print(pmv_ecotype), print(Wilks_pairwise_ecotype))
writeLines(results, con = file("output_PERMANOVA_year_ecotype.txt"))



###---------- Alpha diversity of all samples between years ---------

AD<- diversity(myckel_sp, index = "shannon")
AD

# Shannon's H'
H <- diversity(myckel_sp)

# Observed Richness
richness <- specnumber(myckel_sp)  

# Pielou's Evenness
evenness <- H/log(richness)

# Create alpha diversity dataframe including environmental data
alpha <- cbind(shannon = H, richness = richness, pielou = evenness, myckel_env)
head(alpha)


alpha_summary <- capture.output(print(alpha))
writeLines(alpha_summary, con = file("alpha_diversity_2017_2018.txt"))


#plot
AD <- ggplot(alpha, aes(x = Year, y = shannon)) + 
  geom_boxplot(aes(color=as.factor(Year)), size = 0.7) +
  scale_color_manual(values = c("steelblue4", "orangered4")) +
  ylab("Shannon's Alpha Diversity") + 
  xlab("") +
  theme_bw()
AD

#save plot
ggsave("AD_year.tiff", units="in", width=10, height=8, dpi=1000, compression = 'lzw')

#Statistics

#Summary statistics  
result_AD <- alpha %>%
  group_by(Year) %>%
  get_summary_stats(shannon, type = "full")
result_AD

Wilks_pairwise_ecotype_AD <- pairwise.perm.manova(dist(alpha$shannon,"euclidean"),alpha$Year,nperm=10000, test = "Wilks")
Wilks_pairwise_ecotype_AD


#Save the details
results <- capture.output(print(result_AD), print(Wilks_pairwise_ecotype_AD))
writeLines(results, con = file("AD_statistics.txt"))



###-------------Beta Diversity-----------------##########################


### ------------ Distance based dispersion test -------------------------------
bd <- betadisper(dist_myckel, myckel_env$Ecotype, bias.adjust = TRUE)

#Save output to create a datamatrix usable by ggplot
distance<- bd$distances
group <- bd$group

#Save the details
results <- capture.output(print(distance), print(group))
writeLines(results, con = file("output_bd_plotting.txt"))

# Once saved, open the file and re-order within excel

#once re-ordered import new data
#import data set beta
my_data_beta <- read.delim("~/R/Skogaryd/Taxanomic/PERMANOVA/bd_ggplot.txt")
head(my_data_beta)

#Summary statistics  
result_bd <- my_data_beta %>%
  group_by(Ecotype) %>%
  get_summary_stats(Beta_Diversity, type = "full")
result_bd

#Save the details
results <- capture.output(print(result_bd))
writeLines(results, con = file("output_bd_summary.txt"))

#plot
bd_plot <- ggplot(my_data_beta, aes(x = Ecotype, y = Beta_Diversity)) + 
  geom_boxplot(aes(color=as.factor(Year)), size = .7) +
  scale_color_manual(values = c("steelblue4", "orangered4")) +
  ylab(expression("Beta Diversity (Distance to centroid)")) +
  theme_bw()
bd_plot

#save plot
ggsave("bd_ecotype_year.tiff", units="in", width=10, height=8, dpi=1000, compression = 'lzw')

Wilks_pairwise_ecotype_bd <- pairwise.perm.manova(dist(my_data_beta,"euclidean"),my_data$Ecotype,nperm=10000, test = "Wilks")
Wilks_pairwise_ecotype_bd

#Save the details
results <- capture.output(print(Wilks_pairwise_ecotype_bd))
writeLines(results, con = file("output_wilks_pairwise_ecotype_bd.txt"))


### ------------ SIMPER --------------------------------------------------------
sim <- simper(myckel_sp, group = myckel_env$Ecotype, permutations = 999)
summary <- summary(sim)
# contr :   contribution to dissimilarity between upstream and downstream
# sd    :   standard deviation of contribution (is the species response consitent?)
# ratio :   ratio between contr and sd (high ratio = high, consisten contribution)
# av.   :   average abundance per groups
# cumsum:   cumulative contribution (rule of thumb : species till 70% are investigated)


#Save the details
gene_summary <- capture.output(print(summary))
writeLines(gene_summary, con = file("taxa_summary_ecotype.txt"))
