##################################################################################
# 		  created by Eduard Szöcs,	 modified by Joel D. White, 	2020				  #
# 						        PERMANOVA					  		  	  		  #
###################################################################################

###set working directory
setwd()


### ------------ Load data and package -----------------------------------------
### Load data and package
require(vegan)
library(ggplot2)

myckel_sp <-read.csv2()
head()
str()
class()

myckel_env <- read.table()
head()
str()



### ------------ Distance matrix ----------------------------------------------
### Compute distance matrix using Bray-Curtis on double-root transformed abundances

#use this to ensure transformed data lies between 0 - 10. This makes sure highly abundant species dont over influence the results
range(myckel_sp)
range(myckel_sp ^ 0.5)
range(myckel_sp ^ 0.25)
range(myckel_sp ^ 0.22)
range(myckel_sp ^0.20)

#I have chosen to use 0.25 convertion

dist_myckel <- vegdist(myckel_sp ^ 0.25, method = "bray")
dist_myckel

# otherwise highly abundant genes would dominate the distance measures

### ------------ NMDS ----------------------------------------------------------
### Run NMDS
nmds <- metaMDS(dist_myckel, distance = "bray", k=2, trymax = 1000)


#Export "metaMDS" values for ggplot format
scrs <- scores(nmds, display = 'sites')
scrs <- cbind(as.data.frame(scrs), Management = myckel_env$Group)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Management, data = scrs, FUN = mean)
segs <- merge(scrs, setNames(cent, c('Management','oNMDS1','oNMDS2')),
              by = 'Management', sort = FALSE)

# Plot NMDS using ggplot
p <- ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = Management, alpha = 0.2)) +
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
  geom_point(data = cent, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed() +
  theme_bw()
p


ggsave("nmds_group.tiff", units="in", width=10, height=8, dpi=300, compression = 'lzw')


### ------------ PERMANOVA -----------------------------------------------------

#PERMANOVA Group
pmv_Group <- adonis(
  myckel_sp ^ 0.25 ~ Group, data = myckel_env,
  permutations = 999,
  method = "bray")
pmv_Group


#Save the details
results <- capture.output(print(pmv_Group), print(pmv_Group))
writeLines(results, con = file("output_PERMANOVA.txt"))

#Post Hoc test - pairwise PERMANOVA
library(RVAideMemoire)
Wilks_pairwise_Group <- pairwise.perm.manova(dist(myckel_sp,"euclidean"),myckel_env$Group,nperm=999, test = "Wilks")
Wilks_pairwise_Group

#Save the details
results <- capture.output(print(Wilks_pairwise_Group))
writeLines(results, con = file("output_pairwise_Group.txt"))


### ------------ Distance based dispersion test -------------------------------
bd <- betadisper(dist_myckel, myckel_env$Group)
bd
# also an eigenvalue based method

boxplot(bd)
# boxplot of Average distance to median shows also that there might be lower dispersion
# for upstream sites

# F-Test
anova(bd)

# permutaion test
permutest(bd)
# We cannot find a statistically significantly different dispersion
# -> assumption of homogeneity is met


### ------------ SIMPER --------------------------------------------------------
sim <- simper(myckel_sp, group = myckel_env$Group, permutations = 999)
summary(sim)
# contr :   contribution to dissimilarity between upstream and downstream
# sd    :   standard deviation of contribution (is the species response consitent?)
# ratio :   ratio between contr and sd (high ratio = high, consisten contribution)
# av.   :   average abundance per groups
# cumsum:   cumulative contribution (rule of thumb : species till 70% are investigated)


#Save the details
gene_summary <- capture.output(print(sim))
writeLines(gene_summary, con = file("gene_summary_group.txt"))

