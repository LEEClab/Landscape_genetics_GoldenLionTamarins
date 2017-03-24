###########################################################################
#
# Exploting the distribution of allele distance displacements 
#   for Golden-Lion Tamarins (Leontopithecus rosalia)
# 
# Bernardo Niebuhr <bernardo_brandaum@yahoo.com.br>
# Andreia Moraes <andreia_magro@hotmail.com>
#
# Citation: Moraes et al. Landscape connectivity influences gene flow of 
#   endangered golden lion tamarins (Leontopithecus rosalia) from Atlantic Forest.
#   Submitted to Molecular Ecology.
#
# Feel free to use, modify, and share
# No rights in this world are reserved
###########################################################################

# Loading packages
library(TeachingDemos)

# OPTIONS
# Plot figures?
plot_figures <- FALSE

# Directories
datadir <- "/home/leecb/Dropbox/Moraes et al_connectivity/Allele_movement_final_code/Data"
resultsdir <- "/home/leecb/Dropbox/Moraes et al_connectivity/Allele_movement_final_code/Results"

# Working directory - to load data
setwd(datadir)

# Loading data
dist <- read.table("Distance_between_shared_alleles.csv", header = T, sep = "\t")
head(dist)

#############################
# Ploting Euclidean and Resistance-based distances

# Setting to results folder
setwd(resultsdir)

# Average Euclidean distance
hist(dist$Euclidean_distance/1000, xlab = "Euclidean Distance (km)", ylab = "Frequency", 
     main = "", breaks = 100)

# Average Resistance-based distances
hist(dist$Resistance_based_distance/1000, xlab = "Resistance-based Distance (km)", ylab = "Frequency", 
     main = "", breaks = 100)

# Cutting by the 95th percentile of Euclidean distance
quantile(dist$Euclidean_distance, 0.95) # around 14 km

# Data correspondent to Euclidean distance only up ate 14 km
d <- subset(dist, Euclidean_distance <= 14000) 

# One type of plot
if(plot_figures) {
  tiff("euclid_resist_distance_v1.tif", width = 20, height = 15, units = "cm", res = 300)
  hist(d$Euclidean_distance/1000, xlab = "Euclidean Distance (km)", ylab = "Number of shared alleles", 
       main = "", breaks = seq(0,14,0.5), xlim = c(0, 14), cex.lab = 1.3)
  graphics::text(7, 140, "Resistance-based distance (km)", cex = 1)
  hefdis <- function()
    hist(d$Resistance_based_distance/1000, breaks = seq(0,25,0.5), main = "",
         xlab = "", ylab = "")
  subplot(hefdis(), c(4,10), c(200, 400))
  dev.off()
}

# Another type of plot
if(plot_figures) {
  tiff("euclid_resist_distance_v2.tif", width = 169, height = 100, units = "mm", res = 300)
  par(mfrow = c(1,2), oma = c(0,3,1,0) + 0.1, mar = c(4,2,1,0) + 0.1)
  hist(d$Euclidean_distance/1000, xlab = "Euclidean distance (km)", ylab = "", font.lab=2,
       main = "", breaks = seq(0,14,0.5), xlim = c(0, 14), cex.lab = 1.2, family="serif")
  legend("topright", legend = "(a)", bty = "n", cex = 1.1)
  hist(d$Resistance_based_distance/1000, breaks = seq(0,25,0.5), main = "",
      xlab = "Resistance-based distance (km)", ylab = "", cex.lab = 1.2, font.lab=2, axes = F, family="serif")
  axis(2)
  axis(1, at = seq(0,25,5), labels = seq(0,25,5))
  legend("topright", legend = "(b)", bty = "n", cex = 1.1)
  mtext("Number of shared alleles", side = 2, outer = T, line = 1, font = 2, cex = 1.2, family="serif")
  dev.off()
}

############################
# Calculating statistics

# Euclidean distance
mean(dist$Euclidean_distance/1000); sd(dist$Euclidean_distance/1000)
median(dist$Euclidean_distance/1000); quantile(dist$Euclidean_distance/1000, probs = c(0.05, 0.91, 0.95))

# Resistance-based distance
mean(dist$Resistance_based_distance/1000, na.rm = T); sd(dist$Resistance_based_distance/1000, na.rm = T)
median(dist$Resistance_based_distance/1000, na.rm = T); quantile(dist$Resistance_based_distance/1000, na.rm = T, probs = c(0.05, 0.95))
max(d$Resistance_based_distance/1000, na.rm = T)

median(dist$Resistance_based_distance/1000, na.rm = T)/median(dist$Euclidean_distance/1000)