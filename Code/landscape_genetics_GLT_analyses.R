###########################################################################
#
# Landscape genetics of Golden-Lion Tamarins (Leontopithecus rosalia)
# 
# Andreia Moraes <andreia_magro@hotmail.com>
# Bernardo Niebuhr <bernardo_brandaum@yahoo.com.br>
# Renata Muylaert <renatamuy@gmail.com>
#
# Citation: Moraes et al. Landscape connectivity influences gene flow of 
#   endangered golden lion tamarins (Leontopithecus rosalia) from Atlantic Forest.
#   Submitted to Molecular Ecology.
#
# The data presents data on the kinship between pairs of individuals and has
# the following columns:
#   - Name.i: name of the individual i
#   - Name.j: name of the individual j
#   - Source_patch: ID of the forest patch where individual i was captured
#   - Target_patch: ID og the forest patch where individual j was captured
#   - Sex.i: sex of the individual i
#   - Sex.j: sex of the individual j
#   - Fij: Pairwise kinship, Fij
#   - Management: categories of comparison between pairs of individuals i and j, based on the management 
#       categories native (1), reintroducted (2), translocated (3), and unknown (4)
#       The combinations are e.g. 12 (native-reintroducted) and 24 (reintroducted-unknown)
#   - Roads: binary variable, represents whether there are (1) or not (0) roads between the capture
#       locations of individuals i and j
#   - Resistance_based_distance: Length of the corridor between the Source_patch and the Target_patch, simulated
#       using LSCorridors package and considering the least cost path and matrix resistances
#   - Euclidean_distance: Euclidean distance between the capturelocations of individuals i and j
#   - Landscape_resistance: Landscae resistance, corresponfing to the total cost of the corridor 
#       simulated between the Source_patch and the Target_patch
#   - Landscape_connectivity: Connectivity index calculated as Euclidean_distance/Landscape_resistance
#
# Feel free to use, modify, and share
# No rights in this world are reserved
###########################################################################

#######
# Loading packages
require(bbmle) # for fitting GLMs and applying AIC

#######
# Setting working directory
datadir <- "/home/leecb/Github/Landscape_genetics_GoldenLionTamarins/Data"
setwd(datadir)

#######
# Loading data
gendata <- read.table("Moraes_etal_data_landscape_genetics_Leontopithecus_rosalia.csv", header = T, sep = "\t", dec = ".")
View(gendata) 
str(gendata)

########
# Re-organizing data

# Changing the type of the variables "Management" and "Road" to factor
gendata$Management <- as.factor(gendata$Management)
gendata$Roads <- as.factor(gendata$Roads)

##########################
# Exploratory Analysis

###########
## Outliers

# Response variable: Pairwise kinship (Fij)

# Maximum and minimum values
max(gendata$Fij)
min(gendata$Fij)

# Plotting
plot(gendata$Fij ~ 1) 
# Removing outliers: Fij < -0.4 (only 1 observation) and Fij > 0.599 (4 observations)
gendata <- gendata[gendata$Fij > -0.4 & gendata$Fij <= 0.599,]

# Maximum and minimum values
max(gendata$Fij)
min(gendata$Fij)

# Exploring the relation between "Management" and other variables
out1 <- par(mfrow = c(3, 2))
plot(gendata$Fij ~ gendata$Management, ylab = "Fij", xlab = "Management")
plot(gendata$Resistance_based_distance ~ gendata$Management, ylab = "Corridor length", xlab = "Management")
plot(gendata$Euclidean_distance ~ gendata$Management, ylab = "Euclidean distance", xlab = "Management")
plot(gendata$Landscape_resistance ~ gendata$Management, ylab = "Landscape resistance", xlab = "Management")
plot(gendata$Landscape_connectivity ~ gendata$Management, ylab = "Connectivity index", xlab = "Management")
par(out1)

########
# Collinearity of the explanatory variables.
# We used Spearman correlation, since it does not assume linearity betwee varibales (Zar, 1996)
# We consider a correlation positive if r > 0.7 (Zuur et al 2009)
# In case r > 0.7, we fitted univariate GLMs, compared them through AIC and selected the most plausible variables
# As a result, we kept the varibles: Management, DE, Roads, and Landscape_resistance

# Plotting and calculating correlations

# Function to plot correlation values in the upper part of a pair plot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x,y,  method ="spearman")$estimate
  p <- cor.test(x,y,  method ="spearman")$p.value
  #text(0.5, 0.5, txt, cex = cex.cor * r * 1.1)
  text(0.5, 0.5, paste("r = ", round(r, 2)), cex = 1.1, font = ifelse(r > 0.7, 2, 1))
  text(0.5, 0.2, paste("p = ", round(p, 3)), cex = 1.1)
}

# Function to plot histograms in the diagonal of a pair plot
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}

# Pair plot
pairs(~ Management + Roads + Landscape_resistance + Euclidean_distance + Resistance_based_distance + Landscape_connectivity,
      data=gendata,lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist)

##########
# Selecting the most plausible between the correlated variables: Euclidean_distance, Resistance_based_distance, and Landscape_connectivity

# Euclidean_distance
M1 <- glm(Fij ~ Euclidean_distance, data = gendata, family = gaussian)
summary(M1)

# Resistance_based_distance
M2 <- glm(Fij ~ Resistance_based_distance, data = gendata, family = gaussian)
summary(M2)

# Landscape_connectivity
M3 <- glm(Fij ~ Landscape_connectivity, data = gendata, family = gaussian)
summary(M3)

# Comparing models and selecting variables -> only Euclidean distance was kept!
AICtab(M1, M2, M3, base = T, delta = T, weights = T)

#########
## Relationships between the response variable and the explanatory variables

par(mfrow=c(2,2))
plot(gendata$Fij ~ gendata$Management, ylab = "Fij", xlab = "Management")  
plot(gendata$Fij ~ gendata$Roads, ylab = "Fij", xlab = "Presence of roads")
plot(gendata$Fij ~ gendata$Euclidean_distance, ylab = "Fij", xlab = "Euclidean distance")
plot(gendata$Fij ~ gendata$Landscape_resistance, ylab = "Fij", xlab = "Landscape resistance")
par(mfrow=c(1,1))

# Distribution of the response variable
# Fij is nearly Gaussian - we will use the Gaussian family when fitting GLMs
hist(gendata$Fij)

############
# Exploring residuals
M.full <- glm(Fij ~ Management-1 + Euclidean_distance + Roads-1 + Landscape_resistance, data=gendata, family=gaussian)
summary(M.full)

# Residuals
op <- par(mfrow = c(2, 2))
plot(M1) 
par(op)

# Normality 
E1 <- rstandard(M1)
hist(E1) #normal 
qqnorm(E1)

# Check for independence and homogeneity: residuals
op <- par(mfrow = c(2, 2))
# Management
plot(y = E1, x = gendata$Management, xlab = "Management", ylab = "Residuals")
abline(0, 0, col = 2)

# Roads
plot(y = E1, x = gendata$Roads, xlab = "Roads", ylab = "Residuals")
abline(0, 0, col = 2)

# Euclidean_distance
plot(y = E1, x = gendata$Euclidean_distance, xlab = "Euclidean_distance", ylab = "Residuals") # high variation
abline(0, 0, col = 2)

# Landscape resistance
plot(E1~gendata$Landscape_resistance, xlab = "Landscape_resistance", ylab = "Residuals") #alta variacao
abline(0, 0, col = 2)
par(op)

############################################################
# Analyzing data 
# format: response ~ explanantory variables
# -1 only for categorical variables, to easy the interpretation: Management and Roads
# Response variable: Pairwise kinship, Fij

#############
# Univariable models

# Null model
M0F <- glm(Fij ~ 1, data = gendata, family = gaussian)
M0F
summary(M0F)

# Management
M1.man <- glm(Fij ~ Management - 1, data = gendata, family = gaussian) # coefficients show group means
M1.man
summary(M1.man)

# Roads
M1.road <- glm(Fij ~ Roads - 1, data = gendata, family = gaussian) # coefficients show group means
M1.road
summary(M1.road)

# DE
M1.DE <- glm(Fij ~ Euclidean_distance, data = gendata, family = gaussian)
M1.DE
summary(M1.DE)

# Landscape_resistance
M1.res <- glm(Fij ~ Landscape_resistance, data = gendata, family=gaussian)
M1.res
summary(M1.res)

#############
# Comparing univariate models  
AIC.uni <- AICtab(M0F, M1.man, M1.DE, M1.road, M1.res, base = T, delta = T, weights = T)
AIC.uni

#############
# Multivariate models

# Euclidean_distance + Management 
M1.DE_man <- glm(Fij ~ Euclidean_distance + Management - 1, data = gendata, family = gaussian) 
M1.DE_man
summary(M1.DE_man)

# Roads + management 
M1.road_man <- glm(Fij ~ Roads - 1 + Management - 1, data = gendata, family = gaussian)
M1.road_man
summary(M1.road_man)

# Landscape_resistance + management
M1.res_man <- glm(Fij ~ Landscape_resistance + Management - 1, data = gendata, family = gaussian) 
M1.res_man
summary(M1.res_man)

# Comparing multivariate models for Fij 
AIC_multi <- AICtab(M0F, M1.man, M1.DE, M1.road, M1.res, M1.DE_man, M1.road_man, M1.res_man, base = T, delta = T, weights = T)
AIC_multi

# Model results
summary(M1.res_man)
summary(M1.road_man)
summary(M1.DE_man)
summary(M1.man)

################################
# Representing the results in figures

tiff("fig_landscape_genetics.tif", width = 112, height = 120, units = "mm", res = 300)

plot <- plot(Fij ~ Landscape_resistance, pch=19, col="grey", cex=0.3, ylab= "Pairwise Kinship", xlab= "Landscape Resistance", cex.lab=1, font.lab=2, family="serif", data = gendata)
# cols <- palette(gray(0:10 / 10))
# cols <- rainbow(10)
cols <- c("red", "red", "red", "red", "green", "green", "green", "blue", "blue", "purple")
# line_type < -c(1:10)
line_type <- c(1, 2, 3, 4, 1, 3, 4, 1, 4, 1)

for(i in 1:(length(coef(M1.res_man))-1))
{
  abline(coef(M1.res_man)[i+1], coef(M1.res_man)[1], col = cols[i], lwd = 2, lty = line_type[i])
}

legend("topright", legend = c("native-native", "native-reintro", "native-translo", "native-unknown", 
                              "reintro-reintro", "reintro-translo", "reintro-unknown", 
                              "translo-translo", "translo-unknown", "unknown-unknown"), 
       col = cols, lwd = 1.5, ncol = 2,
       cex = 0.6, title = "MANAGEMENT CATEGORIES", bty = "n",
       lty = line_type)
dev.off()

