# Boral for 97% clustered OTUs
rm(list = ls())
library("boral")
library("GUniFrac")
library(ggplot2)
library(mvabund)
library(corrplot) 
library(reshape2)
setwd("/path")  # set working directory

#######################################################
## DATA ##

## Predictors ##

sampledata <- read.table(file.choose(), header=TRUE) # "samplingdata_181006.txt"
predictors <- droplevels(sampledata)  

# Transform categorial variable bleaching in to numerical 
# The distances between little (2), no (1) and yes (3) bleaching can be represented by equally spaced values 1,2,3, so can be regarded as continuous 
predictors$bleachingNum <- factor(predictors$bleaching)
levels(predictors$bleachingNum) <- c(2,1,3)
predictors$bleachingNum <- as.numeric(levels(predictors$bleachingNum)[predictors$bleachingNum])
#str(predictors)


## Community matrix ##

otus <- read.table(file.choose(), header=TRUE)  # "97_individ_cleaned_190730.txt", purged 1/10000
data <- otus[1:47,-c(1)] #remove sample column, select samples (1-47 or all (67))
rownames(data) <- predictors$sample[1:47]

# Check taxa prevalance: Identify taxa with too few species (OTUs) to be removed 
colSums(data > 0)
table(colSums(data > 0))
# Identify OTUs present in too few (1) specimen
rowSums(data > 0)

#######################################################
## FIT MODELS WITHOUT predictors X ##

## Fit the pure latent variable model to perform model-based unconstrained ordination and perform an exploratory analysis 
## Widen the hypparameters to use flatter prior distributions due to large counts

# make BINOMIAL
sub2 <- data[,-which(colSums(data > 1) <= 1)]  # Remove taxa present in only 1 specimen
names(sub2)
sub2pa <- as.matrix((sub2 > 0)+0) # transform to presence-absence
head(sub2pa)


fit.boral.bin <- boral(y=sub2pa, family = "binomial", num.lv = 2, row.eff="none", save.model=TRUE)
#save(fit.boral.n, fit.boral.bin, file="purelvmfit.67samples.181204.RData")
#load(file = "purelvmfit.RData")

summary(fit.boral.bin)
par(mfrow=c(2,2))
plot(fit.boral.bin)
par(mfrow=c(1,1))
fit.boral.bin$geweke.diag

# ordination with latent variables
getlvs <- lvsplot(fit.boral.bin, alpha = 0.6, return.vals = TRUE, ind.spp = 10)


###########################################################
## Explore predictors and create model matrix X  ###

## Explore predictors
table(factor(predictors$SSTmean), predictors$location)
table(factor(predictors$TSitu), predictors$location)
# Drop location from the model due to complete confounding with SSTmean and TSitu

## for 67 samples:
X <- model.matrix(~ depth + SSTmean + clade1, data = predictors)[,-1] 
## for 47 samples:
# X <- model.matrix(~  size + bleachingNum + depth + SSTmean +TSitu+ clade1, data = predictors[1:47,])[,-1]

## Scale continuousv variables 
X[,1:5] <- scale(X[,1:5]) 
head(X)


#######################################################
## FIT MODELS WITH X ##

## Fit model with each covariate and LVs
## Widen the hypparameters to use flatter prior distributions due to large counts
## Use binomial distribution and presence-absence data (sub2pa)

fit.boral.bin.X <- boral(y=sub2pa, X = X, family = "binomial", num.lv = 2, save.model=TRUE)
save(fit.boral.bin.X, file="fit.boral.bin.X.RData")
#load(file="fit.X.RData") 

summary(fit.boral.bin.X)
par(mfrow=c(2,2))
plot(fit.boral.bin.X)
par(mfrow=c(1,1))
fit.boral.bin.X$geweke.diag

## CORRELATIONS ##

envcors <- get.enviro.cor(fit.boral.bin.X) 
rescors <- get.residual.cor(fit.boral.bin.X) 

par(mfrow=c(1,2))
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, title = "Co-occurrence due to covariates", 
         mar = c(1,0.5,1,0.5), tl.srt = 45, order="hclust", tl.cex=0.6, tl.col=c(1)) #,addgrid.col=NA )
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, title = "Co-occurrence due to latent variables",
         mar = c(1,0.5,1,0.5), tl.srt = 45, order="hclust", tl.cex=0.6,tl.col=c(1)) #, addgrid.col=NA)
par(mfrow=c(1,1))


## VARIATION PARTITIONING ##

colors <- c("black","lightgrey")

varpart <- calc.varpart(fit.boral.bin.X, groupX= c(1,1,1,1,1))   # add plus1 as intercept in beginning
varpart_matrix <- as.data.frame(varpart$varpart.X)
varpart_matrix[2,] <- t(as.vector(varpart$varpart.lv))
rownames(varpart_matrix) <- c("Predictors","LV") 
varpart_matrix_sort <- varpart_matrix[,order(colnames(varpart_matrix))]

dat.m <- melt(t(varpart_matrix_sort), id.vars= c( "Predictors","LV")) 
colnames(dat.m) <- c("Nodes", "Effects", "Variance")

ggplot(dat.m, aes(x=Nodes, y=Variance, fill=Effects)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle= 90, hjust=1)) +
  scale_fill_manual(values= colors) +
  geom_bar(stat='identity')

#######################################################
## FIT MODELS WITH X EACH INDIVIDUALLY ##

## Choose individual predictors for separate fitting
Xs <- X[,1]  # depth
Xs <- X[,2]  # SSTmean 
Xs <- X[,3:4] # clade1

fit.boral.bin.Xi <- boral(y=sub2pa, X = Xs, family = "binomial", num.lv = 2, save.model=TRUE, x)
save(fit.boral.bin.Xi, file="fit.bin.67samples.depth-190731.RData")
#load(file="fit.bin.67samples.depth-190731.RData")

## VARIATION PARTITIONING ##

# repeat for each each predictor:
varpart <- calc.varpart(fit.boral.bin.Xi, groupX= c(1,1,1))   # predictors plus 1 for intercept at the beginning
varpart_matrix <- as.data.frame(varpart$varpart.X)
varpart_matrix[2,] <- t(as.vector(varpart$varpart.lv))
rownames(varpart_matrix) <- c("Predictor","LV") 
varpart_matrix_sort <- varpart_matrix[,order(colnames(varpart_matrix))]

t(varpart_matrix_sort)


#######################################################
# Trace plots
# Reorder OTUs alphabetically for plotting (if not done so before)
fit.boral.bin.X_reordered <- fit.boral.bin.X
hpintervals <- fit.boral.bin.X$hpdintervals$X.coefs
medians <- fit.boral.bin.X$X.coefs.median
fit.boral.bin.X_reordered$hpdintervals$X.coefs <- fit.boral.bin.X$hpdintervals$X.coefs[order(rownames(hpintervals)),,]
fit.boral.bin.X_reordered$X.coefs.median <- fit.boral.bin.X$X.coefs.median[order(rownames(medians)),]

fittedPar <- fitted.boral(fit.boral.bin.X, est="median")

par(mfrow=c(3,3))
coefsplot("SSTmean",fit.boral.bin.X_reordered, cex.axis=0.7)
coefsplot("TSitu",fit.boral.bin.X_reordered, cex.axis=0.7)
coefsplot("bleachingNum",fit.boral.bin.X_reordered, cex.axis=0.7)
coefsplot("size", fit.boral.bin.X_reordered, cex.axis=0.7)
coefsplot("depth", fit.boral.bin.X_reordered, cex.axis=0.7)
coefsplot("clade1L+",fit.boral.bin.X_reordered, cex.axis=0.7)
coefsplot("clade1S",fit.boral.bin.X_reordered, cex.axis=0.7)
par(mfrow=c(1,1))


