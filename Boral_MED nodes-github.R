# Boral for 97% clustered OTUs
rm(list = ls())
library("boral")
library("GUniFrac")
library(pheatmap)
library(ggplot2)
library(mvabund)
library(corrplot) 
library(reshape2)
setwd("/path")  # set working directory

#######################################################
## DATA ##

## Predictors ##

sampledata <- read.table(file.choose(), header=TRUE) # "sampledata_final_190730.txt"
predictors <- droplevels(sampledata)  

## The distances between little (2), no (1) and yes (3) bleaching can be represented by equally spaced values 1,2,3, so can be regarded as continuous 
predictors$bleachingNum <- factor(predictors$bleaching)
levels(predictors$bleachingNum) <- c(2,1,3)
predictors$bleachingNum <- as.numeric(levels(predictors$bleachingNum)[predictors$bleachingNum])
#str(predictors)


## Community matrix ##

otus <- read.table("MATRIX-COUNT_sorted.txt", header=TRUE) 
rownames(otus) <- otus$sample

data <- otus[1:47,-c(1)] #remove sample column, select samples (1-47 or 1-67)
data <- data[,-c(21)]  # remove C10a, because is not Symbiodinium

taxatable <- read.table("taxatable_171027") 
taxatable <- taxatable[-c(21),]  # remove "C10a" (is not Symbiodinium)
colnames(data) <-taxatable$label

## for 47 samples, remove nodes that are only present in Zanpa or Iriomote:
#data <- within(data, rm("C21a_n142" ,"C21a_n578","C21a_n2638" ,"C3_n2875" ,"D1_n1059","D1_n1128","D1_n1310","D1_n2248","D1_n2459" ,"D1_n2592","D1_n31","D1a_n2451",
#                        "D1_n1060","D1_n1127","D1_n1329", "D1_n1331","D1_n1932","D1_n1896", "D1_n1899","D1_n2286","D1_n2466","D1_n2598","D2_n1167","D2_n1168"))



#######################################################
## FIT MODELS WITHOUT predictors X ##

## Fit the pure latent variable model (but include log sequencing depth cause it acts as sampling effort) 
## to perform model-based unconstrained ordination and perform an exploratory analysis 
## Widen the hypparameters to use flatter prior distributions due to large counts

fit.boral.n <- boral(y=sub2, X = log(predictors$seqDepth), family = "negative.binomial", num.lv = 2, row.eff="none", prior.control = list(hypparams = c(100,100,100,30)), lv.control = list(num.lv = 2))
save(fit.boral.n, file="purelvmfit.n.67samples-MED.RData")
#load(file = "fit.nb.67samples-MED.RData")

summary(fit.boral.n)
par(mfrow=c(2,2))
plot(fit.boral.n)
par(mfrow=c(1,1))
fit.boral.n$geweke.diag    


###########################################################
## Create model matrix X  ###

## for 67 samples:
X67 <- model.matrix(~   depth + SSTmean + clade1 +seqDepth, data = predictors)[,-1] 
## Scale continuousv variables except for seqdepth, logarithmically transform seqdepth
X67[,1:2] <- scale(X[,1:2]) 
X67[,5] <- log(X[,5]) 

## for 47 samples:
#X47 <- model.matrix(~  size + bleachingNum + depth + SSTmean +TSitu+ clade1+seqDepth, data = predictors[1:47,])[,-1]
#X47[,1:5] <- scale(X[,1:5]) 
#X47[,8] <- log(X[,8]) 

head(X67)


#######################################################
## FIT MODELS WITH X ##

## Fit model with predictors and LVs
## Widen the hypparameters to use flatter prior distributions due to large counts
## Use negative-binomial distribution, sub2

## 67 samples with 4 predictors (repeat for case of 47 samples with X47) 
fit.boral.nb.X <- boral(y=data, X = X67, family = "negative.binomial", num.lv = 2, save.model=TRUE, prior.control = list(hypparams = c(100,100,100,30)), lv.control = list(num.lv = 2))
save(fit.boral.bin.X, file="fit.nb.67samples.MED.RData")
#load(file="fit.nb.67samples.MED.RData") 

summary(fit.boral.nb.X)
par(mfrow=c(2,2))
plot(fit.boral.nb.X)
par(mfrow=c(1,1))
fit.boral.nb.X$geweke.diag


## CORRELATIONS ##

envcorSig <- envcors$sig.cor
rescorSig <- rescors$sig.cor
envcorSigPh <- envcorSig  
envcorSigPh[1,1] <- -1  # -1 needs to be present in matrix so that colorscale is fromed properly in pheatmap

mat = rescorSig
mat = envcorSigPh
hc = hclust(as.dist(1-mat), method="centroid") 
mat = mat[hc$order, hc$order]
mat[upper.tri(mat,diag=F)] = NA

colourlist <- list( Type = c(C1 = "#3333FF", C1002 = "#0080FF", C1005="#66B2FF", C1013="#004C99", C1060="#00CCCC",
                             C1085 ="#003366", C1148 ="#99CCFF", C1226="#000099", C1228="#99FFFF",  C10a = "#0066CC",
                             C1230 ="#6600CC", C1234 ="#9933FF", "C1b/C1e" ="#FF33FF", "C1c/C45" ="#CC00CC",
                             C1h="#330066", C1p="#4C0099", "C21/3d/C3k"="#009900", C21a="#006600", "C27/C30" = "#00CC00",
                             C3="#FFFF66",C3b="#FFFFCC", C7="#990000",  
                             D1="#FF8000", D106="#CC0000", D1a="#FF3333", D2="#FF6666")) 

# Set annotation colours
# MED colums ('species'):
colAno <- as.data.frame(taxatable$type) #get column for annotation by pheatmap
rownames(colAno) <- taxatable$label
colnames(colAno) <- "Type"

pheatmap(mat, cluster_col = F, cluster_row = F, show_rownames=F,
         annotation_row=colAno, # must not have variable specification, but contain identical column name as coulourlist
         annotation_colors= colourlist, border_color = NA,
         color = rev(colorRampPalette(c("#004C99","#3399FF","white","#FF8000","#CC0000"))(99)))


#######################################################
## TRACE PLOTS ## 
# (Repeat for 67 and 47-samples case)
# Reorder OTUs alphabetically for plotting 
fit.boral.nb.X_reordered <- fit.boral.nb.X
hpintervals <- fit.boral.nb.X$hpdintervals$X.coefs
medians <- fit.boral.nb.X$X.coefs.median
fit.boral.nb.X_reordered$hpdintervals$X.coefs <- fit.boral.nb.X$hpdintervals$X.coefs[order(rownames(hpintervals)),,]
fit.boral.nb.X_reordered$X.coefs.median <- fit.boral.nb.X$X.coefs.median[order(rownames(medians)),]

fittedPar <- fitted.boral(fit.boral.nb.X, est="median")

par(mfrow=c(3,3))
coefplotSST <- coefsplot("SSTmean",fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("TSitu",fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("bleachingNum",fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("size", fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("depth", fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("clade1L+",fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("clade1S",fit.boral.nb.X_reordered, cex.axis=0.7)
coefsplot("seqDepth", fit.boral.nb.X_reordered, cex.axis=0.7)
par(mfrow=c(1,1))