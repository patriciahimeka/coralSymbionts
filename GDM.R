### Generlaized Dissimilarity Modeling GDM ###

rm(list = ls())
install.packages("gdm")
library(gdm)
library(pcaMethods)
library(GUniFrac)
setwd("/path")

#################

## DATA ## 

taxatable <- read.table("taxatable_171027")  # MAIN
taxatable <- taxatable[-c(21),]  # remove "C10a" (not Symbdionium)


## Predictors 

predictors <- read.table(file.choose(), header=TRUE) # "sampledata_final_190730.txt"

# Transform categorial variable bleaching in to numerical 
predictors$bleachingNum <- factor(predictors$bleaching)
levels(predictors$bleachingNum) <- c(2,1,3)
predictors$bleachingNum <- as.numeric(levels(predictors$bleachingNum)[predictors$bleachingNum])

# Transform categorical variables clade into pc coordinates
modelmatrix = model.matrix( ~ clade1, data =predictors) [,-1]  
cov.mat.t=cov(t(modelmatrix[,1:2]),use='pairwise')  
dat.t <- prep(cov.mat.t, scale="none", center=TRUE)
resPCA.t <- pca(dat.t, method="svd", center=FALSE, nPcs=2) 
pointsPC <-loadings(resPCA.t)
colnames(pointsPC) <- c("pct1","pct2")


##  Community data 

otus <- read.table("MATRIX-COUNT_sorted.txt", header=T)
otus <- otus[, -c(16)] #remove C10a
rownames(otus) <- otus[,1]
otus <- otus[,-c(1)] #remove sample column
colnames(otus) <- taxatable$label

# For analysis on 47 samples delete nodes not present in samples 1:47:
#colnames(data) <-taxatable$label
#data <- as.data.frame(data[,order(colnames(data))])
# Remove otus of count <880 
#data <- within(data, rm("C21a_n142" ,"C21a_n578","C21a_n2638" ,"C3_n2875" ,"D1_n1059","D1_n1128","D1_n1310","D1_n2248","D1_n2459" ,"D1_n2592","D1_n31","D1a_n2451",
#                        "D1_n1060","D1_n1127","D1_n1329", "D1_n1331","D1_n1932","D1_n1896", "D1_n1899","D1_n2286","D1_n2466","D1_n2598","D2_n1167","D2_n1168"))

# rarefy counts to the minimum
matrix <- data.matrix(otus[1:67,1:ncol(otus)], rownames.force = T)   # select samples 1:47 or 1:67
matrix_raref <- Rarefy(matrix, depth = min(rowSums(matrix)))
matrix_rare <- matrix_raref$otu.tab.rff


############################################################################################################
##  GDM

inputBiodata <- as.data.frame(matrix_rare[1:67,])  # select samples 1:67 or 1:47
inputBiodata[,length(inputBiodata)+1] <- sampledata[1:67,"sample"]  #235 for all, 211 for only 47 samples  
names(inputBiodata)[length(inputBiodata)] <- "sample"   
head(inputBiodata)

# for 67 samples:
inputPredata <- cbind(droplevels(sampledata[1:67,c("sample", "lat", "lon", "SSTmean","depth")]),pointsPC[1:67,c("pct1","pct2")]) #need to be numeric
rownames(inputPredata) <- sampledata[1:67, "sample"]

# for 47 samples:
#inputPredata <- cbind(droplevels(sampledata[1:47,c("sample", "lat", "lon", "SSTmean","TSitu", "size","bleachingNum","depth")]),pointsPC[1:47,c("pct1","pct2")]) #need to be numeric
#rownames(inputPredata) <- sampledata[1:47, "sample"]

head(inputPredata)

###################
sitepairs <-formatsitepair(inputBiodata, 1, abundance=T, siteColumn = "sample",
                XColumn="lon", YColumn="lat", predData=inputPredata)

gdm.1 <- gdm(sitepairs, geo=T)
gdm1Summary <-summary.gdm(gdm.1)

# significance 
gdm.varImp(sitepairs, geo=T, nPerm=100, fullModelOnly = F)

# spline plots
length(gdm.1$predictors)
dev.off()
plot(gdm.1, plot.layout  =c(3,3))



