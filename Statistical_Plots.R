#Ibrahim Adel 20178001
#Haydey Ali 20188056

# ----install_packages---------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install()
BiocManager::install("Biobase")

BiocManager::install("antiProfilesData")

data=antiProfilesData::apColonData

class(data)
edata<-exprs(data)
pdata<-pData(data)
fdata<-fData(data)
# dim(data)
# dim(edata)
# dim(pdata)
# dim(fdata)

#################################################
#Q1_A

# str(data)
str(edata)
sapply(as.data.frame(edata),class)

str(pdata)
sapply(pdata,class)

str(fdata)
sapply(fdata,class)

#################################################
#Q1_B

colnames(edata)
row.names(edata)

colnames(pdata)
row.names(pdata)

colnames(fdata)
row.names(fdata)

#################################################
#Q1_C

sapply(as.data.frame(edata), summary)

sapply(pdata, summary)

sapply(fdata, summary)

#################################################
#Q1_D

# table(complete.cases(edata))
# as.data.frame(table(edata,exclude=NULL))

#run table if you want to see without freq 
table(complete.cases(pdata))
as.data.frame(table(pdata$ClinicalGroup,exclude=NULL))
as.data.frame(table(pdata$Tissue,exclude=NULL))
as.data.frame(table(pdata$ExperimentID,exclude=NULL))
as.data.frame(table(pdata$SubType,exclude=NULL))
as.data.frame(table(pdata$Status,exclude=NULL))

# table(complete.cases(fdata))
# as.data.frame(table(fdata,exclude=NULL))

#################################################
#Q1_E
install.packages("devtools")

# load the libraries
library(devtools)
library(Biobase)
library(broom)
library(corrplot)

cov<-cov(edata[,1:10],use = "everything", method = c("pearson", "kendall", "spearman"))
cov

cor<-cor(edata[,1:10],use = "everything", method = c("pearson", "kendall", "spearman"))
cor
jpeg()
corrplot(cor)
dev.off()

#################################################
#Q1_F

# cor2<-cor.test(edata[,"GSM95478"],edata[,"GSM95473"])
g1<-edata[,"GSM95478"]
g2<-edata[,"GSM95473"]
newdata<-as.data.frame(edata)

dev.copy(device = jpeg, file = "Relation plot")
plot(g1, g2, main = "Relation plot",
     xlab = "GSM95478 ", ylab = "GSM95473 ",
     pch = 19, frame = FALSE)

abline(lm(g2 ~ g1, data = newdata), col = "red",lw=7)
dev.off()

#################################################
#Q2

norm_edata = t(t(edata) - colMeans(edata))

svd <- svd(norm_edata)
svd
names(svd)

pca <- prcomp(norm_edata)
pca
names(pca)

par(mfrow=c(1,2))

jpeg()

plot(svd$v[,1],svd$v[,2], main = "SVD1", xlab = "v1", ylab = "v2")
plot(pca$rotation[,1],pca$rotation[,2], main = "PCA1", xlab = "rotation1", ylab = "rotation2")

plot(svd$u[,1],svd$u[,2],main = "SVD2", xlab = "U1", ylab = "U2")
plot(pca$x[,1],pca$x[,2], main = "PCA2", xlab = "x1", ylab = "x2")

plot(pca$rotation[, 1], svd$v[, 1], pch = 19, xlab = "pca$rotation", ylab = "svd$v")
abline(c(0, 1),col="red",lw=3)

dev.off()

#################################################
#Q3

category <- c("Aries","Taurus","Gemini","Cancer","Leo","Virgo","Libra","Scorpio","Sagittarius","Capricorn",
              "Aquarius","Pisces")

observed <- c(29,24,22,19,21,18,19,20,23,18,20,23)

expected <- c(256/12)

residual <- c(observed - expected)

obs_exp<- (residual^2)

component <- c(obs_exp/expected)

zodiac_sign <- (data.frame(category,observed,expected,residual,obs_exp,component))
zodiac_sign

H0<-chisq.test(zodiac_sign$component, correct = FALSE)
H0

H1<-chisq.test(zodiac_sign$observed, correct = FALSE)
H1

#p-value=0.9265(very large),which p-value>0.05 so we should not reject H0 and reject H1.

#################################################
#Q4

dist1<- dist(t(edata[,1:10]))
hclust1<-hclust(dist1)

dev.copy(device = jpeg, file = "hclust1 plot")
plot(hclust1,hang = -1)
dev.off()

set.seed(10)

kmeans<- kmeans(edata,centers=2) #Elbow Method
kmeans
table(kmeans$cluster)
kmeans$centers

install.packages("FactoMineR")
install.packages("factoextra")
library("factoextra")
library("FactoMineR")

dev.copy(device = jpeg, file = "kmeans plot")
fviz_cluster(kmeans, data = edata, xlab = "x-axis label", ylab = "y-axis label")
dev.off()