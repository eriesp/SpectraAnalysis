################################################################################
#######                Libraries and dataset loading                     #######
################################################################################

# Very important color initialization
library("RColorBrewer")

# Import the datasets
library(readxl)
library(fda)
library(KernSmooth)
library(fields)
library(fdakma)
library(roahd)

srp <- read_excel("Dataset/Srp_Emissivity.xlsx")
head(srp)
View(srp)
summary(srp)

vul <- read_excel("Dataset/Vulcano_Emissivity.xlsx")
head(vul)
View(vul)
summary(vul)

srp2 <- read_excel("Dataset/Srp_Reflectance.xlsx")
head(srp2)
View(srp2)
summary(srp2)

vul2 <- read_excel("Dataset/Vulcano_Reflectance.xlsx")
head(vul2)
View(vul2)
summary(vul2)

###########
### FDA ###
###########

srp.fda=srp[417:1, -c(1)]
srp2.fda=srp2[298:714, -1]
vul.fda=vul[417:1,-1]
vul2.fda=vul2[298:714,-1]
emis.x=srp[417:1,1]$`Wavelength (micron)`
refl.x=srp2[298:714,1]$wl # si può tagliare ancora un pochino volendo

data_emis=data.frame(srp.fda, vul.fda)
data_refl=data.frame(srp2.fda, vul2.fda)

x11()
par(mfrow=c(2,2))
matplot(emis.x, srp.fda, type='l')
matplot(emis.x, vul.fda, type='l')
matplot(refl.x, srp2.fda, type='l')
matplot(refl.x, vul2.fda, type='l')

x11()
par(mfrow=c(1,2))
matplot(emis.x, data_emis, type='l')
matplot(refl.x, data_refl, type='l')

m=5

######################
##### EMISSIVITY #####

nbasis <- 6:30
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(emis.x[1],emis.x[length(emis.x)]), nbasis[i], m)
  gcv[i] <- smooth.basis(emis.x, as.matrix(data_emis), basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

nbasis=25
basis <- create.bspline.basis(c(emis.x[1],emis.x[length(emis.x)]), nbasis, m)

basismat <- eval.basis(emis.x, basis)
basismat1<- eval.basis(emis.x, basis, Lfdobj=1)
basismat2<- eval.basis(emis.x, basis, Lfdobj=2)


emis.smooth <- basismat %*% lsfit(basismat, data_emis, intercept=FALSE)$coef
emis.der1 <- basismat1 %*% lsfit(basismat, data_emis, intercept=FALSE)$coef
emis.der2 <- basismat2 %*% lsfit(basismat, data_emis, intercept=FALSE)$coef

#srp.CF=numeric(23)
#srp.TF=numeric(23)
#srp.TFval=numeric(23)
#srp.last = numeric(23)
#for(i in 1:23) {
#  srp.CF[i]=which.max(srp.smooth[1:130,i])
#  srp.CF[i]=emis.x[srp.CF[i]]
#  srp.TF[i]=which.min(srp.smooth[130:250,i]) + 130
#  srp.TFval[i]=srp.smooth[srp.TF[i],i]
#  srp.TF[i]=emis.x[srp.TF[i]]
#  srp.last[i] = srp.smooth[length(emis.x),i]
#}

par(mfrow=c(2,2))
matplot(emis.x,data_emis,xlab="t",ylab="observed data", type='l')
matplot(emis.x,emis.smooth ,type="l",lwd=1)
matplot(emis.x,emis.der1 ,type="l",lwd=1)
matplot(emis.x,emis.der2 ,type="l",lwd=1)

dev.off()

"PCA"

data2fd_emis <- Data2fd(y = emis.smooth,argvals = emis.x, basisobj = basis)

pca_emis <- pca.fd(data2fd_emis,nharm=5,centerfns=TRUE)

# scree plot
plot(pca_emis$values[1:nbasis],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_emis$values)[1:nbasis]/sum(pca_emis$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

par(mfrow=c(1,2))
plot.pca.fd(pca_emis, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)
# la componente principale 1 ci suggerisce che la maggior parte della variabilità
# è concentrata nella tranparency feature e nella parte di funzione che segue

par(mfrow=c(1,2))
plot(pca_emis$scores[,1],pca_emis$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
##points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)
plot(pca_emis$scores[,1],pca_emis$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_emis$scores[,1],pca_emis$scores[,2],dimnames(data_emis)[[2]], cex=1)

layout(1)
matplot(emis.x, emis.smooth,type='l')
lines(emis.x,emis.smooth[,24],lwd=4, col=2)



#######################
##### REFLECTANCE #####

nbasis <- 6:50
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), nbasis[i], m)
  gcv[i] <- smooth.basis(refl.x, as.matrix(data_refl), basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

nbasis=22
basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), nbasis, m)

basismat <- eval.basis(refl.x, basis)
basismat1<- eval.basis(refl.x, basis, Lfdobj=1)
basismat2<- eval.basis(refl.x, basis, Lfdobj=2)


refl.smooth <- basismat %*% lsfit(basismat,data_refl, intercept=FALSE)$coef
refl.der1 <- basismat1 %*% lsfit(basismat, data_refl, intercept=FALSE)$coef
refl.der2 <- basismat2 %*% lsfit(basismat, data_refl, intercept=FALSE)$coef

#refl.CF=numeric(12)
#refl.TF=numeric(12)
#refl.TFval=numeric(12)
#refl.last=numeric(12)
#for(i in 1:12) {
#  refl.CF[i]=which.max(refl.smooth[1:130,i]) 
#  refl.CF[i]=refl.x[refl.CF[i]]
#  refl.TF[i]=which.min(refl.smooth[130:250,i]) + 130
#  refl.TFval[i]=refl.smooth[refl.TF[i],i]
#  refl.TF[i]=refl.x[refl.TF[i]]
#  refl.last[i]=refl.smooth[length(refl.x),i]
#}

x11()
par(mfrow=c(2,2))
matplot(refl.x,data_refl,xlab="t",ylab="observed data", type='l')
matplot(refl.x,refl.smooth ,type="l",lwd=1)
matplot(refl.x,refl.der1 ,type="l",lwd=1)
matplot(refl.x,refl.der2 ,type="l",lwd=1)

dev.off()

"PCA"

data2fd_refl <- Data2fd(y = refl.smooth,argvals = refl.x, basisobj = basis)

pca_refl <- pca.fd(data2fd_refl,nharm=5,centerfns=TRUE)

# scree plot
plot(pca_refl$values[1:12],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_refl$values)[1:12]/sum(pca_refl$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

par(mfrow=c(1,2))
plot.pca.fd(pca_refl, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)

par(mfrow=c(1,2))
plot(pca_refl$scores[,1],pca_refl$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
##points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)
plot(pca_refl$scores[,1],pca_refl$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_refl$scores[,1],pca_refl$scores[,2],dimnames(data_refl)[[2]], cex=1)

layout(1)
matplot(refl.x, refl.smooth,type='l')
lines(refl.x,refl.smooth[,1],lwd=4, col=2)
lines(refl.x,refl.smooth[,11],lwd=4, col=3)


##########################
#### OUTLIER DETECTION ###

f_data_emis = fData(emis.x, t(data_emis))

x11()
invisible(fbplot(f_data_emis, main="Magnitude outliers")) ### FUNCTIONAL BOXPLOT
invisible(outliergram(f_data_emis))                       ### OUTLIERGRAM

outliers=fbplot(f_data_emis)
outliers$ID_outliers                         ### !!!!






######################
### SRP EMISSIVITY ###

nbasis <- 6:100
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(emis.x[1],emis.x[length(emis.x)]), nbasis[i], m)
  gcv[i] <- smooth.basis(emis.x, as.matrix(srp.fda), basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

nbasis=25
basis <- create.bspline.basis(c(emis.x[1],emis.x[length(emis.x)]), nbasis, m)

basismat <- eval.basis(emis.x, basis)
basismat1<- eval.basis(emis.x, basis, Lfdobj=1)
basismat2<- eval.basis(emis.x, basis, Lfdobj=2)


srp.smooth <- basismat %*% lsfit(basismat, srp.fda, intercept=FALSE)$coef
srp.der1 <- basismat1 %*% lsfit(basismat, srp.fda, intercept=FALSE)$coef
srp.der2 <- basismat2 %*% lsfit(basismat, srp.fda, intercept=FALSE)$coef

srp.CF=numeric(23)
srp.TF=numeric(23)
srp.TFval=numeric(23)
srp.last = numeric(23)
for(i in 1:23) {
  srp.CF[i]=which.max(srp.smooth[1:130,i])
  srp.CF[i]=emis.x[srp.CF[i]]
  srp.TF[i]=which.min(srp.smooth[130:250,i]) + 130
  srp.TFval[i]=srp.smooth[srp.TF[i],i]
  srp.TF[i]=emis.x[srp.TF[i]]
  srp.last[i] = srp.smooth[length(emis.x),i]
}

par(mfrow=c(2,2))
matplot(emis.x,srp.fda,xlab="t",ylab="observed data", type='l')
matplot(emis.x,srp.smooth ,type="l",lwd=2)
matplot(emis.x,srp.der1 ,type="l",lwd=2)
matplot(emis.x,srp.der2 ,type="l",lwd=2)

dev.off()

"PCA"

data_srp <- Data2fd(y = srp.smooth,argvals = emis.x, basisobj = basis)

pca_srp <- pca.fd(data_srp,nharm=5,centerfns=TRUE)

# scree plot
plot(pca_srp$values[1:nbasis],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_srp$values)[1:nbasis]/sum(pca_srp$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

par(mfrow=c(1,2))
plot.pca.fd(pca_srp, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)
# la componente principale 1 ci suggerisce che la maggior parte della variabilità
# è concentrata nella tranparency feature e nella parte di funzione che segue

par(mfrow=c(1,2))
plot(pca_srp$scores[,1],pca_srp$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
##points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)
plot(pca_srp$scores[,1],pca_srp$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_srp$scores[,1],pca_srp$scores[,2],dimnames(srp.fda)[[2]], cex=1)

layout(1)
matplot(emis.x, srp.smooth,type='l')
#lines(emis.x,srp.smooth[,24],lwd=4, col=2)


"ALLIGNMENT"

kma.compare_srp <- kma.compare (
  x=emis.x, y0=t(srp.smooth), y1=t(srp.der1), 
  n.clust = 1:6, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means',
  plot.graph=1)

fdakma_srp <- kma(
  x=emis.x, y0=t(srp.smooth), y1=t(srp.der1), n.clust = 5, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means')          

kma.show.results(fdakma_srp)

fdakma_srp$labels # 1 1 2 2 2 5 | 1 1 2 2 2 5 | 1 3 1 1 3 5 | 4 4 4 4 3 1 (5 clust) x         600° fucks everything up
                
# new labels:

# 1 1 1 5 5 2 | 1 1 1 1 5 2 | 3 1 1 1 5 2 | 3 3 4 4 4 5 (5 clust)
# 4 5 5 1 1 1 | 4 4 5 5 1 1 | 4 5 2 2 1 1 | 3 3 3 3 2 3
# 5 5 5 3 4 4 | 5 5 5 5 4 4 | 1 2 5 5 3 4 | 1 2 1 1 3                in ogni caso continuano a cambiare...

fdakma_srp$shift
fdakma_srp$dilation

plot(emis.x, type="n", xlim=c(min(emis.x),max(emis.x)), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:24)(
  abline(a=fdakma_srp$shift[i],b=fdakma_srp$dilation[i], 
         col=fdakma_srp$labels[i])
)  

matplot(emis.x, srp.smooth[,which(fdakma_srp$labels==1)], type='l')
matplot(emis.x, srp.smooth[,which(fdakma_srp$labels==2)], type='l')
matplot(emis.x, srp.smooth[,which(fdakma_srp$labels==3)], type='l')
matplot(emis.x, srp.smooth[,which(fdakma_srp$labels==4)], type='l')
matplot(emis.x, srp.smooth[,which(fdakma_srp$labels==5)], type='l')
matplot(emis.x, srp.smooth[,which(fdakma_srp$labels==6)], type='l')

##########################
### VULCANO EMISSIVITY ###

nbasis <- 6:100
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(emis.x[1],emis.x[length(emis.x)]), nbasis[i], m)
  gcv[i] <- smooth.basis(emis.x, as.matrix(vul.fda), basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

nbasis=21
basis <- create.bspline.basis(c(emis.x[1],emis.x[length(emis.x)]), nbasis, m)

basismat <- eval.basis(emis.x, basis)
basismat1<- eval.basis(emis.x, basis, Lfdobj=1)
basismat2<- eval.basis(emis.x, basis, Lfdobj=2)


vul.smooth <- basismat %*% lsfit(basismat, vul.fda, intercept=FALSE)$coef
vul.der1 <- basismat1 %*% lsfit(basismat, vul.fda, intercept=FALSE)$coef
vul.der2 <- basismat2 %*% lsfit(basismat, vul.fda, intercept=FALSE)$coef

vul.CF=numeric(20)
vul.TF=numeric(20)
vul.TFval=numeric(20)
vul.last=numeric(20)
for(i in 1:20) {
  vul.CF[i]=which.max(vul.smooth[1:130,i])
  vul.CF[i]=emis.x[vul.CF[i]]
  vul.TF[i]=which.min(vul.smooth[130:250,i]) + 130
  vul.TFval[i]=vul.smooth[vul.TF[i],i]
  vul.TF[i]=emis.x[vul.TF[i]]
  vul.last[i]=vul.smooth[length(emis.x),i]
}

par(mfrow=c(2,2))
matplot(emis.x,vul.fda,xlab="t",ylab="observed data", type='l')
matplot(emis.x,vul.smooth ,type="l",lwd=2)
matplot(emis.x,vul.der1 ,type="l",lwd=2)
matplot(emis.x,vul.der2 ,type="l",lwd=2)

dev.off()

"PCA"

data_vul <- Data2fd(y = vul.smooth,argvals = emis.x, basisobj = basis)

pca_vul <- pca.fd(data_vul,nharm=5,centerfns=TRUE)

# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first 
# N-1=34 are non-null
plot(pca_vul$values[1:20],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_vul$values)[1:20]/sum(pca_vul$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

par(mfrow=c(1,2))
plot.pca.fd(pca_vul, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)

par(mfrow=c(1,2))
plot(pca_vul$scores[,1],pca_vul$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
##points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)
plot(pca_vul$scores[,1],pca_vul$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_vul$scores[,1],pca_vul$scores[,2],dimnames(vul.fda)[[2]], cex=1)

#layout(1)
#matplot(emis.x, srp.smooth,type='l')
#lines(emis.x,srp.smooth[,24],lwd=4, col=2)

"ALLIGNMENT"

kma.compare_vul <- kma.compare (
  x=emis.x, y0=t(vul.smooth), y1=t(vul.der1), 
  n.clust = 1:6, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means',
  plot.graph=1)

fdakma_vul <- kma(
  x=emis.x, y0=t(vul.smooth), y1=t(vul.der1), n.clust = 4, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means')

kma.show.results(fdakma_vul)

fdakma_vul$labels 
# 1 1 4 1 4 | 1 1 4 4 4 | 1 1 4 2 4 | 2 2 2 3 2 (4 clust) x           600° fucks everything up

# new labels
# 2 2 4 4 4 | 2 2 4 4 1 | 1 2 1 1 3 | 1 3 3 3 3 (4 clust)


fdakma_vul$shift
fdakma_vul$dilation

plot(emis.x, type="n", xlim=c(min(emis.x),max(emis.x)), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:20)(
  abline(a=fdakma_vul$shift[i],b=fdakma_vul$dilation[i], 
         col=fdakma_vul$labels[i])
)  

matplot(emis.x, vul.smooth[,which(fdakma_vul$labels==1)], type='l')
matplot(emis.x, vul.smooth[,which(fdakma_vul$labels==2)], type='l')
matplot(emis.x, vul.smooth[,which(fdakma_vul$labels==3)], type='l')
matplot(emis.x, vul.smooth[,which(fdakma_vul$labels==4)], type='l')


#######################
### SRP REFLECTANCE ###

nbasis <- 6:100
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), nbasis[i], m)
  gcv[i] <- smooth.basis(refl.x, as.matrix(srp2.fda), basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

nbasis=22
basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), nbasis, m)

basismat <- eval.basis(refl.x, basis)
basismat1<- eval.basis(refl.x, basis, Lfdobj=1)
basismat2<- eval.basis(refl.x, basis, Lfdobj=2)


srp2.smooth <- basismat %*% lsfit(basismat, srp2.fda, intercept=FALSE)$coef
srp2.der1 <- basismat1 %*% lsfit(basismat, srp2.fda, intercept=FALSE)$coef
srp2.der2 <- basismat2 %*% lsfit(basismat, srp2.fda, intercept=FALSE)$coef

srp2.CF=numeric(12)
srp2.TF=numeric(12)
srp2.TFval=numeric(12)
srp2.last=numeric(12)
for(i in 1:12) {
  srp2.CF[i]=which.max(srp2.smooth[1:130,i]) 
  srp2.CF[i]=refl.x[srp2.CF[i]]
  srp2.TF[i]=which.min(srp2.smooth[130:250,i]) + 130
  srp2.TFval[i]=srp2.smooth[srp2.TF[i],i]
  srp2.TF[i]=refl.x[srp2.TF[i]]
  srp2.last[i]=srp2.smooth[length(refl.x),i]
}

par(mfrow=c(2,2))
matplot(refl.x,srp2.fda,xlab="t",ylab="observed data", type='l')
matplot(refl.x,srp2.smooth ,type="l",lwd=2)
matplot(refl.x,srp2.der1 ,type="l",lwd=2)
matplot(refl.x,srp2.der2 ,type="l",lwd=2)

dev.off()

"PCA"

data_srp2 <- Data2fd(y = srp2.smooth,argvals = refl.x, basisobj = basis)

pca_srp2 <- pca.fd(data_srp2,nharm=5,centerfns=TRUE)

# scree plot
plot(pca_srp2$values[1:12],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_srp2$values)[1:12]/sum(pca_srp2$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

par(mfrow=c(1,3))
plot.pca.fd(pca_srp2, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)

par(mfrow=c(1,2))
plot(pca_srp2$scores[,1],pca_srp2$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
##points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)
plot(pca_srp2$scores[,1],pca_srp2$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_srp2$scores[,1],pca_srp2$scores[,2],dimnames(srp2.fda)[[2]], cex=1)

layout(1)
matplot(refl.x, srp2.smooth,type='l')
lines(refl.x,srp2.smooth[,1],lwd=4, col=2)
lines(refl.x,srp2.smooth[,11],lwd=4, col=3)

"ALLIGNMENT"

kma.compare_srp2 <- kma.compare (
  x=refl.x, y0=t(srp2.smooth), y1=t(srp2.der1), 
  n.clust = 1:6, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means',
  plot.graph=1)

fdakma_srp2 <- kma(
  x=refl.x, y0=t(srp2.smooth), y1=t(srp2.der1), n.clust = 5, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means')

kma.show.results(fdakma_srp2)

fdakma_srp2$labels  # 1 1 5 5 5 5 2 2 4 4 3 4 (5 clust) x            similar to clustering page 3


# new labels
# 2 4 5 5 5 5 1 1 1 1 3 1   1 1 4 4 4 4 3 3 3 3 5 2      1 3 3 3 2 2 2 5 5 5 4 5

# labels 5 clust shift, d0.L2.centered, k-medoids
# 1 4 | 4 4 | 5 5 | 5 3 | 5 2 | 4 2
# labels 5 clust shift, d0.L2, k-medoids
# 4 3 | 2 2 | 2 2 | 2 5 | 1 2 | 1 5            2 3 | 5 1 | 5 5 | 4 1 | 4 4 | 4 4          5 1 | 4 4 | 4 4 | 3 2 | 3 4 | 3 2
# 2 2 | 4 4 | 4 4 | 1 5 | 1 5 | 3 5
 
fdakma_srp2$shift
fdakma_srp2$dilation

plot(refl.x, type="n", xlim=c(min(refl.x),max(refl.x)), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:20)(
  abline(a=fdakma_srp2$shift[i],b=fdakma_srp2$dilation[i], 
         col=fdakma_srp2$labels[i])
)  

matplot(refl.x, srp2.smooth[,which(fdakma_srp2$labels==1)], type='l')
matplot(refl.x, srp2.smooth[,which(fdakma_srp2$labels==2)], type='l')
matplot(refl.x, srp2.smooth[,which(fdakma_srp2$labels==3)], type='l')
matplot(refl.x, srp2.smooth[,which(fdakma_srp2$labels==4)], type='l')
matplot(refl.x, srp2.smooth[,which(fdakma_srp2$labels==5)], type='l')


###########################
### VULCANO REFLECTANCE ###

nbasis <- 6:100
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), nbasis[i], m)
  gcv[i] <- smooth.basis(refl.x, as.matrix(vul2.fda), basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

nbasis=21
basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), nbasis, m)

basismat <- eval.basis(refl.x, basis)
basismat1<- eval.basis(refl.x, basis, Lfdobj=1)
basismat2<- eval.basis(refl.x, basis, Lfdobj=2)


vul2.smooth <- basismat %*% lsfit(basismat, vul2.fda, intercept=FALSE)$coef
vul2.der1 <- basismat1 %*% lsfit(basismat, vul2.fda, intercept=FALSE)$coef
vul2.der2 <- basismat2 %*% lsfit(basismat, vul2.fda, intercept=FALSE)$coef

vul2.CF=numeric(10)
vul2.TF=numeric(10)
vul2.TFval=numeric(10)
vul2.last = numeric(10)
for(i in 1:10) {
  vul2.CF[i]=which.max(vul2.smooth[1:130,i])
  vul2.CF[i]=refl.x[vul2.CF[i]]
  vul2.TF[i]=which.min(vul2.smooth[130:250,i]) + 130
  vul2.TFval[i]=vul2.smooth[vul2.TF[i],i]
  vul2.TF[i]=refl.x[vul2.TF[i]]
  vul2.last[i]=vul2.smooth[length(refl.x),i]
}

par(mfrow=c(2,2))
matplot(refl.x,vul2.fda,xlab="t",ylab="observed data", type='l')
matplot(refl.x,vul2.smooth ,type="l",lwd=2)
matplot(refl.x,vul2.der1 ,type="l",lwd=2)
matplot(refl.x,vul2.der2 ,type="l",lwd=2)

dev.off()

"PCA"

data_vul2 <- Data2fd(y = vul2.smooth,argvals = refl.x, basisobj = basis)

pca_vul2 <- pca.fd(data_vul2,nharm=5,centerfns=TRUE)

# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first 
# N-1=34 are non-null
plot(pca_vul2$values[1:10],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_vul2$values)[1:10]/sum(pca_vul2$values),xlab='j',ylab='CPV',ylim=c(0.6,1))

par(mfrow=c(1,2))
plot.pca.fd(pca_vul2, nx=100, pointplot=TRUE, harm=c(1,2), expand=0, cycle=FALSE)

par(mfrow=c(1,2))
plot(pca_vul2$scores[,1],pca_vul2$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)
##points(pca_W.1$scores[35,1],pca_W.1$scores[35,2],col=2, lwd=4)
plot(pca_vul2$scores[,1],pca_vul2$scores[,2],type="n",xlab="Scores FPC1",
     ylab="Scores FPC2")
text(pca_vul2$scores[,1],pca_vul2$scores[,2],dimnames(vul2.fda)[[2]], cex=1)

layout(1)
matplot(refl.x, vul2.smooth,type='l')
lines(refl.x,vul2.smooth[,9],lwd=4, col=2)


"ALLIGNMENT"

kma.compare_vul2 <- kma.compare (
  x=refl.x, y0=t(vul2.smooth), y1=t(vul2.der1), 
  n.clust = 1:5, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means',
  plot.graph=1)

fdakma_vul2 <- kma(
  x=refl.x, y0=t(vul2.smooth), y1=t(vul2.der1), n.clust = 4, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means')

kma.show.results(fdakma_vul2)

fdakma_vul2$labels # 4 4 2 2 3 3 1 1 1 1 (4 clust)          similar to clustering page 3


fdakma_vul2$shift
fdakma_vul2$dilation

plot(refl.x, type="n", xlim=c(min(refl.x),max(refl.x)), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:20)(
  abline(a=fdakma_vul2$shift[i],b=fdakma_vul2$dilation[i], 
         col=fdakma_vul2$labels[i])
)  

matplot(refl.x, vul2.smooth[,which(fdakma_vul2$labels==1)], type='l')
matplot(refl.x, vul2.smooth[,which(fdakma_vul2$labels==2)], type='l')
matplot(refl.x, vul2.smooth[,which(fdakma_vul2$labels==3)], type='l')
matplot(refl.x, vul2.smooth[,which(fdakma_vul2$labels==4)], type='l')
matplot(refl.x, vul2.smooth[,which(fdakma_vul2$labels==5)], type='l')


"Some analysis"

silica.srp = c(72.28, 67.47, 62.64, 57.07, 53.03, 48.03)
silica.vul = c(53.33, 59.66, 64.08, 67.96, 73.96)
silica.srp2 = c(72.28, 72.28, 67.47, 67.47, 62.64, 62.64, 57.07, 57.07, 53.03, 53.03, 48.03, 48.03)
silica.vul2 = c(53.33, 53.33, 59.66, 59.66, 64.08, 64.08, 67.96, 67.96, 73.96, 73.96)

# plot SiO2 vs CF

x11()
par(mfrow=c(2,2))
plot(srp.CF, c(rep(silica.srp, 3),silica.srp[-6]))
points(srp.CF[1:6], silica.srp, col='blue', pch=19)
points(srp.CF[7:12], silica.srp, col='green', pch=19)
points(srp.CF[13:18], silica.srp, col='gold', pch=19)
points(srp.CF[19:23], silica.srp[-6], col='red', pch=19)
plot(vul.CF, rep(silica.vul, 4))
points(vul.CF[1:5], silica.vul, col='blue', pch=19)
points(vul.CF[6:10], silica.vul, col='green', pch=19)
points(vul.CF[11:15], silica.vul, col='gold', pch=19)
points(vul.CF[16:20], silica.vul, col='red', pch=19)
plot(srp2.CF, silica.srp2)
points(srp2.CF[seq(1,11,2)], silica.srp2[seq(1,11,2)], col='blue', pch=19)
points(srp2.CF[seq(2,12,2)], silica.srp2[seq(2,12,2)], col='red', pch=19)
plot(vul2.CF, silica.vul2)
points(vul2.CF[seq(1,9,2)], silica.vul2[seq(1,9,2)], col='blue', pch=19)
points(vul2.CF[seq(2,10,2)], silica.vul2[seq(2,10,2)], col='red', pch=19)

# plot SiO2 vs TF

x11()
par(mfrow=c(2,2))
plot(srp.TF, c(rep(silica.srp, 3),silica.srp[-6]))
points(srp.TF[1:6], silica.srp, col='blue', pch=19)
points(srp.TF[7:12], silica.srp, col='green', pch=19)
points(srp.TF[13:18], silica.srp, col='gold', pch=19)
points(srp.TF[19:23], silica.srp[-6], col='red', pch=19)
plot(vul.TF, rep(silica.vul, 4))
points(vul.TF[1:5], silica.vul, col='blue', pch=19)
points(vul.TF[6:10], silica.vul, col='green', pch=19)
points(vul.TF[11:15], silica.vul, col='gold', pch=19)
points(vul.TF[16:20], silica.vul, col='red', pch=19)
plot(srp2.TF, silica.srp2)
points(srp2.TF[seq(1,11,2)], silica.srp2[seq(1,11,2)], col='blue', pch=19)
points(srp2.TF[seq(2,12,2)], silica.srp2[seq(2,12,2)], col='red', pch=19)
plot(vul2.TF, silica.vul2)
points(vul2.TF[seq(1,9,2)], silica.vul2[seq(1,9,2)], col='blue', pch=19)
points(vul2.TF[seq(2,10,2)], silica.vul2[seq(2,10,2)], col='red', pch=19)

# plot SiO2 vs TFval

x11()
par(mfrow=c(2,2))
plot(srp.TFval, c(rep(silica.srp, 3),silica.srp[-6]))
points(srp.TFval[1:6], silica.srp, col='blue', pch=19)
points(srp.TFval[7:12], silica.srp, col='green', pch=19)
points(srp.TFval[13:18], silica.srp, col='gold', pch=19)
points(srp.TFval[19:23], silica.srp[-6], col='red', pch=19)
plot(vul.TFval, rep(silica.vul, 4))
points(vul.TFval[1:5], silica.vul, col='blue', pch=19)
points(vul.TFval[6:10], silica.vul, col='green', pch=19)
points(vul.TFval[11:15], silica.vul, col='gold', pch=19)
points(vul.TFval[16:20], silica.vul, col='red', pch=19)
plot(srp2.TFval, silica.srp2)
points(srp2.TFval[seq(1,11,2)], silica.srp2[seq(1,11,2)], col='blue', pch=19)
points(srp2.TFval[seq(2,12,2)], silica.srp2[seq(2,12,2)], col='red', pch=19)
plot(vul2.TFval, silica.vul2)
points(vul2.TFval[seq(1,9,2)], silica.vul2[seq(1,9,2)], col='blue', pch=19)
points(vul2.TFval[seq(2,10,2)], silica.vul2[seq(2,10,2)], col='red', pch=19)

dummy.srp=numeric(23)
dummy.srp[19:23]=1
dummy.vul=numeric(20)
dummy.vul[16:20]=1
dummy.srp2=numeric(12)
dummy.srp2[seq(2,12,2)]=1
dummy.vul2=numeric(10)
dummy.vul2[seq(2,10,2)]=1


#LM con TF

srp.lm=lm(rep(silica.srp, 4)[-24] ~ dummy.srp + srp.CF + I(srp.CF^2) + srp.TF + 
            dummy.srp:I(srp.CF^2) + dummy.srp:srp.CF + dummy.srp:srp.TF + I(srp.TFval^2))  ### interazioni poco significative
summary(srp.lm)

x11()
par(mfrow=c(2,2))
plot(srp.lm)
shapiro.test(srp.lm$residuals)
x11()
par(mfrow=c(1,2))
plot(srp.CF, c(rep(silica.srp, 3),silica.srp[-6]))
points(srp.CF[1:6], silica.srp, col='blue', pch=19)
points(srp.CF[7:12], silica.srp, col='green', pch=19)
points(srp.CF[13:18], silica.srp, col='gold', pch=19)
points(srp.CF[19:23], silica.srp[-6], col='red', pch=19)
plot(srp.CF, srp.lm$fitted.values)
points(srp.CF[1:6], srp.lm$fitted.values[1:6], col='blue', pch=19)
points(srp.CF[7:12], srp.lm$fitted.values[7:12], col='green', pch=19)
points(srp.CF[13:18], srp.lm$fitted.values[13:18], col='gold', pch=19)
points(srp.CF[19:23], srp.lm$fitted.values[19:23], col='red', pch=19)

vul.lm=lm(rep(silica.vul,4) ~  I(vul.CF^2) + I(vul.CF^4) + vul.CF + vul.TFval)
summary(vul.lm)
x11()
par(mfrow=c(2,2))
plot(vul.lm)

srp2.lm=lm(silica.srp2 ~ srp2.CF )
summary(srp2.lm)
x11()
par(mfrow=c(2,2))
plot(srp2.lm)

vul2.lm=lm(silica.vul2 ~ vul2.CF + dummy.vul2 + dummy.vul2:vul2.CF)
summary(vul2.lm)
x11()
par(mfrow=c(2,2))
plot(vul2.lm)

#LM con valore finale

srp.lm=lm(rep(silica.srp, 4)[-24] ~ srp.CF + srp.last:dummy.srp + I(srp.CF^2) + srp.TF
          + srp.TF:dummy.srp) #-->B600 outlier
summary(srp.lm)
x11()
par(mfrow=c(2,2))
plot(srp.lm)
shapiro.test(srp.lm$residuals)

x11()
plot(srp.lm$fitted.values)
points(c(rep(silica.srp, 3),silica.srp[-6]))

vul.lm=lm(rep(silica.vul,4) ~ 0+vul.CF + vul.last + vul.TF)
summary(vul.lm)
x11()
par(mfrow=c(2,2))
plot(vul.lm)

srp2.lm=lm(silica.srp2 ~ 0+  I(srp2.CF^2) + I(srp2.last^2))
summary(srp2.lm)
x11()
par(mfrow=c(2,2))
plot(srp2.lm)


vul2.lm=lm(silica.vul2 ~ 0+ vul2.last + vul2.CF) #---> 1) Sfresh, 10)RShot outliers ??
summary(vul2.lm)
x11()
par(mfrow=c(2,2))
plot(vul2.lm)


#=???????? si comporta a caso senza 1 e 10
vul2.lm=lm(silica.vul2[-c(1,10)] ~ 0  + vul2.CF[-c(1,10)] + I(vul2.CF[-c(1,10)]^2)) #---> 1) Sfresh, 10)RShot outliers
summary(vul2.lm)
x11()
par(mfrow=c(2,2))
plot(vul2.lm)


# plot fitted values srp
y = c(rep(silica.srp, 3),silica.srp[-6])
fit_val=fitted.values(srp.lm)
labs=c('RB-150','B2-150','B4-150','B6-150','B8-150','B-150','RB-300','B2-300','B4-300','B6-300','B8-300','B-300',
       'RB-450','B2-450','B4-450','B6-450','B8-450','B-450','RB-600','B2-600','B4-600','B6-600','B8-600')
X11()
par(mfrow = c(2,3))
barplot(height=rbind(y[c(1,7,13,19)], fit_val[c(1,7,13,19)]), beside=TRUE, names=labs[c(1,7,13,19)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='RB Blend')
barplot(height=rbind(y[c(2,8,14,20)], fit_val[c(2,8,14,20)]), beside=TRUE, names=labs[c(2,8,14,20)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='B2 Blend')
barplot(height=rbind(y[c(3,9,15,21)], fit_val[c(3,9,15,21)]), beside=TRUE, names=labs[c(3,9,15,21)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='B4 Blend')
barplot(height=rbind(y[c(4,10,16,22)], fit_val[c(4,10,16,22)]), beside=TRUE, names=labs[c(4,10,16,22)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='B6 Blend')
barplot(height=rbind(y[c(5,11,17,23)], fit_val[c(5,11,17,23)]), beside=TRUE, names=labs[c(5,11,17,23)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='B8 Blend')
barplot(height=rbind(y[c(6,12,18)], fit_val[c(6,12,18)]), beside=TRUE, names=labs[c(6,12,18)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='B Blend')

X11()
par(mfrow = c(2,2))
barplot(height=rbind(y[c(1:6)], fit_val[c(1:6)]), beside=TRUE, names=labs[c(1:6)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 150°')
barplot(height=rbind(y[c(7:12)], fit_val[c(7:12)]), beside=TRUE, names=labs[c(7:12)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 300°')
barplot(height=rbind(y[c(13:18)], fit_val[c(13:18)]), beside=TRUE, names=labs[c(13:18)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 450°')
barplot(height=rbind(y[c(19:23)], fit_val[c(19:23)]), beside=TRUE, names=labs[c(19:23)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 600°')


X11()
plot(c(1:23), y[c(1:23)], pch=19, ylim=c(40,80), xlab='Blends', ylab='SiO2 %', xaxt='n', main='Real vs predicted percentage of silica')
text(c(1:23), y[c(1:23)],labels=labs[c(1:23)],pos=3)
points(c(1:23), fit_val[c(1:23)], pch=19, col='red')
for (i in (1:22)){
  segments(i,0,i,y[i],col='black')
  segments(i,fit_val[i],i+1,fit_val[i+1],col='red')
}
segments(23,0,23,y[23],col='black')
legend('bottomleft', legend=c('real values','fitted values'), bg='white', pch=19, col=c('black','red'))


errors.srp = y - fit_val
square.errors.srp=errors.srp^2
errors.srp.df=data.frame(labs,square.errors.srp)

x11()
barplot(errors.srp.df[,2], names.arg = errors.srp.df[,1])

# plot fitted values vul
y = rep(silica.vul,4)
fit_val=fitted.values(vul.lm)
labs=c('S-150','S7-150','S5-150','S3-150','Rs-150','S-300','S7-300','S5-300','S3-300','RS-300',
       'S-450','S7-450','S5-450','S3-450','RS-450','S-600','S7-600','S5-600','S3-600','RS-600')
X11()
par(mfrow = c(2,3))
barplot(height=rbind(y[c(1,6,11,16)], fit_val[c(1,6,11,16)]), beside=TRUE, names=labs[c(1,6,11,16)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='S Blend')
barplot(height=rbind(y[c(2,7,12,17)], fit_val[c(2,7,12,17)]), beside=TRUE, names=labs[c(2,7,12,17)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='S7 Blend')
barplot(height=rbind(y[c(3,8,13,18)], fit_val[c(3,8,13,18)]), beside=TRUE, names=labs[c(3,8,13,18)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='S5 Blend')
barplot(height=rbind(y[c(4,9,14,19)], fit_val[c(4,9,14,19)]), beside=TRUE, names=labs[c(4,9,14,19)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='S3 Blend')
barplot(height=rbind(y[c(5,10,15,20)], fit_val[c(5,10,15,20)]), beside=TRUE, names=labs[c(5,10,15,20)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='RS Blend')

X11()
par(mfrow = c(2,2))
barplot(height=rbind(y[c(1:5)], fit_val[c(1:5)]), beside=TRUE, names=labs[c(1:5)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 150°')
barplot(height=rbind(y[c(6:10)], fit_val[c(6:10)]), beside=TRUE, names=labs[c(6:10)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 300°')
barplot(height=rbind(y[c(11:15)], fit_val[c(11:15)]), beside=TRUE, names=labs[c(11:15)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 450°')
barplot(height=rbind(y[c(16:20)], fit_val[c(16:20)]), beside=TRUE, names=labs[c(16:20)], las=2, ylim=c(0,100), col=c('darkcyan','orange'), main='Blends @ 600°')


X11()
plot(c(1:20), y[c(1:20)], pch=19, ylim=c(40,80), xlab='Blends', ylab='SiO2 %', xaxt='n', main='Real vs predicted percentage of silica')
text(c(1:20), y[c(1:20)],labels=labs[c(1:20)],pos=3)
points(c(1:20), fit_val[c(1:20)], pch=19, col='red')
for (i in (1:19)){
  segments(i,0,i,y[i],col='black')
  segments(i,fit_val[i],i+1,fit_val[i+1],col='red')
}
segments(20,0,20,y[20],col='black')
legend('bottomleft', legend=c('real values','fitted values'), bg='white', pch=19, col=c('black','red'))


errors.srp = y - fit_val
square.errors.srp=errors.srp^2
errors.srp.df=data.frame(labs,square.errors.srp)

x11()
barplot(errors.srp.df[,2], names.arg = errors.srp.df[,1])