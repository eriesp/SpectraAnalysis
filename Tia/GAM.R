################################################################################

### Import the datasets ###

library(readxl)
library(car)
library(ISLR2)
library(mgcv)
library(rgl)
library(splines)
library(pbapply)


srp <- read_excel("Srp.xlsx")
head(srp)
View(srp)
summary(srp)

vul <- read_excel("Vulcano.xlsx")
head(vul)
View(vul)
summary(vul)

srp2 <- read_excel("SrpReflectance.xlsx")
head(srp2)
View(srp2)
summary(srp2)

vul2 <- read_excel("VulReflectance.xlsx")
head(vul2)
View(vul2)
summary(vul2)

################################################################################

### Feature extraction ###

srp.emis=srp[713:1, -c(1,25)] # remove B600, FROM OUTLIER DETECTION (FUNCTIONAL BOXPLOT) 
                              # WE CAN CONSIDER B600 AS A SPECIAL CASE, FOR THE MOMENT I REMOVE IT
srp.refl=srp2[, -1]
vul.emis=vul[713:1,-1]
vul.refl=vul2[,-1]
emis=cbind(srp.emis,vul.emis)
refl=cbind(srp.refl, vul.refl)
emis.x=srp[713:1,1]$`Wavelength (micron)`
refl.x=srp2[,1]$wl

x11()
par(mfrow=c(2,2))
matplot(emis.x, srp.emis, type='l')
matplot(emis.x, vul.emis, type='l')
matplot(refl.x, srp.refl, type='l')
matplot(refl.x, vul.refl, type='l')

x11()
par(mfrow=c(1,2))
matplot(emis.x, emis, type='l')
matplot(refl.x, refl, type='l')

emis.CF=numeric(43)
emis.CFval=numeric(43)
emis.TF=numeric(43)
emis.TFval=numeric(43)
emis.last=numeric(43)
for(i in 1:43) {
  emis.CF[i]=which.max(unlist(emis[213:518,i])) + 213
  emis.CFval[i]=emis[emis.CF[i],i]
  emis.CF[i]=emis.x[emis.CF[i]]
  emis.TF[i]=which.min(unlist(emis[413:523,i])) + 413
  emis.TFval[i]=emis[emis.TF[i],i]
  emis.TF[i]=emis.x[emis.TF[i]]
  emis.last[i]=emis[713,i]
}

srp.emis.CF=numeric(23)
srp.emis.CFval=numeric(23)
srp.emis.TF=numeric(23)
srp.emis.TFval=numeric(23)
srp.emis.last=numeric(23)
for(i in 1:23) {
  srp.emis.CF[i]=which.max(unlist(srp.emis[213:518,i])) + 213
  srp.emis.CFval[i]=srp.emis[srp.emis.CF[i],i]
  srp.emis.CF[i]=emis.x[srp.emis.CF[i]]
  srp.emis.TF[i]=which.min(unlist(srp.emis[413:523,i])) + 413
  srp.emis.TFval[i]=srp.emis[srp.emis.TF[i],i]
  srp.emis.TF[i]=emis.x[srp.emis.TF[i]]
  srp.emis.last[i]=srp.emis[713,i]
}

vul.emis.CF=numeric(20)
vul.emis.CFval=numeric(20)
vul.emis.TF=numeric(20)
vul.emis.TFval=numeric(20)
vul.emis.last=numeric(20)
for(i in 1:20) {
  vul.emis.CF[i]=which.max(unlist(vul.emis[213:518,i])) + 213
  vul.emis.CFval[i]=vul.emis[vul.emis.CF[i],i]
  vul.emis.CF[i]=emis.x[vul.emis.CF[i]]
  vul.emis.TF[i]=which.min(unlist(vul.emis[413:523,i])) + 413
  vul.emis.TFval[i]=vul.emis[vul.emis.TF[i],i]
  vul.emis.TF[i]=emis.x[vul.emis.TF[i]]
  vul.emis.last[i]=vul.emis[713,i]
}

refl.CF=numeric(22)
refl.CFval=numeric(22)
refl.TF=numeric(22)
refl.TFval=numeric(22)
refl.last=numeric(22)
for(i in 1:22) {
  refl.CF[i]=which.max(unlist(refl[213:518,i])) + 213
  refl.CFval[i]=refl[refl.CF[i],i]
  refl.CF[i]=refl.x[refl.CF[i]]
  refl.TF[i]=which.min(unlist(refl[413:550,i])) + 413
  refl.TFval[i]=refl[refl.TF[i],i]
  refl.TF[i]=refl.x[refl.TF[i]]
  refl.last[i]=refl[714,i]
}

srp.refl.CF=numeric(12)
srp.refl.CFval=numeric(12)
srp.refl.TF=numeric(12)
srp.refl.TFval=numeric(12)
srp.refl.last=numeric(12)
for(i in 1:12) {
  srp.refl.CF[i]=which.max(unlist(srp.refl[213:518,i])) + 213
  srp.refl.CFval[i]=srp.refl[srp.refl.CF[i],i]
  srp.refl.CF[i]=refl.x[srp.refl.CF[i]]
  srp.refl.TF[i]=which.min(unlist(srp.refl[413:550,i])) + 413
  srp.refl.TFval[i]=srp.refl[srp.refl.TF[i],i]
  srp.refl.TF[i]=refl.x[srp.refl.TF[i]]
  srp.refl.last[i]=srp.refl[714,i]
}

vul.refl.CF=numeric(10)
vul.refl.CFval=numeric(10)
vul.refl.TF=numeric(10)
vul.refl.TFval=numeric(10)
vul.refl.last=numeric(10)
for(i in 1:10) {
  vul.refl.CF[i]=which.max(unlist(vul.refl[213:518,i])) + 213
  vul.refl.CFval[i]=vul.refl[vul.refl.CF[i],i]
  vul.refl.CF[i]=refl.x[vul.refl.CF[i]]
  vul.refl.TF[i]=which.min(unlist(vul.refl[413:550,i])) + 413
  vul.refl.TFval[i]=vul.refl[vul.refl.TF[i],i]
  vul.refl.TF[i]=refl.x[vul.refl.TF[i]]
  vul.refl.last[i]=vul.refl[714,i]
}

temp.emis.srp = c(rep(150,6), rep(300,6), rep(450,6), rep(600,5))
temp.emis.vul = c(rep(150,5), rep(300,5), rep(450,5), rep(600,5))
temp.emis = c(temp.emis.srp, temp.emis.vul)

temp.refl.srp = rep(c(20,500), 6)
temp.refl.vul = rep(c(20,500), 5)
temp.refl = c(temp.refl.srp, temp.refl.vul)

dummy.srp.emis=numeric(23)
dummy.srp.emis[19:23]=1
dummy.vul.emis=numeric(20)
dummy.vul.emis[16:20]=1
dummy.srp.refl=numeric(12)
dummy.srp.refl[seq(2,12,2)]=1
dummy.vul.refl=numeric(10)
dummy.vul.refl[seq(2,10,2)]=1
dummy.emis=c(dummy.srp.emis, dummy.vul.emis)
dummy.refl=c(dummy.srp.refl, dummy.vul.refl)

silica.srp = c(72.28, 67.47, 62.64, 57.07, 53.03, 48.03)
silica.vul = c(53.33, 59.66, 64.08, 67.96, 73.96)
silica.srp.emis = rep(silica.srp,4)[-24]
silica.vul.emis = rep(silica.vul,4)
silica.srp.refl = c(72.28, 72.28, 67.47, 67.47, 62.64, 62.64, 57.07, 57.07, 53.03, 53.03, 48.03, 48.03)
silica.vul.refl = c(53.33, 53.33, 59.66, 59.66, 64.08, 64.08, 67.96, 67.96, 73.96, 73.96)
silica.emis=c(silica.srp.emis, silica.vul.emis)
silica.refl=c(silica.srp.refl, silica.vul.refl)

alcali.srp = c(7.51, 6.62, 5.72, 4.78, 3.72, 2.71)
alcali.vul = c(8.12, 8.40, 8.47, 8.45, 8.75)
alcali.srp.emis = rep(alcali.srp,4)[-24]
alcali.vul.emis = rep(alcali.vul,4)
alcali.srp.refl = c(7.51, 7.51, 6.62, 6.62, 5.72, 5.72, 4.78, 4.78, 3.72, 3.72, 2.71, 2.71)
alcali.vul.refl = c(8.12, 8.12, 8.40, 8.40, 8.47, 8.47, 8.45, 8.45, 8.75, 8.75)
alcali.emis=c(alcali.srp.emis, alcali.vul.emis)
alcali.refl=c(alcali.srp.refl, alcali.vul.refl)

################################################################################

### Some plots of our features ###

x11()
plot(silica.emis, alcali.emis)
points(silica.srp, alcali.srp, col='red', pch=19)
points(silica.vul, alcali.vul, col='blue', pch=19)


x11()
par(mfrow=c(2,3))
plot(emis.CF, silica.emis)
points(srp.emis.CF, silica.srp.emis, col='red', pch=19)
points(vul.emis.CF, silica.vul.emis, col='blue', pch=19)
plot(emis.CFval, silica.emis)
points(srp.emis.CFval, silica.srp.emis, col='red', pch=19)
points(vul.emis.CFval, silica.vul.emis, col='blue', pch=19)
plot(emis.TF, silica.emis)
points(srp.emis.TF, silica.srp.emis, col='red', pch=19)
points(vul.emis.TF, silica.vul.emis, col='blue', pch=19)
plot(emis.TFval, silica.emis)
points(srp.emis.TFval, silica.srp.emis, col='red', pch=19)
points(vul.emis.TFval, silica.vul.emis, col='blue', pch=19)
plot(emis.last, silica.emis)
points(srp.emis.last, silica.srp.emis, col='red', pch=19)
points(vul.emis.last, silica.vul.emis, col='blue', pch=19)
plot(dummy.emis, silica.emis)
points(dummy.srp.emis, silica.srp.emis, col='red', pch=19)
points(dummy.vul.emis, silica.vul.emis, col='blue', pch=19)

scatterplotMatrix(data.frame(silica.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.last,temp.emis))


x11()
par(mfrow=c(2,3))
plot(refl.CF, silica.refl)
points(srp.refl.CF, silica.srp.refl, col='red', pch=19)
points(vul.refl.CF, silica.vul.refl, col='blue', pch=19)
plot(refl.CFval, silica.refl)
points(srp.refl.CFval, silica.srp.refl, col='red', pch=19)
points(vul.refl.CFval, silica.vul.refl, col='blue', pch=19)
plot(refl.TF, silica.refl)
points(srp.refl.TF, silica.srp.refl, col='red', pch=19)
points(vul.refl.TF, silica.vul.refl, col='blue', pch=19)
plot(refl.TFval, silica.refl)
points(srp.refl.TFval, silica.srp.refl, col='red', pch=19)
points(vul.refl.TFval, silica.vul.refl, col='blue', pch=19)
plot(refl.last, silica.refl)
points(srp.refl.last, silica.srp.refl, col='red', pch=19)
points(vul.refl.last, silica.vul.refl, col='blue', pch=19)
plot(dummy.refl, silica.refl)
points(dummy.srp.refl, silica.srp.refl, col='red', pch=19)
points(dummy.vul.refl, silica.vul.refl, col='blue', pch=19)

scatterplotMatrix(data.frame(silica.refl,refl.CF,refl.CFval,refl.TF,refl.TFval,refl.last,temp.refl))


x11()
par(mfrow=c(2,3))
plot(emis.CF, alcali.emis)
points(srp.emis.CF, alcali.srp.emis, col='red', pch=19)
points(vul.emis.CF, alcali.vul.emis, col='blue', pch=19)
plot(emis.CFval, alcali.emis)
points(srp.emis.CFval, alcali.srp.emis, col='red', pch=19)
points(vul.emis.CFval, alcali.vul.emis, col='blue', pch=19)
plot(emis.TF, alcali.emis)
points(srp.emis.TF, alcali.srp.emis, col='red', pch=19)
points(vul.emis.TF, alcali.vul.emis, col='blue', pch=19)
plot(emis.TFval, alcali.emis)
points(srp.emis.TFval, alcali.srp.emis, col='red', pch=19)
points(vul.emis.TFval, alcali.vul.emis, col='blue', pch=19)
plot(emis.last, alcali.emis)
points(srp.emis.last, alcali.srp.emis, col='red', pch=19)
points(vul.emis.last, alcali.vul.emis, col='blue', pch=19)
plot(dummy.emis, alcali.emis)
points(dummy.srp.emis, alcali.srp.emis, col='red', pch=19)
points(dummy.vul.emis, alcali.vul.emis, col='blue', pch=19)

scatterplotMatrix(data.frame(alcali.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.last,temp.emis))


x11()
par(mfrow=c(2,3))
plot(refl.CF, alcali.refl)
points(srp.refl.CF, alcali.srp.refl, col='red', pch=19)
points(vul.refl.CF, alcali.vul.refl, col='blue', pch=19)
plot(refl.CFval, alcali.refl)
points(srp.refl.CFval, alcali.srp.refl, col='red', pch=19)
points(vul.refl.CFval, alcali.vul.refl, col='blue', pch=19)
plot(refl.TF, alcali.refl)
points(srp.refl.TF, alcali.srp.refl, col='red', pch=19)
points(vul.refl.TF, alcali.vul.refl, col='blue', pch=19)
plot(refl.TFval, alcali.refl)
points(srp.refl.TFval, alcali.srp.refl, col='red', pch=19)
points(vul.refl.TFval, alcali.vul.refl, col='blue', pch=19)
plot(refl.last, alcali.refl)
points(srp.refl.last, alcali.srp.refl, col='red', pch=19)
points(vul.refl.last, alcali.vul.refl, col='blue', pch=19)
plot(dummy.refl, alcali.refl)
points(dummy.srp.refl, alcali.srp.refl, col='red', pch=19)
points(dummy.vul.refl, alcali.vul.refl, col='blue', pch=19)

scatterplotMatrix(data.frame(alcali.refl,refl.CF,refl.CFval,refl.TF,refl.TFval,refl.last,temp.refl))

################################################################################

### GAM for emissivity explaining silica ### SATISFYING RESULTS

### reduced model
silica.emis.gam.reduced = gam(silica.emis ~ emis.CF + s(emis.CFval, bs='cr') + s(emis.TF, bs='cr') + s(emis.TFval, bs='cr'))
summary(silica.emis.gam.reduced)

hist(silica.emis.gam.reduced$residuals)
qqnorm(silica.emis.gam.reduced$residuals)
shapiro.test(silica.emis.gam.reduced$residuals)

silica.emis.gam = gam(silica.emis ~ s(emis.CF, bs='cr') + s(emis.CFval, bs='cr') + s(emis.TF, bs='cr') + s(emis.TFval, bs='cr'))
summary(silica.emis.gam)

hist(silica.emis.gam$residuals)
qqnorm(silica.emis.gam$residuals)
shapiro.test(silica.emis.gam$residuals)

plot(silica.emis.gam.reduced)
plot(silica.emis.gam)

anova(silica.emis.gam.reduced, silica.emis.gam, test = "F") ### no need to smooth also the CF --> we stick with reduced model

x11()
plot(silica.emis, silica.emis.gam.reduced$fitted.values)
abline(a=0, b=1)
plot(silica.emis.gam.reduced$residuals)
abline(a=0, b=0)
abline(v=23.5)

### GAM for emissivity explaining alcali ### NOT GOOD ENOUGH 
 
alcali.emis.gam = gam(alcali.emis ~ s(emis.CF, bs='cr') + s(emis.CFval, bs='cr') + s(emis.TF, bs='cr') + s(emis.TFval, bs='cr'))
summary(alcali.emis.gam)

hist(alcali.emis.gam$residuals)
qqnorm(alcali.emis.gam$residuals)
shapiro.test(alcali.emis.gam$residuals)

plot(alcali.emis.gam)

x11()
plot(alcali.emis, alcali.emis.gam$fitted.values)
abline(a=0, b=1)
plot(alcali.emis.gam$residuals)
abline(a=0, b=0)
abline(v=23.5)



### GAM for reflectance explaining silica ### PRETTY GOOD

### reduced model
silica.refl.gam.reduced = gam(silica.refl ~ refl.CF + s(refl.TFval, bs='cr'))
summary(silica.refl.gam.reduced)

hist(silica.refl.gam.reduced$residuals)
qqnorm(silica.refl.gam.reduced$residuals)
shapiro.test(silica.refl.gam.reduced$residuals)

silica.refl.gam = gam(silica.refl ~ s(refl.CF, bs='cr') + s(refl.TFval, bs='cr'))
summary(silica.refl.gam)

hist(silica.refl.gam$residuals)
qqnorm(silica.refl.gam$residuals)
shapiro.test(silica.refl.gam$residuals)

plot(silica.refl.gam.reduced)
plot(silica.refl.gam)

anova(silica.refl.gam.reduced, silica.refl.gam, test = "F") ### smoothing also the CF seems useful

x11()
plot(silica.refl, silica.refl.gam$fitted.values)
abline(a=0, b=1)
plot(silica.refl.gam$residuals)
abline(a=0, b=0)
abline(v=12.5)

### GAM for reflectance explaining alcali ### THIS ONE LOOKS GOOD ASWELL

alcali.refl.gam = gam(alcali.refl ~ s(refl.CF, bs='cr') + s(refl.TFval, bs='cr'))
summary(alcali.refl.gam)

hist(alcali.refl.gam$residuals)
qqnorm(alcali.refl.gam$residuals)
shapiro.test(alcali.refl.gam$residuals)

plot(alcali.refl.gam)

x11()
plot(alcali.refl, alcali.refl.gam$fitted.values)
abline(a=0, b=1)
plot(alcali.refl.gam$residuals)
abline(a=0, b=0)
abline(v=12.5)

################################################################################

### TAS classification ###

### real values vs estimated values emissivity
x11()
plot(silica.emis, alcali.emis, col='green', pch=19)
points(silica.emis.gam.reduced$fitted.values, alcali.emis.gam$fitted.values, col='red', pch=4)

### real values vs estimated values reflectance
x11()
plot(silica.refl, alcali.refl, col='green', pch=19)
points(silica.refl.gam$fitted.values, alcali.refl.gam$fitted.values, col='red', pch=4)

################################################################################

### More diagnostic on the model: effect of temperature on features ###

i150.emis=c(seq(1,6), seq(24,28))
i300.emis=c(seq(7,12), seq(29,33))
i450.emis=c(seq(13,18), seq(34,38))
i600.emis=c(seq(19,23), seq(39,43))

i20.refl=seq(1,22,2)
i500.refl=seq(2,22,2)

### Christiansen Feature

x11()
par(mfrow=c(1,2))
plot(temp.emis, emis.CF)
plot(temp.refl, refl.CF)

box.CF=boxplot(emis.CF[i150.emis], emis.CF[i300.emis], emis.CF[i450.emis], emis.CF[i600.emis]) ### clear differences
boxplot(refl.CF[i20.refl], refl.CF[i500.refl]) ### same


### TEST 1: JUST USE MEAN
m150.CF=mean(emis.CF[i150.emis])
m300.CF=mean(emis.CF[i300.emis])
m450.CF=mean(emis.CF[i450.emis])
m600.CF=mean(emis.CF[i600.emis])
T01.CF=abs(m150.CF-m300.CF)
T02.CF=abs(m150.CF-m450.CF)
T03.CF=abs(m150.CF-m600.CF)

set.seed(54)

B=10000
T1.CF=numeric(B)
T2.CF=numeric(B)
T3.CF=numeric(B)
n=length(emis.CF)

for(i in 1:B){
  perm = sample(1:n)
  CF.perm = emis.CF[perm]
  m150.CF.perm = mean(CF.perm[i150.emis])
  m300.CF.perm = mean(CF.perm[i300.emis])
  m450.CF.perm = mean(CF.perm[i450.emis])
  m600.CF.perm = mean(CF.perm[i600.emis])
  
  # Test statistic:
  T1.CF[i] = abs(m150.CF.perm-m300.CF.perm)
  T2.CF[i] = abs(m150.CF.perm-m450.CF.perm)
  T3.CF[i] = abs(m150.CF.perm-m600.CF.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.CF, breaks = 30)
abline(v=T01.CF,col=3,lwd=2)
hist(T2.CF, breaks=30)
abline(v=T02.CF,col=3,lwd=2)
hist(T3.CF, breaks=30)
abline(v=T03.CF,col=3,lwd=2)


# p-value
p1 <- sum(T1.CF>=T01.CF)/B
p1
p2 <- sum(T2.CF>=T02.CF)/B
p2
p3 <- sum(T3.CF>=T03.CF)/B
p3                         ### we have evidence that 600 degrees have a significant impact on the Christiansen Feature

### TEST 2 --> EXPLOIT BOXPLOT
T01.CF=as.double((box.CF$conf[,1] - box.CF$conf[,2])%*%(box.CF$conf[,1] - box.CF$conf[,2]))
T02.CF=as.double((box.CF$conf[,1] - box.CF$conf[,3])%*%(box.CF$conf[,1] - box.CF$conf[,3]))
T03.CF=as.double((box.CF$conf[,1] - box.CF$conf[,4])%*%(box.CF$conf[,1] - box.CF$conf[,4]))

set.seed(54)

B=10000
T1.CF=numeric(B)
T2.CF=numeric(B)
T3.CF=numeric(B)
n=length(emis.CF)

for(i in 1:B){
  perm = sample(1:n)
  CF.perm = emis.CF[perm]
  box.CF.perm=boxplot(CF.perm[i150.emis], CF.perm[i300.emis], CF.perm[i450.emis], CF.perm[i600.emis], plot=FALSE)
  
  # Test statistic:
  T1.CF[i]=(box.CF.perm$conf[,1] - box.CF.perm$conf[,2])%*%(box.CF.perm$conf[,1] - box.CF.perm$conf[,2])
  T2.CF[i]=(box.CF.perm$conf[,1] - box.CF.perm$conf[,3])%*%(box.CF.perm$conf[,1] - box.CF.perm$conf[,3])
  T3.CF[i]=(box.CF.perm$conf[,1] - box.CF.perm$conf[,4])%*%(box.CF.perm$conf[,1] - box.CF.perm$conf[,4])
}

x11()
par(mfrow=c(1,3))
hist(T1.CF, breaks = 30)
abline(v=T01.CF,col=3,lwd=2)
hist(T2.CF, breaks=30)
abline(v=T02.CF,col=3,lwd=2)
hist(T3.CF, breaks=30)
abline(v=T03.CF,col=3,lwd=2)


# p-value
p1 <- sum(T1.CF>=T01.CF)/B
p1
p2 <- sum(T2.CF>=T02.CF)/B
p2
p3 <- sum(T3.CF>=T03.CF)/B
p3                            ### not anymore evidence

### TEST 3: USE MEDIAN
med150.CF=median(emis.CF[i150.emis])
med300.CF=median(emis.CF[i300.emis])
med450.CF=median(emis.CF[i450.emis])
med600.CF=median(emis.CF[i600.emis])
T01.CF=abs(med150.CF-med300.CF)
T02.CF=abs(med150.CF-med450.CF)
T03.CF=abs(med150.CF-med600.CF)

set.seed(54)

B=10000
T1.CF=numeric(B)
T2.CF=numeric(B)
T3.CF=numeric(B)
n=length(emis.CF)

for(i in 1:B){
  perm = sample(1:n)
  CF.perm = emis.CF[perm]
  med150.CF.perm = median(CF.perm[i150.emis])
  med300.CF.perm = median(CF.perm[i300.emis])
  med450.CF.perm = median(CF.perm[i450.emis])
  med600.CF.perm = median(CF.perm[i600.emis])
  
  # Test statistic:
  T1.CF[i] <- abs(med150.CF.perm-med300.CF.perm)
  T2.CF[i] <- abs(med150.CF.perm-med450.CF.perm)
  T3.CF[i] <- abs(med150.CF.perm-med600.CF.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.CF, breaks = 30)
abline(v=T01.CF,col=3,lwd=2)
hist(T2.CF, breaks=30)
abline(v=T02.CF,col=3,lwd=2)
hist(T3.CF, breaks=30)
abline(v=T03.CF,col=3,lwd=2)


# p-value
p1 <- sum(T1.CF>=T01.CF)/B
p1
p2 <- sum(T2.CF>=T02.CF)/B
p2
p3 <- sum(T3.CF>=T03.CF)/B
p3                         ### again no strong evidence

### Christiansen Feature value

x11()
par(mfrow=c(1,2))
plot(temp.emis, emis.CFval)
plot(temp.refl, refl.CFval)

box.CFval=boxplot(emis.CFval[i150.emis], emis.CFval[i300.emis], emis.CFval[i450.emis], emis.CFval[i600.emis]) ### differences
boxplot(refl.CFval[i20.refl], refl.CFval[i500.refl]) ### same

### TEST 1: JUST USE MEAN
m150.CFval=mean(emis.CFval[i150.emis])
m300.CFval=mean(emis.CFval[i300.emis])
m450.CFval=mean(emis.CFval[i450.emis])
m600.CFval=mean(emis.CFval[i600.emis])
T01.CFval=abs(m150.CFval-m300.CFval)
T02.CFval=abs(m150.CFval-m450.CFval)
T03.CFval=abs(m150.CFval-m600.CFval)

set.seed(54)

B=10000
T1.CFval=numeric(B)
T2.CFval=numeric(B)
T3.CFval=numeric(B)
n=length(emis.CFval)

for(i in 1:B){
  perm = sample(1:n)
  CFval.perm = emis.CFval[perm]
  m150.CFval.perm = mean(CFval.perm[i150.emis])
  m300.CFval.perm = mean(CFval.perm[i300.emis])
  m450.CFval.perm = mean(CFval.perm[i450.emis])
  m600.CFval.perm = mean(CFval.perm[i600.emis])
  
  # Test statistic:
  T1.CFval[i] = abs(m150.CFval.perm-m300.CFval.perm)
  T2.CFval[i] = abs(m150.CFval.perm-m450.CFval.perm)
  T3.CFval[i] = abs(m150.CFval.perm-m600.CFval.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.CFval, breaks = 30)
abline(v=T01.CFval,col=3,lwd=2)
hist(T2.CFval, breaks=30)
abline(v=T02.CFval,col=3,lwd=2)
hist(T3.CFval, breaks=30)
abline(v=T03.CFval,col=3,lwd=2)


# p-value
p1 <- sum(T1.CFval>=T01.CFval)/B
p1
p2 <- sum(T2.CFval>=T02.CFval)/B
p2
p3 <- sum(T3.CFval>=T03.CFval)/B
p3                                  ### nothing

### TEST 2 --> EXPLOIT BOXPLOT
T01.CFval=as.double((box.CFval$conf[,1] - box.CFval$conf[,2])%*%(box.CFval$conf[,1] - box.CFval$conf[,2]))
T02.CFval=as.double((box.CFval$conf[,1] - box.CFval$conf[,3])%*%(box.CFval$conf[,1] - box.CFval$conf[,3]))
T03.CFval=as.double((box.CFval$conf[,1] - box.CFval$conf[,4])%*%(box.CFval$conf[,1] - box.CFval$conf[,4]))

set.seed(54)

B=10000
T1.CFval=numeric(B)
T2.CFval=numeric(B)
T3.CFval=numeric(B)
n=length(emis.CFval)

for(i in 1:B){
  perm = sample(1:n)
  CFval.perm = emis.CFval[perm]
  box.CFval.perm=boxplot(CFval.perm[i150.emis], CFval.perm[i300.emis], CFval.perm[i450.emis], CFval.perm[i600.emis], plot=FALSE)
  
  # Test statistic:
  T1.CFval[i]=(box.CFval.perm$conf[,1] - box.CFval.perm$conf[,2])%*%(box.CFval.perm$conf[,1] - box.CFval.perm$conf[,2])
  T2.CFval[i]=(box.CFval.perm$conf[,1] - box.CFval.perm$conf[,3])%*%(box.CFval.perm$conf[,1] - box.CFval.perm$conf[,3])
  T3.CFval[i]=(box.CFval.perm$conf[,1] - box.CFval.perm$conf[,4])%*%(box.CFval.perm$conf[,1] - box.CFval.perm$conf[,4])
}

x11()
par(mfrow=c(1,3))
hist(T1.CFval, breaks = 30)
abline(v=T01.CFval,col=3,lwd=2)
hist(T2.CFval, breaks=30)
abline(v=T02.CFval,col=3,lwd=2)
hist(T3.CFval, breaks=30)
abline(v=T03.CFval,col=3,lwd=2)


# p-value
p1 <- sum(T1.CFval>=T01.CFval)/B
p1
p2 <- sum(T2.CFval>=T02.CFval)/B
p2
p3 <- sum(T3.CFval>=T03.CFval)/B
p3                                  ### nothing

### TEST 3: USE MEDIAN
med150.CFval=median(emis.CFval[i150.emis])
med300.CFval=median(emis.CFval[i300.emis])
med450.CFval=median(emis.CFval[i450.emis])
med600.CFval=median(emis.CFval[i600.emis])
T01.CFval=abs(med150.CFval-med300.CFval)
T02.CFval=abs(med150.CFval-med450.CFval)
T03.CFval=abs(med150.CFval-med600.CFval)

set.seed(54)

B=10000
T1.CFval=numeric(B)
T2.CFval=numeric(B)
T3.CFval=numeric(B)
n=length(emis.CFval)

for(i in 1:B){
  perm = sample(1:n)
  CFval.perm = emis.CFval[perm]
  med150.CFval.perm = median(CFval.perm[i150.emis])
  med300.CFval.perm = median(CFval.perm[i300.emis])
  med450.CFval.perm = median(CFval.perm[i450.emis])
  med600.CFval.perm = median(CFval.perm[i600.emis])
  
  # Test statistic:
  T1.CFval[i] <- abs(med150.CFval.perm-med300.CFval.perm)
  T2.CFval[i] <- abs(med150.CFval.perm-med450.CFval.perm)
  T3.CFval[i] <- abs(med150.CFval.perm-med600.CFval.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.CFval, breaks = 30)
abline(v=T01.CFval,col=3,lwd=2)
hist(T2.CFval, breaks=30)
abline(v=T02.CFval,col=3,lwd=2)
hist(T3.CFval, breaks=30)
abline(v=T03.CFval,col=3,lwd=2)


# p-value
p1 <- sum(T1.CFval>=T01.CFval)/B
p1
p2 <- sum(T2.CFval>=T02.CFval)/B
p2
p3 <- sum(T3.CFval>=T03.CFval)/B
p3                                 ### nothing

### Transparency Feature

x11()
par(mfrow=c(1,2))
plot(temp.emis, emis.TF)
plot(temp.refl, refl.TF)

box.TF=boxplot(emis.TF[i150.emis], emis.TF[i300.emis], emis.TF[i450.emis], emis.TF[i600.emis]) ### apparently same
boxplot(refl.TF[i20.refl], refl.TF[i500.refl]) ### same

### TEST 1: JUST USE MEAN
m150.TF=mean(emis.TF[i150.emis])
m300.TF=mean(emis.TF[i300.emis])
m450.TF=mean(emis.TF[i450.emis])
m600.TF=mean(emis.TF[i600.emis])
T01.TF=abs(m150.TF-m300.TF)
T02.TF=abs(m150.TF-m450.TF)
T03.TF=abs(m150.TF-m600.TF)

set.seed(54)

B=10000
T1.TF=numeric(B)
T2.TF=numeric(B)
T3.TF=numeric(B)
n=length(emis.TF)

for(i in 1:B){
  perm = sample(1:n)
  TF.perm = emis.TF[perm]
  m150.TF.perm = mean(TF.perm[i150.emis])
  m300.TF.perm = mean(TF.perm[i300.emis])
  m450.TF.perm = mean(TF.perm[i450.emis])
  m600.TF.perm = mean(TF.perm[i600.emis])
  
  # Test statistic:
  T1.TF[i] = abs(m150.TF.perm-m300.TF.perm)
  T2.TF[i] = abs(m150.TF.perm-m450.TF.perm)
  T3.TF[i] = abs(m150.TF.perm-m600.TF.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.TF, breaks = 30)
abline(v=T01.TF,col=3,lwd=2)
hist(T2.TF, breaks=30)
abline(v=T02.TF,col=3,lwd=2)
hist(T3.TF, breaks=30)
abline(v=T03.TF,col=3,lwd=2)


# p-value
p1 <- sum(T1.TF>=T01.TF)/B
p1
p2 <- sum(T2.TF>=T02.TF)/B
p2
p3 <- sum(T3.TF>=T03.TF)/B
p3                         ### nothing

### TEST 2 --> EXPLOIT BOXPLOT
T01.TF=as.double((box.TF$conf[,1] - box.TF$conf[,2])%*%(box.TF$conf[,1] - box.TF$conf[,2]))
T02.TF=as.double((box.TF$conf[,1] - box.TF$conf[,3])%*%(box.TF$conf[,1] - box.TF$conf[,3]))
T03.TF=as.double((box.TF$conf[,1] - box.TF$conf[,4])%*%(box.TF$conf[,1] - box.TF$conf[,4]))

set.seed(54)

B=10000
T1.TF=numeric(B)
T2.TF=numeric(B)
T3.TF=numeric(B)
n=length(emis.TF)

for(i in 1:B){
  perm = sample(1:n)
  TF.perm = emis.TF[perm]
  box.TF.perm=boxplot(TF.perm[i150.emis], TF.perm[i300.emis], TF.perm[i450.emis], TF.perm[i600.emis], plot=FALSE)
  
  # Test statistic:
  T1.TF[i]=(box.TF.perm$conf[,1] - box.TF.perm$conf[,2])%*%(box.TF.perm$conf[,1] - box.TF.perm$conf[,2])
  T2.TF[i]=(box.TF.perm$conf[,1] - box.TF.perm$conf[,3])%*%(box.TF.perm$conf[,1] - box.TF.perm$conf[,3])
  T3.TF[i]=(box.TF.perm$conf[,1] - box.TF.perm$conf[,4])%*%(box.TF.perm$conf[,1] - box.TF.perm$conf[,4])
}

x11()
par(mfrow=c(1,3))
hist(T1.TF, breaks = 30)
abline(v=T01.TF,col=3,lwd=2)
hist(T2.TF, breaks=30)
abline(v=T02.TF,col=3,lwd=2)
hist(T3.TF, breaks=30)
abline(v=T03.TF,col=3,lwd=2)


# p-value
p1 <- sum(T1.TF>=T01.TF)/B
p1
p2 <- sum(T2.TF>=T02.TF)/B
p2
p3 <- sum(T3.TF>=T03.TF)/B
p3                            ### nothing
### TEST 3: USE MEDIAN
med150.TF=median(emis.TF[i150.emis])
med300.TF=median(emis.TF[i300.emis])
med450.TF=median(emis.TF[i450.emis])
med600.TF=median(emis.TF[i600.emis])
T01.TF=abs(med150.TF-med300.TF)
T02.TF=abs(med150.TF-med450.TF)
T03.TF=abs(med150.TF-med600.TF)

set.seed(54)

B=10000
T1.TF=numeric(B)
T2.TF=numeric(B)
T3.TF=numeric(B)
n=length(emis.TF)

for(i in 1:B){
  perm = sample(1:n)
  TF.perm = emis.TF[perm]
  med150.TF.perm = median(TF.perm[i150.emis])
  med300.TF.perm = median(TF.perm[i300.emis])
  med450.TF.perm = median(TF.perm[i450.emis])
  med600.TF.perm = median(TF.perm[i600.emis])
  
  # Test statistic:
  T1.TF[i] <- abs(med150.TF.perm-med300.TF.perm)
  T2.TF[i] <- abs(med150.TF.perm-med450.TF.perm)
  T3.TF[i] <- abs(med150.TF.perm-med600.TF.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.TF, breaks = 30)
abline(v=T01.TF,col=3,lwd=2)
hist(T2.TF, breaks=30)
abline(v=T02.TF,col=3,lwd=2)
hist(T3.TF, breaks=30)
abline(v=T03.TF,col=3,lwd=2)


# p-value
p1 <- sum(T1.TF>=T01.TF)/B
p1
p2 <- sum(T2.TF>=T02.TF)/B
p2
p3 <- sum(T3.TF>=T03.TF)/B
p3                         ### nothing


### Transparency Feature value

x11()
par(mfrow=c(1,2))
plot(temp.emis, emis.TFval)
plot(temp.refl, refl.TFval)

box.TFval=boxplot(emis.TFval[i150.emis], emis.TFval[i300.emis], emis.TFval[i450.emis], emis.TFval[i600.emis]) ### very different
box2.TFval=boxplot(refl.TFval[i20.refl], refl.TFval[i500.refl]) ### same

### TEST 1: JUST USE MEAN
m150.TFval=mean(emis.TFval[i150.emis])
m300.TFval=mean(emis.TFval[i300.emis])
m450.TFval=mean(emis.TFval[i450.emis])
m600.TFval=mean(emis.TFval[i600.emis])
T01.TFval=abs(m150.TFval-m300.TFval)
T02.TFval=abs(m150.TFval-m450.TFval)
T03.TFval=abs(m150.TFval-m600.TFval)

set.seed(54)

B=10000
T1.TFval=numeric(B)
T2.TFval=numeric(B)
T3.TFval=numeric(B)
n=length(emis.TFval)

for(i in 1:B){
  perm = sample(1:n)
  TFval.perm = emis.TFval[perm]
  m150.TFval.perm = mean(TFval.perm[i150.emis])
  m300.TFval.perm = mean(TFval.perm[i300.emis])
  m450.TFval.perm = mean(TFval.perm[i450.emis])
  m600.TFval.perm = mean(TFval.perm[i600.emis])
  
  # Test statistic:
  T1.TFval[i] = abs(m150.TFval.perm-m300.TFval.perm)
  T2.TFval[i] = abs(m150.TFval.perm-m450.TFval.perm)
  T3.TFval[i] = abs(m150.TFval.perm-m600.TFval.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.TFval, breaks = 30)
abline(v=T01.TFval,col=3,lwd=2)
hist(T2.TFval, breaks=30)
abline(v=T02.TFval,col=3,lwd=2)
hist(T3.TFval, breaks=30)
abline(v=T03.TFval,col=3,lwd=2)


# p-value
p1 <- sum(T1.TFval>=T01.TFval)/B
p1
p2 <- sum(T2.TFval>=T02.TFval)/B
p2
p3 <- sum(T3.TFval>=T03.TFval)/B
p3                         ### all cases seems pretty different in particular 2 and 3!!

### TEST 2 --> EXPLOIT BOXPLOT
T01.TFval=as.double((box.TFval$conf[,1] - box.TFval$conf[,2])%*%(box.TFval$conf[,1] - box.TFval$conf[,2]))
T02.TFval=as.double((box.TFval$conf[,1] - box.TFval$conf[,3])%*%(box.TFval$conf[,1] - box.TFval$conf[,3]))
T03.TFval=as.double((box.TFval$conf[,1] - box.TFval$conf[,4])%*%(box.TFval$conf[,1] - box.TFval$conf[,4]))

set.seed(54)

B=10000
T1.TFval=numeric(B)
T2.TFval=numeric(B)
T3.TFval=numeric(B)
n=length(emis.TFval)

for(i in 1:B){
  perm = sample(1:n)
  TFval.perm = emis.TFval[perm]
  box.TFval.perm=boxplot(TFval.perm[i150.emis], TFval.perm[i300.emis], TFval.perm[i450.emis], TFval.perm[i600.emis], plot=FALSE)
  
  # Test statistic:
  T1.TFval[i]=(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,2])%*%(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,2])
  T2.TFval[i]=(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,3])%*%(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,3])
  T3.TFval[i]=(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,4])%*%(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,4])
}

x11()
par(mfrow=c(1,3))
hist(T1.TFval, breaks = 30)
abline(v=T01.TFval,col=3,lwd=2)
hist(T2.TFval, breaks=30)
abline(v=T02.TFval,col=3,lwd=2)
hist(T3.TFval, breaks=30)
abline(v=T03.TFval,col=3,lwd=2)


# p-value
p1 <- sum(T1.TFval>=T01.TFval)/B
p1
p2 <- sum(T2.TFval>=T02.TFval)/B
p2
p3 <- sum(T3.TFval>=T03.TFval)/B
p3                            ### case 2 and 3!!

### TEST 3: USE MEDIAN
med150.TFval=median(emis.TFval[i150.emis])
med300.TFval=median(emis.TFval[i300.emis])
med450.TFval=median(emis.TFval[i450.emis])
med600.TFval=median(emis.TFval[i600.emis])
T01.TFval=abs(med150.TFval-med300.TFval)
T02.TFval=abs(med150.TFval-med450.TFval)
T03.TFval=abs(med150.TFval-med600.TFval)

set.seed(54)

B=10000
T1.TFval=numeric(B)
T2.TFval=numeric(B)
T3.TFval=numeric(B)
n=length(emis.TFval)

for(i in 1:B){
  perm = sample(1:n)
  TFval.perm = emis.TFval[perm]
  med150.TFval.perm = median(TFval.perm[i150.emis])
  med300.TFval.perm = median(TFval.perm[i300.emis])
  med450.TFval.perm = median(TFval.perm[i450.emis])
  med600.TFval.perm = median(TFval.perm[i600.emis])
  
  # Test statistic:
  T1.TFval[i] <- abs(med150.TFval.perm-med300.TFval.perm)
  T2.TFval[i] <- abs(med150.TFval.perm-med450.TFval.perm)
  T3.TFval[i] <- abs(med150.TFval.perm-med600.TFval.perm)
}

x11()
par(mfrow=c(1,3))
hist(T1.TFval, breaks = 30)
abline(v=T01.TFval,col=3,lwd=2)
hist(T2.TFval, breaks=30)
abline(v=T02.TFval,col=3,lwd=2)
hist(T3.TFval, breaks=30)
abline(v=T03.TFval,col=3,lwd=2)


# p-value
p1 <- sum(T1.TFval>=T01.TFval)/B
p1
p2 <- sum(T2.TFval>=T02.TFval)/B
p2
p3 <- sum(T3.TFval>=T03.TFval)/B
p3                         ### again 


### Based on this many tests we can conclude that we have a strong evidence that high temperatures have an impact
### on the value of the transparency feature in emissivity. So maybe we should keep this in mind when building the model. 

### explore for this case also the reflectance:

### TEST 1: JUST USE MEAN
m20.TFval=mean(refl.TFval[i20.refl])
m500.TFval=mean(refl.TFval[i500.refl])

T0.TFval=abs(m20.TFval-m500.TFval)


set.seed(54)

B=10000
T.TFval=numeric(B)
n=length(refl.TFval)

for(i in 1:B){
  perm = sample(1:n)
  TFval.perm = refl.TFval[perm]
  m20.TFval.perm = mean(TFval.perm[i20.refl])
  m500.TFval.perm = mean(TFval.perm[i500.refl])
  
  # Test statistic:
  T.TFval[i] = abs(m20.TFval.perm-m500.TFval.perm)
}

x11()
hist(T.TFval, breaks = 30)
abline(v=T0.TFval,col=3,lwd=2)


# p-value
p = sum(T.TFval>=T0.TFval)/B
p                            ### nope

### TEST 2 --> EXPLOIT BOXPLOT
T0.TFval=as.double((box2.TFval$conf[,1] - box2.TFval$conf[,2])%*%(box2.TFval$conf[,1] - box2.TFval$conf[,2]))


set.seed(54)

B=10000
T.TFval=numeric(B)

n=length(refl.TFval)

for(i in 1:B){
  perm = sample(1:n)
  TFval.perm = refl.TFval[perm]
  box.TFval.perm=boxplot(TFval.perm[i20.refl], TFval.perm[i500.refl], plot=FALSE)
  
  # Test statistic:
  T.TFval[i]=(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,2])%*%(box.TFval.perm$conf[,1] - box.TFval.perm$conf[,2])
}

x11()
hist(T.TFval, breaks = 30)
abline(v=T0.TFval,col=3,lwd=2)


# p-value
p1 <- sum(T.TFval>=T0.TFval)/B
p1                       ### nope

### TEST 3: USE MEDIAN
med20.TFval=median(refl.TFval[i20.refl])
med500.TFval=median(refl.TFval[i500.refl])

T0.TFval=abs(med20.TFval-med500.TFval)


set.seed(54)

B=10000
T.TFval=numeric(B)

n=length(refl.TFval)

for(i in 1:B){
  perm = sample(1:n)
  TFval.perm = refl.TFval[perm]
  med20.TFval.perm = median(TFval.perm[i20.refl])
  med500.TFval.perm = median(TFval.perm[i500.refl])
  
  # Test statistic:
  T.TFval[i] <- abs(med20.TFval.perm-med500.TFval.perm)
}

x11()
hist(T.TFval, breaks = 30)
abline(v=T0.TFval,col=3,lwd=2)



# p-value
p1 <- sum(T.TFval>=T0.TFval)/B
p1          ### nope again



