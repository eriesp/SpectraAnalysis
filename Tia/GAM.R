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

temp.emis=

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

scatterplotMatrix(data.frame(silica.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.last))


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

scatterplotMatrix(data.frame(silica.refl,refl.CF,refl.CFval,refl.TF,refl.TFval,refl.last))


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

scatterplotMatrix(data.frame(alcali.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.last))


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

scatterplotMatrix(data.frame(alcali.refl,refl.CF,refl.CFval,refl.TF,refl.TFval,refl.last))

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
