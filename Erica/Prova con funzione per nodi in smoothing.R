## provo a distribuire meglio i nodi per la spline in modo che non ci sia "overfitting"

library(fda)


srp <- read_excel("Dataset/Srp_Emissivity.xlsx")
head(srp)
summary(srp)

vul <- read_excel("Dataset/Vulcano_Emissivity.xlsx")
head(vul)
summary(vul)

srp2 <- read_excel("Dataset/Srp_Reflectance.xlsx")
head(srp2)
summary(srp2)

vul2 <- read_excel("Dataset/Vulcano_Reflectance.xlsx")
head(vul2)
summary(vul2)

# taglio gli estemi in quanto non utili
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


#### fda ####

find_knots <- function(n,assex){
  
  nodes1 <- seq(from = assex[1], to = assex[334],length.out = n*3/4)
  nodes2 <- seq(from = assex[334], to = assex[length(assex)], length.out = n*1/4)
  nodes2 <- nodes2[2:length(nodes2)]
  return(c(nodes1,nodes2))
}

m = 5

##### REFLEXIVITY #####

# prova di come posiziona i nodi la funzione

n = 26
nodi = find_knots(n,refl.x)

x11()
par(mfrow = c(1,1))
matplot(refl.x,data_refl, type = 'l')
points(nodi,rep(0.94,length(nodi)), pch = 4)
abline(h = 0.94, col = 'gray')
abline(v = 12.75)

nknots <- 6:50
gcv <- numeric(length(nknots))
for (i in 1:length(nknots)){
  nodi <- find_knots(n = nknots[i], assex = refl.x)
  basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), breaks = nodi, norder = m)
  gcv[i] <- smooth.basis(refl.x, as.matrix(data_refl), basis)$gcv
}
par(mfrow=c(1,1))
plot(nknots,gcv)
nknots[which.min(gcv)]

n.opt.knots=24
nodi <- find_knots(n = n.opt.knots, assex = refl.x)
basis <- create.bspline.basis(c(refl.x[1],refl.x[length(refl.x)]), breaks = nodi, norder = m)

basismat <- smooth.basis(argvals=refl.x, y=as.matrix(data_refl), fdParobj=basis)
refl.smooth <- eval.fd(refl.x, basismat$fd)

# first derivative
refl.der1 <- eval.fd(refl.x, basismat$fd, Lfd=1)
# second derivative
refl.der2 <- eval.fd(refl.x, basismat$fd, Lfd=2)

x11()
par(mfrow=c(2,2))
matplot(refl.x,data_refl,xlab="t",ylab="observed data", type='l')
matplot(refl.x,refl.smooth ,type="l",lwd=1)
matplot(refl.x,refl.der1 ,type="l",lwd=1)
matplot(refl.x,refl.der2 ,type="l",lwd=1)


# emissivity

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


#### l2 norm delle derivate dei segnali di emissività ####

emis.der.l2norm = numeric(44)
for (i in 1:44){
  emis.der.l2norm[i] = norm(emis.der1[i,], type = '2')
}

par(mfrow=c(1,1))
matplot(emis.x,emis.der1[,c(1,44)], type = 'l')

plot(emis.der.l2norm,tupe='l')

plot(emis.der.l2norm, silica.emis)

#### l2 norm delle derivate dei segnali di riflettanza ####

n.refl = dim(data_refl)[2]

refl.der.l2norm = numeric(n.refl)
for (i in 1:n.refl){
  refl.der.l2norm[i] = norm(refl.der1[i,], type = '2')
}

# par(mfrow=c(1,1))
# matplot(refl.x,refl.der1[,c(1,n.refl)], type = 'l')

plot(refl.der.l2norm, silica.refl)
plot(refl.der.l2norm, alcali.refl)

scatterplotMatrix(data.frame(silica.refl,refl.CF,refl.CFval,refl.TF,refl.TFval,refl.der.l2norm))


#### GAM silica ~ f(refl) ####
silica.refl.gam = gam(silica.refl ~ refl.CF + s(refl.TFval, bs='cr') + refl.TF + 
                        I(refl.TF^2) + s(refl.der.l2norm, bs = 'cr'))
summary(silica.refl.gam)

hist(silica.refl.gam$residuals)
qqnorm(silica.refl.gam$residuals)
shapiro.test(silica.refl.gam$residuals)

silica.refl.gam.reduced = gam(silica.refl ~ refl.CF + refl.TFval + s(refl.der.l2norm, bs = 'cr') )
summary(silica.refl.gam.reduced)

hist(silica.refl.gam.reduced$residuals)
qqnorm(silica.refl.gam.reduced$residuals)
shapiro.test(silica.refl.gam.reduced$residuals)

# potrebbe funzionare



#### GAM alcali ~ f(refl) ####

scatterplotMatrix(data.frame(alcali.refl,refl.CF,refl.CFval,refl.TF,refl.TFval,refl.der.l2norm,temp.refl))

alcali.refl.gam = gam(alcali.refl ~ refl.CF + s(refl.TFval, bs='cr') +  s(refl.der.l2norm, bs = 'cr'))
summary(alcali.refl.gam)

hist(alcali.refl.gam$residuals)
qqnorm(alcali.refl.gam$residuals)
shapiro.test(alcali.refl.gam$residuals)

# no funziona meglio quello senza derivate


#### GAM silica ~ f(emis) ####
# dati i plot penso che la norma l2 dei grafici di emissività non possa contribuire 


scatterplotMatrix(data.frame(silica.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.der.l2norm))
# useless


scatterplotMatrix(data.frame(alcali.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.der.l2norm))
# useless



