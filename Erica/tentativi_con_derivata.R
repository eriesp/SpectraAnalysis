# provo a lavorare con le derivate:
# cosa posso estrarre dai miei dati? 
# la parte che cambia significativamente è quella tra il CF in poi 

# ho le derivate calcolate già nel file FDA di Tia


#### l2 norm delle derivate dei segnali di emissività ####

emis.der.l2norm = numeric(44)
for (i in 1:44){
  emis.der.l2norm[i] = norm(emis.der1[i,], type = '2')
}

par(mfrow=c(1,1))
matplot(emis.x,emis.der1[,c(1,44)], type = 'l')

plot(emis.der.l2norm,tupe='l')

plot(emis.der.l2norm, silica.emis)
plot(emis.der.l2norm, alcali.emis)

# vedo se vi è differenza della derivata tra una temperatura e l'altra (?)
library(fdatest)

mean.line.150 = rowMeans(emis.der1[,c(seq(1,6), seq(25,29))])
mean.line.300 = rowMeans(emis.der1[,c(seq(7,12), seq(30,34))])
mean.line.450 = rowMeans(emis.der1[,c(seq(13,18), seq(35,39))])
mean.line.600 = rowMeans(emis.der1[,c(seq(19,24), seq(40,44))])

x11()
matplot(emis.x, mean.line.150, type='l', col='blue', ylim = c(-0.1,0.1))
matlines(emis.x, mean.line.300, type='l', col='green')
matlines(emis.x, mean.line.450, type='l', col='gold')
matlines(emis.x, mean.line.600, type='l', col='red')

#150 vs 300
test1=IWT2(t(data.matrix(emis.der1[seq(1,417,10),c(seq(1,6), seq(25,29))])),
           t(data.matrix(emis.der1[seq(1,417,10),c(seq(7,12), seq(30,34))])))
x11()
par(mfrow=c(1,2))
plot(test1)

#150 vs 450
test2=IWT2(t(data.matrix(emis.der1[seq(1,417,10),c(seq(1,6), seq(25,29))])),
           t(data.matrix(emis.der1[seq(1,417,10),c(seq(13,18), seq(35,39))])))
x11()
par(mfrow=c(1,2))
plot(test2)

#150 vs 600
test3=IWT2(t(data.matrix(emis.der1[seq(1,417,10),c(seq(1,6), seq(25,29))])),
           t(data.matrix(emis.der1[seq(1,417,10),c(seq(19,24), seq(40,44))])))
x11()
par(mfrow=c(1,2))
plot(test3)

#300 vs 600
test3bis=IWT2(t(data.matrix(emis.der1[seq(1,417,10),c(seq(7,12), seq(30,34))])),
              t(data.matrix(emis.der1[seq(1,417,10),c(seq(19,24), seq(40,44))])))
x11()
par(mfrow=c(1,2))
plot(test3bis)

#450 vs 600
test3tris=IWT2(t(data.matrix(emis.der1[seq(1,417,10),c(seq(13,18), seq(35,39))])),
               t(data.matrix(emis.der1[seq(1,417,10),c(seq(19,24), seq(40,44))])))
x11()
par(mfrow=c(1,2))
plot(test3tris)

# stessi risultati che nelle funzioni non derivate (d'altronde ...)

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
# ok useless


scatterplotMatrix(data.frame(alcali.emis,emis.CF,emis.CFval,emis.TF,emis.TFval,emis.der.l2norm))
# ok useless
