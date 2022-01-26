################################################################################
#######                Libraries and dataset loading                     #######
################################################################################

# Very important color initialization
library("RColorBrewer")
x11()
display.brewer.all()

# Import the datasets
library(readxl)

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

################################################################################
#######                        Exploratory plots                         #######
################################################################################

# Graphs wavelength vs emissivity same blend different temperatures in Srp
X11()
par(mfrow=c(2,3))

matplot(srp[,1], srp[,c(2,8,14,20)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='Spectral')[4:1], xlab='Wavelength (micron)', ylab='emissivity', main='RB Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='Spectral')[4:1])

matplot(srp[,1], srp[,c(3,9,15,21)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='Spectral')[4:1], xlab='Wavelength (micron)', ylab='emissivity', main='B2 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='Spectral')[4:1])

matplot(srp[,1], srp[,c(4,10,16,22)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='Spectral')[4:1], xlab='Wavelength (micron)', ylab='emissivity', main='B4 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='Spectral')[4:1])

matplot(srp[,1], srp[,c(5,11,17,23)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='Spectral')[4:1], xlab='Wavelength (micron)', ylab='emissivity', main='B6 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='Spectral')[4:1])

matplot(srp[,1], srp[,c(6,12,18,24)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='Spectral')[4:1], xlab='Wavelength (micron)', ylab='emissivity', main='B8 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='Spectral')[4:1])

matplot(srp[,1], srp[,c(7,3,19,25)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='Spectral')[4:1], xlab='Wavelength (micron)', ylab='emissivity', main='B Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='Spectral')[4:1])


# Graphs wavelength vs emissivity same blend different temperatures in Srp
X11()
par(mfrow=c(2,2))

matplot(srp[,1], srp[,c(2,3,4,5,6,7)], type='l', lty=1, lwd=2, col=brewer.pal(n=6,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='150°')
legend('bottomright', legend=c('RB','B2','B4','B6','B8','B'), lty=1, col=brewer.pal(n=6,name='BrBG'))

matplot(srp[,1], srp[,c(8,9,10,11,12,13)], type='l', lty=1, lwd=2, col=brewer.pal(n=6,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='300°')
legend('bottomright', legend=c('RB','B2','B4','B6','B8','B'), lty=1, col=brewer.pal(n=6,name='BrBG'))

matplot(srp[,1], srp[,c(14,15,16,17,18,19)], type='l', lty=1, lwd=2, col=brewer.pal(n=6,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='450°')
legend('bottomright', legend=c('RB','B2','B4','B6','B8','B'), lty=1, col=brewer.pal(n=6,name='BrBG'))

matplot(srp[,1], srp[,c(20,21,22,23,24,25)], type='l', lty=1, lwd=2, col=brewer.pal(n=6,name='BrBG'),
        xlab='Wavelength (micron)', ylab='emissivity', main='Emissivity - 600°C')
legend('bottomright', legend=c('RB','B2','B4','B6','B8','B'), lty=1, col=brewer.pal(n=6,name='BrBG'))


# Graphs wavelength vs emissivity same blend different temperatures in Vulcano
X11()
par(mfrow=c(2,3))

matplot(vul[,1], vul[,c(2,7,12,17)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='YlOrRd'), xlab='Wavelength (micron)', ylab='emissivity', main='S Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='YlOrRd'))

matplot(vul[,1], vul[,c(3,8,13,18)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='YlOrRd'), xlab='Wavelength (micron)', ylab='emissivity', main='S7 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='YlOrRd'))

matplot(vul[,1], vul[,c(4,9,14,19)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='YlOrRd'), xlab='Wavelength (micron)', ylab='emissivity', main='S5 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='YlOrRd'))

matplot(vul[,1], vul[,c(5,10,15,20)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='YlOrRd'), xlab='Wavelength (micron)', ylab='emissivity', main='S3 Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='YlOrRd'))

matplot(vul[,1], vul[,c(6,11,16,21)], type='l', lty=1, lwd=2, col=brewer.pal(n=4,name='YlOrRd'), xlab='Wavelength (micron)', ylab='emissivity', main='RS Blend')
legend('bottomright', legend=c('150','300','450','600'), lty=1, col=brewer.pal(n=4,name='YlOrRd'))

# Graphs wavelength vs emissivity same blend different temperatures in Vulcano
X11()
par(mfrow=c(2,2))

matplot(vul[,1], vul[,c(2,3,4,5,6)], type='l', lty=1, lwd=2, col=brewer.pal(n=5,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='150°')
legend('bottomright', legend=c('S','S7','S5','S3','RS'), lty=1, col=brewer.pal(n=5,name='BrBG'))

matplot(vul[,1], vul[,c(7,8,9,10,11)], type='l', lty=1, lwd=2, col=brewer.pal(n=5,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='300°')
legend('bottomright', legend=c('S','S7','S5','S3','RS'), lty=1, col=brewer.pal(n=5,name='BrBG'))

matplot(vul[,1], vul[,c(12,13,14,15,16)], type='l', lty=1, lwd=2, col=brewer.pal(n=5,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='450°')
legend('bottomright', legend=c('S','S7','S5','S3','RS'), lty=1, col=brewer.pal(n=5,name='BrBG'))

matplot(vul[,1], vul[,c(17,18,19,20,21)], type='l', lty=1, lwd=2, col=brewer.pal(n=5,name='BrBG'), xlab='Wavelength (micron)', ylab='emissivity', main='600°')
legend('bottomright', legend=c('S','S7','S5','S3','RS'), lty=1, col=brewer.pal(n=5,name='BrBG'))


# Graphs wavelength vs 1-Reflectance same blend different temperatures in Srp
X11()
par(mfrow=c(2,3))

matplot(srp2[,1], srp2[,c(2,3)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='RB Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(srp2[,1], srp2[,c(4,5)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='B2 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(srp2[,1], srp2[,c(6,7)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='B4 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(srp2[,1], srp2[,c(8,9)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='B6 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(srp2[,1], srp2[,c(10,11)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='B8 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(srp2[,1], srp2[,c(12,13)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='B Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

# Graphs wavelength vs 1-Reflectance same temperature different blends in Srp
X11()
par(mfrow=c(1,2))

matplot(srp2[,1], srp2[,c(2,4,6,8,10,12)], type='l', lty=1, lwd=2, col=brewer.pal(n=6,name='BrBG'), xlab='Wavelength (micron)', ylab='1-reflectance', main='room temp')
legend('bottomright', legend=c('RB','B2','B4','B6','B8','B'), lty=1, col=brewer.pal(n=6,name='BrBG'))

matplot(srp2[,1], srp2[,c(3,5,7,9,11,13)], type='l', lty=1, lwd=2, col=brewer.pal(n=6,name='BrBG'), xlab='Wavelength (micron)', ylab='1-reflectance', main='500°')
legend('bottomright', legend=c('RB','B2','B4','B6','B8','B'), lty=1, col=brewer.pal(n=6,name='BrBG'))


# Graphs wavelength vs 1-Reflectance same blend different temperatures in Vulcano
X11()
par(mfrow=c(2,3))

matplot(vul2[,1], vul2[,c(2,3)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='S Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(vul2[,1], vul2[,c(4,5)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='S7 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(vul2[,1], vul2[,c(6,7)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='S5 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(vul2[,1], vul2[,c(8,9)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='S3 Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

matplot(vul2[,1], vul2[,c(10,11)], type='l', lty=1, lwd=2, col=c('lightblue','red'), xlab='Wavelength (micron)', ylab='1-reflectance', main='RS Blend')
legend('bottomright', legend=c('room temp','500°'), lty=1, col=c('lightblue','red'))

# Graphs wavelength vs 1-Reflectance same temperature different blends in Vulcano
X11()
par(mfrow=c(1,2))

matplot(vul2[,1], vul2[,c(2,4,6,8,10)], type='l', lty=1, lwd=2, col=brewer.pal(n=5,name='BrBG'), xlab='Wavelength (micron)', ylab='1-reflectance', main='room temp')
legend('bottomright', legend=c('S','S7','S5','S3','RS'), lty=1, col=brewer.pal(n=5,name='BrBG'))

matplot(vul2[,1], vul2[,c(3,5,7,9,11)], type='l', lty=1, lwd=2, col=brewer.pal(n=5,name='BrBG'), xlab='Wavelength (micron)', ylab='1-reflectance', main='500°')
legend('bottomright', legend=c('S','S7','S5','S3','RS'), lty=1, col=brewer.pal(n=5,name='BrBG'))