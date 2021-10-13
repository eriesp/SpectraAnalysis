library(readxl)

Base=read_excel("Baseline.xlsx")
Tier2=read_excel("Tier2.xlsx")
Tier3=read_excel("Tier3.xlsx")

area=unlist(Base[,3], use.names = FALSE)
country=unlist(Base[,4], use.names = FALSE)
name=unlist(Base[,7], use.names = FALSE)
type=unlist(Base[,9], use.names = FALSE)
coordinates=Base[,c(10,11)]

### Plot of geographical coordinates
x11()
plot(coordinates)
points(coordinates[area=='East',], col='red', pch=16)
points(coordinates[area=='West',], col='Blue', pch=16)
points(coordinates[area=='Central',], col='Green', pch=16)
points(coordinates[area=='Southern',], col='Gold', pch=16)

east_n=length(which(area=='East'))
west_n=length(which(area=='West'))
center_n=length(which(area=='Central'))
south_n=length(which(area=='Southern'))

### Distance from national grid and national border
dis_bord=unlist(Base[,12], use.names = FALSE)
dis_grid=unlist(Base[,13], use.names = FALSE)

bp=boxplot(dis_grid)
unlist(Base[which(dis_grid>200),3])
unlist(Base[which(dis_grid>200),4])
bp$out

boxplot(dis_bord)
unlist(Base[which(dis_bord>150),3])
unlist(Base[which(dis_bord>150),4])

x11()
plot(dis_bord, dis_grid)
points(dis_bord[area=='East'],dis_grid[area=='East'], col='red', pch=16)
points(dis_bord[area=='West'],dis_grid[area=='West'], col='Blue', pch=16)
points(dis_bord[area=='Central'],dis_grid[area=='Central'], col='Green', pch=16)
points(dis_bord[area=='Southern'],dis_grid[area=='Southern'], col='Gold', pch=16)

### PLOT TOTAL COST VS POPULATION

cost=unlist(Base[,36], use.names = FALSE)
pop=unlist(Base[,14], use.names = FALSE)

x11()
plot(pop,cost)
points(pop[area=='East'],cost[area=='East'], col='red', pch=16)
points(pop[area=='West'],cost[area=='West'], col='Blue', pch=16)
points(pop[area=='Central'],cost[area=='Central'], col='Green', pch=16)
points(pop[area=='Southern'],cost[area=='Southern'], col='Gold', pch=16)

unlist(Base[which(cost>20000000),4])

### PLOT TOTAL COST VS DISTANCE FROM NATIONAL GRID

x11()
plot(dis_grid, cost)
points(dis_grid[area=='East'],cost[area=='East'], col='red', pch=16)
points(dis_grid[area=='West'],cost[area=='West'], col='Blue', pch=16)
points(dis_grid[area=='Central'],cost[area=='Central'], col='Green', pch=16)
points(dis_grid[area=='Southern'],cost[area=='Southern'], col='Gold', pch=16)

### PLOT COST VS DEMAND

demand=unlist(Base[,21], use.names = FALSE)

x11()
plot(demand, cost)
points(demand[area=='East'],cost[area=='East'], col='red', pch=16)
points(demand[area=='West'],cost[area=='West'], col='Blue', pch=16)
points(demand[area=='Central'],cost[area=='Central'], col='Green', pch=16)
points(demand[area=='Southern'],cost[area=='Southern'], col='Gold', pch=16)
abline(a=0, b=4000)


### Classification of refugee settlement --> difference in population mainly and distance from border
x11()
plot(dis_bord, pop)
points(dis_bord[type=='Refugee Camp'], pop[type=='Refugee Camp'], col='red', pch=16)
points(dis_bord[type=='Refugee Settlement'], pop[type=='Refugee Settlement'], col='Blue', pch=16)
points(dis_bord[type=='Refugee Location'], pop[type=='Refugee Location'], col='Green', pch=16)
points(dis_bord[type=='Refugee Urban Location'], pop[type=='Refugee Urban Location'], col='Gold', pch=16)

m_camp=mean(pop[type=='Refugee Camp'])
m_settlement=mean(pop[type=='Refugee Settlement'])
m_location=mean(pop[type=='Refugee Location'])
m_urban=mean(pop[type=='Refugee Urban Location'])


library(plyr)
n_of_camps=count(Base$country)





### DEPTH MEASURE ###

library(MASS)
library(rgl)
library(DepthProc)
library(hexbin)
library(packagefinder)
library(aplpack)
library(robustbase)












