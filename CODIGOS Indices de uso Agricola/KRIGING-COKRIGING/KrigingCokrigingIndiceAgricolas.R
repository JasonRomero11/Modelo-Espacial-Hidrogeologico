#########################################################################################
#######"MODELADO ESPACIAL HIDROGEOLÓGICO PARA DETERMINAR ÍNDICES DE CALIDAD##############
#######Y VULNERABILIDAD DE LAS AGUAS SUBTERRÁNEAS EN LA ZONA CENTRO DE BOYACÁ"###########

#################################################################
####CÁLCULO DE INDICE DE CALIDAD PARA USO AGRICOLA -KRIGING ####
#################################################################

####LIBRERIAS A UTILIZAR

library(sp)
library(spgwr)
library(raster)
library(ggplot2)
library(sf)
library(adespatial)
library(foreign)
library(spatial)
library(spData)
library(rgdal)
library(spatialreg)
library(spdep)
library(tmap)
library(RColorBrewer)
library(intervals)
library(classInt)
library(MASS)
library(gstat)
library(geoR)
library(sgeostat)
library(geospt)
library(scatterplot3d)
library(ggplot2)
library(car)
library(xtable)
library(stargazer)
library(RGeostats)
library(nortest)
library(intamap)
library(ggmap)
library(maptools)
library(maps)
library(ggspatial)
library(corrplot)
library(gridExtra)
MunicipiosZonaEstudio<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO_B/MunicipiosZonaEstudio.shp")
DatosCampo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/MuestroBalanceIonico.shp")

UsoDelSuelo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/BOYACA_CAPACIDAD_VF/UsoSueloZonaEstudio.shp")
ShapeZona<-readOGR("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO/unidad_crono_Boyaca_zona_centro.shp")


Coordenadas<-read.table("C:/PROYECTODEINVESTIGACION/SHAPEFILE/Depuracion.csv",header=TRUE,sep=",")


epsg3116<-"+init=epsg:3116"
ShapeZona<-spTransform(ShapeZona,CRS(epsg3116))

PuntosInterpolar<- spsample(ShapeZona,n=70000,type="regular")

##############################################################################
################  CALCULO DEL % DE SODIO  #############################
##############################################################################
##############################################################################





##########################################################
#######CALCULO DEL %Na ###################################
##########################################################


#####################################
### CONVERSION DE MG/L A meq/L ######
#####################################

####SE TRABAJA CON LOS CATIONES CA,MG,NA,K###########
#####################################################


IonesNa<-data.frame(DatosCampo$CA,DatosCampo$MG,DatosCampo$NA.,DatosCampo$K)

PesosMolecularIonesNa<-c(40.078 ,24.305,22.989,39.1)

ValenciaIonesNa<-c(2,2,1,1)

PesoEquivalenteIones<-PesosMolecularIonesNa/ValenciaIonesNa

IonesNa$Cameq<-((IonesNa$DatosCampo.CA)/PesoEquivalenteIones[1])
IonesNa$Mgmeq<-((IonesNa$DatosCampo.MG)/PesoEquivalenteIones[2])
IonesNa$Nameq<-((IonesNa$DatosCampo.NA.)/PesoEquivalenteIones[3])
IonesNa$Kmeq<-((IonesNa$DatosCampo.K)/PesoEquivalenteIones[4])



########################################
## CALCULO DEL PORCENTAJE DE SODIO #####
########################################

IonesNa$PorcentajeSodio<-(((IonesNa$Nameq+IonesNa$Kmeq)/(IonesNa$Cameq+IonesNa$Mgmeq+IonesNa$Nameq)))*100
DatosCampo$PorceSodio<-IonesNa$PorcentajeSodio
summary(DatosCampo$PorceSodio)


p1<-qplot(DatosCampo$CE__25C_,geom="histogram",binwidth = 50, xlab = "Conductividad eléctrica", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p1<-p1+scale_x_continuous(breaks =c(0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800))

p2<-qplot(DatosCampo$PorceSodio,geom="histogram", xlab = "%Na", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))

p2<-p2+scale_x_continuous (breaks = c(10,20,30,40,50,60,70,80,90,100))

grid.arrange(p2,p1,nrow = 1, widths=c(2,2))

C1S1<-DatosCampo[which(DatosCampo$PorceSodio<20 & DatosCampo$CE__25C_<250),]  
C2S1<-DatosCampo[which(DatosCampo$PorceSodio>20 & DatosCampo$PorceSodio<40),] #42
C2S2<-DatosCampo[which(DatosCampo$CE__25C_>250 & DatosCampo$CE__25C_<750),]
C3S1<-DatosCampo[which(DatosCampo$PorceSodio>40 & DatosCampo$PorceSodio<60 & DatosCampo$CE__25C_>750 & DatosCampo$CE__25C_<2000 ),] #42
C4S2<-DatosCampo[which(DatosCampo$PorceSodio>80 &  DatosCampo$CE__25C_>3000),] #42
C4S1<-RAS[which(RAS$RAS<4 & RAS$DatosCampo.CE__25C_>2250),]

###################################################################3
###CLASIFICACION DEL INDICE SEGÚN WILCOX ##########################
Excelente<-DatosCampo[which(DatosCampo$PorceSodio<20 & DatosCampo$CE__25C_<250),]  #38
Buena<-DatosCampo[which(DatosCampo$PorceSodio>=20 & DatosCampo$PorceSodio<=40 & DatosCampo$CE__25C_>=250& DatosCampo$CE__25C_<750),] #17
Permisible<-DatosCampo[which(DatosCampo$PorceSodio>40 & DatosCampo$PorceSodio<=60 & DatosCampo$CE__25C_>750& DatosCampo$CE__25C_<=2000),] #2
Deficiente<-DatosCampo[which(DatosCampo$PorceSodio>60 & DatosCampo$PorceSodio<=80 & DatosCampo$CE__25C_>2000& DatosCampo$CE__25C_<=3000),] #2
Inadecuada<-DatosCampo[which(DatosCampo$PorceSodio>80 & DatosCampo$CE__25C_>3000),] #2

#################################################
#################################################
##### ANÁLISIS DESCRIPTIVO %Na ##################
#################################################

##ANALSIS DE NORMALIDAD
##METODO GRÁFICO

X11()
qplot(DatosCampo$PorceSodio,geom = 'blank', xlab = "(%Na)") +   
  geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
  theme(legend.position = c(0.5, 0.5)) 

shapiro.test(DatosCampo$PorceSodio) #p-value < 8.564e-08 no existe normalidad
ks.test(as.numeric(scale(sort(DatosCampo$PorceSodio))),pnorm) #p-value < 0.006684 no existe normalidad

###################################################################################################
###################################################################################################
###############################CORRECCION DE NORMALIDAD ###########################################
###################################################################################################
###################################################################################################


###########################################################################
########CORRECCION DE NORMALIDAD POR ANAMORFOSIS GAUSSIANA#################
###########################################################################

PorceSodioselect = DatosCampo[,c("X","Y","PorceSodio")] 
Na.rgdb <- db.create(PorceSodioselect,ndim=2,autoname=F)
Na.herm <- anam.fit(Na.rgdb,name="PorceSodio",type="gaus")
Na.hermtrans <- anam.z2y(Na.rgdb,names="PorceSodio",anam=Na.herm)
Na.trans <- Na.hermtrans@items$Gaussian.PorceSodio


DatosCampo$NaTransGauss<-Na.trans

shapiro.test(DatosCampo$NaTransGauss) #p-value 0.9154 existe normalidad
ks.test(as.numeric(scale(sort(DatosCampo$NaTransGauss))),pnorm) #p-value 0.9999 existe normalidad



############################################################################
#######CORRECCION DE NORMALIDAD POR BOX-COX ################################
############################################################################


NABoxCox<-boxcox(DatosCampo$PorceSodio~1)
lambdaNa <- NABoxCox$x[which.max(NABoxCox$y)]
BoxNA<- (((DatosCampo$PorceSodio)^lambdaNa)-1)/lambdaNa
shapiro.test(BoxNA) #p-value 0.09327  > 0.05 % Ya existe normalidad
ks.test(as.numeric(scale(sort(BoxNA))),pnorm) #p-value  0.8859 ya existe normalidad

DatosCampo$NaTrasBox<-BoxNA

p1<-qplot(DatosCampo$NaTrasBox,geom="histogram",main = "", xlab = "%Na Transformado Box-Cox", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p2<-qplot(DatosCampo$NaTrasBox,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="%Na Transformado Box-Cox",ylab="Densidad")
p3<-qplot(DatosCampo$NaTransGauss,geom="histogram",main = "", xlab = "%Na Transformado Anamorfosis Gaussiana", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p4<-qplot(DatosCampo$NaTransGauss,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="%Na Transformado Anamorfosis Gaussiana",ylab="Densidad")

grid.arrange(p1,p2,p3,p4,ncol = 2)
######NOTA###
#SE DECIDE TRABAJAR CON ANAMORFOSIS GAUSSIANA YA NORMALIZA LA VARIABLE %Na


####################################################
##Análisis de tendencia y elección del mejor modelo#
####################################################

ModeloGeneral<-lm(DatosCampo$NaTrasBox~DatosCampo$T+DatosCampo$SDT+DatosCampo$PH+DatosCampo$COLOR+DatosCampo$DUREZA+DatosCampo$CE__25C_+DatosCampo$ALCALINIDA+DatosCampo$HCO3+DatosCampo$SO4)
summary(ModeloGeneral)
shapiro.test(ModeloGeneral$residuals)# 

ModeloGenera2<-lm(DatosCampo$NaTrasBox~DatosCampo$T+DatosCampo$SDT+DatosCampo$PH+DatosCampo$COLOR+DatosCampo$DUREZA)
summary(ModeloGenera2)
shapiro.test(ModeloGeneral$residuals)#

shapiro.test(ModeloGenera2$residuals)#
ks.test(as.numeric(scale(sort(ModeloGenera2$residuals))),pnorm) #p-value0.1473 existe normalidad

Modelo2<-lm(DatosCampo$NaTrasBox~DatosCampo$T+DatosCampo$COLOR+DatosCampo$DUREZA+DatosCampo$SO4+DatosCampo$ALCALINIDA)  
summary(Modelo2)
shapiro.test(Modelo2$residuals)##P-value 0.3282
ks.test(as.numeric(scale(sort(Modelo2$residuals))),pnorm) #p-value0.1473 existe normalidad


AICModeloG<-stepAIC(ModeloGenera2,direction = "both") #AIC=



p1<-qplot(ModeloGenera2$residuals,geom="histogram",main = "", xlab = "Histograma residuos modelo lineal ", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p2<-qplot(ModeloGenera2$residuals,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="Densidad residuos modelo lineal ",ylab="Densidad")
grid.arrange(p1,p2, ncol = 2)
#########NOTA
#########NINGUN MODELO PRESENTA TENDENCIA, POR LO TANTO SE TRABAJA CON K.O Y K.S
######################################################################
################ANALISIS DE ANISOTROPIA###############################
######################################################################
TablaPorcSodio<-data.frame(cbind(DatosCampo$X,DatosCampo$Y,DatosCampo$NaTransGauss))

geodata.anisotropia<-as.geodata(TablaPorcSodio,coords.col = 1:2,data.col =3)
Dmax<-sqrt((max(TablaPorcSodio[, 1])-min(TablaPorcSodio[, 1]))^2+(max(TablaPorcSodio[, 2])-min(TablaPorcSodio[, 2]))^2)
Dusar<-Dmax/2  ###Distancia a usar para modelar el variograma 29617.13
x11()
plot(variog4(geodata.anisotropia,max.dist = Dusar), xlab = "Distancia (m)", ylab = "Semivariograma estimado", legend = F)
legend(locator(1), legend = c(expression(0*degree), expression(45*degree), expression(90*degree), expression(135*degree)), col = 1:4, lty = 1:4)
title("Semivariogramas experimentales (estimador clasico)")
##métodos "estadisticos" para la anisotropia

library(intamap)
Coorxy<- geodata.anisotropia$coords
ZNa<-TablaPorcSodio$X3
Puntos= data.frame(Coorxy,DatosCampo$Z, z=ZNa)

##Asignamos coordenadas a los puntos de muestreo
coordinates(Puntos) <- ~X1+X2
SpatialPoints(coordinates(Puntos))
class(Puntos)
estimateAnisotropy(Puntos)
####NO EXISTE ANISOTROPÍA

#######################################################################
#######################################################################
######################CALCULO VARIOGRAMA EXPERIMENTAL##################
#######################################################################

Base<-data.frame(geodata.anisotropia)  
Point <- point(Base, x = "X1", y = "X2") 
Pair <- pair(Point,num.lags=20,maxdist=Dusar)
EstimaSV<- est.variograms(Point,Pair,"data",trim=0.1) 
x11() 

layout(matrix(c(1:4), nrow=2, byrow=FALSE))

layout.show(4) 
x11()
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$med,type="b",ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental a partir de la mediana",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$classic,type="b",ylim=c(0,3), col =1,pch= 16,main = "Modelo experimental Clásico",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,3), col =1,pch= 16,main = "Modeloxperimental media recortada",xlab="Distancia", ylab="Semivarianza") 

##NOTA
##SE SELECCIONA EL ROBUSTO 
#####################################################
############VARIOGRAMA SIMPLE########################
#####################################################

x11()
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia", ylab="Semivarianza") 

###MODELOS EXPERIMENTALES 
Exponencial<-likfit(geodata = geodata.anisotropia,nugget = 0.15,ini = c(1.0,6500),fix.nug=T)
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,3), col =1,pch= 16, main = "", xlab="Distancia", ylab="Semivarianza") 

lines(Exponencial,max.dist=Dusar,lwd=3,col="blue") 



Esferico<-likfit(geodata = geodata.anisotropia,nugget = 0.16, ini = c(1.0,8500),cov.model="sph",fix.nug=T) #Definitivo

lines(Esferico,max.dist=Dusar,lwd=3,col="red") 

Matern<-likfit(geodata = geodata.anisotropia,nugget = 0.18, ini = c(1.5,8900),cov.model="mat",kappa=0.41,fix.kappa = T,fix.nug=T) 

lines(Matern,max.dist=Dusar,lwd=3,col="Black")

Circular<-likfit(geodata = geodata.anisotropia,nugget = 0.2, ini = c(1.2,8500),cov.model="cir",fix.nug=T) 


lines(Circular,max.dist=Dusar,lwd=3,col="yellow")


PowerExponencial<-likfit(geodata = geodata.anisotropia, nugget = 0.2,ini = c(1.3,9300),
                         cov.model="powered.exponential",kappa=0.80,fix.nug=T) 

lines(PowerExponencial,max.dist=Dusar,lwd=3,col="grey") 


lines(Exponencial,max.dist=Dusar,lwd=3,col="black") 
lines(Esferico,max.dist=Dusar,lwd=3,col="red") 
lines(Matern,max.dist=Dusar,lwd=3,col="yellow") 
lines(Circular,max.dist=Dusar,lwd=3,col="blue") 
lines(PowerExponencial,max.dist=Dusar,lwd=3,col="orange") 
legend(locator(1),c('Exponencial','Esferico','Matern','Circular','Potencial'),col=c("black","red","yellow","blue", "orange"), lty=c(1,1,1,1,1,1,1))

#AIC
Exponencial$AIC #594.6754
Esferico$AIC #595.1347
Matern$AIC #588.2187
Circular$AIC #594.797
PowerExponencial$AIC #586.0402


##############################################
#########vARIOGRAMA ROBUSTO###################
##############################################

##CON MATHERN ############

PowerExponencial<-likfit(geodata = geodata.anisotropia, nugget = 0.2,ini = c(1.3,9300),
                         cov.model="powered.exponential",kappa=0.80,fix.nug=T) 

VariogramaRobusto<-variog( geodata.anisotropia,trend="cte",max.dist=Dusar, 
                           option = "cloud",estimator.type="modulus")

Maternwls<-variofit(vario = VariogramaRobusto, nugget = 0.85,ini = c(1.8, 10000),
                    kappa = 0.58,fix.nugget= T,weights="npairs",cov.model = "stable")

Maternml<-likfit(geodata = geodata.anisotropia, nugget = 0.08,ini = c(1.6, 10000),
                 kappa = 0.6,fix.nugget= T,cov.model = "stable")

Maternrml<-likfit(geodata = geodata.anisotropia, nugget = 0.08,ini = c(1.5,9800),
                  kappa = 0.7,fix.nugget= T,method='RML',cov.model = "stable")

x11()##Eleccion del método de estimación por medios gráficos
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,3), col =1,pch= 16, main = "", xlab="Distancia", ylab="Semivarianza") 

lines(Maternwls,max.dist=Dusar,lwd=2,lty=1,col="coral")##mejor ajuste
lines(Maternml,max.dist=Dusar,lwd=2,lty=1,col="orange")
lines(Maternrml,max.dist=Dusar,lwd=3,lty=1, col="yellow")
legend(locator(1),c('ML','RML','MCP'),lty=c(1,1,1),col=c("orange","yellow","coral"))

#############################################
#############################################
#############################################
##Métodos de interpolación geoestadísticos###
#############################################
#############################################



##############################################
##### KRIGING SIMPLE #########################
zNa<-DatosCampo$NaTransGauss
Media  = mean(zNa)                                                         
So <- c(1118001,1123900)
So <- as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So) <- c("x", "y")

xy <- TablaPorcSodio[,1:2]


pts = data.frame(xy, z=zNa) ##COORDENADAS %Na
coordinates(pts) <- c("X1", "X2")
proj4string(pts)<-CRS("+init=epsg:3116")


gridded(PuntosInterpolar) <- TRUE



Prediccion.ks <- krige(z~1, pts, PuntosInterpolar, vgm(1.8,"Exc",10000,0.85,kappa = 0.58,beta=Media)) ####kRIGING SIMPLE
Prediccion.ko <- krige(z~1, pts, PuntosInterpolar, vgm(1.8,"Exc",10000,0.85,kappa = 0.58)) ####kRIGING ORDINARIO



####################################################
####################################################
##KRIGING ORDINARIO ANTITRANSFORMACION PREDICCIONES
####################################################

coords1 <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords1, PuntosInterpolar)
Na.datos = cbind(coordinates(DataCoords),Prediccion.ko$var1.pred)
Na.db = db.create(Na.datos,autoname = F)
Napred = anam.y2z(Na.db,names="V5",anam = Na.herm)

#Prediction map
Prediccion.ko$var1.pred<- Napred@items$Raw.V5
summary(Prediccion.ko$var1.pred)

pred_meanNaKo<-Prediccion.ko$var1.pred
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_meanNaKo, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
koNa<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()

#VARIANZA

varianzaKoNa<-Prediccion.ko$var1.var
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = varianzaKoNa, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
VkoNa<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()


#################################################
##KRIGING SIMPLE ANTITRANSFORMACION PREDICCIONES
#################################################

coords1 <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords1, PuntosInterpolar)
Na.datos = cbind(coordinates(DataCoords),Prediccion.ks$var1.pred)
Na.db = db.create(Na.datos,autoname = F)
NapredKS = anam.y2z(Na.db,names="V5",anam = Na.herm)


#Prediction map

Prediccion.ks$var1.pred<- NapredKS@items$Raw.V5
summary(Prediccion.ko$var1.pred)
summary(Prediccion.ks$var1.pred)

pred_meanKsNa<-Prediccion.ks$var1.pred
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_meanKsNa, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
ksNa<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "%Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()
grid.arrange(ko,ks, ncol = 2)

#VARIANZA

varianzaKsNa<-Prediccion.ks$var1.var
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = varianzaKsNa, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
Vks<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()



stable<-vgm(1.8,"Exc",10000,0.85,kappa = 0.58)
ValidacionKs<-krige.cv(z~1,pts,TablaPorcSodio,model=stable,beta=Media,maxdist=3948.951)
ValidacionKo<-krige.cv(z~1,pts,TablaPorcSodio,model=stable,maxdist=3948.951)

class(Prediccion.ks)
ValidacionFinalKs<- criterio.cv(ValidacionKs)
#MPE              ASEPE     RMSPE         MSPE
#-0.002274689   1.22267 0.91116246 -0.001509596
#0.001135893 1.230936 0.9084671 0.001079823 0.7278262 2.302668 0.4242891
criterio.cv(ValidacionKo)
#0.004818776 1.237228 0.912048 0.003909972 0.7254257 2.572995 0.4216885
library(SpecsVerification)


GaussCrps(mean(ValidacionKs$var1.pred),sd(ValidacionKs$var1.pred),mean(ValidacionKs$observed))
#0.1011479


GaussCrps(mean(ValidacionKo$var1.pred),sd(ValidacionKo$var1.pred),mean(ValidacionKo$observed))
#0.111

dataa<-data.frame(ValidacionKo$var1.pred,ValidacionKo$var1.var,ValidacionKo$observed,ValidacionKo$residual)
Export(dataa,"ResultadosFinaleskrigingOrdinarioPorSodio.csv")

dataa<-data.frame(ValidacionKs$var1.pred,ValidacionKs$var1.var,ValidacionKs$observed,ValidacionKs$residual)
Export(dataa,"ResultadosFinaleskrigingSimplePorSodio.csv")


####################################################
#### PREDICCION COKRIGING  ORDINARIO################
####################################################


################################################
#####MATRIZ DE CORRELACION #####################
Dureza<-DatosCampo$DUREZA
Color<-DatosCampo$COLOR
CE<-DatosCampo$CE__25C_
SDT<-DatosCampo$SDT
HCO3<-DatosCampo$HCO3
PH<-DatosCampo$PH
NO3<-DatosCampo$NO3
CL<-DatosCampo$CL
Tem<-DatosCampo$T
SO4<-DatosCampo$SO4
Alca<-DatosCampo$ALCALINIDA
NaTransformado<-DatosCampo$NaTransGauss
VariablesInteres<-data.frame(Tem,SDT,PH,Color,Dureza,CE,HCO3,CL,SO4,Alca,NaTransformado)
library(corrplot)
Correlacion<-cor(VariablesInteres,method="pearson")
round(Correlacion,digits = 2)
x11()
corrplot(Correlacion,type="upper",method="shade", tl.col = "black", tl.srt = 45,addCoef.col = "Black", diag=F, addshade = "all") 

##LAS VARIABLES MÁS CORRELACIONAS SON: 
#PH, DUREZA,CL PROBAR CON SDT Y CE

##CON DUREZA
DUREZAselect = DatosCampo[,c("X","Y","DUREZA")] 
DUREZA.rgdb <- db.create(DUREZAselect,ndim=2,autoname=F)
DUREZA.herm <- anam.fit(DUREZA.rgdb,name="DUREZA",type="gaus")
DUREZA.hermtrans <- anam.z2y(DUREZA.rgdb,names="DUREZA",anam=DUREZA.herm)
DUREZA.trans <- DUREZA.hermtrans@items$Gaussian.DUREZA
shapiro.test(DUREZA.trans) #normalizo 
ks.test(as.numeric(scale(sort(DUREZA.trans))),pnorm) #

##CON PH
PHselect = DatosCampo[,c("X","Y","PH")] 
PH.rgdb <- db.create(PHselect,ndim=2,autoname=F)
PH.herm <- anam.fit(PH.rgdb,name="PH",type="gaus")
PH.hermtrans <- anam.z2y(PH.rgdb,names="PH",anam=PH.herm)
PH.trans <- PH.hermtrans@items$Gaussian.PH
shapiro.test(PH.trans) #normalizo
ks.test(as.numeric(scale(sort(PH.trans))),pnorm) #

#hco3
HCO3select = DatosCampo[,c("X","Y","HCO3")] 
HCO3.rgdb <- db.create(HCO3select,ndim=2,autoname=F)
HCO3.herm <- anam.fit(HCO3.rgdb,name="HCO3",type="gaus")
HCO3.hermtrans <- anam.z2y(HCO3.rgdb,names="HCO3",anam=HCO3.herm)
HCO3.trans <- HCO3.hermtrans@items$Gaussian.HCO3
shapiro.test(HCO3.trans) ##
ks.test(as.numeric(scale(sort(HCO3.trans))),pnorm) #

##Color
CLOROselect = DatosCampo[,c("X","Y","COLOR")] 
CLORO.rgdb <- db.create(CLOROselect,ndim=2,autoname=F)
CLORO.herm <- anam.fit(CLORO.rgdb,name="COLOR",type="gaus")
CLORO.hermtrans <- anam.z2y(CLORO.rgdb,names="COLOR",anam=CLORO.herm)
CLORO.trans <- CLORO.hermtrans@items$Gaussian.COLOR

shapiro.test(CLORO.trans) #
ks.test(as.numeric(scale(sort(CLORO.trans))),pnorm) #


##CON CL
CLselect = DatosCampo[,c("X","Y","CL")] 
CL.rgdb <- db.create(CLselect,ndim=2,autoname=F)
CL.herm <- anam.fit(CL.rgdb,name="CL",type="gaus")
CL.hermtrans <- anam.z2y(CL.rgdb,names="CL",anam=CL.herm)
CL.trans <- CL.hermtrans@items$Gaussian.CL
shapiro.test(CL.trans) # 
ks.test(as.numeric(scale(sort(CL.trans))),pnorm) #


##Conductividad electrica
CEselect = DatosCampo[,c("X","Y","CE__25C_")] 
CE.rgdb <- db.create(CEselect,ndim=2,autoname=F)
CE.herm <- anam.fit(CE.rgdb,name="CE__25C_",type="gaus")
CE.hermtrans <- anam.z2y(CE.rgdb,names="CE__25C_",anam=CE.herm)
CE.trans <- CE.hermtrans@items$Gaussian.CE__25C_
shapiro.test(CE.trans) ####normalizo
ks.test(as.numeric(scale(sort(CE.trans))),pnorm) #


##SDT
SDTselect = DatosCampo[,c("X","Y","SDT")] 
SDT.rgdb <- db.create(SDTselect,ndim=2,autoname=F)
SDT.herm <- anam.fit(SDT.rgdb,name="SDT",type="gaus")
SDT.hermtrans <- anam.z2y(SDT.rgdb,names="SDT",anam=SDT.herm)
SDT.trans <- SDT.hermtrans@items$Gaussian.SDT
shapiro.test(SDT.trans) ####normalizo
ks.test(as.numeric(scale(sort(SDT.trans))),pnorm) 

##TEM 
TEMselect = DatosCampo[,c("X","Y","T")] 
TEM.rgdb <- db.create(TEMselect,ndim=2,autoname=F)
TEM.herm <- anam.fit(TEM.rgdb,name="T",type="gaus")
TEM.hermtrans <- anam.z2y(TEM.rgdb,names="T",anam=TEM.herm)
TEM.trans <- TEM.hermtrans@items$Gaussian.T
shapiro.test(TEM.trans) ##Normalizo 
ks.test(as.numeric(scale(sort(TEM.trans))),pnorm) #

Coorx<-DatosCampo$X
Coory<-DatosCampo$Y
DataCokriging<-data.frame(Coorx,Coory,DatosCampo$NaTransGauss,PH.trans,DUREZA.trans,SDT.trans,HCO3.trans,TEM.trans)
coordinates(DataCokriging)<-~Coorx+Coory
proj4string(DataCokriging)<-CRS("+init=epsg:3116")

CoKo=gstat(id = "NAPorce", formula = DatosCampo.NaTransGauss~1, data = DataCokriging)
CoKo=gstat(CoKo,id = "Tem", formula = TEM.trans~1, data = DataCokriging)
CoKo=gstat(CoKo,id = "PH", formula = PH.trans~1, data = DataCokriging)
CoKo=gstat(CoKo,id = "DUREZA", formula = DUREZA.trans~1, data = DataCokriging)
CoKo=gstat(CoKo,id = "CL", formula = CL.trans~1, data = DataCokriging)
CoKo=gstat(CoKo,id = "SDT", formula = SDT.trans~1, data = DataCokriging)


CoKo=gstat(CoKo,id = "NAPorce", model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58))  
CoKo=gstat(CoKo,id = "Tem", model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58))
CoKo=gstat(CoKo,id = "PH", model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58))
CoKo=gstat(CoKo,id = "DUREZA", model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58))
CoKo=gstat(CoKo,id = "CL", model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58))
CoKo=gstat(CoKo,id = "SDT", model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58))

Vcros<-variogram(CoKo,cutoff=Dusar)

CoKo<-gstat(CoKo,id="NAPorce",model=vgm(1.8,"Exc",10000,0.85,kappa = 0.58), fill.all=T)

###SE REALIZA EL AJUSTE DE LOS MODELOS CORRELACIONADOS
AjusteCokrigingNa<-fit.lmc(Vcros,CoKo)



plot(variogram(CoKo), model=AjusteCokrigingNa$model,pch=20)

ValidacionCoKoNa=gstat.cv(AjusteCokrigingNa,maxdist=3948.951)


ValidacionCokoNa<-criterio.cv(ValidacionCoKoNa)
#FINAL !!

#           MPE    ASEPE     RMSPE         MSPE   RMSSPE    MAPPE      CCPE        R2
#0.001559102 0.5134206 0.5473662 0.00150978 1.055729 4.140783 0.8384197 0.7018609 0.7029476
#0.001559102 0.5134206 0.5473662 0.00150978 1.055729 4.140783 0.838419

dataa<-data.frame(ValidacionCoKoNa$NAPorce.pred,ValidacionCoKoNa$NAPorce.var,ValidacionCoKoNa$observed,ValidacionCoKoNa$residual)
Export(dataa,"ResultadosTESISCokrigingNa.csv") 



mean(GaussCrps(mean(ValidacionCoKoNa$NAPorce.pred),sd(ValidacionCoKoNa$NAPorce.pred),mean(ValidacionCoKoNa$observed)))
#0.19
#0.1891528

#########################PREDICCIONES COKRIGING ###############

PrediccionCoksNa<-predict(AjusteCokrigingNa,PuntosInterpolar)



pred_meanNaCo <-PrediccionCoksNa$NAPorce.pred
coords <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords, PuntosInterpolar)
Na.datos = cbind(coordinates(DataCoords),pred_meanNaCo)
NA.db = db.create(Na.datos,autoname = F)
Napred = anam.y2z(NA.db,names="pred_meanNaCo",anam = Na.herm)

summary(Napred$Raw.pred_meanNaCo)

pred_meanNaCo<-Napred$Raw.pred_meanNaCo
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_meanNaCo, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable)
ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "%Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()










##############################################################################
################  CALCULO DEL RAS Y RAS AJUSTADO #############################
##############################################################################
##############################################################################

DatosCampo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/MuestroBalanceIonico.shp")





######################################################
#########CONCENTRACION MOLAR DE LOS IONES ############
######IONES A TENER EN CUENTA 
Iones<-data.frame(DatosCampo$CA,DatosCampo$MG,DatosCampo$NA.,DatosCampo$K,DatosCampo$HCO3,DatosCampo$CL,DatosCampo$SO4)

PesosMolecularIones<-c(40.078 ,24.305,22.989,39.1,61.0168,35.453,96.06)

ValenciaIones<-c(2,2,1,1,1,1,2)

PesoEquivalenteIones<-PesosMolecularIones/ValenciaIones

Iones$Cameq<-((Iones$DatosCampo.CA)/PesoEquivalenteIones[1])
Iones$Mgmeq<-((Iones$DatosCampo.MG)/PesoEquivalenteIones[2])
Iones$Nameq<-((Iones$DatosCampo.NA.)/PesoEquivalenteIones[3])
Iones$Kmeq<-((Iones$DatosCampo.K)/PesoEquivalenteIones[4])
Iones$HCO3meq<-((Iones$DatosCampo.HCO3)/PesoEquivalenteIones[5])
Iones$Clmeq<-((Iones$DatosCampo.CL)/PesoEquivalenteIones[6])
Iones$SO4meq<-((Iones$DatosCampo.SO4)/PesoEquivalenteIones[7])


###################################################
#############CALCULO DEL RAS ######################

RAS<-((Iones$Nameq))/((sqrt(((Iones$Cameq)+(Iones$Mgmeq))/2)))

RAS<-data.frame(RAS,DatosCampo$CE__25C_)
summary(RAS)
##CLASIFICACION DEL RAS Y LA CONDUCTIVIDAD ELECTRICA SEGÚN WILCOX 
C1S1<-RAS[which(RAS$RAS<10 & RAS$DatosCampo.CE__25C_<250),]  #150
C2S1<-RAS[which(RAS$RAS<10 & RAS$DatosCampo.CE__25C_>250 & RAS$DatosCampo.CE__25C_<750),] #42
C3S1<-RAS[which(RAS$RAS<10 & RAS$DatosCampo.CE__25C_>750 & RAS$DatosCampo.CE__25C_<2250),] #16
C1S2<-RAS[which(RAS$RAS>10 &RAS$DatosCampo.CE__25C_<250),]
C4S1<-RAS[which(RAS$RAS<4 & RAS$DatosCampo.CE__25C_>2250),]

######################################################
######CALCULO DEL RAS AJUSTADO #######################
######################################################



 I<-((Iones$Cameq*(ValenciaIones[1]^2))+(Iones$Mgmeq*(ValenciaIones[2]^2))+(Iones$Nameq*(ValenciaIones[3]^2))+(Iones$Kmeq*(ValenciaIones[4]^2))+(Iones$HCO3meq*(ValenciaIones[5]^2))+(Iones$Clmeq*(ValenciaIones[6]^2))+(Iones$SO4meq*(ValenciaIones[7]^2)))/2
# 
 IntensidadIones<-data.frame(I)


#########CALCULO DEL PK2-PKC

 pk2pkc<-(2.0269*(0.5092*(((4*(sqrt(IntensidadIones)))/(1+2*(sqrt(IntensidadIones))))+((sqrt(IntensidadIones))/(1+1.45*(sqrt(IntensidadIones)))))))
# 
# 
# ########CALCULO DEL P(CA+MG)
# 
 PCAMG<- -1*log10((Iones$Cameq/1000)+(Iones$Mgmeq/1000))
# 
# ########CALCULO DEL P(Alk)
# 
 meqAlc<-DatosCampo$ALCALINIDA/50.04345
 pAlk<- -1*log10(meqAlc/1000)
# 
# ############################################################
# #############CALCULO DEL RAS AJUSTADO ######################
# ############################################################
 PHc=pk2pkc+PCAMG+pAlk
 RASajustado<-RAS$RAS*(1+(8.4-PHc))
 RASajustado<-data.frame(RASajustado)
 summary(RASajustado)
DatosCampo$RASajustado<-RASajustado$I
 p<-qplot(RASajustado$I,geom="histogram",main = "", xlab = "RASajustado", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
 p<-p+scale_x_continuous (breaks = c(0,5,10,15,20,25,30))


###################################
#######HISTOGRAMA DEL RAS ########
#################################
DatosCampo$RAS<-RAS$RAS
 

p<-qplot(RAS$RAS,geom="histogram", xlab = "RAS (meq/l)", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p<-p+scale_x_continuous (breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))


p1<-qplot(DatosCampo$CE__25C_,geom="histogram",binwidth = 50, xlab = "Conductividad eléctrica", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p1<-p1+scale_x_continuous(breaks =c(0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800))

grid.arrange(p,p1,nrow = 1, widths=c(2,2))


 p3<-ggplot(data = DatosCampo) + geom_sf(aes(size=RAS), color="red", alpha=0.7) +
   geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +ggtitle(" ")+ annotation_north_arrow(location='tr')+annotation_scale()
 p3<-p3+scale_size(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12))
# 
 p3<-ggplot(data = DatosCampo) + geom_sf(aes(size=RASajustado), color="red", alpha=0.7) +
   geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +ggtitle(" ")+ annotation_north_arrow(location='tr')+annotation_scale()
 p3<-p3+scale_size(breaks = c(0,4,8,12,16,20,24,28,32))
 
p4<-ggplot(data = DatosCampo) + geom_sf(aes(size=CE__25C_), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud")+ ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()
p4<-p4+scale_size(breaks = c(100,250,400,550,700,850,1000,1150,1300,1450,1600,1750,1900,2050,2200))
grid.arrange(p3,p4,nrow = 1, widths=c(2,2))





##########################################################
##########################################################
##### ANÁLISIS DESCRIPTIVO RAS AJUSTADO ##################
##########################################################

##ANALSIS DE NORMALIDAD
##METODO GRÁFICO

# X11()
# qplot(DatosCampo$Rasajustado,geom = 'blank', xlab = "ICA") +   
#   geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
#   stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
#   geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
#   scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
#   theme(legend.position = c(0.5, 0.5)) 
# 

##########################################################
##########################################################
##### ANÁLISIS DESCRIPTIVO RAS ###########################
##########################################################

##ANALSIS DE NORMALIDAD
##METODO GRÁFICO

qplot(DatosCampo$RAS,geom = 'blank', xlab = "RAS") +   
  geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
  theme(legend.position = c(0.5, 0.5)) 
qplot(DatosCampo$RAS,geom="histogram",main = "", xlab = "Temperatura", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))


shapiro.test(DatosCampo$RAS) #p-value < 2.2e-16 no existe normalidad
ks.test(as.numeric(scale(sort(DatosCampo$RAS))),pnorm) #p-value < 6.661e-16 no existe normalidad


###################################################################################################
###################################################################################################
###############################CORRECCION DE NORMALIDAD ###########################################
###################################################################################################
###################################################################################################


###########################################################################
########CORRECCION DE NORMALIDAD POR ANAMORFOSIS GAUSSIANA#################
###########################################################################

RASselect = DatosCampo[,c("X","Y","RAS")] 
RAS.rgdb <- db.create(RASselect,ndim=2,autoname=F)
RAS.herm <- anam.fit(RAS.rgdb,name="RAS",type="gaus")
RAS.hermtrans <- anam.z2y(RAS.rgdb,names="RAS",anam=RAS.herm)
RAS.trans <- RAS.hermtrans@items$Gaussian.RAS


DatosCampo$RASTransGauss<-RAS.trans


shapiro.test(DatosCampo$RASTransGauss) #p-value 0.0247 no existe normalidad
ks.test(as.numeric(scale(sort(DatosCampo$RASTransGauss))),pnorm) #p-value 0.3081 existe normalidad


############################################################################
#######CORRECCION DE NORMALIDAD POR BOX-COX ################################
############################################################################


RASBoxCox<-boxcox(DatosCampo$RAS~1)
lambdaRAS <- RASBoxCox$x[which.max(RASBoxCox$y)]
BoxRAS<- (((DatosCampo$RAS)^lambdaRAS)-1)/lambdaRAS
shapiro.test(BoxRAS) #p-value 0.9305 > 0.05 % Ya existe normalidad
ks.test(as.numeric(scale(sort(BoxRAS))),pnorm) #p-value  0.8804 ya existe normalidad
DatosCampo$RAS.TRANS<-BoxRAS
X11()
theme_set(theme_update())
p1<-qplot(DatosCampo$RAS.TRANS,geom="histogram",main = "", xlab = "RAS Transformado Box-Cox", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p2<-qplot(DatosCampo$RAS.TRANS,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="RAS Transformado Box-Cox",ylab="Densidad")
p3<-qplot(DatosCampo$RASTransGauss,geom="histogram",main = "", xlab = "RAS Transformado Anamorfosis Gaussiana", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p4<-qplot(DatosCampo$RASTransGauss,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="RAS Transformado Anamorfosis Gaussiana",ylab="Densidad")


grid.arrange(p1,p2,p3,p4,ncol = 2)
############
#AGREGAMOS LA VARIABLE TRANSFORMADA AL ARCHIVO SP Y DATA FRAME
#POR UN LADO EL MANEJO DE MAPAS Y OTRO EL USO DEL DATA.FRAME  

####################################################
##Análisis de tendencia y elección del mejor modelo#
####################################################

ModeloGeneral<-lm(DatosCampo$RAS.TRANS~DatosCampo$X+DatosCampo$Y+DatosCampo$SDT+DatosCampo$T+DatosCampo$DUREZA+DatosCampo$CE__25C_+DatosCampo$ALCALINIDA)
summary(ModeloGeneral)
shapiro.test(ModeloGeneral$residuals)##P-value 0.0003715

Modelo1<-lm(DatosCampo$RAS.TRANS~DatosCampo$SDT+DatosCampo$T+DatosCampo$DUREZA)  
summary(Modelo1)
shapiro.test(Modelo1$residuals)##P-value 0.3282
ks.test(as.numeric(scale(sort(Modelo1$residuals))),pnorm) #p-value0.1473 existe normalidad

p1<-qplot(Modelo1$residuals,geom="histogram",main = "", xlab = "Histograma residuos modelo lineal ", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p2<-qplot(Modelo1$residuals,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="Densidad residuos modelo lineal ",ylab="Densidad")
grid.arrange(p1,p2, ncol = 2)

AICModeloG<-stepAIC(ModeloGeneral,direction = "both") #AIC=-26.52
AICModelo1<-stepAIC(Modelo1,direction = "both") #AIC=-30.85


##################################################################
######### ANÁLISIS KRIGING #######################################

######################################################################
################ANALISIS DE ANISOTROPIA###############################
######################################################################
Datos<-data.frame(DatosCampo$X,DatosCampo$Y,DatosCampo$RAS.TRANS)

geodata.anisotropia<-as.geodata(Datos,coords.col = 1:2,data.col =3)
Dmax<-sqrt((max(Datos[, 1])-min(Datos[, 1]))^2+(max(Datos[, 2])-min(Datos[, 2]))^2)
Dusar<-Dmax/2  ###Distancia a usar para modelar el variograma 29617.13
x11()
plot(variog4(geodata.anisotropia,max.dist = Dusar), xlab = "Distancia (m)", ylab = "Semivariograma estimado", legend = F)
legend(locator(1), legend = c(expression(0*degree), expression(45*degree), expression(90*degree), expression(135*degree)), col = 1:4, lty = 1:4)
title("Semivariogramas experimentales (estimador clasico)")

Coorxy<- geodata.anisotropia$coords
ZRAS<-Datos$DatosCampo.RAS.TRANS
Puntos= data.frame(Coorxy,DatosCampo$Z, z=ZRAS)
coordinates(Puntos) <- ~DatosCampo.X+DatosCampo.Y
SpatialPoints(coordinates(Puntos))
class(Puntos)
estimateAnisotropy(Puntos)

#######################################################################
#######################################################################
######################CALCULO VARIOGRAMA EXPERIMENTAL##################
#######################################################################

Base<-data.frame(geodata.anisotropia)  
Point <- point(Base, x = "DatosCampo.X", y = "DatosCampo.Y") 
Pair <- pair(Point,num.lags=20,maxdist=Dusar)
EstimaSV<- est.variograms(Point,Pair,"data",trim=0.1) 
x11() 

layout(matrix(c(1:4), nrow=2, byrow=FALSE))

layout.show(4) 
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,3.9), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$med,type="b",ylim=c(0,3.9), col =1,pch= 16, main = "Modelo experimental a partir de la mediana",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$classic,type="b",ylim=c(0,3.9), col =1,pch= 16,main = "Modelo experimental Clásico",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,3.9), col =1,pch= 16,main = "Modeloxperimental media recortada",xlab="Distancia", ylab="Semivarianza") 


######### SE SELECCIONA MEDIA RECORTADA

#####################################################
############VARIOGRAMA SIMPLE########################
#####################################################


plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,4), col =1,pch= 16,main = "",xlab="Distancia", ylab="Semivarianza") 


Exponencial<-likfit(geodata = geodata.anisotropia,nugget = 0.18,ini = c(2.2,6900),fix.nug=T)

lines(Exponencial,max.dist=Dusar,lwd=3,col="blue") 

Esferico<-likfit(geodata = geodata.anisotropia,nugget = 0.19, ini = c(2.4,8200),cov.model="sph",fix.nug=T) #Definitivo

lines(Esferico,max.dist=Dusar,lwd=3,col="red")

Matern<-likfit(geodata = geodata.anisotropia,nugget = 0.20, ini = c(2.4,9500),cov.model="mat",kappa=0.38,fix.kappa = T,fix.nug=T) 

lines(Matern,max.dist=Dusar,lwd=3,col="Black")

Circular<-likfit(geodata = geodata.anisotropia,nugget = 0.25, ini = c(2.4,9800),cov.model="cir",fix.nug=T) 

lines(Circular,max.dist=Dusar,lwd=3,col="yellow")

PowerExponencial<-likfit(geodata = geodata.anisotropia, nugget = 0.2,ini = c(2.5,5500),cov.model="powered.exponential",kappa=1.08,fix.nug=T) 

lines(PowerExponencial,max.dist=Dusar,lwd=3,col="grey") 



plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,3.9), col =1,pch= 16,main = "",xlab="Distancia", ylab="Semivarianza") 

lines(Exponencial,max.dist=Dusar,lwd=3,col="black") 
lines(Esferico,max.dist=Dusar,lwd=3,col="red") 
lines(Matern,max.dist=Dusar,lwd=3,col="yellow") 
lines(Circular,max.dist=Dusar,lwd=3,col="blue") 
lines(PowerExponencial,max.dist=Dusar,lwd=3,col="orange") 
legend(locator(1),c('Exponencial','Esferico','Matern','Circular','Potencial'),col=c("black","red","yellow","blue", "orange"), lty=c(1,1,1,1,1,1,1))


#AIC
Exponencial$AIC #657.2371
Esferico$AIC #656.8485
Matern$AIC #649.5498
Circular$AIC #653.3807
PowerExponencial$AIC #656.7995


##############################################
#########AJUSTE SEMIVARIOGRAMA###################
##############################################

Matern<-likfit(geodata = geodata.anisotropia,nugget = 0.20, ini = c(2.4,9500),cov.model="mat",kappa=0.38,fix.kappa = T,fix.nug=T) 

VariogramaRobusto<-variog(geodata.anisotropia,trend="cte",max.dist=Dusar, option = "cloud",estimator.type="modulus")


Maternwls<-variofit(vario = VariogramaRobusto,nugget = 0.9, ini = c(2.35,9800),kappa=0.20,fix.nugget= T,weights="npairs",cov.model = "matern")

plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,3.9), col =1,pch= 16,main = "Modeloxperimental media recortada",xlab="Distancia", ylab="Semivarianza") 

lines(Maternwls,max.dist=Dusar,lwd=2,lty=1,col="blue")##mejor ajuste


Maternml<-likfit(geodata = geodata.anisotropia,nugget = 0.15, ini = c(2.4,9200),kappa=0.36,fix.nugget= T,cov.model = "matern")


Maternrml<-likfit(geodata = geodata.anisotropia,nugget = 0.22, ini = c(2.2,9300),kappa=0.39,fix.nugget= T,method='RML',cov.model = "matern")

x11()##Eleccion del método de estimación por medios gráficos
plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,3.9), col =1,pch= 16,main = "",xlab="Distancia", ylab="Semivarianza") 
lines(Maternrml,max.dist=Dusar,lwd=3,lty=1, col="yellow")
lines(Maternwls,max.dist=Dusar,lwd=2,lty=1,col="blue")##mejor ajuste
lines(Maternml,max.dist=Dusar,lwd=2,lty=1,col="black")
legend(locator(1),c('ML','RML','MCP'),lty=c(1,1,1),col=c("Black","yellow","blue"))


#############################################
#############################################
#############################################
##Métodos de interpolación geoestadísticos###
#############################################
#############################################



##############################################
##### KRIGING SIMPLE #########################
z<-DatosCampo$RAS.TRANS
Media  = mean(z)                                                         
So <- c(1118001,1123900)
So <- as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So) <- c("x", "y")

TablaRAS<-data.frame(DatosCampo$X,DatosCampo$Y, DatosCampo$RAS.TRANS)
xy <- TablaRAS[,1:2]

pts = data.frame(xy, z=z) #

coordinates(pts) <- c("DatosCampo.X", "DatosCampo.Y")
proj4string(pts)<-CRS("+init=epsg:3116")


gridded(PuntosInterpolar) <- TRUE

Prediccion.ks <- krige(z~1, pts, PuntosInterpolar, vgm(2.35,"Mat",9800,0.9,kappa = 0.22,beta=Media)) ####kRIGING SIMPLE
Prediccion.ko <- krige(z~1, pts, PuntosInterpolar, vgm(2.35,"Mat",9800,0.9,kappa = 0.22)) ####kRIGING ORDINARIO

#####VALIDACION CURZADA ###################
Mater<-vgm(2.35,"Mat",9800,0.9,kappa = 0.22)

ValidacionKo<-krige.cv(z~1,pts,TablaRAS,model=Mate,maxdist=3948.951)
ValidacionFinalKo<- criterio.cv(ValidacionKo)
#          MPE    ASEPE   RMSPE         MSPE    RMSSPE    MAPPE      CCPE
#-0.007928203 1.434003 1.03184 -0.002713085 0.7135717 2.919254 0.5531964

ValidacionKS<-krige.cv(z~1,pts,TablaRAS,model=Mater,beta=Media,maxdist=3948.951)
ValidacionFinalKs<- criterio.cv(ValidacionKS)
# -0.02205275 1.433699 1.033115 -0.01236388 0.7145439 2.872352 0.5522514 0.3014815

######################################
##########MAPA DE PREDICCIONES #######
######################################

##############################
##PARA KRIGING ORDINARIO #####
##############################  
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
Prediccion.ko$Prediccion<-InversaBoxcox(Prediccion.ko$var1.pred,lambdaRAS)

summary(Prediccion.ko$Prediccion)
summary(Prediccion.ko$var1.var)
summary(DatosCampo$RAS)

pred_mean<-Prediccion.ko$Prediccion
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
ko<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()

##VARIANZA
varianza<-Prediccion.ko$var1.var
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = varianza, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
Varko<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()


##############################
##PARA KRIGING SIMPLE #####
##############################  

Prediccion.ks$Prediccion<-InversaBoxcox(Prediccion.ks$var1.pred,lambdaRAS)


pred_mean<-Prediccion.ks$Prediccion
dpm1 <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean, variable = ""
  ))
dpm1$variable <- as.factor(dpm1$variable) 
ks<-ggplot(dpm1) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()

grid.arrange(ko,ks, ncol = 2)

##VARIANZA
varianza<-Prediccion.ks$var1.var
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = varianza, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
Varks<-ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()
grid.arrange(Varko,Varks, ncol = 2)

library(SpecsVerification)

GaussCrps(mean(ValidacionKo$var1.pred),sd(ValidacionKo$var1.pred),mean(ValidacionKo$observed))
#0.1449561

GaussCrps(mean(ValidacionKS$var1.pred),sd(ValidacionKS$var1.pred),mean(ValidacionKS$observed))
#0.143

dataa<-data.frame(ValidacionKo$var1.pred,ValidacionKo$var1.var,ValidacionKo$observed,ValidacionKo$residual)
Export(dataa,"ResultadosTESISKORAS.csv") 
dataa<-data.frame(ValidacionKS$var1.pred,ValidacionKS$var1.var,ValidacionKS$observed,ValidacionKS$residual)
Export(dataa,"ResultadosTESISKsRAS.csv") 


#######################################################
###### COKRIGING ORDINARIO  ###########################
#######################################################



RASTRANS<-BoxRAS
RAS.Transformado<-TablaRAS$DatosCampo.RAS.TRANS

Dframe2<-data.frame(Tem,SDT,PH,Color,Dureza,CE,Alca,HCO3,SO4,CL,RASTRANS)
CorrelacionRAS<-cor(Dframe2,method="pearson")
round(CorrelacionRAS,digits = 2)
library(corrplot)
x11()
corrplot(CorrelacionRAS,type="upper",method="shade", tl.col = "black", tl.srt = 45,addCoef.col = "Black", diag=F, addshade = "all") 

###VARIABLES CORRELACIONADAS SON 

#TEM, SDT, CE, HCO3,SO4, CL 

#CL
CLBoxCox<-boxcox(DatosCampo$CL~1)
lambdaCL<- CLBoxCox$x[which.max(CLBoxCox$y)]
BoxCL<-(((Iones$Clmeq)^lambdaCL)-1 )/lambdaCL
shapiro.test(BoxCL)# p-value = 0.01155
ks.test(as.numeric(scale(sort(BoxCL))),pnorm) 
hist(BoxCL)
#pasa
#SDT
SDTBoxCox<-boxcox(DatosCampo$SDT~1)
lambdaSDT<- SDTBoxCox$x[which.max(SDTBoxCox$y)]
BoxSDT<-(((DatosCampo$SDT)^lambdaSDT)-1 )/lambdaSDT
shapiro.test(BoxSDT)# p-value = 0.2692
ks.test(as.numeric(scale(sort(BoxSDT))),pnorm) #p-value  0.6238 existe normalidad
#pasa
#DUREZA
DureBoxCox<-boxcox(DatosCampo$DUREZA~1)
lambdaDureza<- DureBoxCox$x[which.max(DureBoxCox$y)]
BoxDureza<-(((DatosCampo$DUREZA)^lambdaDureza)-1 )/lambdaDureza
shapiro.test(BoxDureza)# p-value = 0.09604 
ks.test(as.numeric(scale(sort(BoxDureza))),pnorm) #p-value  0.674 existe normalidad

#CE
CEBoxCox<-boxcox(DatosCampo$CE__25C_~1)
lambdaCE<- CEBoxCox$x[which.max(CEBoxCox$y)]
BoxCE<-(((DatosCampo$CE__25C_)^lambdaCE)-1 )/lambdaCE
shapiro.test(BoxCE)# p-value = 0.4138 
ks.test(as.numeric(scale(sort(BoxCE))),pnorm) #p-va existe normalidad


#HC03
HCO3BoxCox<-boxcox(DatosCampo$HCO3~1)
lambdaHCO3<- HCO3BoxCox$x[which.max(HCO3BoxCox$y)]
BoxHCO3<-(((DatosCampo$HCO3)^lambdaHCO3)-1 )/lambdaHCO3
shapiro.test(BoxHCO3)# p-value = 0.2429
ks.test(as.numeric(scale(sort(BoxHCO3))),pnorm) #p-va existe normalidad
#Tem
TemBoxCox<-boxcox(DatosCampo$T~1)
lambdaTem<- TemBoxCox$x[which.max(TemBoxCox$y)]
BoxTem<-(((DatosCampo$T)^lambdaTem)-1 )/lambdaTem
shapiro.test(BoxTem)# p-value = 0.2429
ks.test(as.numeric(scale(sort(BoxHCO3))),pnorm) #p-va existe normalidad

#SO4

####se estandariza ya que tiene valores negativos ...
Estandarizacion <-function(x) {(x-mean(x))/sd(x)}

SO4Estan<-Estandarizacion(DatosCampo$SO4)




Coorx<-DatosCampo$X
Coory<-DatosCampo$Y
Dframe.Trans<-data.frame(Coorx,Coory,BoxSDT,BoxDureza,BoxCE,BoxHCO3,BoxCL,BoxTem,BoxRAS,SO4Estan)
coordinates(Dframe.Trans)<-~Coorx+Coory
proj4string(Dframe.Trans)<-CRS("+init=epsg:3116")

CoKo=gstat(id = "RAS", formula = BoxRAS~1, data = Dframe.Trans)
CoKo=gstat(CoKo,id = "SDT", formula = BoxSDT~1, data = Dframe.Trans)
#CoKo=gstat(CoKo,id = "Dure", formula = BoxDureza~1, data = Dframe.Trans)
#CoKo=gstat(CoKo,id = "CE", formula = BoxCE~1, data = Dframe.Trans)
#CoKo=gstat(CoKo,id = "HCO3", formula = BoxHCO3~1, data = Dframe.Trans)
CoKo=gstat(CoKo,id = "CL", formula = BoxCL~1, data = Dframe.Trans)
#CoKo=gstat(CoKo,id = "Tem", formula = BoxTem~1, data = Dframe.Trans)
#CoKo=gstat(CoKo,id = "SO4", formula = SO4Estan~1, data = Dframe.Trans)

CoKo=gstat(CoKo,id = "RAS", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))  
CoKo=gstat(CoKo,id = "SDT", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))
#CoKo=gstat(CoKo,id = "Dure", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))
#CoKo=gstat(CoKo,id = "CE", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))
#CoKo=gstat(CoKo,id = "HCO3", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))
CoKo=gstat(CoKo,id = "CL", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))
#CoKo=gstat(CoKo,id = "Tem", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))
#CoKo=gstat(CoKo,id = "SO4", model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22))

Vcros<-variogram(CoKo,Dusar)

CoKo<-gstat(CoKo,id="RAS",model=vgm(2.35,"Mat",9800,0.9,kappa = 0.22), fill.all=T)


###SE REALIZA EL AJUSTE DE LOS MODELOS CORRELACIONADOS
AjusteCokriging<-fit.lmc(Vcros,CoKo)

X11()
plot(variogram(CoKo), model=AjusteCokriging$model,pch=20)


#####################################################
###Validación cruzcada de CoKriging###

hist(DatosCampo$RAS)

ValidacionCoKo=gstat.cv(AjusteCokriging,maxdist=3948.951)
criterio.cv(ValidacionCoKo)
summary(ValidacionCoKo$observed)
dataa<-data.frame(ValidacionCoKo$RAS.pred,ValidacionCoKo$RAS.var,ValidacionCoKo$observed,ValidacionCoKo$residual)
Export(dataa,"ResultadosTESISCokrigingRAS1.csv")

PrediccionCok<-predict(AjusteCokriging,PuntosInterpolar)
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
PrediccionCok$RAS.pred<-InversaBoxcox(PrediccionCok$RAS.pred,lambdaRAS)

summary(InversaBoxcox(PrediccionCok$RAS.pred,lambdaRAS))
summary(ValidacionCoKo$observed)
summary(PrediccionCok$RAS.pred)
mean(GaussCrps(mean(ValidacionCoKo$RAS.pred),sd(ValidacionCoKo$RAS.pred),mean(ValidacionCoKo$observed)))
#0.2120644



pred_mean<-PrediccionCok$RAS.pred
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean, variable = ""
  ))
dpm$variable <- as.factor(dpm$variable) 
ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()




