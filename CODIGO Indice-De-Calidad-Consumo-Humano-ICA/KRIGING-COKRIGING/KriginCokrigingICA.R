#########################################################################################
#######"MODELADO ESPACIAL HIDROGEOLÓGICO PARA DETERMINAR ÍNDICES DE CALIDAD##############
#######Y VULNERABILIDAD DE LAS AGUAS SUBTERRÁNEAS EN LA ZONA CENTRO DE BOYACÁ"###########

#################################################################
####CÁLCULO DE INDICE DE CALIDAD PARA CONSUMO HUMANO-KRIGING ####
#################################################################
####kLIBRERIAS A UTILIZAR


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

##########################################################################
######################LECTURA DE LA BASE DE DATOS#########################
##########################################################################


###SHAPE DE MUNICIPIOS DE LA ZONA ##
MunicipiosZonaEstudio<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO_B/MunicipiosZonaEstudio.shp")
#BASE DE DATOS LUEGO DEL BALANCE IÓNICO ##
DatosCampo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/MuestroBalanceIonico.shp")
##SHAPE DE USOS DEL SUELO ZONA CENTRO ##
UsoDelSuelo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/BOYACA_CAPACIDAD_VF/UsoSueloZonaEstudio.shp")
#PERÍMETRO DE LA ZONA
ShapeZona<-readOGR("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO/unidad_crono_Boyaca_zona_centro.shp")



##SISTEMA DE REFERENCIA COORDENDAS DE LOS DATOS DE CAMPO ##
Coordenadas<-read.table("C:/PROYECTODEINVESTIGACION/SHAPEFILE/Depuracion.csv",header=TRUE,sep=",")


epsg3116<-"+init=epsg:3116"
ShapeZona<-spTransform(ShapeZona,CRS(epsg3116))

##GRILLA 
PuntosInterpolar<- spsample(ShapeZona,n=70000,type="regular")

#########################################################################
#### GRÁFICO DE LOS PARÁMETROS DE INTERES PARA EL ICA ###################
#########################################################################

##CONDUCTIVIDAD ELECTRICA

ggplot(data = DatosCampo) + geom_sf(aes(size=CE__25C_), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud")+ ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

##CALCIO
ggplot(data = DatosCampo) + geom_sf(aes(size=CA), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

##MAGNESIO
ggplot(data = DatosCampo) + geom_sf(aes(size=MG), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

##SODIO

ggplot(data = DatosCampo) + geom_sf(aes(size=NA.), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()


##ALCALINIDA
ggplot(data = DatosCampo) + geom_sf(aes(size=ALCALINIDA), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

#NITRATOS
ggplot(data = DatosCampo) + geom_sf(aes(size=NO3), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

##SULFATO
ggplot(data = DatosCampo) + geom_sf(aes(size=SO4), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

#FOSFATO
ggplot(data = DatosCampo) + geom_sf(aes(size=PO4), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

##PH

ggplot(data = DatosCampo) + geom_sf(aes(size=PH), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

#TURBIEDAD
ggplot(data = DatosCampo) + geom_sf(aes(size=TURBIEDAD), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

#MANGANESO
ggplot(data = DatosCampo) + geom_sf(aes(size=MN), color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()

###############
##HISTOGRAMAS 
###############

p1<-qplot(DatosCampo$CE__25C_,geom="histogram",main = "Histograma de la CE (??S/cm) ", xlab = "Conductividad eléctrica", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p2<-qplot(DatosCampo$CA,geom="histogram",main = "Histograma de  calcio", xlab = "Calcio Ca", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p3<-qplot(DatosCampo$MG,geom="histogram",main = "Histograma de  Magnesio", xlab = "Magnesio Mg", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p4<-qplot(DatosCampo$NA.,geom="histogram",main = "Histograma de Sodio", xlab = "Sodio", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))

grid.arrange(p1,p2,p3,p4,nrow = 2, widths=c(2,2))

p1<-qplot(DatosCampo$NO3,geom="histogram",main = "Histograma de Nitratos ", xlab = "Nitratos NO3", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p2<-qplot(DatosCampo$PO4,geom="histogram",main = "Histograma de  Fosfato", xlab = "Fosfato PO4", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p3<-qplot(DatosCampo$SO4,geom="histogram",main = "Histograma de  Sulfato", xlab = "Sulfato SO4", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p4<-qplot(DatosCampo$ALCALINIDA,geom="histogram",main = "Histograma de Alcalinida", xlab = "Alcalinida", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))

grid.arrange(p1,p2,p3,p4,nrow = 2, widths=c(2,2))


p1<-qplot(DatosCampo$PH,geom="histogram",main = "Histograma de PH ", xlab = "Potencial de hidrógeno PH", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p2<-qplot(DatosCampo$TURBIEDAD,geom="histogram",main = "Histograma de  Turbiedad", xlab = "Turbiedad", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))
p3<-qplot(DatosCampo$MN,geom="histogram",main = "Histograma de  Manganeso", xlab = "Manganeso Mn", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))

grid.arrange(p1,p2, ncol = 2, p3)

qplot(Dframe$CE,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="WQI",ylab="Densidad")
qqnorm(TablaSi$WQI, ylab="WQI", xlab="Cuantiles teóricos",main="",pch=20,col=rgb(0.3,0.5,1,0.4))


###################################################################################
##############CALCULO DEL INDICE DE CALIDAD IWQ####################################
###################################################################################

###VARIABLES A UTILIZAR PARA EL IWQ
##Parámetros tenidos en cuenta SE CREA UN DATA.FRAME
BDIWQ<-data.frame(DatosCampo$PH,DatosCampo$CE__25C_,DatosCampo$TURBIEDAD,DatosCampo$ALCALINIDA,DatosCampo$PO4,DatosCampo$NO3,DatosCampo$SO4,DatosCampo$CA,DatosCampo$MG,DatosCampo$MN,DatosCampo$NA.) ###n=11 parámetros
summary(BDIWQ)
##pesos de cada parámetro #SEGUN EL CALCULO DEL ÍNDICE...
Parametros<-c("PH","Conductividad Electrica","Turbiedad","Alcalinidad","Fosfato","Nitratos","Sulfato","Calcio","Magensio","Manganeso","Sodio")
Pesos<-c(1,5,4,4,1,5,4,4,3,3,4) #VECTOR DE PESOS
SumPesos<-sum(Pesos)

Wi<-c(Pesos/SumPesos) 
SumaWi=sum(Wi)##VALOR FINAL DEL ÍNDICE

#####TABLA DE PESOS RELATIVOS

TpesosRelativos<-data.frame (cbind(Parametros,Pesos,Wi))

###Calculo qi= 100 [(Vi - Vid) /(Si - Vid), 
##donde Vi= valor medido, Vid= valor ideal del parámetros, Si=Valor estandar permitido
##Esto solo se usa para el PH, en los demás casos Vid=0. Para el PH Vid=7

BDIWQ$qiPH<-(((DatosCampo$PH-7)/(8.5-7))*100)
BDIWQ$qiCE<-(DatosCampo$CE__25C_/500)*100
BDIWQ$qiTurbi<-(DatosCampo$TURBIEDAD/2)*100
BDIWQ$qiAlca<-(DatosCampo$ALCALINIDA/200)*100
BDIWQ$qiFosfato<-(DatosCampo$PO4/0.5)*100
BDIWQ$qiNitratos<-(DatosCampo$NO3/10)*100
BDIWQ$qiSulfato<-(DatosCampo$SO4/250)*100
BDIWQ$qiCa<-(DatosCampo$CA/60)*100
BDIWQ$qiMg<-(DatosCampo$MG/36)*100
BDIWQ$qiMN<-(DatosCampo$MN/0.1)*100
BDIWQ$qiSodio<-(DatosCampo$NA./200)*100

#######Calculo de Si=Wi*qi####################

SiPH<-BDIWQ$qiPH*Wi[1]
SiCE<-BDIWQ$qiCE*Wi[2]
SiTurbi<-BDIWQ$qiTurbi*Wi[3]
SiAlca<-BDIWQ$qiAlca*Wi[4]
SiFosfato<-BDIWQ$qiFosfato*Wi[5]
SiNitratos<-BDIWQ$qiNitratos*Wi[6]
SiSulfato<-BDIWQ$qiSulfato*Wi[7]
SiCalcio<-BDIWQ$qiCa*Wi[8]
SiMagnesio<-BDIWQ$qiMg*Wi[9]
SiManganeso<-BDIWQ$qiMN*Wi[10]
SiSodio<-BDIWQ$qiSodio*Wi[11]

##################################################
###### ÍNDICE DE CALIDAD WQI #####################

TablaWQI<-data.frame (cbind(DatosCampo$X,DatosCampo$Y,SiPH,SiCE,SiTurbi,SiAlca,SiFosfato,SiNitratos,SiSulfato,SiCalcio,SiMagnesio,SiManganeso,SiSodio))

TablaWQI$WQI<-rowSums(TablaWQI[,3:13])   
DatosCampo$WQI<-rowSums(TablaWQI[,3:13]) 

summary(TablaWQI$WQI) 

###############################################################
###MATRIZ DE CORRELACION PARÁMETRIS DEL IWQ ###################
###############################################################
MatrizCorrelacion<-data.frame(DatosCampo$PH,DatosCampo$CE__25C_,DatosCampo$TURBIEDAD,DatosCampo$ALCALINIDA,DatosCampo$PO4,DatosCampo$NO3,DatosCampo$SO4,DatosCampo$CA,DatosCampo$MG,DatosCampo$MN,DatosCampo$NA.)

MatrizCorrelacion$PH<-MatrizCorrelacion$DatosCampo.PH
MatrizCorrelacion$CE<-MatrizCorrelacion$DatosCampo.CE__25C_
MatrizCorrelacion$Turbieda<-MatrizCorrelacion$DatosCampo.TURBIEDAD
MatrizCorrelacion$Alcalinidad<-MatrizCorrelacion$DatosCampo.ALCALINIDA
MatrizCorrelacion$PO4<-MatrizCorrelacion$DatosCampo.PO4
MatrizCorrelacion$NO3<-MatrizCorrelacion$DatosCampo.NO3
MatrizCorrelacion$SO4<-MatrizCorrelacion$DatosCampo.SO4
MatrizCorrelacion$Ca<-MatrizCorrelacion$DatosCampo.CA
MatrizCorrelacion$Mg<-MatrizCorrelacion$DatosCampo.MG
MatrizCorrelacion$Mn<-MatrizCorrelacion$DatosCampo.MN
MatrizCorrelacion$Na<-MatrizCorrelacion$DatosCampo.NA.

CorrelacionParametros<-cor(MatrizCorrelacion[,12:22],method="pearson")

round(CorrelacionParametros,digits = 2)

corrplot(CorrelacionParametros,type="upper",method="shade", tl.col = "black", tl.srt = 45,addCoef.col = "Black", diag=F, addshade = "all") 

###########################################################
##RESUMEN ESTADISTICO DEL ÍNDICE DE CALIDAD DEL AGUA ######
###########################################################
summary(DatosCampo$WQI) 

################################################################
#########GRÁFICOS DEL INDICE DE CALIDAD DEL AGUA ###############
################################################################

qplot(DatosCampo$WQI,geom="histogram",main = "", xlab = "ICA", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))+scale_x_continuous(limits = c(0, 1200),breaks = c(50,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,12200))+annotate(geom = "text", x = 700, y = 80, label = "Rango", color="firebrick")+annotate(geom = "text", x = 900, y = 80, label = "Tipo de agua", color="firebrick")+annotate(geom = "text", x = 700, y = 70, label = "<50")+annotate(geom = "text", x = 700, y = 65, label = "50-100")+annotate(geom = "text", x = 700, y = 60, label = "100-200")+annotate(geom = "text", x = 700, y = 55, label = "200-300")+annotate(geom = "text", x = 700, y = 50, label = ">300")+annotate(geom = "text", x = 895, y = 70, label = "Agua excelente")+annotate(geom = "text", x = 870, y = 65, label = "Agua buena")+annotate(geom = "text", x = 868, y = 60, label = "Agua pobre")+annotate(geom = "text", x = 900, y = 55, label = "Agua muy pobre")+annotate(geom = "text", x = 960, y = 50, label = "Agua no apta para beber")

qplot(DatosCampo$WQI,geom="histogram",main = "", xlab = "ICA", ylab = "Frecuecnia",fill=I("Red"),col=I("black"), alpha=I(.5))+scale_x_continuous(limits = c(1200,12200),breaks=c(1200,2200,3200,4200,5200,6200,7200,8200,9200,10200,11200,12100))+annotate(geom = "text", x = 8200, y = 80, label = "Rango", color="firebrick")+annotate(geom = "text", x = 9800, y = 80, label = "Tipo de agua", color="firebrick")+annotate(geom = "text", x = 8200, y = 70, label = "<50")+annotate(geom = "text", x = 8200, y = 65, label = "50-100")+annotate(geom = "text", x = 8200, y = 60, label = "100-200")+annotate(geom = "text", x = 8200, y = 55, label = "200-300")+annotate(geom = "text", x =8200, y = 50, label = ">300")+annotate(geom = "text", x = 10000, y = 70, label = "Agua excelente")+annotate(geom = "text", x = 9800, y = 65, label = "Agua buena")+annotate(geom = "text", x = 9780, y = 60, label = "Agua pobre")+annotate(geom = "text", x = 10100, y = 55, label = "Agua muy pobre")+annotate(geom = "text", x = 10600, y = 50, label = "Agua no apta para beber")

DatosCampo$ICA<-DatosCampo$WQI
p<-ggplot(data = DatosCampo) + geom_sf(aes(size=ICA),color="red", alpha=0.7) +
  geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("")+ annotation_north_arrow(location='tr')+annotation_scale()
p+scale_size(breaks =c(50,200,400,600,800,1000,1500,2000,2500,3000,3500,11500))

##########################################################################
##########################################################################
###### ANALISIS DE USO ACTUAL DEL SUELO ##################################
##########################################################################
##########################################################################

p<-ggplot(data =UsoDelSuelo) + geom_sf(aes( fill=USO_ACTUAL))+geom_sf(data =DatosCampo, aes(size=WQI), color="Black", alpha=0.5)+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("Uso del suelo zona centro de Boyacá vs ICA")+ annotation_north_arrow(location='tr')+annotation_scale()
p+scale_size(breaks =c(50,200,400,600,800,1000,1500,2000,2500,3000,3500,11500))+guides(size= guide_legend(nrow = 2))

ggplot(data =UsoDelSuelo) + geom_sf(aes( fill=USO_ACTUAL))+xlab("Longitud")+ ylab("Latitud") +
  ggtitle("Uso del suelo zona centro de Boyacá")+ annotation_north_arrow(location='tr')+annotation_scale()

theme(legend.position = "top")


#################################################
#################################################
##### ANÁLISIS DESCRIPTIVO ICA ##################
#################################################

##ANALSIS DE NORMALIDAD
##METODO GRÁFICO

X11()
qplot(DatosCampo$WQI,geom = 'blank', xlab = "ICA") +   
  geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
  theme(legend.position = c(0.5, 0.5)) 
x11()
theme_set(theme_update())
qplot(DatosCampo$WQI,geom="histogram",main = "", xlab = "Temperatura", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
qplot(DatosCampo$WQI,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="WQI",ylab="Densidad")
qqnorm(DatosCampo$WQI, ylab="WQI", xlab="Cuantiles teóricos",main="",pch=20,col=rgb(0.3,0.5,1,0.4))
qqline(DatosCampo$WQI,lwd=1.9)
boxplot(DatosCampo$WQI, data=TablaSi,col=rgb(0.3,0.5,1,0.4),alpha=0.1,main=" ",horizontal = T, xlab="WQI",pch=20)

x11()

shapiro.test(DatosCampo$WQI) #p-value < 2.2e-16 no existe normalidad
ks.test(as.numeric(scale(sort(DatosCampo$WQI))),pnorm) #p-value < 2.2e-16 no existe normalidad

###################################################################################################
###################################################################################################
###############################CORRECCION DE NORMALIDAD ###########################################
###################################################################################################
###################################################################################################


###########################################################################
########CORRECCION DE NORMALIDAD POR ANAMORFOSIS GAUSSIANA#################
###########################################################################

WQIselect = DatosCampo[,c("X","Y","WQI")] 
WQI.rgdb <- db.create(WQIselect,ndim=2,autoname=F)
WQI.herm <- anam.fit(WQI.rgdb,name="WQI",type="gaus")
WQI.hermtrans <- anam.z2y(WQI.rgdb,names="WQI",anam=WQI.herm)
WQI.trans <- WQI.hermtrans@items$Gaussian.WQI

mean(WQI.trans) 

DatosCampo$WQITransGauss<-WQI.trans

shapiro.test(DatosCampo$WQITransGauss) #p-value < 7.34e-16 no existe normalidad
ks.test(as.numeric(scale(sort(DatosCampo$WQITransGauss))),pnorm) #p-value < 2.656e-08 no existe normalidad



############################################################################
#######CORRECCION DE NORMALIDAD POR BOX-COX ################################
############################################################################


WQIBoxCox<-boxcox(DatosCampo$WQI~1)
lambdaWQI <- WQIBoxCox$x[which.max(WQIBoxCox$y)]
BoxWQI<- (((DatosCampo$WQI)^lambdaWQI)-1)/lambdaWQI 
shapiro.test(BoxWQI) #p-value 0.9573  > 0.05 % Ya existe normalidad
ks.test(as.numeric(scale(sort(BoxWQI))),pnorm) #p-value  0.9891 ya existe normalidad

############
#AGREGAMOS LA VARIABLE TRANSFORMADA AL ARCHIVO SP Y DATA FRAME
#POR UN LADO EL MANEJO DE MAPAS Y OTRO EL USO DEL DATA.FRAME 

DatosCampo$WQI.TRANS<-BoxWQI
TablaWQI$WQI.TRANS<-BoxWQI  ###ASIGNAMOS AL DATA FRAME 

X11()
theme_set(theme_update())
p1<-qplot(DatosCampo$WQI.TRANS,geom="histogram",main = "", xlab = "ICA Transformado Box-Cox", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p2<-qplot(DatosCampo$WQI.TRANS,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="ICA Transformado Box-Cox",ylab="Densidad")
p3<-qplot(DatosCampo$WQITransGauss,geom="histogram",main = "", xlab = "ICA Transformado Anamorfosis Gaussiana", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p4<-qplot(DatosCampo$WQITransGauss,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="ICA Transformado Anamorfosis Gaussiana",ylab="Densidad")


grid.arrange(p1,p2,p3,p4,ncol = 2)
######NOTA###
#CON BOX COX YA NORMALIZA LA VARIABLE WQI 

####################################################
##Análisis de tendencia y elección del mejor modelo#
####################################################

ModeloGeneral<-lm(DatosCampo$WQI.TRANS~DatosCampo$X+DatosCampo$Y+DatosCampo$Z+DatosCampo$T+DatosCampo$ST+DatosCampo$COLOR+DatosCampo$DUREZA)
summary(ModeloGeneral)
shapiro.test(ModeloGeneral$residuals)##P-value 0.1457 

Modelo1<-lm(DatosCampo$WQI.TRANS~DatosCampo$Y+DatosCampo$T+DatosCampo$COLOR+DatosCampo$DUREZA)  
summary(Modelo1)
shapiro.test(Modelo1$residuals)##P-value 0.3282
ks.test(as.numeric(scale(sort(Modelo1$residuals))),pnorm) #p-value0.1473 existe normalidad



AICModeloG<-stepAIC(ModeloGeneral,direction = "both") #
AICModelo1<-stepAIC(Modelo1,direction = "both") #

p1<-qplot(Modelo1$residuals,geom="histogram",main = "", xlab = "Histograma residuos modelo lineal ", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))
p2<-qplot(DatosCampo$WQI.TRANS,geom="density",col=I("blue"),fill=I("blue"), alpha=I(0.2), main="", xlab="Densidad residuos modelo lineal ",ylab="Densidad")
grid.arrange(p1,p2, ncol = 2)

#########NOTA
#########NINGUN MODELO PRESENTA TENDENCIA, POR LO TANTO SE TRABAJA CON K.O Y K.S


######################################################################
################ANALISIS DE ANISOTROPIA###############################
######################################################################

geodata.anisotropia<-as.geodata(TablaWQI,coords.col = 1:2,data.col =15)
Dmax<-sqrt((max(Coordenadas[,9])-min(Coordenadas[,9]))^2+(max(Coordenadas[,10])-min(Coordenadas[,10]))^2)
Dusar<-Dmax/2###Distancia a usar para modelar el variograma 29617.13
x11()
plot(variog4(geodata.anisotropia,max.dist = Dusar), xlab = "Distancia (m)", ylab = "Semivariograma estimado", legend = F)
legend(locator(1), legend = c(expression(0*degree), expression(45*degree), expression(90*degree), expression(135*degree)), col = 1:4, lty = 1:4)
title("Semivariogramas experimentales (estimador clasico)")
##métodos "estadisticos" para la anisotropia
Coorxy<- geodata.anisotropia$coords
ZICA<-DatosCampo$WQI.TRANS
Puntos= data.frame(Coorxy,DatosCampo$Z, z=ZICA)
coordinates(Puntos) <- ~V1+V2
SpatialPoints(coordinates(Puntos))
class(Puntos)
estimateAnisotropy(Puntos)

####NO EXISTE ANISOTROPÍA


#######################################################################
#######################################################################
######################CALCULO VARIOGRAMA EXPERIMENTAL##################
#######################################################################

Base<-data.frame(geodata.anisotropia)  
Point <- point(Base, x = "V1", y = "V2") 
Pair <- pair(Point,num.lags=20,maxdist=Dusar)
EstimaSV<- est.variograms(Point,Pair,"data",trim=0.1) 
x11() 

layout(matrix(c(1:4), nrow=2, byrow=FALSE))

layout.show(4) 
x11()
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,0.6), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$med,type="b",ylim=c(0,0.6), col =1,pch= 16, main = "Modelo experimental a partir de la mediana",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$classic,type="b",ylim=c(0,0.6), col =1,pch= 16,main = "Modelo experimental Clásico",xlab="Distancia", ylab="Semivarianza") 

x11()
plot(EstimaSV$bins,EstimaSV$trimmed.mean,type="b",ylim=c(0,0.6), col =1,pch= 16,main = "Modelo experimental media recortada",xlab="Distancia", ylab="Semivarianza") 


###NOTA
####ROBUSTO Y CLASICO LOS MEJORES, SIN EMBARGO SE SELECCIONA ROBUSTO 
####

#####################################################
############VARIOGRAMA SIMPLE########################
#####################################################

x11()
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,0.6), col =1,pch= 16, main = "Modelo perimental Robusto", xlab="Distancia", ylab="Semivarianza") 

###MODELOS EXPERIMENTALES 
Exponencial<-likfit(geodata = geodata.anisotropia,nugget = 0.055,ini = c(0.35,6500),fix.nug=T)

lines(Exponencial,max.dist=Dusar,lwd=3,col="blue") 

Esferico<-likfit(geodata = geodata.anisotropia,nugget = 0.05, ini = c(0.35,8200),cov.model="sph",fix.nug=T) #Definitivo

lines(Esferico,max.dist=Dusar,lwd=3,col="red") 

plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,0.6), col =1,pch= 16, main = "Modelo perimental Robusto", xlab="Distancia", ylab="Semivarianza") 

Matern<-likfit(geodata = geodata.anisotropia,nugget = 0.09, ini = c(0.35,9000),cov.model="mat",kappa=0.8,fix.kappa = T,fix.nug=T) 

lines(Matern,max.dist=Dusar,lwd=3,col="Black")

Circular<-likfit(geodata = geodata.anisotropia,nugget = 0.07, ini = c(0.36,9800),cov.model="cir",fix.nug=T) 

lines(Circular,max.dist=Dusar,lwd=3,col="yellow")

PowerExponencial<-likfit(geodata = geodata.anisotropia, nugget = 0.055,ini = c(0.30,5300),cov.model="powered.exponential",kappa=1.08,fix.nug=T) 

lines(PowerExponencial,max.dist=Dusar,lwd=3,col="grey") 

plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,0.6), col =1,pch= 16, main = "", xlab="Distancia", ylab="Semivarianza") 


lines(Exponencial,max.dist=Dusar,lwd=3,col="black") 
lines(Esferico,max.dist=Dusar,lwd=3,col="red") 
lines(Matern,max.dist=Dusar,lwd=3,col="yellow") 
lines(Circular,max.dist=Dusar,lwd=3,col="blue") 
lines(PowerExponencial,max.dist=Dusar,lwd=3,col="orange") 
legend(locator(1),c('Exponencial','Esferico','Matern','Circular','Potencial'),col=c("black","red","yellow","blue", "orange"), lty=c(1,1,1,1,1,1,1))

#AIC
Exponencial$AIC #300.6184
Esferico$AIC #306.3185
Matern$AIC #296.2527
Circular$AIC #297.9547
PowerExponencial$AIC #303.0985

###SE ELIGE EL MATERN POR ESTADÍSTICO AIC

##############################################
#########AJUSTE SEMIVARIOGRAMA###################
##############################################

##CON MATHERN ############

Matern<-likfit(geodata = geodata.anisotropia,nugget = 0.09, ini = c(0.35,9000),cov.model="mat",kappa=0.8,fix.kappa = T,fix.nug=T) 
VariogramaRobusto<-variog(geodata.anisotropia,trend="cte",max.dist=Dusar, option = "cloud",estimator.type="modulus")
Maternwls<-variofit(vario = VariogramaRobusto, nugget = 0.24,ini = c(0.30, 8800),kappa = 0.33,fix.nugget= T,weights="npairs",cov.model = "matern")
Maternml<-likfit(geodata = geodata.anisotropia, nugget = 0.10,ini = c(0.28, 7800),kappa = 1.2,fix.nugget= T,cov.model = "matern")
Maternrml<-likfit(geodata = geodata.anisotropia,nugget = 0.10,ini = c(0.25, 9200),kappa = 0.9,fix.nugget= T,method='RML',cov.model = "matern")

x11()##Eleccion del método de estimación por medios gráficos
plot(EstimaSV$bins,EstimaSV$robust,type="b",ylim=c(0,0.6), col =1,pch= 16, main = "", xlab="Distancia", ylab="Semivarianza") 

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
z<-DatosCampo$WQI.TRANS
Media  = mean(z)                                                         
So <- c(1118001,1123900)
So <- as.data.frame(cbind(x=So[1],y=So[2]))
coordinates(So) <- c("x", "y")

xy <- TablaWQI[,1:2]

pts = data.frame(xy, z=z) ##COORDENADAS WQI

coordinates(pts) <- c("V1", "V2")
proj4string(pts)<-CRS("+init=epsg:3116")

gridded(PuntosInterpolar) <- T


MunicipiosZonaEstudioTrans<-st_transform(MunicipiosZonaEstudio,crs=st_crs(3116))

Prediccion.ks <- krige(z~1, pts, PuntosInterpolar, vgm(0.30,"Mat",8800,0.24,kappa = 0.33,beta=Media)) ####kRIGING SIMPLE
Prediccion.ko <- krige(z~1, pts, PuntosInterpolar, vgm(0.30,"Mat",8800,0.24,kappa = 0.33)) ####kRIGING ORDINARIO



#############################################################
#####REESCALAR LOS DATOS DE LAS PREDICCIONES ################
###   ANTI-TRANSFORMACION  DE LAS PREDICCIONES  ############
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
##############################
##PARA KRIGING ORDINARIO #####
##############################  

Prediccion.ko$Prediccion<-InversaBoxcox(Prediccion.ko$var1.pred,lambdaWQI)

summary(Prediccion.ko$Prediccion)

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
    name = "ICA",
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
      name = "ICA",
      low = "blue", high = "orange"
    ) +  
    ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
  theme_bw()
  ##############################
  ##PARA KRIGING SIMPLE #####
  ##############################  
  
  Prediccion.ks$Prediccion<-InversaBoxcox(Prediccion.ks$var1.pred,lambdaWQI)
  summary(Prediccion.ko$var1.var)
  
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
      name = "ICA",
      low = "blue", high = "orange"
    ) +  
    ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
  theme_bw()
  
  grid.arrange(ko,ks, ncol = 2)


  p<-ggplot(data = dpm)+
    geom_sf(data = MunicipiosZonaEstudio, fill = "Brown", alpha=0.2,  color = "black")++xlab("Longitud")+ ylab("Latitud") +
    ggtitle("Distribución del índice de calidad de agua subterránea")+ annotation_north_arrow(location='tr')+annotation_scale()
  p+scale_size(breaks =c(50,200,400,600,800,1000,1500,2000,2500,3000,3500,11500))

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
      name = "ICA",
      low = "blue", high = "orange"
    ) +  
    ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
  theme_bw()
  grid.arrange(Varks,Varko, ncol = 2)
  
######################
##VALIDACIÓN CRUZADA##
######################

Mater<-vgm(0.30,"Mat",8800,0.24,kappa = 0.33)
  
rad <- min(3029.715, max(dist(xy)) /4)
ValidacionKo<-krige.cv(z~1,pts,TablaWQI,model=Mater, maxdist=3948.951)
ValidacionKs<-krige.cv(z~1,pts,TablaWQI,model=Mater,maxdist=3948.951)

##VALIDACION CRUZADA DEJANDO UNO AFUERA 
ValidacionFinalKo<- criterio.cv(ValidacionKo)
#            MPE     ASEPE   RMSPE          MSPE    RMSSPE     MAPPE      CCPE
#1 -0.0002585684 0.6001908 0.43992 -0.0002256916 0.7308515 0.1377737 0.3030553
#-0.003977452 0.60731 0.4455347 -0.00686346 0.7305851 0.1395784 0.2963595
ValidacionFinalKs<- criterio.cv(ValidacionKs)
#           MPE     ASEPE     RMSPE         MSPE    RMSSPE     MAPPE     CCPE
#1 -0.003727993 0.6000559 0.4398119 -0.005895226 0.7308837 0.1380159 0.303031
#-0.003977452 0.60731 0.4455347 -0.00686346 0.7305851 0.1395784 0.2963595
library(SpecsVerification)

GaussCrps(mean(ValidacionKo$var1.pred),sd(ValidacionKo$var1.pred),mean(ValidacionKo$observed))
#0.04359991
#0.0505

GaussCrps(mean(ValidacionKs$var1.pred),sd(ValidacionKs$var1.pred),mean(ValidacionKs$observed))
#0.0505

dataa<-data.frame(ValidacionKo$var1.pred,ValidacionKo$var1.var,ValidacionKo$observed,ValidacionKo$residual)
Export(dataa,"ResultadosICATESISKO.csv") 
dataa<-data.frame(ValidacionKs$var1.pred,ValidacionKs$var1.var,ValidacionKs$observed,ValidacionKs$residual)
Export(dataa,"ResultadosICATESISKs.csv") 

##########################################
#### PREDICCION COKRIGING ################

################################################
#####MATRIZ DE CORRELACION #####################

##CORRELACION CON IWQ 

library(corrplot)

SDT<-DatosCampo$SDT
Tem<-DatosCampo$T
Color<-DatosCampo$COLOR
Dureza<-DatosCampo$DUREZA
HCO3<-DatosCampo$HCO3
ICATransformado<-DatosCampo$WQI.TRANS
Dframe2<-data.frame(SDT,Tem,Color,Dureza,HCO3,ICATransformado)
CorrelacionWQI<-cor(Dframe2,method="pearson")
round(CorrelacionWQI,digits = 2)
X11()
corrplot(CorrelacionWQI,type="upper",method="shade", tl.col = "black", tl.srt = 45,addCoef.col = "Black", diag=F, addshade = "all") 


###VARIABLES A TENER EN CUENTA SON , SDT,  DUREZA, COLOR, CL Y HCO3


####NORMALIZAMOS LAS COVARIABLES ######

#SDT
SDTBoxCox<-boxcox(DatosCampo$SDT~1)
lambdaSDT<- SDTBoxCox$x[which.max(SDTBoxCox$y)]
BoxSDT<-(((DatosCampo$SDT)^lambdaSDT)-1 )/lambdaSDT
shapiro.test(BoxSDT)# p-value = 0.2692
ks.test(as.numeric(scale(sort(BoxSDT))),pnorm) #p-value  0.6238 existe normalidad

#DUREZA
DureBoxCox<-boxcox(DatosCampo$DUREZA~1)
lambdaDureza<- DureBoxCox$x[which.max(DureBoxCox$y)]
BoxDureza<-(((DatosCampo$DUREZA)^lambdaDureza)-1 )/lambdaDureza
shapiro.test(BoxDureza)# p-value = 0.09604 #### NO SE UTILIZA EN EL MODELO, cercano a 0.05
ks.test(as.numeric(scale(sort(BoxDureza))),pnorm) #p-value  0.6238 existe normalidad

#CL
CLBoxCox<-boxcox(DatosCampo$CL~1)
lambdaCL<- CLBoxCox$x[which.max(CLBoxCox$y)]
BoxCL<-(((DatosCampo$CL)^lambdaCL)-1 )/lambdaCL
shapiro.test(BoxCL)# p-value = 0.0015 ####SE UTILIZA EN EL MODELO 
ks.test(as.numeric(scale(sort(BoxCL))),pnorm) #p-value  0.6238 existe normalidad

#HCO3
HCO3BoxCox<-boxcox(DatosCampo$HCO3~1)
lambdaHCO3<- HCO3BoxCox$x[which.max(HCO3BoxCox$y)]
BoxHCO3<-(((DatosCampo$HCO3)^lambdaHCO3)-1 )/lambdaHCO3
shapiro.test(BoxHCO3)# p-value = 0.2429 ####NO SE UTILIZA EN EL MODELO 
ks.test(as.numeric(scale(sort(BoxHCO3))),pnorm) #p-value  0.6238 existe normalidad

#Color

ColorBoxCox<-boxcox(DatosCampo$COLOR~1)
lambdaColor<- ColorBoxCox$x[which.max(ColorBoxCox$y)]
BoxColor<-(((DatosCampo$COLOR)^lambdaColor)-1 )/lambdaColor
shapiro.test(BoxColor)# p-value = 0.01155 
ks.test(as.numeric(scale(sort(BoxColor))),pnorm) #p-value  0.5573 existe normalidad
hist(BoxColor)


##########################################################
#####COKRIGING ORDINARIO  ################################


###COKRIGING ORDINARIO

####Data frame con las covariables transformadas
gridded(PuntosInterpolar)<-TRUE
DatosCampo<-st_transform(DatosCampo,crs=st_crs(3116))
Coorx<-DatosCampo$X
Coory<-DatosCampo$Y
Dframe.Trans<-data.frame(Coorx,Coory,BoxWQI,BoxSDT,BoxDureza,BoxHCO3,BoxColor,BoxCL)
coordinates(Dframe.Trans)<-~Coorx+Coory
proj4string(Dframe.Trans)<-CRS("+init=epsg:3116")

CoKo=gstat(id = "WQI", formula = BoxWQI~1, data = Dframe.Trans)
CoKo=gstat(CoKo,id = "SDT", formula = BoxSDT~1, data = Dframe.Trans)
CoKo=gstat(CoKo,id = "Dure", formula = BoxDureza~1, data = Dframe.Trans)
CoKo=gstat(CoKo,id = "HCO3", formula = BoxHCO3~1, data = Dframe.Trans)
CoKo=gstat(CoKo,id = "Color", formula = BoxColor~1, data = Dframe.Trans)


CoKo=gstat(CoKo,id = "WQI", model=vgm(0.30,"Mat",8800,0.24,kappa = 0.33))  
CoKo=gstat(CoKo,id = "SDT", model=vgm(0.30,"Mat",8800,0.24,kappa = 0.33))
CoKo=gstat(CoKo,id = "Dure", model=vgm(0.30,"Mat",8800,0.24,kappa = 0.33))
CoKo=gstat(CoKo,id = "HCO3", model=vgm(0.30,"Mat",8800,0.24,kappa = 0.33))
CoKo=gstat(CoKo,id = "Color", model=vgm(0.30,"Mat",8800,0.24,kappa = 0.33))

Vcros<-variogram(CoKo,cutoff=Dusar)

CoKo<-gstat(CoKo,id="WQI",model=vgm(0.28, "Mat",8800,0.24,kappa = 0.33), fill.all=T)


###SE REALIZA EL AJUSTE DE LOS MODELOS CORRELACIONADOS
AjusteCokriging<-fit.lmc(Vcros,CoKo)

X11()
plot(variogram(CoKo), model=AjusteCokriging$model,pch=20)


#####################################################
###Validación cruzcada de CoKriging###


ValidacionCoKo=gstat.cv(AjusteCokriging,maxdist=3948.951)
criterio<-criterio.cv(ValidacionCoKo)
#           MPE     ASEPE     RMSPE         MSPE   RMSSPE     MAPPE      CCPE        R2
#1 0.0004657915 0.3388565 0.3455512 0.0006807316 1.017734 0.1047993 0.6581808 0.4331997
#-0.0009241966 0.3718417 0.3691441 -0.001504659 0.9871349 0.1137387 0.5948056
dataa<-data.frame(ValidacionCoKo$WQI.pred,ValidacionCoKo$WQI.var,ValidacionCoKo$observed,ValidacionCoKo$residual)
Export(dataa,"ResultadosTESISCokrigingICA.csv") 

mean(GaussCrps(mean(ValidacionCoKo$WQI.pred),sd(ValidacionCoKo$WQI.var),mean(ValidacionCoKo$observed)))
#0.0006279256



#######################################################
########PREDICCIONES COKRIGING #########################
PrediccionCok<-predict(AjusteCokriging,PuntosInterpolar)



PrediccionCok$Prediccion<-InversaBoxcox(PrediccionCok$WQI.pred,lambdaWQI)

summary(PrediccionCok$Prediccion)





##MAPA DE PREDICCION COKRIGING
PrediccionCok$Prediccion<-InversaBoxcox(PrediccionCok$WQI.pred,lambdaWQI)
summary(PrediccionCok$Prediccion)

pred_mean<-PrediccionCok$Prediccion
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
    name = "ICA",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()




ValidacionCoKs=gstat.cv(AjusteCokriging, maxdist=3948.951)
criterio<-criterio.cv(ValidacionCoKs)
summary(Prediccion.ks$var1.pred)
summary(Prediccion.ko$var1.pred)
summary(PrediccionCok$WQI.pred)
  