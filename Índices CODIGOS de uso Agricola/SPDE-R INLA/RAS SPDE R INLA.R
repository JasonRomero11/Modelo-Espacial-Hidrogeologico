#########################################################################################
#######"MODELADO ESPACIAL HIDROGEOLÓGICO PARA DETERMINAR ÍNDICES DE CALIDAD##############
#######Y VULNERABILIDAD DE LAS AGUAS SUBTERRÁNEAS EN LA ZONA CENTRO DE BOYACÁ"###########

#################################################################
####CÁLCULO DE INDICE DE CALIDAD PARA USO AGRICOLA -SPDE ####
#################################################################






library(inlabru)
library(devtools)
library(INLAutils)
library(INLA)
library(sp)
library(spgwr)
library(raster)
library(ggplot2)
library(sf)
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
library(gridExtra)


###############################################################
#######LECTURA DE ARCHIVOS VECTORIALES ########################
###############################################################


MunicipiosZonaEstudio<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO_B/MunicipiosZonaEstudio.shp")
DatosCampo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/MuestroBalanceIonico.shp")

UsoDelSuelo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/BOYACA_CAPACIDAD_VF/UsoSueloZonaEstudio.shp")
ShapeZona<-readOGR("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO/unidad_crono_Boyaca_zona_centro.shp")


epsg3116<-"+init=epsg:3116"
ShapeZona<-spTransform(ShapeZona,CRS(epsg3116))

PuntosInterpolar<- spsample(ShapeZona,n=70000,type="regular")

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
C1S1<-RAS[which(RAS$RAS<10 & RAS$DatosCampo.CE__25C_<250),]  #150
C2S1<-RAS[which(RAS$RAS<10 & RAS$DatosCampo.CE__25C_>250 & RAS$DatosCampo.CE__25C_<750),] #42
C3S1<-RAS[which(RAS$RAS<10 & RAS$DatosCampo.CE__25C_>750 & RAS$DatosCampo.CE__25C_<2250),] #16
C1S2<-RAS[which(RAS$RAS>10 &RAS$DatosCampo.CE__25C_<250),]
C4S1<-RAS[which(RAS$RAS<4 & RAS$DatosCampo.CE__25C_>2250),]

DatosCampo$RAS<-RAS$RAS


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

qplot(DatosCampo$RASTransGauss,geom="histogram",main = "", xlab = "Temperatura", ylab = "Frecuencia",fill=I("blue"),col=I("black"), alpha=I(.3))

############################################################################
#######CORRECCION DE NORMALIDAD POR BOX-COX ################################
############################################################################


RASBoxCox<-boxcox(DatosCampo$RAS~1)
lambdaRAS <- RASBoxCox$x[which.max(RASBoxCox$y)]
BoxRAS<- (((DatosCampo$RAS)^lambdaRAS)-1)/lambdaRAS
shapiro.test(BoxRAS) #p-value 0.9305 > 0.05 % Ya existe normalidad
ks.test(as.numeric(scale(sort(BoxRAS))),pnorm) #p-value  0.8804 ya existe normalidad
DatosCampo$RAS.TRANS<-BoxRAS



###########CORRELACIONES ALTAS ###################

#TEM, SDT, DUREZA, SO4, CL, CE 

#TEM 

TEMBoxCox<-boxcox(DatosCampo$T~1)
lambdaTEM <- TEMBoxCox$x[which.max(TEMBoxCox$y)]
BoxTEM<- (((DatosCampo$T)^lambdaTEM)-1)/lambdaTEM
shapiro.test(BoxTEM) #No existe normalidad
ks.test(as.numeric(scale(sort(BoxTEM))),pnorm) #p-value  0.8804 ya existe normalidad

#SDT
SDTBoxCox<-boxcox(DatosCampo$SDT~1)
lambdaSDT <- SDTBoxCox$x[which.max(SDTBoxCox$y)]
BoxSDT<- (((DatosCampo$SDT)^lambdaSDT)-1)/lambdaSDT
shapiro.test(BoxSDT) #Existe normalidad
ks.test(as.numeric(scale(sort(BoxSDT))),pnorm) # ya existe normalidad

#Dureza
DUREZABoxCox<-boxcox(DatosCampo$DUREZA~1)
lambdaDUREZA <- DUREZABoxCox$x[which.max(DUREZABoxCox$y)]
BoxDUREZA<- (((DatosCampo$DUREZA)^lambdaDUREZA)-1)/lambdaDUREZA
shapiro.test(BoxDUREZA) #Existe normalidad
ks.test(as.numeric(scale(sort(BoxDUREZA))),pnorm) # ya existe normalidad

#CL
CLBoxCox<-boxcox(DatosCampo$CL~1)
lambdaCL <- CLBoxCox$x[which.max(CLBoxCox$y)]
BoxCL<- (((DatosCampo$CL)^lambdaCL)-1)/lambdaCL
shapiro.test(BoxCL) #Existe normalidad
ks.test(as.numeric(scale(sort(BoxCL))),pnorm) # ya existe normalidad

#CE
CEBoxCox<-boxcox(DatosCampo$CE__25C_~1)
lambdaCE <- CEBoxCox$x[which.max(CEBoxCox$y)]
BoxCE<- (((DatosCampo$CE)^lambdaCE)-1)/lambdaCE
shapiro.test(BoxCE) #Existe normalidad
ks.test(as.numeric(scale(sort(BoxCE))),pnorm) # ya existe normalidad
#HC03
HCO3BoxCox<-boxcox(DatosCampo$HCO3~1)
lambdaHCO3<- HCO3BoxCox$x[which.max(HCO3BoxCox$y)]
BoxHCO3<-(((DatosCampo$HCO3)^lambdaHCO3)-1 )/lambdaHCO3
shapiro.test(BoxHCO3)# p-value = 0.2429
ks.test(as.numeric(scale(sort(BoxHCO3))),pnorm) 

Estandarizacion <-function(x) {(x-mean(x))/sd(x)}

SO4Estan<-Estandarizacion(DatosCampo$SO4)


DframeSPDE<-data.frame(DatosCampo$X,DatosCampo$Y,BoxTEM,BoxSDT,BoxDUREZA,BoxCL,BoxCE,SO4Estan,BoxHCO3)

DframeSPDE$RAS<-BoxRAS


#############################################################
####CONSIDERANDO una cobertura no convexa ##################

Coordenadas<-data.frame(DatosCampo$X,DatosCampo$Y)
coordinates(Coordenadas)=~DatosCampo.X+DatosCampo.Y
proj4string(Coordenadas)<-CRS("+init=epsg:3116")
MeshNoConvex <- inla.nonconvex.hull(Coordenadas)

MeshNoConvex1<-inla.mesh.2d(
  boundary = MeshNoConvex, max.edge = c(800,  3500),offset = c(500,1500), min.angle = c(26,16), cutoff = 0,plot.delay = NULL
)



RASSPDEmaterMesNohConvex1 <- inla.spde2.matern(mesh = MeshNoConvex1, alpha = 2,constr = TRUE) ##Modelo Mathern

#Luego generamos el conjunto de índices para el modelo SPDE
#Especificamos el nombre del efecto ( s), s es GAUSS y el número de vértices en el modelo SPDE

s.indexMeshNoConvex1 <- inla.spde.make.index(name="spatial.field",n.spde=RASSPDEmaterMesNohConvex1 $n.spde)


#DATOS DE PREDICCION

s.proyeccionMeshNoConvex1 <- inla.spde.make.A(mesh = MeshNoConvex1, loc = Coordenadas)

s.prediccionMeshNoConvex1<-inla.spde.make.A(mesh = MeshNoConvex1, loc = PuntosInterpolar)

####sin covariables

stk.MeshNoConvex1 <- inla.stack(
  tag = "est",
  data = list(Y.RAS= DframeSPDE$RAS),
  A = list(1, s.proyeccionMeshNoConvex1 ),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas)))), s = s.indexMeshNoConvex1)
)
formMeshNoConvex1 <- Y.RAS ~ 0 + b0 + f(spatial.field, model = RASSPDEmaterMesNohConvex1)
RAS.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.RAS = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,RAS.stackpredict)


RASINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                             data = inla.stack.data(Union_stack, spde =RASSPDEmaterMesNohConvex1),
                             control.predictor = list(compute=TRUE,A = inla.stack.A(Union_stack)),
                             control.compute = list(dic=TRUE, waic=TRUE))
RASINLAMeshNoConvex1$dic$dic # 556.4317
RASINLAMeshNoConvex1$waic$waic #  558.598

RASINLAMeshNoConvex1$summary.fixed
#b0 -1.214203 0.1650689   -1.56054 -1.207648 -0.9058055


sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 RASINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
# 0.5153048 0.1262194 0.3158146 0.4982938 0.8084384 
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(RASINLAMeshNoConvex1,'spatial.field', RASSPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
#1.182894 0.2603136 0.7562067 1.153363 1.774055
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
#4543.293 1178.315 2720.131 4370.683 7315.713

IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data

pred_mean2 <-RASINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
pred_mean2<-InversaBoxcox(pred_mean2,lambdaRAS)
summary(DatosCampo$RAS)
summary(pred_mean2)



dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean2, variable = "Prediccion RAS"
  ))
dpm$variable <- as.factor(dpm$variable)
ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("Índice de calidad del agua RAS \nPredicciones SPDE sin Covariables")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()


################################################
####CON COVARIABLES  ###########################


stk.MeshNoConvex1<- inla.stack(
  tag = "est",
  data = list(Y.RAS = DframeSPDE$RAS),
  A = list(1,s.proyeccionMeshNoConvex1),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas))),list(DframeSPDE[,3:10])), s = s.indexMeshNoConvex1)
)

formMeshNoConvex1 <- Y.RAS ~ 0+b0+SO4Estan+BoxCL+BoxDUREZA+BoxSDT+BoxTEM+BoxHCO3+BoxCE+f(spatial.field, model = RASSPDEmaterMesNohConvex1)
RAS.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.RAS = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,RAS.stackpredict)


RASINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                             data = inla.stack.data(Union_stack, spde =RASSPDEmaterMesNohConvex1),
                             control.predictor = list(A = inla.stack.A(Union_stack)),
                             control.compute = list(dic=TRUE, waic=TRUE,cpo=TRUE))
RASINLAMeshNoConvex1$summary.fixed
                mean         sd 0.025quant   0.5quant 0.975quant
                b0        -6.3958400 0.77515711 -7.91169160 -6.3977725 -4.8721341
                SO4Estan   0.2409619 0.05847275  0.12654201  0.2408085  0.3560888
                BoxCL      0.4208059 0.05220126  0.31832093  0.4208043  0.5232192
                BoxDUREZA -1.1093627 0.10769589 -1.32063913 -1.1094783 -0.8977101
                BoxSDT     1.3238734 0.21202871  0.90761218  1.3237597  1.7403205
                BoxTEM     0.5935498 0.18282722  0.23265319  0.5945593  0.9492762
                BoxHCO3    0.2215272 0.04937489  0.12474001  0.2214493  0.3186400
                BoxCE      0.2196909 0.15829783 -0.09090703  0.2195553  0.5306369
ggplot_inla_residuals(RASINLAMeshNoConvex1, DframeSPDE$RAS.TRANS, binwidth = 0.1)

RASINLAMeshNoConvex1$dic$dic #  -67.2153
RASINLAMeshNoConvex1$waic$waic # -22.95323
p <- autoplot(RASINLAMeshNoConvex1)
cowplot::plot_grid(plotlist = p)

p[[1]] +
  geom_line(aes(colour = var), size = 1.3) +geom_vline ( xintercept = 0 ) +
  palettetown::scale_colour_poke(pokemon = 'Charmander', spread = 9)
par(mfrow = c(1, 2))
hist(WQIINLAMeshNoConvex1$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(WQIINLAMeshNoConvex1$cpo$pit))), 
       WQIINLAMeshNoConvex1$cpo$pit, main = "Q-Q plot for Unif(0,1)", 
       xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales")

####CON COVARIABLES  SIGNIFICATIVAS


stk.MeshNoConvex1<- inla.stack(
  tag = "est",
  data = list(Y.RAS = DframeSPDE$RAS),
  A = list(1,s.proyeccionMeshNoConvex1),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas))),list(DframeSPDE[,3:10])), s = s.indexMeshNoConvex1)
)

formMeshNoConvex1 <- Y.RAS ~ 0+b0+SO4Estan+BoxCL+BoxDUREZA+BoxSDT+BoxTEM+BoxHCO3+f(spatial.field, model = RASSPDEmaterMesNohConvex1)
RAS.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.RAS = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,RAS.stackpredict)


RASINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                             data = inla.stack.data(Union_stack, spde =RASSPDEmaterMesNohConvex1),
                             control.predictor = list(A = inla.stack.A(Union_stack)),
                             control.compute = list(dic=TRUE, waic=TRUE,cpo=TRUE))
RASINLAMeshNoConvex1$summary.fixed
               mean         sd 0.025quant   0.5quant 0.975quant
b0        -6.976 & 0.725& -8.400 &-6.977& -5.550
SO4Estan   0.267 &0.061 & 0.146 &0.267& 0.388
BoxCL      0.431 &0.050 & 0.332 & 0.431&  0.530
BoxDUREZA -1.105 &0.097 &-1.297& -1.105& -0.912
BoxSDT     1.494 &0.207 & 1.086 & 1.494&  1.901
BoxTEM     0.786 &0.161 & 0.468& 0.786 & 1.103
BoxHCO3    0.212 &0.049 & 0.114 & 0.212 & 0.309

RASINLAMeshNoConvex1$dic$dic # 179.663
RASINLAMeshNoConvex1$waic$waic #169.439
p <- autoplot(RASINLAMeshNoConvex1)
cowplot::plot_grid(plotlist = p)

p[[1]] +
  geom_line(aes(colour = var), size = 1.3) +geom_vline ( xintercept = 0 ) +
  palettetown::scale_colour_poke(pokemon = 'Charmander', spread = 9)
par(mfrow = c(1, 2))
hist(RASINLAMeshNoConvex1$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(RASINLAMeshNoConvex1$cpo$pit))), 
       RASINLAMeshNoConvex1$cpo$pit, main = "Q-Q plot for Unif(0,1)", 
       xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales")



IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data


pred_mean2 <-RASINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
summary(pred_mean2)
###ANTI-TRANSFORMACIÓN DE LOS DATOS 
#Función
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
                              
 #####Anti-transformación de los betas---Se tomo el exponencial, ya que el lambda es cercano a cero. 


Intercepto<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$b0)0.00198013
SO4<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$SO4Esta)1.276657
CL<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$BoxCL)1.540736
Dureza<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$BoxDUREZA)0.331
SDT<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$BoxSDT)4.25534
TEM<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$BoxTEM)2.028296
HCO3<-inla.emarginal(exp, RASINLAMeshNoConvex1$marginals.fixed$BoxHCO3)1.244986



InversaBoxcoxCl<-function(z){
  
  res=exp(log(z*lambdaCL+1)/lambdaCL)
  return(res)
}
InversaBoxcoxCl(0.431)
InversaBoxcoxDureza<-function(z){
  
  res=exp(log(z*lambdaDUREZA+1)/lambdaDUREZA)
  return(res)
}
InversaBoxcoxDureza(-1.105)
InversaBoxcoxSDT<-function(z){
  
  res=exp(log(z*lambdaDUREZA+1)/lambdaDUREZA)
  return(res)
}
summary(DframeSPDE$RAS)


sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 RASINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")

#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(RASINLAMeshNoConvex1,'spatial.field', RASSPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")

range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
# 952.6368 139.654 695.6362 947.5958 1238.746 




summary(pred_mean2)
summary(DatosCampo$RAS.TRANS)

pred_mean2 <-RASINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
pred_mean2<-InversaBoxcox(pred_mean2,lambdaRAS)
summary(DatosCampo$RAS)
summary(pred_mean2)



dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean2, variable = "Prediccion RAS"
  ))
dpm$variable <- as.factor(dpm$variable)
ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "RAS",
    low = "blue", high = "orange"
  ) +  
  ggtitle("Índice de calidad del agua RAS \nPredicciones SPDE con Covariables")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()







#######VALIDACION CRUZADA ENFOQUE SPDE ##################


##con COVARIABLES 

Covariables1<-data.frame(DatosCampo$X,DatosCampo$Y,DframeSPDE[3:8])
coords<-scale(Covariables1[, c('DatosCampo.X', 'DatosCampo.Y')])
Covariables<-data.frame(DatosCampo$X,DatosCampo$Y,BoxCE,BoxDUREZA,BoxCL,BoxTEM,SO4Estan,BoxHCO3,BoxSDT,DframeSPDE$RAS)
colnames(coords) <- c('long', 'lat')
dataf1 <- sp::SpatialPointsDataFrame(coords = coords, data = Covariables[, -c(1:2)])
dataframe <- data.frame(dataf1)

Dmax<-sqrt((max(coords[,1])-min(coords[,1]))^2+(max(coords[,2])-min(coords[,2]))^2)
Dusar<-Dmax/15 # Valor numérico que proporciona el radio de la zona de influencia espacial alrededor de la ubicación del punto excluido para la validacion cruzada
ss <- 211
formMeshNoConvex1 <- Y.RAS ~ 0+y.intercept+SO4Estan+BoxCL+BoxDUREZA+BoxSDT+BoxTEM+BoxHCO3
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'DframeSPDE.RAS', ss = ss,
               rad = Dusar,alpha=0.05,
               modform =0.3779013, formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)



dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDERASCON07474.csv") 



##SIN COVARIABLES 


formMeshNoConvex1 <- Y.RAS ~ 0+y.intercept
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'DframeSPDE.RAS', ss = ss,
               rad = 0.3779013,alpha=0.05,
               modform = formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)
dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDERASSIN12436.csv") 


library(SpecsVerification)
setwd("C:/TESIS/RAS/SPDE")
ValidacionCruzadaSPDECO <-read.table("ResultadosFinalesSPDERASCON07474.csv",header=TRUE,sep=",")

MediaPredicciones<-mean(ValidacionCruzadaSPDECO$Predicciones)
DesviacionPrediccion<-sd(ValidacionCruzadaSPDECO$Predicciones) 
MediaObservaciones<-mean(ValidacionCruzadaSPDECO$Observaciones)
GaussCrps(MediaPredicciones,DesviacionPrediccion,MediaObservaciones)
#0.23703


setwd("C:/TESIS/RAS/SPDE")
ValidacionCruzadaSPDE <-read.table("ResultadosFinalesSPDERASSIN124.csv",header=TRUE,sep=",")

MediaPredicciones<-mean(ValidacionCruzadaSPDE$Predicciones)
DesviacionPrediccion<-sd(ValidacionCruzadaSPDE$Predicciones) 
MediaObservaciones<-mean(ValidacionCruzadaSPDE$Observaciones)
GaussCrps(MediaPredicciones,DesviacionPrediccion,MediaObservaciones)
#0.00223027

