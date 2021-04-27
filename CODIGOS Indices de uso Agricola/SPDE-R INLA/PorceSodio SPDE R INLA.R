#########################################################################################
#######"MODELADO ESPACIAL HIDROGEOLÓGICO PARA DETERMINAR ÍNDICES DE CALIDAD##############
#######Y VULNERABILIDAD DE LAS AGUAS SUBTERRÁNEAS EN LA ZONA CENTRO DE BOYACÁ"###########

#################################################################
####CÁLCULO DEl % DE SODIO-SPDE ######## #########################
#################################################################

####LIBRERIAS A UTILIZAR

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


Coordenadas<-read.table("C:/PROYECTODEINVESTIGACION/SHAPEFILE/Depuracion.csv",header=TRUE,sep=",")


epsg3116<-"+init=epsg:3116"
ShapeZona<-spTransform(ShapeZona,CRS(epsg3116))

PuntosInterpolar<- spsample(ShapeZona,n=70000,type="regular")


##########################################################
#######CALCULO DEL %Na ###################################
##########################################################


#####################################
### CONVERSION DE MG/L A meq/L ######
#####################################

####SE TRABAJA CON LOS CATIONES CA,MG,NA,K###########
#####################################################


Iones<-data.frame(DatosCampo$CA,DatosCampo$MG,DatosCampo$NA.,DatosCampo$K)

PesosMolecularIones<-c(40.078 ,24.305,22.989,39.1)

ValenciaIones<-c(2,2,1,1)

PesoEquivalenteIones<-PesosMolecularIones/ValenciaIones

Iones$Cameq<-((Iones$DatosCampo.CA)/PesoEquivalenteIones[1])
Iones$Mgmeq<-((Iones$DatosCampo.MG)/PesoEquivalenteIones[2])
Iones$Nameq<-((Iones$DatosCampo.NA.)/PesoEquivalenteIones[3])
Iones$Kmeq<-((Iones$DatosCampo.K)/PesoEquivalenteIones[4])



########################################
## CALCULO DEL PORCENTAJE DE SODIO #####
########################################

Iones$PorcentajeSodio<-(((Iones$Nameq+Iones$Kmeq)/(Iones$Cameq+Iones$Mgmeq+Iones$Nameq)))*100
DatosCampo$PorceSodio<-Iones$PorcentajeSodio

summary(DatosCampo$PorceSodio) 
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
shapiro.test(BoxNA) #
ks.test(as.numeric(scale(sort(BoxNA))),pnorm) #p

DatosCampo$NaTrasBox<-BoxNA

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
ks.test(as.numeric(scale(sort(PH.trans))),pnorm) #p-value < 0.006684 no existe normalidad



shapiro.test(CLORO.trans) #No normalizo 
ks.test(as.numeric(scale(sort(CLORO.trans))),pnorm) #p-value < 0.006684 no existe normalidad


##CON CL
CLselect = DatosCampo[,c("X","Y","CL")] 
CL.rgdb <- db.create(CLselect,ndim=2,autoname=F)
CL.herm <- anam.fit(CL.rgdb,name="CL",type="gaus")
CL.hermtrans <- anam.z2y(CL.rgdb,names="CL",anam=CL.herm)
CL.trans <- CL.hermtrans@items$Gaussian.CL
shapiro.test(CL.trans) # W = 0.96701, p-value = 7.701e-05
ks.test(as.numeric(scale(sort(CL.trans))),pnorm) #D = 0.069221, p-value = 0.2642




##SDT
SDTselect = DatosCampo[,c("X","Y","SDT")] 
SDT.rgdb <- db.create(SDTselect,ndim=2,autoname=F)
SDT.herm <- anam.fit(SDT.rgdb,name="SDT",type="gaus")
SDT.hermtrans <- anam.z2y(SDT.rgdb,names="SDT",anam=SDT.herm)
SDT.trans <- SDT.hermtrans@items$Gaussian.SDT
shapiro.test(SDT.trans) ####normalizo
ks.test(as.numeric(scale(sort(SDT.trans))),pnorm) #p-value < 0.006684 no existe normalidad

##TEM 
TEMselect = DatosCampo[,c("X","Y","T")] 
TEM.rgdb <- db.create(TEMselect,ndim=2,autoname=F)
TEM.herm <- anam.fit(TEM.rgdb,name="T",type="gaus")
TEM.hermtrans <- anam.z2y(TEM.rgdb,names="T",anam=TEM.herm)
TEM.trans <- TEM.hermtrans@items$Gaussian.T
shapiro.test(TEM.trans) ##Normalizo  #W = 0.99815, p-value = 0.9979
ks.test(as.numeric(scale(sort(TEM.trans))),pnorm) #p-value D = 0.026029, p-value = 0.9988

#hco3
HCO3select = DatosCampo[,c("X","Y","HCO3")] 
HCO3.rgdb <- db.create(HCO3select,ndim=2,autoname=F)
HCO3.herm <- anam.fit(HCO3.rgdb,name="HCO3",type="gaus")
HCO3.hermtrans <- anam.z2y(HCO3.rgdb,names="HCO3",anam=HCO3.herm)
HCO3.trans <- HCO3.hermtrans@items$Gaussian.HCO3
shapiro.test(HCO3.trans) ##
ks.test(as.numeric(scale(sort(HCO3.trans))),pnorm) #

DframeSPDE<-data.frame(DatosCampo$X,DatosCampo$Y,SDT.trans,TEM.trans,CL.trans,PH.trans,DUREZA.trans,HCO3.trans,Na.trans)


Coordenadas<-data.frame(DatosCampo$X,DatosCampo$Y)
coordinates(Coordenadas)=~DatosCampo.X+DatosCampo.Y
proj4string(Coordenadas)<-CRS("+init=epsg:3116")

#############################################################
####CONSIDERANDO una cobertura no convexa ##################


MeshNoConvex <- inla.nonconvex.hull(Coordenadas)

MeshNoConvex1<-inla.mesh.2d(
  boundary = MeshNoConvex, max.edge = c(800,  3500),offset = c(500,1500), min.angle = c(26,16), cutoff = 0,plot.delay = NULL
)


NaSPDEmaterMesNohConvex1 <- inla.spde2.matern(mesh = MeshNoConvex1, alpha = 2,constr = TRUE) ##Modelo Mathern

#Luego generamos el conjunto de índices para el modelo SPDE
#Especificamos el nombre del efecto ( s), s es GAUSS y el número de vértices en el modelo SPDE

s.indexMeshNoConvex1 <- inla.spde.make.index(name="spatial.field",n.spde=NaSPDEmaterMesNohConvex1 $n.spde)


#DATOS DE PREDICCION

s.proyeccionMeshNoConvex1 <- inla.spde.make.A(mesh = MeshNoConvex1, loc = Coordenadas)

s.prediccionMeshNoConvex1<-inla.spde.make.A(mesh = MeshNoConvex1, loc = PuntosInterpolar)



################################################
####sin covariables ############################

stk.MeshNoConvex1 <- inla.stack(
  tag = "est",
  data = list(Y.NA= DframeSPDE$Na.trans),
  A = list(1, s.proyeccionMeshNoConvex1 ),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas)))), s = s.indexMeshNoConvex1)
)
formMeshNoConvex1 <- Y.NA ~ 0  +b0+ f(spatial.field, model = NaSPDEmaterMesNohConvex1)
NA.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.WQI = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,NA.stackpredict)


NaINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                            data = inla.stack.data(Union_stack, spde =NaSPDEmaterMesNohConvex1),
                            control.predictor = list(compute=TRUE,A = inla.stack.A(Union_stack)),
                            control.compute = list(dic=TRUE, waic=TRUE))
NaINLAMeshNoConvex1$dic$dic # 497.2182
NaINLAMeshNoConvex1$waic$waic # 497.2182
NaINLAMeshNoConvex1$summary.fixed

         mean        sd 0.025quant    0.5quant 0.975quant       
b0 -0.0279672 0.1166517 -0.2609934 -0.02756016  0.2022277
sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 NaINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
#0.4105331 0.1149519 0.2266311 0.396545 0.6747191 
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(NaINLAMeshNoConvex1,'spatial.field', NaSPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
#0.7117654 0.1893834 0.4131165 0.6864809 1.152168
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
3781.998 1198.994 1899.512 3627.555 6562.201  

IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data
pred_mean2 <-NaINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
coords <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords, PuntosInterpolar)
NA.datos = cbind(coordinates(DataCoords),pred_mean2)
NA.db = db.create(NA.datos,autoname = F)
Napred = anam.y2z(NA.db,names="pred_mean2",anam = Na.herm)

summary(DatosCampo$PorceSodio)
summary(Napred$Raw.pred_mean2)



pred_mean<-Napred$Raw.pred_mean2

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
    name = "%Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()

coordinates(dpm)<-~Este+Norte
gridded(dpm) <- TRUE
dfr <- rasterFromXYZ(dpm)
writeRaster(dfr, filename="SPDENa.tif", options=c('TFW=YES'),overwrite=TRUE)
Covariables1<-data.frame(DatosCampo$X,DatosCampo$Y,DframeSPDE$DUREZA.trans,DframeSPDE$TEM.trans,DframeSPDE$CL.trans,DframeSPDE$SDT.trans,DframeSPDE$Na.trans)
coords<-scale(Covariables1[, c('DatosCampo.X', 'DatosCampo.Y')])
Covariables<-data.frame(DatosCampo$X,DatosCampo$Y,DUREZA.trans,SDT.trans,CL.trans,TEM.trans,SDT.trans,DframeSPDE$Na.trans)
colnames(coords) <- c('long', 'lat')
dataf1 <- sp::SpatialPointsDataFrame(coords = coords, data = Covariables[, -c(1:2)])
dataframe <- data.frame(dataf1)

###############VALIDACION CRUZADA ##################

ss <- 211
rad <- min(3781.9981,max(dist(coords)) /40)

formMeshNoConvex1 <- Y.NA ~ 0+y.intercept 
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'DframeSPDE.Na.trans', ss = ss,
               rad = rad,alpha=0.05,
               modform = formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)
dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDENASIN1008.csv") # 

#############################################
###con todas covariables
stk.MeshNoConvex1<- inla.stack(
  tag = "est",
  data = list(Y.NA = DframeSPDE$Na.trans),
  A = list(1,s.proyeccionMeshNoConvex1),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas))),list(DframeSPDE[,3:9])), s = s.indexMeshNoConvex1)
)

formMeshNoConvex1 <- Y.NA ~ 0  +CL.trans+SDT.trans+DUREZA.trans+TEM.trans+HCO3.trans+f(spatial.field, model = NaSPDEmaterMesNohConvex1)

NA.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.NA = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,NA.stackpredict)


NaINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                            data = inla.stack.data(Union_stack, spde =NaSPDEmaterMesNohConvex1),
                            control.predictor = list(compute=TRUE,A = inla.stack.A(Union_stack)),
                            control.compute = list(dic=TRUE, waic=TRUE))
NaINLAMeshNoConvex1$dic$dic # 284.3417
NaINLAMeshNoConvex1$waic$waic # 291.5867
NaINLAMeshNoConvex1$summary.fixed

NaINLAMeshNoConvex1$summary.fixed

CL.trans      0.59845393 0.06993621  0.4610885090  0.59842134  0.7358637
SDT.trans     0.83187500 0.11819763  0.6001459287  0.83166217  1.0645483
DUREZA.trans -1.77531389 0.12062437 -2.0135977101 -1.77491624 -1.5394928
TEM.trans     0.08884218 0.04512550  0.0001029968  0.08886314  0.1773824
HCO3.trans    0.28119490 0.08151400  0.1214386692  0.28101882  0.4417569

sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 NaINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
#0.1454091 0.04399564 0.07848653 0.1387236 0.2497298
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(NaINLAMeshNoConvex1,'spatial.field', NaSPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
# 0.2317229 0.06531768 0.1286143 0.2231887 0.3833243 
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
# 3118.278 897.9902 1716.84 2994.972 5216.812

IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data
pred_mean2 <-NaINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
coords <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords, PuntosInterpolar)
NA.datos = cbind(coordinates(DataCoords),pred_mean2)
NA.db = db.create(NA.datos,autoname = F)
Napred = anam.y2z(NA.db,names="pred_mean2",anam = Na.herm)

summary(DatosCampo$PorceSodio)
summary(Napred$Raw.pred_mean2)



pred_mean<-Napred$Raw.pred_mean2
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean, variable = "Prediccion Na"
  ))
dpm$variable <- as.factor(dpm$variable)
ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "%Na",
    low = "blue", high = "orange"
  ) +  
  ggtitle("Índice de %Na \nPredicciones SPDE con Covariables")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()


###con todas covariables significativas
stk.MeshNoConvex1<- inla.stack(
  tag = "est",
  data = list(Y.NA = DframeSPDE$Na.trans),
  A = list(1,s.proyeccionMeshNoConvex1),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas))),list(DframeSPDE[,3:9])), s = s.indexMeshNoConvex1)
)

formMeshNoConvex1 <- Y.NA ~ 0  +TEM.trans+DUREZA.trans+CL.trans+SDT.trans+f(spatial.field, model = NaSPDEmaterMesNohConvex1)

Na.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.NA = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,Na.stackpredict)

NaINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                            data = inla.stack.data(Union_stack , spde =NaSPDEmaterMesNohConvex1),
                            control.predictor = list(A = inla.stack.A(Union_stack)),
                            control.compute = list(dic=TRUE, waic=TRUE))
NaINLAMeshNoConvex1$dic$dic # 270.823
NaINLAMeshNoConvex1$waic$waic # 281.227

NaINLAMeshNoConvex1$summary.fixed

TEM.trans  &   0.113 & 0.045 &  0.024 &  0.113 &  0.203 
DUREZA.trans &-1.571 & 0.106 & -1.783 & -1.571 & -1.363
CL.trans    &  0.533 & 0.068 &  0.397 &  0.533 &  0.668
SDT.trans    & 0.878 & 0.120 &  0.642 &  0.877 &  1.117


p <- autoplot(NaINLAMeshNoConvex1)
cowplot::plot_grid(plotlist = p)

p[[1]] +
  geom_line(aes(colour = var), size = 1.3) +geom_vline ( xintercept = 0 ) +
  palettetown::scale_colour_poke(pokemon = 'Charmander', spread = 5)
par(mfrow = c(1, 2))
hist(WQIINLAMeshNoConvex1$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(WQIINLAMeshNoConvex1$cpo$pit))), 
       WQIINLAMeshNoConvex1$cpo$pit, main = "Q-Q plot for Unif(0,1)", 
       xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales")


sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 NaINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
#0.137 & 0.047 &  0.066 & 0.129 & 0.251
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(NaINLAMeshNoConvex1,'spatial.field', NaSPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
#0.258 &  0.070 & 0.147 & 0.249 & 0.423
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
# 2601.694 & 749.6521 & 1414.11  & 2506.042  & 4337.058


##ANTITRANSFORMACION BETAS, 


##ANTITRANSFORMACION DE LOS DATOS
IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data
pred_mean2 <-NaINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
coords <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords, PuntosInterpolar)
NA.datos = cbind(coordinates(DataCoords),pred_mean2)
NA.db = db.create(NA.datos,autoname = F)
Napred = anam.y2z(NA.db,names="pred_mean2",anam = Na.herm)


##########################################
###VALIDACION CRUZADA 
##########################################


##sin covariables
Covariables1<-data.frame(DatosCampo$X,DatosCampo$Y,DframeSPDE$DUREZA.trans,DframeSPDE$TEM.trans,DframeSPDE$CL.trans,DframeSPDE$Na.trans)
coords<-scale(Covariables1[, c('DatosCampo.X', 'DatosCampo.Y')])
Covariables<-data.frame(DatosCampo$X,DatosCampo$Y,DUREZA.trans,SDT.trans,CL.trans,TEM.trans,HCO3.trans,DframeSPDE$Na.trans)
colnames(coords) <- c('long', 'lat')
dataf1 <- sp::SpatialPointsDataFrame(coords = coords, data = Covariables[, -c(1:2)])
dataframe <- data.frame(dataf1)
Dmax<-sqrt((max(coords[,1])-min(coords[,1]))^2+(max(coords[,2])-min(coords[,2]))^2)
Dusar<-Dmax/15 # Valor numérico que proporciona el radio de la zona de influencia espacial alrededor de la ubicación del punto excluido para la validacion cruzada
ss <- 211
formMeshNoConvex1 <- Y.NA ~ 0+y.intercept+f(spatial.field, model = NaSPDEmaterMesNohConvex1)
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'DframeSPDE.Na.trans', ss = ss,
               rad = 0.3779013,alpha=0.05,
               modform = formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)
dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDENASIN1016.csv") # 

#1.016109
#mae 0.8077261

#CON COVARIABLES 
Covariables1<-data.frame(DatosCampo$X,DatosCampo$Y,DframeSPDE$DUREZA.trans,DframeSPDE$TEM.trans,DframeSPDE$CL.trans,DframeSPDE$SDT.trans,HCO3.trans,DframeSPDE$Na.trans)
coords<-scale(Covariables1[, c('DatosCampo.X', 'DatosCampo.Y')])
Covariables<-data.frame(DatosCampo$X,DatosCampo$Y,DUREZA.trans,SDT.trans,CL.trans,TEM.trans,SDT.trans,HCO3.trans,DframeSPDE$Na.trans)
colnames(coords) <- c('long', 'lat')
dataf1 <- sp::SpatialPointsDataFrame(coords = coords, data = Covariables[, -c(1:2)])
dataframe <- data.frame(dataf1)

formMeshNoConvex1 <- Y.NA ~ 0  +CL.trans+SDT.trans+DUREZA.trans+TEM.trans+HCO3.trans+f(spatial.field, model = NaSPDEmaterMesNohConvex1)
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'DframeSPDE.Na.trans', ss = ss,
               rad = 0.3779013,alpha=0.05,
               modform = formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)
dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDECON0586871.csv") # 

##ANTITRANSFORMACION DE LOS DATOS
IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data
pred_mean2 <-NaINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
coords <- cbind(PuntosInterpolar$x1,PuntosInterpolar$x2)
DataCoords <- data.frame(coords, PuntosInterpolar)
NA.datos = cbind(coordinates(DataCoords),pred_mean2)
NA.db = db.create(NA.datos,autoname = F)
Napred = anam.y2z(NA.db,names="pred_mean2",anam = Na.herm)


#####Anti-transformación de los betas 

#para los betas por anamorfosis

#####Anti-transformación de los betas---Se tomo el exponencial, ya que el lambda es cercano a cero. 

Tem<-inla.emarginal(exp, NaINLAMeshNoConvex1$marginals.fixed$)
Colorex<-inla.emarginal(exp, WQIINLAMeshNoConvex1$marginals.fixed$BoxColor)

SDTex<-inla.emarginal(exp, WQIINLAMeshNoConvex1$marginals.fixed$BoxSDT)



summary(DatosCampo$PorceSodio)
summary(Napred$Raw.pred_mean2)



Napred$Prediccion<-Napred$Raw.pred_mean2
pred_mean2$Prediccion<-Napred$Raw.pred_mean2
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean2, variable = "pred_mean"
  ))
dpm$variable <- as.factor(dpm$variable)
ggplot(dpm) + geom_tile(aes(Este, Norte, fill = value)) +
  facet_wrap(~variable, nrow = 1) +
  coord_fixed(ratio = 1) +
  scale_fill_gradient(
    name = "ICA",
    low = "blue", high = "orange"
  ) +  
  ggtitle("Índice %Na \nPredicciones SPDE sin Covariables")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()


library(SpecsVerification)
setwd("C:/TESIS/Na/SPDE")
ValidacionCruzadaSPDECO <-read.table("ResultadosFinalesSPDECON0586871.csv",header=TRUE,sep=",")

MediaPredicciones<-mean(ValidacionCruzadaSPDECO$Predicciones)
DesviacionPrediccion<-sd(ValidacionCruzadaSPDECO$Predicciones) 
MediaObservaciones<-mean(ValidacionCruzadaSPDECO$Observaciones)
GaussCrps(MediaPredicciones,DesviacionPrediccion,MediaObservaciones)
#0.1934171

setwd("C:/TESIS/Na/SPDE")
ValidacionCruzadaSPDE <-read.table("ResultadosFinalesSPDENASIN1008.csv",header=TRUE,sep=",")

MediaPredicciones<-mean(ValidacionCruzadaSPDE$Predicciones)
DesviacionPrediccion<-sd(ValidacionCruzadaSPDE$Predicciones) 
MediaObservaciones<-mean(ValidacionCruzadaSPDE$Observaciones)
GaussCrps(MediaPredicciones,DesviacionPrediccion,MediaObservaciones)
#0.001763153

mean(GaussCrps(MediaPredicciones, DesviacionPrediccion, ValidacionCruzadaSPDECO$Observaciones))
