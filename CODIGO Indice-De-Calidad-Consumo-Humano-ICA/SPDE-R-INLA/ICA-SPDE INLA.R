 
#########################################################################################
#######"MODELADO ESPACIAL HIDROGEOLÓGICO PARA DETERMINAR ÍNDICES DE CALIDAD##############
#######Y VULNERABILIDAD DE LAS AGUAS SUBTERRÁNEAS EN LA ZONA CENTRO DE BOYACÁ"###########

#################################################################
####CÁLCULO DE INDICE DE CALIDAD PARA CONSUMO HUMANO-SPDE ####
#################################################################
####kLIBRERIAS A UTILIZAR


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
  
  MunicipiosZonaEstudio<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO/MunicipiosZonaEstudio.shp")
  DatosCampo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/MuestroBalanceIonico.shp")
  
  UsoDelSuelo<-st_read("C:/PROYECTODEINVESTIGACION/SHAPEFILE/BOYACA_CAPACIDAD_VF/UsoSueloZonaEstudio.shp")
  ShapeZona<-readOGR("C:/PROYECTODEINVESTIGACION/SHAPEFILE/ZONA_DE_ESTUDIO/unidad_crono_Boyaca_zona_centro.shp")
  
  
  Coordenadas<-read.table("C:/PROYECTODEINVESTIGACION/SHAPEFILE/Depuracion.csv",header=TRUE,sep=",")
  
  
  epsg3116<-"+init=epsg:3116"
  ShapeZona<-spTransform(ShapeZona,CRS(epsg3116))
  
  PuntosInterpolar<- spsample(ShapeZona,n=70000,type="regular")
  
  
  ###VARIABLES A UTILIZAR PARA EL IWQ
  ##Parámetros tenidos en cuenta
  BDIWQ<-data.frame(DatosCampo$PH,DatosCampo$CE__25C_,DatosCampo$TURBIEDAD,DatosCampo$ALCALINIDA,DatosCampo$PO4,DatosCampo$NO3,DatosCampo$SO4,DatosCampo$CA,DatosCampo$MG,DatosCampo$MN,DatosCampo$NA.) ###n=11 parámetros
  summary(BDIWQ)
  ##pesos de cada parámetro
  Parametros<-c("PH","Conductividad Electrica","Turbiedad","Alcalinidad","Fosfato","Nitratos","Sulfato","Calcio","Magensio","Manganeso","Sodio")
  Pesos<-c(1,5,4,4,1,5,4,4,3,3,4)
  SumPesos<-sum(Pesos)
  
  Wi<-c(Pesos/SumPesos)
  SumaWi=sum(Wi)
  
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
  
  mean(WQI.trans) #0.000227504
  
  DatosCampo$WQITransGauss<-WQI.trans
  
  shapiro.test(DatosCampo$WQITransGauss) #p-value < 7.34e-16 no existe normalidad
  ks.test(as.numeric(scale(sort(DatosCampo$WQITransGauss))),pnorm) #p-value < 2.656e-08 no existe normalidad
  
  hist(WQI.trans)
  
  ############################################################################
  #######CORRECCION DE NORMALIDAD POR BOX-COX ################################
  ############################################################################
  
  
  WQIBoxCox<-boxcox(DatosCampo$WQI~1)
  lambdaWQI <- WQIBoxCox$x[which.max(WQIBoxCox$y)]
  BoxWQI<- (((DatosCampo$WQI)^lambdaWQI)-1)/lambdaWQI 
  shapiro.test(BoxWQI) #p-value 0.7932 > 0.05 % Ya existe normalidad
  ks.test(as.numeric(scale(sort(BoxWQI))),pnorm) #p-value  0.9891 ya existe normalidad
  hist(BoxWQI)
  ############
  #AGREGAMOS LA VARIABLE TRANSFORMADA AL ARCHIVO SP Y DATA FRAME
  #POR UN LADO EL MANEJO DE MAPAS Y OTRO EL USO DEL DATA.FRAME  
  DatosCampo$WQI.TRANS<-BoxWQI
  TablaWQI$WQI.TRANS<-BoxWQI  ###ASIGNAMOS AL DATA FRAME 
  
  
  #############################################################
  #######CONSTRUCCIÓN DE LA MEJOR MALLA #######################
  #############################################################
  Coordenadas<-data.frame(DatosCampo$X,DatosCampo$Y)
  coordinates(Coordenadas)=~DatosCampo.X+DatosCampo.Y
  proj4string(Coordenadas)<-CRS("+init=epsg:3116")
  
  #############################################################
  ####CONSIDERANDO una cobertura no convexa ##################
  
  
  MeshNoConvex <- inla.nonconvex.hull(Coordenadas)
  
  MeshNoConvex1<-inla.mesh.2d(
    boundary = MeshNoConvex, max.edge = c(800,  3500),offset = c(500,1500), min.angle = c(26,16), cutoff = 0,plot.delay = NULL
  )
  plot(MeshNoConvex1)
  
  
  
  ####construccion SPDE COVARIANZA MATHERN
  WQISPDEmaterMesNohConvex1 <- inla.spde2.matern(mesh = MeshNoConvex1, alpha = 2,constr = TRUE) ##Modelo Mathern
  
  #Luego generamos el conjunto de índices para el modelo SPDE
  #Especificamos el nombre del efecto ( s), s es GAUSS y el número de vértices en el modelo SPDE
  
  s.indexMeshNoConvex1 <- inla.spde.make.index(name="spatial.field",n.spde=WQISPDEmaterMesNohConvex1 $n.spde)
  
  #DATOS DE PREDICCION
  
  s.proyeccionMeshNoConvex1 <- inla.spde.make.A(mesh = MeshNoConvex1, loc = Coordenadas)
  
  s.prediccionMeshNoConvex1<-inla.spde.make.A(mesh = MeshNoConvex1, loc = PuntosInterpolar)


  
######################################################
######################################################
####COVARIABLES DE INTERES ###########################
######################################################

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
shapiro.test(BoxDureza)# p-value = 0.09604 #### 
ks.test(as.numeric(scale(sort(BoxDureza))),pnorm) #p-value  existe normalidad

#HCO3
HCO3BoxCox<-boxcox(DatosCampo$HCO3~1)
lambdaHCO3<- HCO3BoxCox$x[which.max(HCO3BoxCox$y)]
BoxHCO3<-(((DatosCampo$HCO3)^lambdaHCO3)-1 )/lambdaHCO3
shapiro.test(BoxHCO3)# p-value = 0.2429 #### 
ks.test(as.numeric(scale(sort(BoxHCO3))),pnorm) #p-value   existe normalidad

#CLORO
CLBoxCox<-boxcox(DatosCampo$CL~1)
lambdaCL<- CLBoxCox$x[which.max(CLBoxCox$y)]
BoxCL<-(((DatosCampo$CL)^lambdaCL)-1 )/lambdaCL
shapiro.test(BoxCL)# p-value = 0.01155 
ks.test(as.numeric(scale(sort(BoxCL))),pnorm) #p-value  0.5573 existe normalidad

ColorBoxCox<-boxcox(DatosCampo$COLOR~1)
lambdaColor<- ColorBoxCox$x[which.max(ColorBoxCox$y)]
BoxColor<-(((DatosCampo$COLOR)^lambdaColor)-1 )/lambdaColor
shapiro.test(BoxColor)# p-value = 0.01155 
ks.test(as.numeric(scale(sort(BoxColor))),pnorm)


DframeSPDE<-data.frame(DatosCampo$X,DatosCampo$Y,BoxSDT,BoxDureza,BoxHCO3,BoxCL,BoxColor)
DframeSPDE$WQITransBox<-BoxWQI 


###########################################
###MODELO SIN COVARIABLES #################

stk.MeshNoConvex1 <- inla.stack(
  tag = "est",
  data = list(Y.WQI = DframeSPDE$WQITransBox),
  A = list(1, s.proyeccionMeshNoConvex1 ),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas)))), s = s.indexMeshNoConvex1)
)
formMeshNoConvex1 <- Y.WQI ~ 0 + b0 + f(spatial.field, model = WQISPDEmaterMesNohConvex1)

WQI.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.WQI = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,WQI.stackpredict)

WQIINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                             data = inla.stack.data(Union_stack, spde =WQISPDEmaterMesNohConvex1),
                             control.predictor = list(A = inla.stack.A(Union_stack)),
                             control.compute = list(dic=TRUE, waic=TRUE,cpo=TRUE))

WQIINLAMeshNoConvex1$dic$dic # 247.6371
WQIINLAMeshNoConvex1$waic$waic # 250.7069

WQIINLAMeshNoConvex1$summary.fixed


#        mean         sd 0.025quant 0.5quant 0.975quant     mode          kld
#b0 2.674074 0.05660289   2.552458 2.677438   2.777104
#-- Extract results for sigma2eps (1/precision)
sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 WQIINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
#0.1553871 0.02470659 0.1154646 0.1522066 0.2119214 
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(WQIINLAMeshNoConvex1,'spatial.field', WQISPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
#0.06706722 0.02831373 0.02555067 0.06260151 0.134808
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
#

#8163.904 5588.365 2301.083 6585.472 23215.68 


IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data

pred_mean2 <-WQIINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
pred_ll2 <- WQIINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "0.025quant"]
pred_ul2 <- WQIINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "0.975quant"]

InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}
pred_mean2<-InversaBoxcox(pred_mean2,lambdaWQI)
pred_ll2<-InversaBoxcox(pred_ll2,lambdaWQI)
pred_ul2<-InversaBoxcox(pred_ul2,lambdaWQI)

summary(pred_mean2)
summary(pred_ll2)
dpm <- rbind(
  data.frame(
    Este = PuntosInterpolar$x1, Norte = PuntosInterpolar$x2,
    value = pred_mean2, variable = ""
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


#################################################
###  MODELO CON COVARIABLES  ####################
#################################################

###CON TODAS

stk.MeshNoConvex1<- inla.stack(
  tag = "est",
  data = list(Y.WQI = DframeSPDE$WQITransBox),
  A = list(1,s.proyeccionMeshNoConvex1),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas))),list(DframeSPDE[,3:7])), s = s.indexMeshNoConvex1)
)

formMeshNoConvex1 <- Y.WQI ~ 0 + b0 +BoxSDT+BoxDureza+BoxColor+BoxHCO3+f(spatial.field, model = WQISPDEmaterMesNohConvex1)


WQIINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                             data = inla.stack.data(stk.MeshNoConvex1 , spde =WQISPDEmaterMesNohConvex1),
                             control.predictor = list(A = inla.stack.A(stk.MeshNoConvex1)),
                             control.compute = list(dic=TRUE, waic=TRUE))
WQIINLAMeshNoConvex1$dic$dic # 127.314
WQIINLAMeshNoConvex1$waic$waic # 140.071

WQIINLAMeshNoConvex1$summary.fixed

b0  $        1.274 & 0.184 &   0.912 &  1.274 &  1.637
BoxSDT $     0.280 & 0.092 &  0.099 &  0.280 & 0.460 
BoxDureza $ -0.021 &  0.053 & -0.125 & -0.021 & 0.083
BoxColor  $ 0.248 & 0.034 &  0.180 & 0.248 & 0.316 
BoxHCO3   $ 0.028 & 0.026 & -0.023 &  0.028 & 0.080
out.field <- inla.spde2.result(WQIINLAMeshNoConvex1,'spatial.field', WQISPDEmaterMesNohConvex1, do.transf = TRUE)
range.out <- inla.emarginal(function(x) x,
                            
                            out.field$marginals.range.nominal[[1]]) #rango 3772.292

sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 WQIINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
#0.1553871 0.02470659 0.1154646 0.1522066 0.2119214 
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(WQIINLAMeshNoConvex1,'spatial.field', WQISPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
#0.06706722 0.02831373 0.02555067 0.06260151 0.134808
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
# 7788.424 4936.08 2328.179 6461.127 20904.9

p <- autoplot(WQIINLAMeshNoConvex1)
cowplot::plot_grid(plotlist = p)

p[[1]] +
  geom_line(aes(colour = var), size = 1.3) +geom_vline ( xintercept = 0 ) +
  palettetown::scale_colour_poke(pokemon = 'Charmander', spread = 5)
par(mfrow = c(1, 2))
hist(WQIINLAMeshNoConvex1$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(WQIINLAMeshNoConvex1$cpo$pit))), 
       WQIINLAMeshNoConvex1$cpo$pit, main = "Q-Q plot for Unif(0,1)", 
       xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales")


########################################################################
### SON SIGNIFICATIVOS COLO SDT Y COLOR ################################


stk.MeshNoConvex1<- inla.stack(
  tag = "est",
  data = list(Y.WQI = DframeSPDE$WQITransBox),
  A = list(1,s.proyeccionMeshNoConvex1),
  effects = list(data.frame(b0 = rep(1, nrow(coordinates(Coordenadas))),list(DframeSPDE[,3:7])), s = s.indexMeshNoConvex1)
)

formMeshNoConvex1 <- Y.WQI ~ 0 + b0 +BoxSDT+BoxColor+f(spatial.field, model = WQISPDEmaterMesNohConvex1)
ICA.stackpredict <- inla.stack(
  tag = "pred",
  data = list(Y.WQI = NA),
  A = list(1,s.prediccionMeshNoConvex1),
  effects = list(data.frame(b0= rep(1, nrow(coordinates(PuntosInterpolar)))), s = s.indexMeshNoConvex1))
Union_stack <- inla.stack(stk.MeshNoConvex1,ICA.stackpredict)


WQIINLAMeshNoConvex1 <- inla(formMeshNoConvex1 , family = 'normal',
                             data = inla.stack.data(Union_stack, spde =WQISPDEmaterMesNohConvex1),
                             control.predictor = list(A = inla.stack.A(Union_stack)),
                             control.compute = list(dic=TRUE, waic=TRUE,cpo=TRUE))
WQIINLAMeshNoConvex1$dic$dic # 140.0035
WQIINLAMeshNoConvex1$waic$waic # 148.7923

WQIINLAMeshNoConvex1$summary.fixed

b0  &     1.248 & 0.142 &  0.969 & 1.248 &  1.530
BoxSDT &  0.294 & 0.037 & 0.219 & 0.294 &  0.367
BoxColor & 0.253 & 0.033 &  0.187 & 0.253 &  0.319

out.field <- inla.spde2.result(WQIINLAMeshNoConvex1,'spatial.field', WQISPDEmaterMesNohConvex1, do.transf = TRUE)
range.out <- inla.emarginal(function(x) x,
                            
                            out.field$marginals.range.nominal[[1]]) #rango 3772.292
#rango es 3063.268
sigma2e_marg =
  inla.tmarginal(function(x) 1/x,
                 WQIINLAMeshNoConvex1$marginals.hyperpar$"Precision for the Gaussian observations")
sigma2e_m1 = inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 = inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev = sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)
cat("-----Results for sigma2eps-----\n",
    c(sigma2e_m1, sigma2e_stdev, sigma2e_quantiles),
    "\n-----------------------------")
#0.092 & 0.023 & 0.054 & 0.089 & 0.147  
#-- Extract results for the covariance function parameters --#
mod.field = inla.spde2.result(WQIINLAMeshNoConvex1,'spatial.field', WQISPDEmaterMesNohConvex1, do.transf = TRUE)

var.nom.marg = mod.field$marginals.variance.nominal[[1]]
var.nom.m1 = inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 = inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev = sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)
cat("-----Results for sigma2omega-----\n",
    c(var.nom.m1, var.nom.stdev, var.nom.quantiles),
    "\n-----------------------------")
#0.049 & 0.034 &  0.011 & 0.040 & 0.141
range.nom.marg = mod.field$marginals.range.nominal[[1]]
range.nom.m1 = inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 = inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev = sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles = inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)
cat("-----Results for the range-----\n",
    c(range.nom.m1, range.nom.stdev, range.nom.quantiles),
    "\n-----------------------------")
#3063.268 1958.305 777.1765 2576.093 8189.847 

#MAPA DE PREDICCIONES
IndexPredic <- inla.stack.index(Union_stack, tag = "pred")$data


pred_mean2 <-WQIINLAMeshNoConvex1$summary.fitted.values[IndexPredic, "mean"]
summary(pred_mean2)

round(WQIINLAMeshNoConvex1$summary.linear.pred[IndexPredic,], 4)
###ANTI-TRANSFORMACIÓN DE LOS DATOS 
#Función
InversaBoxcox<-function(z,lambda){
  
  res=exp(log(z*lambda+1)/lambda)
  return(res)
}

summary(InversaBoxcox(pred_mean2,lambdaWQI))
pred_mean2<-InversaBoxcox(pred_mean2,lambdaWQI)
pred_mean2<-exp(pred_mean2)

summary(pred_mean2)

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
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5))
theme_bw()


#####Anti-transformación de los betas---Se tomo el exponencial, ya que el lambda es cercano a cero. 

Intercepto<-inla.emarginal(exp, WQIINLAMeshNoConvex1$marginals.fixed$b0)
exp(Intercepto)# 34.13769
SDTex<-inla.emarginal(exp, WQIINLAMeshNoConvex1$marginals.fixed$BoxSDT)
exp(SDTex)#3.827826
Colorex<-inla.emarginal(exp, WQIINLAMeshNoConvex1$marginals.fixed$BoxColor)
exp(Colorex)#3.632244


round(WQIINLAMeshNoConvex1$summary.fixed,6)$mean*100

######SIGNIFICANCIA DE LOS BETAS ###############

p <- autoplot(WQIINLAMeshNoConvex1)
cowplot::plot_grid(plotlist = p)

p[[1]] +
  geom_line(aes(colour = var), size = 1.3) +geom_vline ( xintercept = 0 ) +
  palettetown::scale_colour_poke(pokemon = 'Charmander', spread = 5)
par(mfrow = c(1, 2))
hist(WQIINLAMeshNoConvex1$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(WQIINLAMeshNoConvex1$cpo$pit))), 
       WQIINLAMeshNoConvex1$cpo$pit, main = "Q-Q plot for Unif(0,1)", 
       xlab = "Cuantiles teóricos", ylab = "Cuantiles muestrales")





###################################################################
##VALIDACION CRUZADA SPDE CON LA FUCNIÓN INLASLOO #################



##########################
##########################
##SPDE CON COVARIABLES ###
Covariables1<-data.frame(DatosCampo$X,DatosCampo$Y,BoxSDT,BoxColor,DatosCampo$SDT,DatosCampo$COLOR,DframeSPDE$WQITransBox)
coords<-scale(Covariables1[, c('DatosCampo.X', 'DatosCampo.Y')])
colnames(coords) <- c('long', 'lat')
dataf1 <- sp::SpatialPointsDataFrame(coords = coords, data = DframeSPDE[, -c(1:2)])
dataframe <- data.frame(dataf1)
Dmax<-sqrt((max(coords[,1])-min(coords[,1]))^2+(max(coords[,2])-min(coords[,2]))^2)
Dusar<-Dmax/15 # Valor numérico que proporciona el radio de la zona de influencia espacial alrededor de la ubicación del punto excluido para la validacion cruzada
ss <- 211

formMeshNoConvex1 <- Y.WQI ~ 0 + y.intercept+BoxSDT+BoxColor+f(spatial.field, model = WQISPDEmaterMesNohConvex1)
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'WQITransBox', ss = ss,
               rad = 0.3779013,alpha=0.05,
               modform = formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)

dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDEICA036.csv") # RMSE 0.3784728 #MAE 0.3004132
#rmse 0.3600953 mae 0.2926
data(meuse)
#3029.715 /45 rmse 0.3587793 y MAE 0.2915502

library(SpecsVerification)
setwd("C:/TESIS/ICA/SPDE")
ValidacionCruzadaSPDECO <-read.table("ResultadosFinalesSPDEICA03587793.csv",header=TRUE,sep=",")

MediaPredicciones<-mean(ValidacionCruzadaSPDECO$Predicciones)
DesviacionPrediccion<-sd(ValidacionCruzadaSPDECO$Predicciones) 
MediaObservaciones<-mean(ValidacionCruzadaSPDECO$Observaciones)
GaussCrps(MediaPredicciones,DesviacionPrediccion,MediaObservaciones,ValidacionCruzadaSPDECO$Residuos)

#0.0684591

##PARA SPDE SIN COVARIABLES 

formMeshNoConvex1 <- Y.WQI ~ 0 + y.intercept+f(spatial.field, model = WQISPDEmaterMesNohConvex1)
cv <- inlasloo(dataframe = dataframe,
               long = 'long', lat = 'lat',
               y = 'WQITransBox', ss = ss,
               rad = 0.3779013,alpha=0.05,
               modform = formMeshNoConvex1,
               mesh = MeshNoConvex1, family = 'normal',
               mae = TRUE)
dataa <- data.frame (matrix (unlist (cv), byrow = T), stringsAsFactors = FALSE)
Export(dataa,"ResultadosFinalesSPDEICA04645.csv") # 

ValidacionCruzadaSPDE <-read.table("ResultadosFinalesSPDEICA04614.csv",header=TRUE,sep=",")

MediaPredicciones<-mean(ValidacionCruzadaSPDE$Predicciones)
DesviacionPrediccion<-sd(ValidacionCruzadaSPDE$Predicciones) 
MediaObservaciones<-mean(ValidacionCruzadaSPDE$Observaciones)
GaussCrps(MediaPredicciones,DesviacionPrediccion,MediaObservaciones)
#0.0007416698





