

#####VALIDACION CRUZADA ELECCION DEL MEJOR MODELO ICA #################################

#LECTURA DE ARCHIVOS


library(multDM)
setwd("C:/TESIS/ICA/DM.TEST")

ICAkrigingO<-read.table("ResultadosICATESISKO.csv",header=TRUE,sep=",")
ICAkrigingS<-read.table("ResultadosICATESISKs.csv",header=TRUE,sep=",")
ICACokrigin<-read.table("ResultadosTESISCokrigingICA.csv",header=TRUE,sep=",")

DM.test(ICAkrigingO$ValidacionKo.var1.pred,ICAkrigingS$ValidacionKs.var1.pred,ICAkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(ICAkrigingO$ValidacionKo.var1.pred,ICACokrigin$ValidacionCoKo.WQI.pred,ICAkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")




ICASPDESIN <-read.table("ResultadosFinalesSPDEICASINCOVARIABLES04645.csv",header=TRUE,sep=",")
ICASPDECON<-read.table("ResultadosFinalesSPDEICACOVARIABLES036.csv",header=TRUE,sep=",")

DM.test(ICAkrigingO$ValidacionKo.var1.pred,ICASPDESIN$Predicciones,ICAkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(ICAkrigingS$ValidacionKs.var1.pred,ICASPDESIN$Predicciones,ICAkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(ICACokrigin$ValidacionCoKo.WQI.pred,ICASPDESIN$Predicciones,ICAkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(ICACokrigin$ValidacionCoKo.WQI.pred,ICASPDECON$Predicciones,ICAkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(ICASPDESIN$Predicciones,ICASPDECON$Predicciones,ICACokrigin$ValidacionCoKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")


#####VALIDACION CRUZADA ELECCION DEL MEJOR MODELO PORCENTAJE DE SODIO ###################

##dm.test

library(multDM)
setwd("C:/TESIS/Na/DM.TEST")

NakrigingO<-read.table("ResultadosFinaleskrigingOrdinarioPorSodio.csv",header=TRUE,sep=",")
NakrigingS<-read.table("ResultadosFinaleskrigingSimplePorSodio.csv",header=TRUE,sep=",")
NaCokrigin<-read.table("ResultadosTESISCokrigingNa.csv",header=TRUE,sep=",")

DM.test(NakrigingO$ValidacionKo.var1.pred,NakrigingS$ValidacionKs.var1.pred,NakrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(NakrigingS$ValidacionKs.var1.pred,NaCokrigin$ValidacionCoKoNa.NAPorce.pred,NakrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")

setwd("C:/TESIS/Na/DM.TEST")
NaSPDESIN <-read.table("ResultadosFinalesSPDENASIN1016.csv",header=TRUE,sep=",")
NaSPDECON<-read.table("ResultadosFinalesSPDECON058.csv",header=TRUE,sep=",")

DM.test(NakrigingO$ValidacionKo.var1.pred,NaSPDESIN$Predicciones,NaSPDECON$Observaciones,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(NakrigingS$ValidacionKs.var1.pred,NaSPDESIN$Predicciones,NaSPDECON$Observaciones,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(NaCokrigin$ValidacionCoKo.NAPorce.pred,NaSPDESIN$Predicciones,NakrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(NaCokrigin$ValidacionCoKoNa.NAPorce.pred,NaSPDECON$Predicciones,NakrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(NaSPDESIN$Predicciones,NaSPDECON$Predicciones,NaSPDECON$Observaciones,loss.type="SE",h=1,c=FALSE,H1="more")


#####VALIDACION CRUZADA ELECCION DEL MEJOR MODELO RAS ###################
library(multDM)
setwd("C:/TESIS/RAS/KrigingCokriging")

RASkrigingO<-read.table("ResultadosTESISKORAS.csv",header=TRUE,sep=",")
RASkrigingS<-read.table("ResultadosTESISKsRAS.csv",header=TRUE,sep=",")
RASCokrigin<-read.table("ResultadosTESISCokrigingRAS1.csv",header=TRUE,sep=",")

DM.test(RASkrigingO$ValidacionKo.var1.pred,RASkrigingS$ValidacionKs.var1.pred,RASkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(RASkrigingO$ValidacionKo.var1.pred,RASCokrigin$ValidacionCoKo.RAS.pred,RASkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")

setwd("C:/TESIS/RAS/SPDE")
RASSPDESIN <-read.table("ResultadosFinalesSPDERASSIN124.csv",header=TRUE,sep=",")
RASSPDECON<-read.table("ResultadosFinalesSPDERASCON07474.csv",header=TRUE,sep=",")

DM.test(RASkrigingO$ValidacionKo.var1.pred,RASSPDESIN$Predicciones,RASSPDESIN$Observaciones,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(RASkrigingS$ValidacionKs.var1.pred,RASSPDESIN$Predicciones,RASkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(RASCokrigin$ValidacionCoKo.RAS.pred,RASSPDESIN$Predicciones,RASCokrigin$ValidacionCoKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(RASCokrigin$ValidacionCoKo.RAS.pred,RASSPDECON$Predicciones,RASkrigingO$ValidacionKo.observed,loss.type="SE",h=1,c=FALSE,H1="more")
DM.test(RASSPDESIN$Predicciones,RASSPDECON$Predicciones,RASSPDECON$Observaciones,loss.type="SE",h=1,c=FALSE,H1="more")
