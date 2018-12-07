#Preparar paquetes y ambiente
rm(list=ls(all=T))

suppressMessages(suppressWarnings(library(maptools)))
suppressMessages(suppressWarnings(library(rgdal)))
suppressMessages(suppressWarnings(library(spdep)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(lattice)))
suppressMessages(suppressWarnings(library(gstat)))
suppressMessages(suppressWarnings(library(georob)))
suppressMessages(suppressWarnings(library(xtable)))


#Leer archivo
#setwd("C:/Users/rguadamuz/Dropbox/Maestría/Estadística UCR/2018 II/SP1649 - Topicos de Estadistica Espacial Aplicada/Proyectos/Proyecto Final")
setwd("C:/Users/Anibal/Dropbox/Proyectos/Proyecto Final")

poly <- readOGR("Tibas.shp")
poly$VEHIC <- as.numeric(poly$VEHIC)

#Visualización inicial
(spplot(poly, c("INGRESO_HO")))

#Modelo lineal no espacial
mod.lm <- lm(INGRESO_HO~NOMB_UGED+Shape_Area+VIVIENDAS+HABITANTES+EDAD_PROM+DENS_HABIT+Trab_MCant, data=poly)
summary(mod.lm)

#Se dejan sólo las covariables significativas (Encontradas con Selección hacia atrás, hacia adelante con p-value, AIC y BIC)
mod.lm <- lm(INGRESO_HO ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant, data=poly)
summary(mod.lm)
print(xtable(summary(mod.lm)$coefficients, type = "latex"), file = "lm.table.tex")
#Supuestos, outliers e influencia en el otro código

poly$lmresid <- summary(mod.lm)$residuals
(spplot(poly, c("lmresid")))
#Ya no son tan intensas las diferencias, pero todavía se notan patrones espaciales.

#Test de Moran's I
list.queen <- poly2nb(poly, queen=TRUE, snap = 20)
W <- nb2listw(list.queen, style="W", zero.policy=FALSE); W
plot(W,coordinates(poly), col="red")
lm.morantest(mod.lm, W)
#Moran I statistic standard deviate = 8.1942, p-value < 2.2e-16
#A pesar de todas las variables, todavía hay algo espacial que no quedó explicado en el LM

#Otros tipos de W para descartar que no sea solamente por la definición Queen
list.queen <- poly2nb(poly, queen=FALSE, snap = 20)
W <- nb2listw(list.queen, style="W", zero.policy=FALSE)
lm.morantest(mod.lm, W)

W <- nb2listw(list.queen, style="B", zero.policy=FALSE)
lm.morantest(mod.lm, W)

W <- nb2listw(list.queen, style="C", zero.policy=FALSE)
lm.morantest(mod.lm, W)

W <- nb2listw(list.queen, style="U", zero.policy=FALSE)
lm.morantest(mod.lm, W)

list.queen <- poly2nb(poly, queen=TRUE, snap = 20)
W <- nb2listw(list.queen, style="minmax", zero.policy=FALSE)
lm.morantest(mod.lm, W)

W <- nb2listw(list.queen, style="S", zero.policy=FALSE)
lm.morantest(mod.lm, W)
#Todas los estilos de matriz rechazan la independencia espacial.

#Se corre modelo espacial (SAR) para determinar el lambda (rho)
mod.sar <- spautolm(INGRESO_HO ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant, data=poly, listw=W)
mod.sar$lambda
summary(mod.sar)
print(xtable(summary(mod.sar)$Coef, type = "latex"), file = "sar.table.tex")

#Lambda: 0.47647 LR test value: 54.06 p-value: 1.9451e-13 
#AIC: 12948
#hay estructura espacial sifnificativa

poly$trend <- mod.sar$fit$signal_trend
(spplot(poly, c("trend")))
poly$stochastic <- mod.sar$fit$signal_stochastic
(spplot(poly, c("stochastic")))
#Se nota el trend muy plano de color, mientras que el estocástico con diferencias.

#Se puede hacer corrección de pesos en la regresión.
#Corrección por área
mod.sar.w1 <- spautolm(INGRESO_HO ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant, data=poly, listw=W, weights=Shape_Area)
summary(mod.sar.w1)
poly$trend <- mod.sar.w1$fit$signal_trend
(spplot(poly, c("trend")))
poly$stochastic <- mod.sar.w1$fit$signal_stochastic
(spplot(poly, c("stochastic")))

#Corrección por EDAD_PROM
mod.sar.w2 <- spautolm(INGRESO_HO ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant, data=poly, listw=W, weights=EDAD_PROM)
summary(mod.sar.w2)
poly$trend <- mod.sar.w2$fit$signal_trend
(spplot(poly, c("trend")))
poly$stochastic <- mod.sar.w2$fit$signal_stochastic
(spplot(poly, c("stochastic")))
#Sigue habiendo correlación espacial. Los pesos no ayudan en mucho, así que se elige el modelo sin pesos (correrlo de nuevo para sobreescribir variables)

#Modelo CAR
mod.car <- spautolm(INGRESO_HO ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant, data=poly, listw=W, family = "CAR", method = "Matrix")
mod.car$lambda
summary(mod.car)
print(xtable(summary(mod.car)$Coef, type = "latex"), file = "car.table.tex")
#Lambda: 0.039513 LR test value: 2.0458 p-value: 0.15263  
#AIC: 13001
#En este caso ya no existe estructura espacial, y los coeficientes de CAR y SAR son bastante diferentes
poly$trend <- mod.car$fit$signal_trend
(spplot(poly, c("trend")))
poly$stochastic <- mod.car$fit$signal_stochastic
(spplot(poly, c("stochastic")))

#Comparando los AIC el SAR da mejores resultados (12948 vs 13001)

#Geoestadística
data = cbind(coordinates(poly), poly@data)
colnames(data)[1:2] = c("X", "Y")
head(data)
df = data

coordinates(df) <- c("X", "Y")

mod.lm <- lm(log(INGRESO_HO) ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant, data=data)
summary(mod.lm)

# v = variogram(log(INGRESO_HO) ~1, df)
# plot(v)
# fit.variogram(v, vgm(1,"Sph",1200,1))
# 
# plot(variogram(log(INGRESO_HO) ~ 1, df, alpha = c(0, 45, 90, 135)))
# 
# f <- log(INGRESO_HO) ~ NOMB_UGED + VIVIENDAS + HABITANTES + Trab_MCant
# v = variogram(log(INGRESO_HO) ~1, df)
# vt = variogram(f, df)
# plot(v)
# plot(vt)
# v.fit <- fit.variogram(v, vgm(1, "Sph", 1200, 1))
# vt.fit <- fit.variogram(vt, vgm(1, "Exp", 400, 1))
# plot(r.sv <- sample.variogram(mod.lm$residuals, locations=data[, c("X","Y")],
#                               lag.dist.def=12, max.lag=400,
#                               estimator="matheron"), type="l", main="Semivariograma de residuos")
# r.sv.spher <- fit.variogram.model(r.sv, variogram.mode="RMspheric",
#                                   param=c(variance=4.2, nugget=0.01, scale=170), method="BFGS")
# lines(r.sv.spher, col="red")
# 
# plot(variogram(f, df, map = TRUE, cutoff = 1000, width = 100))
# 
# gwls = gstat(NULL, "RES", f, df, model = vt.fit, set = list(gls= 1))

suppressMessages(suppressWarnings(library(georob)))

jpeg(filename = "C:/Users/Anibal/Dropbox/Proyectos/Proyecto Final/Imagen informe/Semivariogramas.jpeg", 
     height = 323, width = 623, units = "px")
par(mfrow = c(1,2), cex = 0.9)
plot(r.sv <- sample.variogram(residuals(mod.lm), locations=data[, c("X","Y")],
                              lag.dist.def=100, max.lag=2000,
                              estimator="matheron"), type="l", , ylab = "Semivariancia", xlab = "Distancia de rezago",
     main="Semivariograma de residuos")

plot(sample.variogram(residuals(mod.lm), locations=data[, c("X","Y")],
                      lag.dist.def=100, max.lag=2000, xy.angle.def=c(0, 22.5, 67.5, 112.5, 157.5, 180),
                      estimator="matheron"), type="l",legend.pos = "bottomright",
     main="Semivariograma de residuos \n con dirección", ylab = "Semivariancia", xlab = "Distancia de rezago")
par(mfrow = c(1,1))
dev.off()

plot(r.sv <- sample.variogram(residuals(mod.lm), locations=data[, c("X","Y")],
                              lag.dist.def=100, max.lag=2000,
                              estimator="matheron"), type="l", , ylab = "Semivariancia", xlab = "Distancia de rezago",
     main="Semivariograma de residuos")

lines(r.sv.spher <- fit.variogram.model(r.sv, variogram.mode="RMspheric",
                                        param=c(variance=0.1, nugget=0.03, scale=1000), method="BFGS"), col = "red")

summary(r.sv.spher)


r.georob.m0.spher.reml <- georob(f, data, locations=~X+Y,
                                 variogram.model="RMspheric", param=c(variance=0.009757, nugget=0.0395, scale=1493),
                                 tuning.psi=1000)
summary(r.georob.m0.spher.reml)
extractAIC(r.georob.m0.spher.reml)

waldtest(r.georob.m0.spher.reml, .~.-NOMB_UGED)

#sí hay diferencia estadística entre los modelos, se prefiere el que tiene todas las variables.
ajustados = exp(r.georob.m0.spher.reml$fitted.values)
plot(ajustados)

poly$ajustados_georob = ajustados
poly$resid_georob = data$INGRESO_HO-ajustados

jpeg(filename = "C:/Users/Anibal/Dropbox/Proyectos/Proyecto Final/Imagen informe/result_georob.jpeg", 
     height = 323, width = 623, units = "px")
par(mfrow = c(1,2), cex = 0.9)

g1 = (spplot(poly, c("ajustados_georob"), xlab = "Valores estimados"))
g2 = (spplot(poly, c("resid_georob"), xlab = "Residuos"))

dev.off()
                        
library(gridExtra)
grid.arrange(g1,g2, ncol  = 2)
