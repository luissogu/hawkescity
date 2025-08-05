# ----------------------- Prediction ------------------

##-------------------------------planar--------------------------
#par(mfrow=c(1,4),  mai = c(0.5, 0.5, 0.3, 0.1))

### Cumulative frequencies of crime vs original
plot(c(0,events$days), 0:nrow(events), type="s",xlab="Time (days)", ylab="Sum. Freq.", lwd=2, xaxs="i", yaxs="i",xlim=c(0,TT), ylim=c(0, nrow(events)),main="(a)",cex.main=1.5)
abline(0, nrow(events)/TT,lwd=1.5, lty=2)  


### Cumulative frequencies of crime vs average tranformed (95% confidence)
bgrates.at.all.locations.no.mu <- ((trend.fun(time.marks)*weekly.fun(time.marks))*
                                     mean(background.spatial.fun(background.basex,background.basey)*background.marks)*
                                     (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))

triggers.at.all.locations.no.A <- 0

for(i in 1:nrow(events)){
  if(as.integer(i/100)*100==i)print(i)
  
  temp <- excite.spatial.fun(background.basex-events$coorx[i],
                             background.basey-events$coory[i])
  
  triggers.at.all.locations.no.A <- (triggers.at.all.locations.no.A+ (excite.temporal.fun(time.marks-events$days[i]))*
                                       mean((background.marks>0)*temp)* (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))
}

lambda.at.all.locations <-  mu* bgrates.at.all.locations.no.mu + A * triggers.at.all.locations.no.A


plot(cumsum(lambda.at.all.locations)*0.005, stepfun(events$days, (0:nrow(events)))(time.marks), type="s",lwd=2,
     xlab="Transformed time", ylab="Cum. Freq.", xaxs="i", yaxs="i",xlim=c(0,nrow(events)), ylim=c(0, nrow(events)),main="Transformed time series",
     cex.main=1.5)
abline(0, 1,lwd=1.5, lty=2) 
points(0:nrow(events), nrow(events)* qbeta(0.975, 0:nrow(events)+1, nrow(events)-(0:nrow(events))+1),type="l", lty=2)
points(0:nrow(events), nrow(events)* qbeta(0.025, 0:nrow(events)+1, nrow(events)-(0:nrow(events))+1),type="l", lty=2)

##########################################################################################################################


library(spatstat)
library(spatstat.utils)
library(spatstat.geom)
library(spatstat.core)
library(akima)       # para `interp`
library(fields)      # para `interp.surface`
library(RColorBrewer)

# --- Paso 1: Crear imagen de intensidad a partir de los puntos -------------------
# events: data frame con columnas coorx, coory, background.total.net

stuff <- list(x = events$coorx, 
              y = events$coory, 
              z = mu*bg.at.events.no.mu)

# Convertimos la intensidad a imagen tipo `im` de spatstat
Z <- as.im(stuff)                 # resolución por defecto
Z <- as.im(Z, dimyx = 256)        # ajustamos la resolución si queremos mayor detalle

# --- Paso 2: Crear máscara del polígono donde se definirá la intensidad ----------
# city: objeto sf con el polígono de la ciudad

# Convertimos el polígono `sf` a un objeto `owin`
city_owin <- as.owin(city)

# Lo transformamos en imagen binaria (TRUE dentro del polígono, FALSE fuera)
city_mask <- as.im(city_owin, dimyx = 128)

# --- Paso 3: Extraer coordenadas de los píxeles válidos dentro del polígono -----
xx <- raster.x(city_mask)
yy <- raster.y(city_mask)
xx <- as.vector(xx)
yy <- as.vector(yy)

# Creamos un patrón de puntos con los centros de píxel
pixelcentres <- ppp(xx, yy, window = as.rectangle(city_mask))
pixelcentres <- unique.ppp(pixelcentres)  # por si acaso hay duplicados

# --- Paso 4: Interpolación en una grilla fina ------------------------------------
# Creamos interpolación en grilla (aquí 10000 x 10000 puede ser excesivo para algunos casos)
interp_grid <- interp(x = events$coorx, 
                      y = events$coory, 
                      z = lambda.at.events,
                      xo = seq(min(xx), max(xx), length = 200), 
                      yo = seq(min(yy), max(yy), length = 200),
                      extrap = TRUE,
                      duplicate = "mean")

# Función que interpola valores a partir de la grilla
interpuchi <- function(x, y) {
  interp.surface(obj = list(x = interp_grid$x, 
                            y = interp_grid$y, 
                            z = interp_grid$z),
                 loc = cbind(x = x, y = y))
}

# Evaluamos los valores interpolados en cada centro de píxel válido
values <- interpuchi(xx, yy)

# --- Paso 5: Construir imagen final con valores interpolados ---------------------
# Inicializamos la imagen final con ceros fuera del polígono
Z <- city_mask
Z[pixelcentres] <- values  # asignamos los valores interpolados a los píxeles válidos

# --- Paso 6: Estimar intensidad observada y calcular residuos ---------------------
# ppp0: patrón de puntos observado (tipo `ppp`) sobre el mismo dominio

lambda <- density(ppp0, sigma = bw.diggle(ppp0), dimyx = 128)  # intensidad suavizada observada

# Normalizamos por duración o cantidad total si es necesario
residuos <- eval.im(lambda / TT - Z)

# --- Paso 7: Visualizar los residuos ---------------------------------------------
par(mfrow = c(1, 2), mai = c(0.5, 0.5, 0.4, 0.2))

# Mapa de residuos
plot(residuos, 
     main = expression('Smoothed raw residuals'), 
     box = FALSE, 
     ribside = 'bottom', 
     col = brewer.pal(n = 11, name = "RdYlBu"), 
     axes = TRUE, 
     cex.axis = 0.7)

# Histograma de residuos
cols_cov <- colorRampPalette(brewer.pal(10, 'RdYlBu'))(20)
hist(residuos, 
     main = '', 
     xlab = 'Residuals', 
     ylab = 'Frequency', 
     col = cols_cov, 
     border = cols_cov, 
     breaks = 40, 
     xlim = range(residuos, finite = TRUE))
par(mfrow = c(1,1))

# Rojo: valores cercano a 0
# Azul: Mas eventos de los estimados

# --- Paso 8: Diagnóstico final (opcional) -----------------------------------------
range(residuos)