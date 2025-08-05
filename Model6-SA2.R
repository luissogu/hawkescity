################################################################################
#             IMPLEMENTATION OF THE EM ALGORITHM                               #
################################################################################

# The objective is the estimation of the background functions, the
# triggering functions and the parameters A and mu_0

# The algorithm consists of the following:
#   1. Estimation of the probabilities
#   2. Estimation of the non-parametric functions with the gaussian kernel
#   3. Estimation of the parameters A and mu_0 using MLE

# In order to achieve this, we need the initial probabilities.
# In this script we calculate the initial probs and also the first iteration of the algorithm

# The code BigLoop does the rest of the iterations

library(osmdata)
library(sf)
library(dplyr)
library(polyCub)
library(spatstat)
library(stopp)
library(viridis)
library(foreach)
library(doParallel)
library('GEOmap')
library('geometry')
library(ggplot2)
library(fields)
library(parallel)
library(progress)



setwd("~/EcuadorOpenData/City")
source("functionsSemiHawkes.R") # load functions

load("guayaquil.RData") # load geometry
load("crime_gyq.RData") # load geocodes + times


##################################################################################
# We make sure our points are within the polygon, and that they have the same CRS

pt <- sf_puntos

pt <- st_transform(pt, crs = 32717)
city <- st_transform(city, st_crs(pt))



plot(st_geometry(city))
plot(st_geometry(pt), add = TRUE, col = "red", pch = 16)

###############################################################################
# Object events will be used for all the code.
# Date "2014-01-01" can be changed depending of the data

events <- data.frame(
  coorx = st_coordinates(pt)[,1],
  coory = st_coordinates(pt)[,2],
  days = as.numeric(as.Date(pt$fecha_infraccion) - as.Date("2014-01-01")),
  month = as.numeric(format(pt$fecha_infraccion, "%m")) - 1
)

events <- events[order(events$days),]

TT<- range(events$days)[2] + 1 # final day
time.marks <- seq (0 , TT ,0.005) # time.marks 

###############################################################################
# Plot of the original Zhuang and Mateu (2019) paper


colores <- viridis(80)[cut(events$days, breaks = 80, labels = FALSE)]

par(mfrow = c(2, 2), mar = c(5, 5, 2, 2), cex = 1.2)

plot(st_geometry(city), xlab = "X", ylab = "Y", main = "(a)", border = "grey")
points(events$coorx, events$coory, col = colores, pch = 16, cex = 0.5)

plot(events$days, events$coory, xlab = "Days", ylab = "Y", col = colores, cex = 0.6, main = "(b)")

plot(events$coorx, events$days, xlab = "X", ylab = "Days", col = colores, cex = 0.6, main = "(c)")

plot(c(0, events$days, TT), c(0:(nrow(events)+1)), type = "s", lwd = 2,
     xlab = "Days", ylab = "Cum. Freq.", main = "(d)")

par(mfrow = c(1, 1))


##############################################################################
# Integral for the spatial background

dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}

events$bandwidth <- rep(2000,nrow(events)) # 100 metros

# Choose of the bandwidth for everypoint: 
# This code assigns the kernel to be =  if the 5th more close point is too close to the point

for(i in 1:nrow(events)){
  temp <- dist.squared(events$coorx[i], events$coory[i], events$coorx[-i], events$coory[-i])
  temp2 <- sqrt(sort(temp[temp>=1000^2])[5])
  if(events$bandwidth[i] <= temp2) events$bandwidth[i] = temp2 
}

events$bg.integral <- rep(0, nrow(events)) # Value of integral initial = 0

w = as.owin(city)
plot(w)
points(events$coorx, events$coory)


gaussian_density <- function(coords, mean, sigma) {
  mvtnorm::dmvnorm(coords, mean = mean, sigma = sigma)
}

for (i in 1:nrow(events)) {
  mean_i <- c(events$coorx[i], events$coory[i])
  sigma_i <- diag(events$bandwidth[i]^2, 2)
  
  events$bg.integral[i] <- polyCub.SV(
    polyregion = w,
    f = function(coords) gaussian_density(coords, mean = mean_i, sigma = sigma_i) # Edge correction, this should be = 1. 
  )
  
  print(paste(i, 'of', nrow(events), ":", events$bg.integral[i]))
}

events$bg.integral<-ifelse(events$bg.integral > 0.9999 , 1,events$bg.integral)

save(events,file = "~/EcuadorOpenData/City/Model5/events.RData")

#####################################################################################################
##                                  BASE VALUES                                                    ##
#####################################################################################################

### monthly base value (?)

monthly<- seq(0,12,0.05)
mew.marks <- events$month

# Count full years
n_total_days <- TT
n_full_years <- floor(n_total_days / 365.25)

# Check for additional months in the final year
last_day <- as.Date("2014-01-01") + n_total_days
last_month <- as.numeric(format(last_day, "%m"))
#  n_full_years veces cada mes
#  +1 vez los meses entre enero y last_month

# Count number of observed months
month_counts <- rep(n_full_years, 12)
month_counts[1:last_month] <- month_counts[1:last_month] + 1

# weights for the histogram
mweights <- 1 / month_counts[events$month + 1]  # +1 porque R indexa desde 1

temp <- hist.weighted(mew.marks, mweights, breaks=monthly)
tband = 1
monthly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)

monthly.basevalue <- monthly.basevalue/mean(monthly.basevalue)

plot(monthly, monthly.basevalue, type="l",xaxt="n",main = "Monthly Initial City")
axis(1, at=0:11, labels=month.abb)


### weekly base value

weekly.base <- seq(0, 7, 0.005)
new.marks <- as.POSIXlt(as.Date("2014-01-01") + events$days)$wday # Changed this to have  0: Sunday

# Total de días
n_total_days <- TT
  
# full weeks
n_full_weeks <- floor(n_total_days / 7)
  
# additonal days
n_extra_days <- n_total_days %% 7
start_wday <- as.POSIXlt(as.Date("2014-01-01"))$wday  # 0: Sunday, ..., 6: Saturday
extra_days_indices <- (start_wday + 0:(n_extra_days - 1)) %% 7  # Esto da los wday de los días extra
  
# Start with same weight
wday_counts <- rep(n_full_weeks, 7)
  
# Aumentar 1 a los días extras
wday_counts[extra_days_indices + 1] <- wday_counts[extra_days_indices + 1] + 1  # +1 por índice R
  
new.marks <- as.POSIXlt(as.Date("2014-01-01") + events$days)$wday
  
weights <- 1 / wday_counts[new.marks + 1]

temp <- hist.weighted(new.marks, weights, breaks=weekly.base)
tband = 0.5 #original
weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)

weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)
plot(weekly.base, weekly.basevalue, type="l",main = "Weekly",sub = "0 : Sunday")


plot(weekly.base, weekly.basevalue, type="l", main="Weekly Initial City", xaxt="n")
axis(1, at=seq(0, 7, by=1), labels=c("Sunday","Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", ""))

############################  Trend     ###################################################

cl <- makeCluster(4) # WARNING : number of cores used 
registerDoParallel(cl)

trend.base <- seq(0, TT, 0.005)
trend.basevalue <- rep(0, length(time.marks))
bandwidth <- 90

clusterExport(cl, varlist = c("events", "time.marks", "TT", "bandwidth"))

partial_trends <- foreach(i = 1:nrow(events), .combine = '+') %dopar% {
  contrib <- dnorm(events$days[i] - time.marks, 0, bandwidth) /
    (pnorm(TT, events$days[i], bandwidth) - pnorm(0, events$days[i], bandwidth))
  contrib
}

trend.basevalue <- partial_trends / mean(partial_trends)

stopCluster(cl)

plot(time.marks, trend.basevalue, type="l",main = "Trend Initial City")


save(trend.basevalue,weekly.basevalue,monthly.basevalue, file = "~/EcuadorOpenData/City/Model6/Loop0/Basevaluetemporal0.RData")
load("~/EcuadorOpenData/City/Model6/Loop0/Basevaluetemporal0.RData")

############################ Spatial base value ###############################################

bbox <- st_bbox(city)

Xrange <- c(bbox["xmin"], bbox["xmax"])
Yrange <- c(bbox["ymin"], bbox["ymax"])

# This will depend on the size of the grid we want to use
# In this case we use 100m because the data is in UTM (m)
background.base <- list(x = seq(Xrange[1], Xrange[2], 100), 
                        y = seq(Yrange[1], Yrange[2], 100))

background.basex <- background.base$x %o% rep(1, length(background.base$y))

background.basey <- rep(1, length(background.base$x))%o% background.base$y

background.basevalue <- matrix(0, nrow=length(background.base$x), 
                               ncol=length(background.base$y))

background.base <- list(
  x = seq(Xrange[1], Xrange[2], 100), 
  y = seq(Yrange[1], Yrange[2], 100)
)

x_coords <- as.vector(background.base$x %o% rep(1, length(background.base$y)))
y_coords <- as.vector(rep(1, length(background.base$x)) %o% background.base$y)

# Matriz de inclusión inicial vacía
inside_total <- rep(FALSE, length(x_coords))

# Procesar cada geometría (puede haber múltiples polígonos)
# This will check if the cell falls within the polygon
# actually, this is not necessary and can be changed for a simple inpoly function in case the geometry is simply connected
for (i in seq_along(st_geometry(city))) {
  coords <- st_coordinates(st_geometry(city)[[i]])
  outer_rings <- unique(coords[coords[, "L2"] == 1, "L1"])
  
  for (ring in outer_rings) {
    ring_coords <- coords[coords[, "L1"] == ring & coords[, "L2"] == 1, ]
    result <- inpoly(x_coords, y_coords, list(x = ring_coords[,1], y = ring_coords[,2]))
    inside_total <- inside_total | result
  }
}

background.marks <- matrix(as.integer(inside_total), 
                           ncol = length(background.base$y))

save(background.marks, file = "~/EcuadorOpenData/City/Model6/background.marks.RData")
load("~/EcuadorOpenData/City/Model6/background.marks.RData")

###################################################################################################
dir.create("~/EcuadorOpenData/City/Model5/Background.Smoothers.City", showWarnings = FALSE, recursive = TRUE)

for(i in 1:nrow(events)){
  fn <- file.path("~/EcuadorOpenData/City/Model5/Background.Smoothers.City/", paste0("bgsmoother-", i, ".val"))
  if(!file.exists(fn)){
    bgsmoother <- dnorm(background.basex, events$coorx[i], events$bandwidth[i])* dnorm(background.basey,events$coory[i], events$bandwidth[i])/events$bg.integral[i]
    save(bgsmoother,file=fn)
  } else{
    load(fn)
  }   
  background.basevalue <- background.basevalue + bgsmoother
}

background.basevalue <-  background.basevalue/mean(background.basevalue[background.marks>0])


save(background.basevalue,trend.basevalue, weekly.basevalue, events,monthly.basevalue, file = "~/EcuadorOpenData/City/Model6/Loop0/CityInitialBaseValue.RData")
load("~/EcuadorOpenData/City/Model6/Loop0/CityInitialBaseValue.RData")


############## Nice Plot for the spatial background ######################

grid_data <- expand.grid(
  y = background.base$y,
  x = background.base$x
)

grid_data$basevalue <- as.vector(t(background.basevalue))
grid_data$marks <- as.vector(t(background.marks))

valid_grid <- grid_data[grid_data$marks == 1, ]

ggplot() +
  geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +         # Polígono
  geom_tile(data = valid_grid, aes(x = x, y = y, fill = log10(basevalue))) +  # Suavizadores
  # geom_point(data = events, aes(x = coorx, y = coory), colour = "lightblue",
  #            shape = 16, size = 0.5, alpha = 0.8) +                               # Eventos
  scale_fill_viridis(option = "magma", name = "(Background rate)") +
  coord_sf() +
  labs(
    title = "Events + Initial background",
    x = "X",
    y = "Y"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )


############################################################################################
# Self exciting components

# Spatial and temporal guesses

excite.temporal.basevalue <- (0.1 + seq(0,90,0.05)/10)^(-1.03) ## why -1.03?
excite.temporal.basevalue <- excite.temporal.basevalue/simpson(excite.temporal.basevalue,0.05)

plot(seq(0,90,0.05), excite.temporal.basevalue,type="l",main = "Temporal triggering function",xlab = "Time [days]")

excite.spatial.base.x <- seq(-45000, 45000, 450)
excite.spatial.base.y <- seq(-45000, 45000, 450)
excite.spatial.basex= excite.spatial.base.x%o% rep(1, length(excite.spatial.base.y))
excite.spatial.basey= t(excite.spatial.base.y%o% rep(1, length(excite.spatial.base.x)))

excite.spatial.basevalue = matrix(1/((abs(excite.spatial.base.x %o% rep(1, length(excite.spatial.base.y)))^2)/(45000^2)
                                     + (abs(rep(1, length(excite.spatial.base.x))%o% excite.spatial.base.y)^2)/(45000^2)+ 0.01), ncol=length(excite.spatial.base.x), nrow=length(excite.spatial.base.y))

excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue, dx = 450, dy = 450)

filled.contour(excite.spatial.base.x,excite.spatial.base.y,
               (excite.spatial.basevalue))

#################################################################################################
# Evaluated lambda 

trend.fun <-  approxfun(time.marks, trend.basevalue, yleft=0, yright=0)

weekly.fun <- function(x){
  approxfun(weekly.base, weekly.basevalue,             
            yleft=0, yright=0)(x- as.integer(x/7)*7)
}

monthly.fun <- function(x){
  approxfun(monthly,monthly.basevalue, yleft=0, yright=0)(x- as.integer(x/12)*12)
}

background.spatial.fun <- function(x,y) (interp.surface(obj=list(x=background.base$x, 
                                                                 y=background.base$y, z=background.basevalue),
                                                        loc=cbind(x=c(x), y=c(y))))

excite.temporal.fun <- approxfun(seq(0,90,0.05)+0.1e-12, excite.temporal.basevalue, 
                                 yleft=0, yright=0)

excite.spatial.fun <- function(x,y){
  temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
                                  y=excite.spatial.base.y, z=excite.spatial.basevalue),
                         loc=cbind(x=c(x), y=c(y))) 
  temp[is.na(temp)]<- 0
  temp
}

mub.events <- background.spatial.fun(events$coorx,events$coory)
bg.at.events.no.mu <- (trend.fun(events$days)*weekly.fun(events$days)*monthly.fun(events$days)*mub.events) ## mu()/mu0 for the observed points


mub <- background.spatial.fun(background.basex,background.basey) # spatial for all the grid

bg.at.all.no.mu <- as.numeric(mean(trend.fun(time.marks)*weekly.fun(time.marks) * monthly.fun(time.marks))*TT*
                                mean(mub*background.marks)*
                                (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))

save(trend.fun,monthly.fun, weekly.fun, background.spatial.fun,excite.spatial.fun,excite.temporal.fun,mub.events,bg.at.events.no.mu,mub, bg.at.all.no.mu, file = "~/EcuadorOpenData/City/Model6-SA2/Loop0/CityFunctionsLoop0.RData")
load("~/EcuadorOpenData/City/Model6-SA2/Loop0/CityFunctionsLoop0.RData")

triggers.at.events.no.events <- rep(0, nrow(events))
triggers.at.all.no.events <- 0

#################################################################################################
nprocs <- 4
cl <- makeCluster(nprocs)

clusterExport(cl, varlist = c(
  "events", "background.basex", "background.basey", "background.marks",
  "excite.temporal.basevalue", "excite.spatial.base.x", "excite.spatial.base.y", "excite.spatial.basevalue",
  "time.marks", "TT", "Xrange", "Yrange"
), envir = environment())

clusterEvalQ(cl, {
  library(fields)
  
  excite.temporal.fun <- approxfun(seq(0,90,0.05) + 0.1e-12, excite.temporal.basevalue, yleft = 0, yright = 0)
  
  excite.spatial.fun <- function(x, y) {
    temp <- interp.surface(
      obj = list(x = excite.spatial.base.x, y = excite.spatial.base.y, z = excite.spatial.basevalue),
      loc = cbind(x = c(x), y = c(y))
    )
    temp[is.na(temp)] <- 0
    temp
  }
  
  calculate_excite <- function(i) {
    mytriggers.at.events.no.A <- excite.temporal.fun(events$days - events$days[i]) *
      excite.spatial.fun(events$coorx - events$coorx[i], events$coory - events$coory[i])
    
    temp <- excite.spatial.fun(background.basex - events$coorx[i],
                               background.basey - events$coory[i])
    
    mytriggers.at.all.no.A <- mean(excite.temporal.fun(time.marks - events$days[i])) * TT *
      mean((background.marks > 0) * temp) *
      (Xrange[2] - Xrange[1]) * (Yrange[2] - Yrange[1])
    
    return(list(event = mytriggers.at.events.no.A, all = mytriggers.at.all.no.A))
  }
})

results <- parLapply(cl, 1:nrow(events), function(i) calculate_excite(i))
stopCluster(cl)

mytriggers.at.events.no.A <- Reduce("+", lapply(results, function(res) res$event))

mytriggers.at.all.no.A <- sum(sapply(results, function(res) res$all))

rm(list = c("results")) # this is normally a big object, better to remove it

save(mytriggers.at.all.no.A, mytriggers.at.events.no.A, file = "~/EcuadorOpenData/City/Model6-SA2/Loop0/CityInitialTriggersLoop0.RData")

load("~/EcuadorOpenData/City/Model6-SA2/Loop0/CityInitialTriggersLoop0.RData")

################################################################################################

#########################################################################################

A = 0.5  # Initial guesses for the expected number of triggered events

mu = 0.00000000075

lambda.at.events_Loop0 <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
lambda.at.all_Loop0 <-  mu* bg.at.all.no.mu + A * mytriggers.at.all.no.A

bgprobs <- mu * bg.at.events.no.mu / lambda.at.events_Loop0

NegLogLikehood <- function(x){ # this is a function of (mu0,A)
  mu <- x[1]^2
  A<- x[2]^2
  
  lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
  lambda.at.all <-  mu* bg.at.all.no.mu + A * mytriggers.at.all.no.A
  
  - sum(log(lambda.at.events)) + lambda.at.all # Approx? of the -log
} # 

A <- (length(events$bandwidth) - sum(bgprobs))/mytriggers.at.all.no.A
mu<- (length(events$bandwidth) - A*mytriggers.at.all.no.A)/bg.at.all.no.mu

print(paste("mu=",mu, "A=", A, "at Loop", 0))

lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
lambda.at.all <-  mu* bg.at.all.no.mu + A * mytriggers.at.all.no.A

bgprobs <- mu * bg.at.events.no.mu / lambda.at.events

llik<- -(- sum(log(lambda.at.events)) + lambda.at.all)


save(mu,A,llik, lambda.at.all,lambda.at.events, file = "~/EcuadorOpenData/City/Model6-SA2/Loop0/Results.RData")
load("~/EcuadorOpenData/City/Model6-SA2/Loop0/Results.RData")


###########################################################################################
# This is the first iteration of the EM Algoritm
#########################################################################################


# next step: 
# Get the weights:  w_T, w_w, w_m

#### lambda
lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A

# probs baackground
bgprobs <- mu * bg.at.events.no.mu / lambda.at.events

# weekly term
weW<- weekly.fun(events$days)* mub.events / lambda.at.events

# monthly term
weM <- monthly.fun(events$days)* mub.events / lambda.at.events

# trend tern
weT <- trend.fun(events$days)* mub.events / lambda.at.events

#weights<- data.frame(weW = weW , weM = weM , weT = weT)

## Weekly term!

weekly.base <- seq(0, 7, 0.005)
new.marks <- as.POSIXlt(as.Date("2014-01-01") + events$days)$wday


temp <- hist.weighted(new.marks, weW*weights, breaks=weekly.base)
tband = 0.5 #original
weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)

weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)
plot(weekly.base, weekly.basevalue, type="l",main = "Weekly at Loop: 1",sub = "0 : Sunday")


# monthly term!

monthly<- seq(0,12,0.05)
mew.marks <- events$month

temp <- hist.weighted(mew.marks, weM*mweights, breaks=monthly)
tband = 1

monthly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)

monthly.basevalue <- monthly.basevalue/mean(monthly.basevalue)

plot(monthly, monthly.basevalue, type="l",xaxt="n",main = "Monthly at Loop : 1")
axis(1, at=0:11, labels=month.abb)

###################################################################################
#### Trend
# without parallel computing ~ 15 minutes

wghs.trend <- trend.fun(events$days) / lambda.at.events

trend.base <- time.marks
trend.basevalue <- rep(0, length(time.marks))
bandwidth <- 90

for(i in 1:nrow(events)){
  wi <- wghs.trend[i]
  trend.basevalue <- trend.basevalue + wi * dnorm(events$days[i] - time.marks, 0, bandwidth) /
    (pnorm(TT, events$days[i], bandwidth) - pnorm(0, events$days[i], bandwidth))
  cat("it: ",i)
  cat("\n")
}
trend.basevalue <- trend.basevalue / mean(trend.basevalue)


png(paste0("~/EcuadorOpenData/City/Model6-SA2/Loop", 1, "/trend.png"))
plot(time.marks, trend.basevalue, type = "l", main = paste0("Trend at Loop : ", 1))
dev.off()

#####################################################################################

# Spatial trend:

background.basevalue <- matrix(0, nrow=length(background.base$x), 
                               ncol=length(background.base$y))

pb <- progress_bar$new(
  format = "  Progreso [:bar] :percent | Iteración :current de :total",
  total = nrow(events), clear = FALSE, width = 60
)

for(i in 1:nrow(events)){
  bgsmoother <- dnorm(background.basex, events$coorx[i], events$bandwidth[i])* dnorm(background.basey,events$coory[i], events$bandwidth[i])/events$bg.integral[i]
  background.basevalue <- background.basevalue + bgprobs[i]*bgsmoother
  pb$tick()  # Actualiza la barra de progreso
}

background.basevalue <-  background.basevalue/mean(background.basevalue[background.marks>0])



save(background.basevalue,trend.basevalue,monthly.basevalue,weekly.basevalue,file= "~/EcuadorOpenData/City/Model6-SA2/Loop1/CityBasevalue.RData")
load("~/EcuadorOpenData/City/Model6-SA2/Loop1/CityBasevalue.RData")


# Triggering functions


# probs excited

# This code calculate the marks for the (s - s_i ) of the excite spatial base
# Asumiendo que city y pt ya tienen el CRS 32717 y están alineados


#Procesar anillos exteriores una vez
city_polygons <- list()

for (i in seq_along(st_geometry(city))) {
  coords <- st_coordinates(st_geometry(city)[[i]])
  outer_rings <- unique(coords[coords[, "L2"] == 1, "L1"])
  
  for (ring in outer_rings) {
    ring_coords <- coords[coords[, "L1"] == ring & coords[, "L2"] == 1, ]
    city_polygons[[length(city_polygons) + 1]] <- list(
      x = ring_coords[, 1],
      y = ring_coords[, 2]
    )
  }
}

# dir.create("~/EcuadorOpenData/City/Model6-SA2/City.Excite.Spatial.Marks", showWarnings = FALSE)
# 
# # Loop por evento
# for (i in 1:nrow(events)) {
#   fn <- file.path("~/EcuadorOpenData/City/Model6-SA2/City.Excite.Spatial.Marks",
#                   paste0("crime1-", substr(100000 + i, 2, 6), ".mark"))
#   
#   if (!file.exists(fn)) {
#     # Coordenadas desplazadas
#     grid_x <- excite.spatial.basex + events$coorx[i]
#     grid_y <- excite.spatial.basey + events$coory[i]
#     
#     # Crear coordenadas combinadas
#     pts_x <- as.vector(grid_x)
#     pts_y <- as.vector(grid_y)
#     
#     # Inicializar inclusión
#     inside_logical <- rep(FALSE, length(pts_x))
#     
#     # Verificar inclusión en cada polígono de la ciudad
#     for (poly in city_polygons) {
#       inside_logical <- inside_logical | inpoly(pts_x, pts_y, poly)
#     }
#     
#     # Convertir a matriz (según dimensiones de la grilla base)
#     mark_temp <- matrix(inside_logical,
#                         ncol = length(excite.spatial.basex),
#                         byrow = TRUE)
#     
#     save(mark_temp, file = fn)
#   }
# }

#####################################################################################################

excite.temporal.base <- seq(0,90,0.05)
excite.spatial.basex <- excite.spatial.base.x %o% rep(1, length(excite.spatial.base.y))
excite.spatial.basey<- rep(1, length(excite.spatial.base.x)) %o% excite.spatial.base.y

temporal.repetance <- 0 * excite.temporal.base
spatial.repetance <- matrix(0, ncol=length(excite.spatial.base.x),
                            nrow=length(excite.spatial.base.y))

excite.temporal.edge.correction <- rep(0, nrow(events))
excite.spatial.edge.correction <- rep(0, nrow(events))


pb <- progress_bar$new(
  format = "  Progreso [:bar] :percent | Iteración :current de :total",
  total = nrow(events), clear = FALSE, width = 60
)

for(i in 1:nrow(events)){
  excite.temporal.edge.correction[i] <-sum(excite.temporal.fun(seq(0, TT-events$days[i], 0.05)+0.6e-5))
  
  temporal.repetance [excite.temporal.base < TT -events$days[i]] <- temporal.repetance [excite.temporal.base < TT -events$days[i]]+1
  
  load(paste("~/EcuadorOpenData/City/Model6-SA1/City.Excite.Spatial.Marks/crime1-",substr(100000+i,2,6), ".mark", sep="")) # load mark_temp
  
  mark_temp <- matrix(mark_temp, nrow = nrow(spatial.repetance), ncol = ncol(spatial.repetance))
  
  spatial.repetance <- spatial.repetance + mark_temp
  
  excite.spatial.edge.correction [i] <- simpson.2D(mark_temp*excite.spatial.basevalue, dx = 450, dy = 450)
  
  pb$tick()
}

temporal.repetance.fun <- approxfun(excite.temporal.base, temporal.repetance, 
                                    yleft=1, yright=1)

spatial.repetance.fun <- function(x,y){
  temp <- interp.surface(obj=list(x=excite.spatial.base.x, y=excite.spatial.base.y, z= spatial.repetance), loc=cbind(x=c(x), y=c(y)))
  temp[is.na(temp)]<-0
  temp
}


temp.mat <- (1:nrow(events))%o% rep(1, nrow(events))
ij.mat <- cbind(c(t(temp.mat)), c(temp.mat))
rm(list = c("temp.mat"))


ij.mat <- ij.mat[
  events$days[ij.mat[,1]] > events$days[ij.mat[,2]] &  # events j after i
    events$days[ij.mat[,1]] <= events$days[ij.mat[,2]] + 90.0 &  # Maximum 15 days after
    abs(events$coorx[ij.mat[,1]] - events$coorx[ij.mat[,2]]) <= 45000 &  # |x_i - x_j| <= 60000
    abs(events$coory[ij.mat[,1]] - events$coory[ij.mat[,2]]) <= 45000,   # |y_i - y_j| <= 60000
]

excite.wghs <- A*(excite.temporal.fun(events$days[ij.mat[,1]]-events$days[ij.mat[,2]])
                  *excite.spatial.fun(events$coorx[ij.mat[,1]]-events$coorx[ij.mat[,2]],events$coory[ij.mat[,1]]-events$coory [ij.mat[,2]])/ lambda.at.events[ij.mat[,1]])


excite.temporal.series <- hist.weighted(events$days[ij.mat[,1]]-events$days[ij.mat[,2]], # days difference
                                        excite.wghs/(excite.temporal.edge.correction[ij.mat[,2]]*temporal.repetance.fun(events$days[ij.mat[,1]]-events$days[ij.mat[,2]])),
                                        breaks=excite.temporal.base)

excite.temporal.basevalue<- ker.smooth.conv(excite.temporal.series$mids, excite.temporal.series$density[1:1800], bandwidth=5)

excite.temporal.basevalue <- excite.temporal.basevalue/simpson(excite.temporal.basevalue, 0.05) 
plot(excite.temporal.base, excite.temporal.basevalue,type="l", lwd=2)


#################################        Spatial          ####################################################
dis.mat <- cbind(events$coorx[ij.mat[,1]] - events$coorx[ij.mat[,2]], # x_i - x_j
                 events$coory[ij.mat[,1]] - events$coory[ij.mat[,2]]) #  y_i - y_j

excite.spatial.series <-  hist.weighted.2D(dis.mat[,1],
                                           dis.mat[,2],
                                           excite.wghs/(excite.spatial.edge.correction[ij.mat[,2]]),
                                           x.breaks= excite.spatial.base.x,
                                           y.breaks= excite.spatial.base.y)

excite.spatial.mark2 <- (excite.spatial.basex ^2 + excite.spatial.basey^2 < 30500^2) # (?)

temp<- spatial.repetance.fun(excite.spatial.series$x.mids%o%rep(1, length(excite.spatial.series$y.mids)),
                             rep(1, length(excite.spatial.series$y.mids))%o%excite.spatial.series$x.mids)


excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,
                                              excite.spatial.series$y.mids,
                                              #excite.spatial.series$density/temp,
                                              excite.spatial.series$density,
                                              x.bandwidth=2000,
                                              y.bandwidth=2000)*excite.spatial.mark2

excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue, dx = 450, dy = 450)

filled.contour(excite.spatial.base.x,excite.spatial.base.y,
               (excite.spatial.basevalue), main='Spatial triggering Loop 1')


excite.temporal.fun <- approxfun(seq(0, 90, 0.05)+0.6e-12, excite.temporal.basevalue, 
                                 yleft=0, yright=0)

excite.spatial.fun <- function(x,y){
  temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
                                  y=excite.spatial.base.y, z=excite.spatial.basevalue),
                         loc=cbind(x=c(x), y=c(y))) 
  temp[is.na(temp)]<- 0
  temp
}

save(excite.temporal.basevalue,excite.spatial.fun,excite.temporal.fun,
     excite.spatial.basevalue, temporal.repetance.fun, spatial.repetance.fun,
     spatial.repetance,temporal.repetance,file = "~/EcuadorOpenData/City/Model6-SA2/Loop1/ExciteLoop1.RData")


