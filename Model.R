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
library(units)

setwd("~/hawkescity")
source("functionsSemiHawkes.R") # load functions
load("guayaquil.RData") # load geometry
load("crime_gyq.RData") # load geocodes + times
load("events.RData")

##################################################################################
# We make sure our points are within the polygon, and that they have the same CRS

pt <- new_sf
pt <- st_transform(pt, crs = 32717)
city <- st_transform(city, st_crs(pt))

plot(st_geometry(city))
plot(st_geometry(pt), add = TRUE, col = "red", pch = 16,cex = 0.1)

###############################################################################
# Object events will be used for all the code.
# Date "2014-01-01" can be changed depending of the data
pt <- pt |> filter(as.Date(fecha_infraccion) >= as.Date("2014-01-01"))

events <- data.frame(
  coorx = st_coordinates(pt)[,1],
  coory = st_coordinates(pt)[,2],
  days = as.numeric(as.Date(pt$fecha_infraccion) - as.Date("2014-01-01")),
  month = as.numeric(format(pt$fecha_infraccion, "%m")) - 1
)

w = as.owin(city)
inside <- spatstat.geom::inside.owin(x = events$coorx,
                                     y = events$coory,
                                     w = w)
plot(w)
points(events[!inside, ]$coorx, events[!inside, ]$coory,
       col = "red", pch = 16,cex = 0.5)

events <- events[inside,]

events <- events |> filter(days < as.numeric(difftime(as.Date("2025-01-01"), as.Date("2014-01-01"), units = "days")))

events <- events[order(events$days),]

TT<- range(events$days)[2] + 1 # final day
time.marks <- seq (0 , TT ,0.05) # time.marks

# Integral for the spatial background

dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}

events$bandwidth <- rep(250,nrow(events)) # 100 metros

# Choose of the bandwidth for everypoint:
# This code assigns the kernel to be =  if the 5th more close point is too close to the point

for(i in 1:nrow(events)){
  temp <- dist.squared(events$coorx[i], events$coory[i], events$coorx[-i], events$coory[-i])
  temp2 <- sqrt(sort(temp[temp>=150^2])[5])
  if(events$bandwidth[i] <= temp2) events$bandwidth[i] = temp2
}

events$bg.integral <- rep(0, nrow(events)) # Value of integral initial = 0

w = as.owin(city)
plot(w)
points(events$coorx, events$coory)


gaussian_density <- function(coords, mean, sigma) {
  detS <- det(sigma)
  invS <- solve(sigma)
  diffs <- t(t(coords) - mean)
  exp_term <- rowSums((diffs %*% invS) * diffs)
  (1 / (2 * pi * sqrt(detS))) * exp(-0.5 * exp_term)
}

pts <- spatstat.random::rpoint(3000, win = w)
coords <- cbind(pts$x, pts$y)
areaW <- spatstat.geom::area.owin(w)

events$bg.integral <- sapply(1:nrow(events), function(i) {
  mean_i <- c(events$coorx[i], events$coory[i])
  sigma_i <- diag(events$bandwidth[i]^2, 2)
  mean(gaussian_density(coords, mean = mean_i, sigma = sigma_i)) * areaW
})


events$bg.integral<-ifelse(events$bg.integral > 0.99 , 1,events$bg.integral)


inside <- spatstat.geom::inside.owin(x = events$coorx,
                                     y = events$coory,
                                     w = w)
plot(w)
points(events[!inside, ]$coorx, events[!inside, ]$coory,
       col = "red", pch = 16,cex = 0.2)

events <- events[inside,]

save(events,file = "events.RData")


#####################################################################################################
##                                  BASE VALUES                                                    ##
#####################################################################################################

### weekly base value

weekly.base <- seq(0, 7, 0.05)
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
tband = 1#original
weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)

weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)
plot(weekly.base, weekly.basevalue, type="l",main = "Weekly",sub = "0 : Sunday")


plot(weekly.base, weekly.basevalue, type="l", main="Weekly Initial City", xaxt="n")
axis(1, at=seq(0, 7, by=1), labels=c("Sunday","Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", ""))


############################  Monthly    ###################################################

monthly.base <- seq(0, 12, 0.05)
mew.marks <- events$days - as.integer(events$days/12)*12 
mweights <- 1/ (as.integer(TT/12) + (events$days - as.integer(events$days/12)*12 > TT - as.integer(TT/12)*12)) 
temp <- hist.weighted(mew.marks, mweights, breaks=monthly.base)
monthly.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.9)
monthly.basevalue <- monthly.basevalue/mean(monthly.basevalue)
plot(monthly.base, monthly.basevalue, type="l")

plot(monthly.base, monthly.basevalue, type="l",xaxt="n",main = "Monthly Initial City")
axis(1, at=0:11, labels=month.abb)

############################  Trend     ###################################################

trend.base <- seq(0, TT, 0.05)
trend.basevalue <- rep(1, length(time.marks))
bandwidth <- 150

# partial_trends <- rep(0, length(time.marks))
# 
# for (i in 1:nrow(events)) {
#   contrib <- dnorm(events$days[i] - time.marks, 0, bandwidth) /
#     (pnorm(TT, events$days[i], bandwidth) - pnorm(0, events$days[i], bandwidth))
#   partial_trends <- partial_trends + contrib
# }
# 
# trend.basevalue <- partial_trends / mean(partial_trends)

plot(trend.basevalue,type="l")
save(trend.basevalue,weekly.basevalue,monthly.basevalue,file = "Loop/Loop0/Basevaluetemporal0.RData")

############################ Spatial base value ###############################################
bbox <- st_bbox(city)

Xrange <- c(bbox["xmin"], bbox["xmax"])
Yrange <- c(bbox["ymin"], bbox["ymax"])

# This will depend on the size of the grid we want to use
# In this case we use 100m because the data is in UTM (m)
background.base <- list(x = seq(Xrange[1], Xrange[2], 150), 
                        y = seq(Yrange[1], Yrange[2], 150))

background.basex <- background.base$x %o% rep(1, length(background.base$y))

background.basey <- rep(1, length(background.base$x))%o% background.base$y

background.basevalue <- matrix(0, nrow=length(background.base$x), 
                               ncol=length(background.base$y))

background.base <- list(
  x = seq(Xrange[1], Xrange[2], 150), 
  y = seq(Yrange[1], Yrange[2], 150)
)

# Crear el grid como puntos sf directamente
grid_points <- expand.grid(
  x = seq(Xrange[1], Xrange[2], 150),
  y = seq(Yrange[1], Yrange[2], 150)
)

grid_sf <- st_as_sf(grid_points, coords = c("x", "y"), crs = st_crs(city))

inside_total <- lengths(st_intersects(grid_sf, city)) > 0

background.marks <- matrix(as.integer(inside_total),
                           ncol = length(seq(Yrange[1], Yrange[2], 150)))


save(background.marks, file = "background.marks.RData")
###################################################################################################
for(i in 1:nrow(events)){
  bgsmoother <- dnorm(background.basex, events$coorx[i], events$bandwidth[i])* dnorm(background.basey,events$coory[i], events$bandwidth[i])/events$bg.integral[i]
  background.basevalue <- background.basevalue + bgsmoother
  cat(i,"\n")
}

background.basevalue <-  background.basevalue/mean(background.basevalue[background.marks>0])


save(background.basevalue,trend.basevalue, weekly.basevalue,monthly.basevalue ,events, file = "Loop/Loop0/CityInitialBaseValue.RData")

############## Nice Plot for the spatial background ######################

grid_data <- expand.grid(
  y = background.base$y,
  x = background.base$x
)

grid_data$basevalue <- as.vector(t(background.basevalue))
grid_data$basevalue <- as.vector(t(background.basevalue))

grid_data$marks <- as.vector(t(background.marks))

valid_grid <- grid_data[grid_data$marks == 1, ]

ggplot() +
  geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +         # Polígono
  geom_tile(data = valid_grid, aes(x = x, y = y, fill = log(basevalue))) +  # Suavizadores
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

excite.spatial.base.x <- seq(-4000, 4000, 200)
excite.spatial.base.y <- seq(-4000, 4000, 200)
excite.spatial.basex= excite.spatial.base.x%o% rep(1, length(excite.spatial.base.y))
excite.spatial.basey= t(excite.spatial.base.y%o% rep(1, length(excite.spatial.base.x)))

excite.spatial.basevalue = matrix(1/((abs(excite.spatial.base.x %o% rep(1, length(excite.spatial.base.y)))^2)/(40000^2)
                                     + (abs(rep(1, length(excite.spatial.base.x))%o% excite.spatial.base.y)^2)/(40000^2)+ 0.01), ncol=length(excite.spatial.base.x), nrow=length(excite.spatial.base.y))

excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue, dx = 200, dy = 200)

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
  approxfun(monthly.base,monthly.basevalue, yleft=0, yright=0)(x- as.integer(x/12)*12)
}

background.spatial.fun <- function(x,y) (interp.surface(obj=list(x=background.base$x, 
                                                                 y=background.base$y, z=background.basevalue),
                                                        loc=cbind(x=c(x), y=c(y))))

excite.temporal.fun <- approxfun(seq(0,90,0.05)+0.1e-12, excite.temporal.basevalue, 
                                 yleft=0, yright=0)

excite.spatial.fun <- function(x, y) {
  temp <- interp.surface(
    obj = list(x = excite.spatial.base.x, y = excite.spatial.base.y, z = excite.spatial.basevalue),
    loc = cbind(x = c(x), y = c(y))
  )
  temp[is.na(temp)] <- 0
  temp
}


mub.events <- background.spatial.fun(events$coorx,events$coory)
bg.at.events.no.mu <- (trend.fun(events$days)*weekly.fun(events$days)*monthly.fun(events$days)*mub.events) ## mu()/mu0 for the observed points


mub <- background.spatial.fun(background.basex,background.basey) # spatial for all the grid

bg.at.all.no.mu <- as.numeric(mean(trend.fun(time.marks)*weekly.fun(time.marks)*monthly.fun(time.marks))*TT*
                                mean(mub*background.marks)*
                                (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))

save(trend.fun, weekly.fun,monthly.fun,background.spatial.fun,excite.temporal.fun,mub.events,bg.at.events.no.mu,mub, bg.at.all.no.mu, file = "Loop/Loop0/CityFunctionsLoop0.RData")

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

save(mytriggers.at.all.no.A, mytriggers.at.events.no.A, file = "Loop/Loop0/CityInitialTriggersLoop0.RData")

################################################################################################

#########################################################################################

A = 0.6  # Initial guesses for the expected number of triggered events

mu = (dim(events)[1] - A*mytriggers.at.all.no.A)/bg.at.all.no.mu

lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
lambda.at.all <-  mu* bg.at.all.no.mu + A * mytriggers.at.all.no.A

bgprobs <- mu * bg.at.events.no.mu / lambda.at.events
llik<- -(- sum(log(lambda.at.events)) + lambda.at.all)


save(mu,A,llik, lambda.at.all,lambda.at.events, file = "Loop/Loop0/Results.RData")

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
###################################################################################
## Weekly term

weekly.base <- seq(0, 7, 0.05)
new.marks <- as.POSIXlt(as.Date("2014-01-01") + events$days)$wday
temp <- hist.weighted(new.marks, weW*weights, breaks=weekly.base)
tband = 1 #original
weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)

weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)
plot(weekly.base, weekly.basevalue, type="l",main = "Weekly at Loop: 1",sub = "0 : Sunday")


# Monthly term

monthly.base<- seq(0,12,0.05)
mew.marks <- events$month

temp <- hist.weighted(mew.marks, weM*mweights, breaks=monthly.base)
tband = 1
monthly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)
monthly.basevalue <- monthly.basevalue/mean(monthly.basevalue)
plot(monthly.base, monthly.basevalue, type="l",xaxt="n",main = "Monthly at Loop : 1")
axis(1, at=0:11, labels=month.abb)


###################################################################################
#### Trend
# without parallel computing ~ 15 minutes

wghs.trend <- trend.fun(events$days)* mub.events/ lambda.at.events

trend.base <- time.marks
trend.basevalue <- rep(0, length(time.marks))
bandwidth <-trend_band

for(i in 1:nrow(events)){
  wi <- wghs.trend[i]
  trend.basevalue <- trend.basevalue + wi * dnorm(events$days[i] - time.marks, 0, bandwidth) /
    (pnorm(TT, events$days[i], bandwidth) - pnorm(0, events$days[i], bandwidth))
  cat("it: ",i*100/nrow(events))
  cat("\n")
}
trend.basevalue <- trend.basevalue / mean(trend.basevalue)


plot(time.marks, trend.basevalue, type = "l", main = paste0("Trend at Loop : ", 1))

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


save(background.basevalue,trend.basevalue,weekly.basevalue,monthly.basevalue,file= "Loop/Loop1/CityBasevalue.RData")
#####################################################################################

trend.fun <-  approxfun(time.marks, trend.basevalue, yleft=0, yright=0)

weekly.fun <- function(x){
  approxfun(weekly.base, weekly.basevalue,             
            yleft=0, yright=0)(x- as.integer(x/7)*7)
}

monthly.fun <- function(x){
  approxfun(monthly.base,monthly.basevalue, yleft=0, yright=0)(x- as.integer(x/12)*12)
}

background.spatial.fun <- function(x,y) (interp.surface(obj=list(x=background.base$x, 
                                                                 y=background.base$y, z=background.basevalue),
                                                        loc=cbind(x=c(x), y=c(y))))

# excite.temporal.fun <- approxfun(seq(0,90,0.05)+0.1e-12, excite.temporal.basevalue, 
#                                  yleft=0, yright=0)
# 
# excite.spatial.fun <- function(x,y){
#   temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
#                                   y=excite.spatial.base.y, z=excite.spatial.basevalue),
#                          loc=cbind(x=c(x), y=c(y))) 
#   temp[is.na(temp)]<- 0
#   temp
# }

mub.events <- background.spatial.fun(events$coorx,events$coory)
bg.at.events.no.mu <- (trend.fun(events$days)*weekly.fun(events$days)*monthly.fun(events$days)*mub.events) ## mu()/mu0 for the observed points


mub <- background.spatial.fun(background.basex,background.basey) # spatial for all the grid

bg.at.all.no.mu <- as.numeric(mean(trend.fun(time.marks)*weekly.fun(time.marks)*monthly.fun(time.marks))*TT*
                                mean(mub*background.marks)*
                                (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))


###############################################################################
city_polygons <- list()

for (i in seq_along(st_geometry(city))) {
  geom <- st_geometry(city)[[i]]  # MULTIPOLYGON
  # geom es una lista de POLYGONs
  for (polygon in geom) {
    # polygon es una lista de anillos (matrices)
    for (ring_coords in polygon) {
      city_polygons[[length(city_polygons) + 1]] <- list(
        x = ring_coords[,1],
        y = ring_coords[,2]
      )
    }
  }
}

dir.create("City.Excite.Spatial.Marks", showWarnings = FALSE)

# Loop por evento
for (i in 1:nrow(events)) {
  fn <- file.path("City.Excite.Spatial.Marks",
                  paste0("crime1-", substr(100000 + i, 2, 6), ".mark"))
  cat("i: ",i," out of ",nrow(events),"\n")
  if (!file.exists(fn)) {
    # Coordenadas desplazadas
    grid_x <- excite.spatial.basex + events$coorx[i]
    grid_y <- excite.spatial.basey + events$coory[i]

    # Crear coordenadas combinadas
    pts_x <- as.vector(grid_x)
    pts_y <- as.vector(grid_y)

    # Inicializar inclusión
    inside_logical <- rep(FALSE, length(pts_x))

    # Verificar inclusión en cada polígono de la ciudad
    for (poly in city_polygons) {
      inside_logical <- inside_logical | inpoly(pts_x, pts_y, poly)
    }

    # Convertir a matriz (según dimensiones de la grilla base)
    mark_temp <- matrix(inside_logical,
                        ncol = length(excite.spatial.basex),
                        byrow = TRUE)

    save(mark_temp, file = fn)
  }
}

# 
# #####################################################################################################
# 
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
  
  load(paste("City.Excite.Spatial.Marks/crime1-",substr(100000+i,2,6), ".mark", sep="")) # load mark_temp
  
  mark_temp <- matrix(mark_temp, nrow = nrow(spatial.repetance), ncol = ncol(spatial.repetance))
  
  spatial.repetance <- spatial.repetance + mark_temp
  
  excite.spatial.edge.correction [i] <- simpson.2D(mark_temp*excite.spatial.basevalue, dx = 200, dy = 200)
  
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
    events$days[ij.mat[,1]] <= events$days[ij.mat[,2]] + 90.0 &  # Maximum 90 days after
    abs(events$coorx[ij.mat[,1]] - events$coorx[ij.mat[,2]]) <= 4000 &  # |x_i - x_j| <= 4000
    abs(events$coory[ij.mat[,1]] - events$coory[ij.mat[,2]]) <= 4000,   # |y_i - y_j| <= 4000
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

excite.spatial.mark2 <- (excite.spatial.basex ^2 + excite.spatial.basey^2 < 4200^2) # (?)

temp<- spatial.repetance.fun(excite.spatial.series$x.mids%o%rep(1, length(excite.spatial.series$y.mids)),
                             rep(1, length(excite.spatial.series$y.mids))%o%excite.spatial.series$x.mids)


excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,
                                              excite.spatial.series$y.mids,
                                              excite.spatial.series$density/temp,
                                              #excite.spatial.series$density,
                                              x.bandwidth=250,
                                              y.bandwidth=250)*excite.spatial.mark2

excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue, dx = 200, dy = 200)

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
     spatial.repetance,temporal.repetance,file = "Loop/Loop1/ExciteLoop1.RData")

library(parallel)

nprocs <- 4
cl <- makeCluster(nprocs)
clusterExport(cl, varlist = c(
  "events", "background.basex", "background.basey", "background.marks",
  "excite.temporal.basevalue", "excite.spatial.base.x", "excite.spatial.base.y", "excite.spatial.basevalue",
  "time.marks", "TT", "Xrange", "Yrange"
), envir = environment())

clusterEvalQ(cl, {
  library(fields)
  
  excite.temporal.fun <- approxfun(seq(0, 90, 0.05) + 0.6e-12, excite.temporal.basevalue, yleft = 0, yright = 0)
  
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

rm(list = c("results"))

save(mytriggers.at.all.no.A, mytriggers.at.events.no.A, file = paste0("Loop/Loop",1,"/TriggersLoop",1,".RData"))

#########################################################################################
A <- (length(events$bandwidth) - sum(bgprobs))/mytriggers.at.all.no.A
mu<- (length(events$bandwidth) - A*mytriggers.at.all.no.A)/bg.at.all.no.mu

lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
lambda.at.all <-  mu* bg.at.all.no.mu + A * mytriggers.at.all.no.A


llik<- -(- sum(log(lambda.at.events)) + lambda.at.all)
bgprobs <- mu * bg.at.events.no.mu / lambda.at.events # \varphi_i :  prob of event i to be events background event

save(mu,A,llik, lambda.at.all,lambda.at.events, file = paste0("Loop/Loop",1,"/Results.RData"))



