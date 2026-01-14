################################################################################
#             IMPLEMENTATION OF THE EM ALGORITHM                               #
################################################################################

# The objective is the estimation of the background functions, the
# triggering functions and the parameters A and mu_0

# The algorithm consists of the following:
#   1. Estimation of the probabilities
#   2. Estimation of the non-parametric functions with the gaussian kernel
#   3. Estimation of the parameters A and mu_0 using MLE

# This code will run the rest of the iterations

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


setwd("~/EcuadorOpenData/MF_m")
#load("Loop2.RData")
maxit = 50
tolA = 0.001

llik_prev <- llik  # para comparar con la iteración anterior

for (kk in 2:maxit) {
  print("###########################################")
  
  print(paste0("Starting Loop: ",kk))
  
  print("###########################################")
  
  dir.create(paste0("Loop/Loop", kk), showWarnings = FALSE)
  
  lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
  
  bgprobs <- mu * bg.at.events.no.mu / lambda.at.events
  
  weW <- weekly.fun(events$days) * mub.events / lambda.at.events
  weT <- trend.fun(events$days) * mub.events / lambda.at.events
  weM <- monthly.fun(events$days) * mub.events / lambda.at.events
  
  
  ### Weekly
  weekly.base <- seq(0, 7, 0.05)
  new.marks <- as.POSIXlt(as.Date("2014-01-01") + events$days)$wday
  temp <- hist.weighted(new.marks, weW*weights, breaks=weekly.base)
  tband = 1
  weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)
  weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)
  
  png(paste0("Loop/Loop", kk, "/weekly.png"))
  plot(weekly.base, weekly.basevalue, type="l", main=paste0("Weekly at Loop: ", kk), xaxt="n")
  axis(1, at=seq(0, 7, by=1), labels=c("Sunday","Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", ""))
  
  dev.off()
  
  
  ### Monthly
  monthly.base<- seq(0,12,0.05)
  mew.marks <- events$month
  
  temp <- hist.weighted(mew.marks, weM*mweights, breaks=monthly.base)
  tband = 1
  monthly.basevalue <- ker.smooth.fft(temp$mids, temp$density, tband)
  monthly.basevalue <- monthly.basevalue/mean(monthly.basevalue)
  
  png(paste0("Loop/Loop", kk, "/weekly.png"))
  plot(monthly.base, monthly.basevalue, type="l", main=paste0("Weekly at Loop: ", kk), xaxt="n")
  axis(1, at=0:11, labels=month.abb)
  
  dev.off()
  
  
  ### Trend
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  wghs.trend <- trend.fun(events$days)*mub.events/ lambda.at.events
  trend.base <- time.marks
  trend.basevalue <- rep(0, length(time.marks))
  bandwidth <- 150
  
  clusterExport(cl, varlist = c("wghs.trend", "events", "time.marks", "TT", "bandwidth"))
  
  partial_trends <- foreach(i = 1:nrow(events), .combine = '+') %dopar% {
    wi <- wghs.trend[i]
    contrib <- wi * dnorm(events$days[i] - time.marks, 0, bandwidth) /
      (pnorm(TT, events$days[i], bandwidth) - pnorm(0, events$days[i], bandwidth))
    contrib
  }
  
  trend.basevalue <- partial_trends / mean(partial_trends)
  rm(list = c("partial_trends"))
  
  stopCluster(cl)
  
  
  png(paste0("Loop/Loop", kk, "/trend.png"))
  plot(time.marks, trend.basevalue, type = "l", main = paste0("Trend at Loop : ", kk))
  dev.off()
  
  print("Temporal background: Complete")
  
  
  print("Starting spatial background: ")
  
  
  background.basevalue <- matrix(0, nrow=length(background.base$x), 
                                 ncol=length(background.base$y))
  
  for(i in 1:nrow(events)){
    bgsmoother <- dnorm(background.basex, events$coorx[i], events$bandwidth[i])* dnorm(background.basey,events$coory[i], events$bandwidth[i])/events$bg.integral[i]
    background.basevalue <- background.basevalue + bgprobs[i]*bgsmoother
  }
  
  background.basevalue <-  background.basevalue/mean(background.basevalue[background.marks>0])
  print("Spatial background: Complete")
  
  
  # Plot background
  library(viridis)
  
  grid_data <- expand.grid(
    y = background.base$y,
    x = background.base$x
  )
  
  grid_data$basevalue <- as.vector(t(background.basevalue))
  grid_data$marks <- as.vector(t(background.marks))
  
  valid_grid <- grid_data[grid_data$marks == 1, ]
  
  png(paste0("Loop/Loop", kk, "/spatialbackground.png"))
  
  print(
    ggplot() +
      geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +
      geom_tile(data = valid_grid, aes(x = x, y = y, fill = (basevalue))) +
      scale_fill_gradientn(
        colours = c("darkblue", "yellow", "orange", "red"),
        na.value = "white",
        name = expression(μ[b])  # 👈 mu con subíndice b
      )+
      coord_sf() +
      labs(
        title = paste0("Spatial background at Loop: ", kk),
        x = "X", y = "Y"
      ) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14,
          family = "sans"
        )
      )
  )
  
  dev.off()
  
  
  save(background.basevalue, trend.basevalue, weekly.basevalue,
       file = paste0("Loop/Loop", kk, "/BasevalueLoop", kk, ".RData"))
  
  print(paste0("Full background complete Loop: ",kk))
  
  ############# Functions
  
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
  
  
  mub.events <- background.spatial.fun(events$coorx,events$coory)
  bg.at.events.no.mu <- (trend.fun(events$days)*weekly.fun(events$days)*monthly.fun(events$days)*mub.events) ## mu()/mu0 for the observed points
  
  
  mub <- background.spatial.fun(background.basex,background.basey) # spatial for all the grid
  
  bg.at.all.no.mu <- (mean(trend.fun(time.marks)*weekly.fun(time.marks)*monthly.fun(time.marks))*TT*
                        mean(mub*background.marks)*
                        (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))
  
  #save(trend.fun, weekly.fun, background.spatial.fun,mub.events,bg.at.events.no.mu,mub, bg.at.all.no.mu, file = paste0("Loop/Loop", kk, "/FunctionsLoop", kk, ".RData"))
  
  
  print("Starting Exciting: ")
  
  excite.temporal.base <- seq(0,90,0.05)
  excite.spatial.basex <- excite.spatial.base.x %o% rep(1, length(excite.spatial.base.y))
  excite.spatial.basey<- rep(1, length(excite.spatial.base.x)) %o% excite.spatial.base.y
  
  excite.temporal.edge.correction <- rep(0, nrow(events))
  excite.spatial.edge.correction <- rep(0, nrow(events))
  
  pb <- progress_bar$new(
    format = "  Progreso [:bar] :percent | Iteración :current de :total",
    total = nrow(events), clear = FALSE, width = 60
  )
  
  for(i in 1:nrow(events)){
    excite.temporal.edge.correction[i] <-sum(excite.temporal.fun(seq(0, TT-events$days[i], 0.05)+0.6e-5))
    load(paste("~/EcuadorOpenData/MF/City.Excite.Spatial.Marks/crime1-",substr(100000+i,2,6), ".mark", sep="")) # load mark_temp
    mark_temp <- matrix(mark_temp,ncol=length(excite.spatial.base.x),
                        nrow=length(excite.spatial.base.y))
    excite.spatial.edge.correction [i] <- simpson.2D(mark_temp*excite.spatial.basevalue, dx = 200, dy = 200)
    pb$tick()
  }
  
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
                                                x.bandwidth=300,
                                                y.bandwidth=300)*excite.spatial.mark2
  
  excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue, dx = 200, dy = 200)
  
  
  excite.temporal.fun <- approxfun(seq(0, 90, 0.05)+0.6e-12, excite.temporal.basevalue, 
                                   yleft=0, yright=0)
  
  
  excite.spatial.fun <- function(x,y){
    temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
                                    y=excite.spatial.base.y, z=excite.spatial.basevalue),
                           loc=cbind(x=c(x), y=c(y))) 
    temp[is.na(temp)]<- 0
    temp
  }
  
  png(paste0("Loop/Loop", kk, "/excite_temporal.png"))
  plot(excite.temporal.base, excite.temporal.basevalue, type="l", lwd=2, main = paste0("Excite Temporal Loop: ", kk))
  dev.off()
  
  png(paste0("Loop/Loop", kk, "/excite_spatial.png"))
  filled.contour(excite.spatial.base.x,excite.spatial.base.y,
                 (excite.spatial.basevalue),
                 main = paste0("Spatial triggering Loop: ", kk))
  dev.off()
  
  # save(excite.temporal.basevalue, excite.spatial.basevalue, temporal.repetance.fun, spatial.repetance.fun,
  #      excite.temporal.fun, excite.spatial.fun,
  #      file = paste0("Loop/Loop", kk, "/ExciteLoop", kk, ".RData"))
  
  print(paste0("Triggering effect complete Loop: ",kk))
  
  
  # Calculate A and mu
  
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
  
  #save(mytriggers.at.all.no.A, mytriggers.at.events.no.A, file = paste0("Loop/Loop",kk,"/TriggersLoop",kk,".RData"))
  
  ###############################################################################################
  ###                             Loglikelihood                                               ###
  ###############################################################################################
  #bgprobs <- mu * bg.at.events.no.mu / lambda.at.events # \varphi_i :  prob of event i to be events background event
  
  # res.optim <- optim(par=sqrt(c(A, mu)), NegLogLikehood, control=list(trace=6))
  # 
  # 
  # mu <- res.optim$par[1]^2
  # A <- res.optim$par[2]^2
  
  A <- (length(events$bandwidth) - sum(bgprobs))/mytriggers.at.all.no.A
  mu<- (length(events$bandwidth) - A*mytriggers.at.all.no.A)/bg.at.all.no.mu
  
  print(paste("mu=",mu, "A=", A, "at Loop", kk))
  
  lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
  lambda.at.all <-  mu* bg.at.all.no.mu + A * mytriggers.at.all.no.A
  
  llik<- -(- sum(log(lambda.at.events)) + lambda.at.all)
  #ll <- c(ll,llik)
  bgprobs <- mu * bg.at.events.no.mu / lambda.at.events # \varphi_i :  prob of event i to be events background event
  
  save(mu, A, llik, lambda.at.all, lambda.at.events,
       file = paste0("Loop/Loop", kk, "/Results.RData"))
  
  print(paste0("End of Loop: ", kk))
  
  # --- Condiciones de parada ---
  stop_condition <- FALSE
  if (A < 0.001) {
    print("A < 0.001")
    stop_condition <- TRUE
  } else if (!is.na(llik_prev) && abs(llik - llik_prev) < tolA) {
    print(paste0("⚠️ |ΔA| = ", round(abs(llik- llik_prev), 6), " < ", tolA, " → terminando iteraciones."))
    stop_condition <- TRUE
  }
  
  # Guardar todo si se cumple alguna condición de parada
  if (stop_condition || kk == maxit) {
    save.image(file = paste0("Loop","/FullEnvironment_Loop", kk, ".RData"))
    print(paste0("💾 Estado completo guardado en Loop/Loop", kk, "/FullEnvironment_Loop", kk, ".RData"))
    break
  }
  
  llik_prev <- llik
}

# save.image(paste0("Loop",kk,".RData"))
#
# setwd("~/EcuadorOpenData/UsingSubcircuitos2")
# 
# mu_values <- NULL
# A_values <- NULL
# llik_values <- NULL
# iteration <- NULL
# lambda.all <- NULL
# 
# for (i in 1:kk) {
#   file_path <- paste0("Loop/Loop", i, "/Results.RData")
#   env <- new.env()
#   load(file_path, envir = env)
# 
#   mu_values <- c(mu_values, env$mu)
#   A_values <- c(A_values, env$A)
#   llik_values <- c(llik_values, env$llik)
#   iteration <- c(iteration, i)
#   lambda.all <- c(lambda.all,env$lambda.at.all)
# }
# 
# # Summary : for the iterations + 1
# 
# results_summary <- data.frame( Iteration = iteration, mu = mu_values, A =
# A_values, llik = llik_values, lambda = lambda.all )
# 
# print(results_summary)
# 
# plot(llik_values,type = "l",main = "Loglikelihood no month")
# plot((mu_values),type = "l",main = "Mu_0 no month") 
# plot((A_values),type = "l",main = "A no month")
# plot(lambda.all,type = "l",main = "Integral lambda no month")
# 
# 
# 
# ##############################################################
# library(ggplot2)
# library(patchwork)
# 
# df1 <- data.frame(x = seq_along(llik_values), y = llik_values)
# df2 <- data.frame(x = seq_along(mu_values), y = mu_values)
# df3 <- data.frame(x = seq_along(A_values), y = A_values)
# df4 <- data.frame(x = seq_along(lambda.all), y = lambda.all)
# 
# p1 <- ggplot(df1, aes(x, y)) + geom_line(color = "darkblue", size = 1) +
# labs(title = "Loglikelihood", x = "Index", y = "LogLik") + theme_minimal()
# 
# p2 <- ggplot(df2, aes(x, y)) + geom_line(color = "darkgreen", size = 1) +
# labs(title = "Mu_0", x = "Index", y = "Mu_0") + theme_minimal()
# 
# p3 <- ggplot(df3, aes(x, y)) + geom_line(color = "tomato", size = 1) +
# labs(title = "A ", x = "Index", y = "A") + theme_minimal()
# 
# # Combinar con patchwork
# combined_plot <- (p1) / (p2 | p3) +
# plot_annotation(title = paste0("Model Report 3"))
# # 
# combined_plot
# # # ggsave("model_iterations.png", combined_plot, width = 12, height= 8, dpi = 300)