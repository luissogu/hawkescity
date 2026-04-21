# ---- EM ALGORITHM — MAIN LOOP ----
# Iterates the following steps until convergence:
#   1. E-step: compute background probabilities (bgprobs) and weights
#   2. M-step: re-estimate background (weekly, monthly, trend, spatial)
#              and triggering (temporal, spatial) functions
#   3. Update A and mu via closed-form MLE
#   4. Check convergence: stop if |delta(llik)| < tolA or A < 0.001

library(osmdata); library(sf); library(dplyr); library(spatstat)
library(fields);  library(parallel); library(foreach); library(doParallel)
library(progress); library(GEOmap); library(geometry); library(ggplot2)

maxit <- 50
tolA  <- 0.001
llik_prev <- llik

for (kk in 2:maxit) {
  
  message(sprintf("===== Loop %d =====", kk))
  dir.create(here("output", "Loop", paste0("Loop", kk)), showWarnings = FALSE)
  
  # ---- E-STEP ----
  
  lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
  
  bgprobs <- mu * bg.at.events.no.mu / lambda.at.events
  weW     <- weekly.fun(events$days)  * mub.events / lambda.at.events
  weM     <- monthly.fun(events$days) * mub.events / lambda.at.events
  weT     <- trend.fun(events$days)   * mub.events / lambda.at.events
  
  # ---- M-STEP: TEMPORAL BACKGROUND ----
  
  weekly.base      <- seq(0, 7,  0.05)
  monthly.base     <- seq(0, 12, 0.05)
  new.marks        <- as.POSIXlt(as.Date(START_DATE) + events$days)$wday
  
  temp             <- hist.weighted(new.marks, weW * weights, breaks = weekly.base)
  weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, bw = 1)
  weekly.basevalue <- weekly.basevalue / mean(weekly.basevalue)
  
  temp              <- hist.weighted(events$month, weM * mweights, breaks = monthly.base)
  monthly.basevalue <- ker.smooth.fft(temp$mids, temp$density, bw = 1)
  monthly.basevalue <- monthly.basevalue / mean(monthly.basevalue)
  
  # Trend: parallelised kernel smoothing with boundary correction
  wghs.trend  <- weT
  cl          <- makeCluster(N_CORES)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("wghs.trend", "events", "time.marks", "TT", "TREND_BAND"))
  
  trend.basevalue <- foreach(i = seq_len(nrow(events)), .combine = "+") %dopar% {
    wghs.trend[i] * dnorm(events$days[i] - time.marks, 0, TREND_BAND) /
      (pnorm(TT, events$days[i], TREND_BAND) - pnorm(0, events$days[i], TREND_BAND))
  }
  
  stopCluster(cl)
  rm(cl)
  trend.basevalue <- trend.basevalue / mean(trend.basevalue)
  message("Temporal background: complete")
  
  # ---- M-STEP: SPATIAL BACKGROUND ----
  
  background.basevalue <- matrix(0, nrow = length(background.base$x),
                                 ncol = length(background.base$y))
  
  for (i in seq_len(nrow(events))) {
    bgsmoother <- dnorm(background.basex, events$coorx[i], events$bandwidth[i]) *
      dnorm(background.basey, events$coory[i], events$bandwidth[i]) /
      events$bg.integral[i]
    background.basevalue <- background.basevalue + bgprobs[i] * bgsmoother
  }
  
  background.basevalue <- background.basevalue /
    mean(background.basevalue[background.marks > 0])
  message("Spatial background: complete")
  
  # ---- UPDATE INTERPOLATION FUNCTIONS ----
  
  trend.fun  <- approxfun(time.marks, trend.basevalue, yleft = 0, yright = 0)
  
  weekly.fun <- function(x)
    approxfun(weekly.base, weekly.basevalue, yleft = 0, yright = 0)(x - as.integer(x / 7) * 7)
  
  monthly.fun <- function(x)
    approxfun(monthly.base, monthly.basevalue, yleft = 0, yright = 0)(x - as.integer(x / 12) * 12)
  
  background.spatial.fun <- function(x, y)
    interp.surface(obj = list(x = background.base$x,
                              y = background.base$y,
                              z = background.basevalue),
                   loc = cbind(x = c(x), y = c(y)))
  
  mub.events       <- background.spatial.fun(events$coorx, events$coory)
  bg.at.events.no.mu <- trend.fun(events$days) * weekly.fun(events$days) *
    monthly.fun(events$days) * mub.events
  
  mub            <- background.spatial.fun(background.basex, background.basey)
  bg.at.all.no.mu <- mean(trend.fun(time.marks) * weekly.fun(time.marks) *
                            monthly.fun(time.marks)) * TT *
    mean(mub * background.marks) *
    (Xrange[2] - Xrange[1]) * (Yrange[2] - Yrange[1])
  
  # ---- M-STEP: TRIGGERING FUNCTIONS ----
  
  excite.temporal.base <- seq(0, 90, 0.05)
  excite.spatial.basex <- excite.spatial.base.x %o% rep(1, length(excite.spatial.base.y))
  excite.spatial.basey <- rep(1, length(excite.spatial.base.x)) %o% excite.spatial.base.y
  
  excite.temporal.edge.correction <- rep(0, nrow(events))
  excite.spatial.edge.correction  <- rep(0, nrow(events))
  
  for (i in seq_len(nrow(events))) {
    excite.temporal.edge.correction[i] <- sum(
      excite.temporal.fun(seq(0, TT - events$days[i], 0.05) + 0.6e-5)
    )
    fn        <- file.path(path_marks, paste0("crime1-", substr(100000 + i, 2, 6), ".mark"))
    load(fn)  # loads mark_temp
    mark_temp <- matrix(mark_temp, nrow = length(excite.spatial.base.y),
                        ncol = length(excite.spatial.base.x))
    excite.spatial.edge.correction[i] <- simpson.2D(mark_temp * excite.spatial.basevalue,
                                                    dx = 200, dy = 200)
  }
  
  excite.wghs <- A * (
    excite.temporal.fun(events$days[ij.mat[, 1]] - events$days[ij.mat[, 2]]) *
      excite.spatial.fun(events$coorx[ij.mat[, 1]] - events$coorx[ij.mat[, 2]],
                         events$coory[ij.mat[, 1]] - events$coory[ij.mat[, 2]]) /
      lambda.at.events[ij.mat[, 1]]
  )
  
  # Temporal triggering kernel
  excite.temporal.series <- hist.weighted(
    events$days[ij.mat[, 1]] - events$days[ij.mat[, 2]],
    excite.wghs / (excite.temporal.edge.correction[ij.mat[, 2]] *
                     temporal.repetance.fun(events$days[ij.mat[, 1]] - events$days[ij.mat[, 2]])),
    breaks = excite.temporal.base
  )
  
  excite.temporal.basevalue <- ker.smooth.conv(excite.temporal.series$mids,
                                               excite.temporal.series$density[1:1800],
                                               bandwidth = 5)
  excite.temporal.basevalue <- excite.temporal.basevalue /
    simpson(excite.temporal.basevalue, 0.05)
  
  # Spatial triggering kernel
  dis.mat <- cbind(events$coorx[ij.mat[, 1]] - events$coorx[ij.mat[, 2]],
                   events$coory[ij.mat[, 1]] - events$coory[ij.mat[, 2]])
  
  excite.spatial.series <- hist.weighted.2D(
    dis.mat[, 1], dis.mat[, 2],
    excite.wghs / excite.spatial.edge.correction[ij.mat[, 2]],
    x.breaks = excite.spatial.base.x,
    y.breaks  = excite.spatial.base.y
  )
  
  # Circular mask: zero out kernel beyond EXCITE_S_MAX metres
  excite.spatial.mask <- (excite.spatial.basex^2 + excite.spatial.basey^2) < EXCITE_S_MAX^2
  
  temp <- spatial.repetance.fun(
    excite.spatial.series$x.mids %o% rep(1, length(excite.spatial.series$y.mids)),
    rep(1, length(excite.spatial.series$x.mids)) %o% excite.spatial.series$y.mids
  )
  
  excite.spatial.basevalue <- ker.smooth.2D.fft(
    excite.spatial.series$x.mids,
    excite.spatial.series$y.mids,
    excite.spatial.series$density / temp,
    x.bandwidth = 300,
    y.bandwidth = 300
  ) * excite.spatial.mask
  
  excite.spatial.basevalue <- excite.spatial.basevalue /
    simpson.2D(excite.spatial.basevalue, dx = 200, dy = 200)
  
  excite.temporal.fun <- approxfun(seq(0, 90, 0.05) + 0.6e-12,
                                   excite.temporal.basevalue, yleft = 0, yright = 0)
  
  excite.spatial.fun <- function(x, y) {
    temp <- interp.surface(
      obj = list(x = excite.spatial.base.x,
                 y = excite.spatial.base.y,
                 z = excite.spatial.basevalue),
      loc = cbind(x = c(x), y = c(y))
    )
    temp[is.na(temp)] <- 0
    temp
  }
  message("Triggering functions: complete")
  
  # ---- PARALLEL COMPUTATION OF TRIGGER INTEGRALS ----
  
  cl <- makeCluster(N_CORES)
  clusterExport(cl, varlist = c(
    "events", "background.basex", "background.basey", "background.marks",
    "excite.temporal.basevalue", "excite.spatial.base.x",
    "excite.spatial.base.y", "excite.spatial.basevalue",
    "time.marks", "TT", "Xrange", "Yrange"
  ), envir = environment())
  
  clusterEvalQ(cl, {
    library(fields)
    excite.temporal.fun <- approxfun(seq(0, 90, 0.05) + 0.6e-12,
                                     excite.temporal.basevalue, yleft = 0, yright = 0)
    excite.spatial.fun <- function(x, y) {
      temp <- interp.surface(
        obj = list(x = excite.spatial.base.x,
                   y = excite.spatial.base.y,
                   z = excite.spatial.basevalue),
        loc = cbind(x = c(x), y = c(y))
      )
      temp[is.na(temp)] <- 0
      temp
    }
    calculate_excite <- function(i) {
      at.events <- excite.temporal.fun(events$days - events$days[i]) *
        excite.spatial.fun(events$coorx - events$coorx[i],
                           events$coory - events$coory[i])
      temp    <- excite.spatial.fun(background.basex - events$coorx[i],
                                    background.basey - events$coory[i])
      at.all  <- mean(excite.temporal.fun(time.marks - events$days[i])) * TT *
        mean((background.marks > 0) * temp) *
        (Xrange[2] - Xrange[1]) * (Yrange[2] - Yrange[1])
      list(event = at.events, all = at.all)
    }
  })
  
  results <- parLapply(cl, seq_len(nrow(events)), function(i) calculate_excite(i))
  stopCluster(cl)
  
  mytriggers.at.events.no.A <- Reduce("+", lapply(results, `[[`, "event"))
  mytriggers.at.all.no.A    <- sum(sapply(results, `[[`, "all"))
  rm(results)
  
  # ---- PARAMETER UPDATE (CLOSED-FORM MLE) ----
  
  A  <- (nrow(events) - sum(bgprobs)) / mytriggers.at.all.no.A
  mu <- (nrow(events) - A * mytriggers.at.all.no.A) / bg.at.all.no.mu
  
  lambda.at.events <- mu * bg.at.events.no.mu + A * mytriggers.at.events.no.A
  lambda.at.all    <- mu * bg.at.all.no.mu    + A * mytriggers.at.all.no.A
  llik             <- -(- sum(log(lambda.at.events)) + lambda.at.all)
  bgprobs          <- mu * bg.at.events.no.mu / lambda.at.events
  
  message(sprintf("Loop %d complete — mu = %.4f, A = %.4f, llik = %.4f", kk, mu, A, llik))
  
  save(mu, A, llik, lambda.at.all, lambda.at.events,
       file = here("output", "Loop", paste0("Loop", kk), "Results.RData"))
  
  # ---- CONVERGENCE CHECK ----
  
  stop_condition <- A < 0.001 || (!is.na(llik_prev) && abs(llik - llik_prev) < tolA)
  
  if (stop_condition || kk == maxit) {
    save.image(file = here("output", "Loop", paste0("FullEnvironment_Loop", kk, ".RData")))
    message(sprintf("Converged at loop %d. Full environment saved.", kk))
    break
  }
  
  llik_prev <- llik
}