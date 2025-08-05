library(tidyr)
library(dplyr)

all_days <- data.frame(days = 0:TT)

tau <- events |>
  group_by(days) |>
  summarise(count = n(), .groups = "drop") |>
  right_join(all_days, by = "days") |>
  replace_na(list(count = 0)) |>
  arrange(days) |>
  mutate(cumulative_count = cumsum(count))

plot(tau$days,tau$cumulative_count,type = "l")
abline(a = 0, b = 10742/4108,lty = "dashed")

##################################

# int sum(lambda)_{ti :  } (Area)

for (t in seq_along(times)) {
  current_time <- times[t]
    bg <- grid_data |>
    mutate(
      bg = mu * trend.fun(current_time) *
        weekly.fun(current_time) *
        monthly.fun(current_time) * basevalue
    )
  history <- events[events$days < current_time & current_time - events$days <= 90, ]
  
  # Excitation components
  temp <- numeric(nrow(bg))
  if (nrow(history) > 0) {
    temporal_component <- excite.temporal.fun(current_time - history$days)
    
    for (i in seq_len(nrow(bg))) {
      spatial_component <- excite.spatial.fun(
        bg$x[i] - history$coorx,
        bg$y[i] - history$coory
      )
      temp[i] <- sum(spatial_component * temporal_component)
    }
  }
  
  bg$ex <- A * temp
  bg$lambda <- bg$bg + bg$ex
  bg$year <- 2013 + t  # Para luego filtrar
  
  results_list[[t]] <- bg
  max_lambda <- max(max_lambda, max(bg$lambda, na.rm = TRUE))
  min_lambda <- min(min_lambda, min(bg$lambda, na.rm = TRUE))
}



#########################################################################################################################

# Supuestos
ti <- unique(events$days)
ti <- sort(ti)
n <- length(events$days)
MU <- c()

for(i in seq_along(ti)){
  t <- ti[i]
  
  indices <- which(time.marks <= t)
  
  val <- mu * as.numeric(
    mean(trend.fun(time.marks[indices]) *
           weekly.fun(time.marks[indices]) *
           monthly.fun(time.marks[indices]))
  ) * t *
    mean(mub * background.marks) *
    (Xrange[2] - Xrange[1]) * (Yrange[2] - Yrange[1])
  
  MU <- c(MU, val)
}

alpha <- 0.05
lower <- upper <- numeric(length(MU))

for (i in seq_along(ti)) {
  if (i <= n) {
    lower[i] <- i*qbeta(alpha / 2, i + 1, n - i + 1)
    upper[i] <- i*qbeta(1 - alpha / 2, i + 1, n - i + 1)
  } else {
    lower[i] <- qgamma(alpha / 2, shape = i - n, scale = 1)
    upper[i] <- qgamma(1 - alpha / 2, shape = i - n, scale = 1)
  }
}

# Plot
plot(ti, MU, type = "l", lwd = 2, col = "blue",
     ylim = range(c(lower, upper, MU)),
     xlab = "Days since 2022-01-01", ylab = "Cumulative Intensity")

lines(ti, lower, col = "red", lty = 2)
lines(ti, upper, col = "red", lty = 2)
abline(a = 0, b = 7975 / 1186, lty = "dashed")  # Línea de referencia
legend("topleft", legend = c("Estimate", "95% CI", "Reference"),
       col = c("blue", "red", "black"), lty = c(1, 2, 2))



#################################################################################################
observed_cumsum <- cumsum(table(factor(events$days, levels = ti)))

plot(ti, MU, type = "l", col = "blue", lwd = 2,
     xlab = "Days since 2022-01-01", ylab = "Cumulative Intensity",
     ylim = range(c(MU, observed_cumsum)))
lines(ti, observed_cumsum, col = "black", lwd = 2)
legend("topleft", legend = c("Estimated (MU)", "Observed"),
       col = c("blue", "black"), lty = 1)
