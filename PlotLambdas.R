# Plot triggers and lambda for all the region.

# This has to be done today lol

# I'll do it for some times
times<- as.numeric(as.Date(paste0(2014:2024, "-12-31")) - as.Date("2014-01-01")) + 1
t = 1

bg <- grid_data |> mutate(bg = mu*trend.fun(times[t])*weekly.fun(times[t]) * monthly.fun(times[t])*basevalue) # mu(t)

history <- events[events$days < times[t] & times[t] - events$days <= 90,] # H_t

excite.temporal.fun(times[t] - history$days)

temp <- numeric(length(grid_data$basevalue))

for(i in 1:length(grid_data$basevalue)){
  temp[i] <- sum(excite.spatial.fun(bg$x[i] - history$coorx ,bg$y[i] - history$coory)*excite.temporal.fun(times[t] - history$days))
}
bg$ex <- A*temp

bg$lambda <- bg$bg + bg$ex


valid <-bg[bg$marks == 1, ]

ggplot() +
  geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +
  geom_tile(data = valid, aes(x = x, y = y, fill = lambda)) +
  scale_fill_gradientn(
    colours = c("blue", "lightblue", "yellow", "orange", "red"),
    name = "Intensity"
  ) +
  coord_sf() +
  labs(
    title = paste0("Intensity 2014" ),
    x = "X",
    y = "Y"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

######################################################################################
# LOOP

times <- as.numeric(as.Date(paste0(2014:2024, "-12-31")) - as.Date("2014-01-01")) + 1
results_list <- list()
max_lambda <- -Inf
min_lambda <- Inf

for (t in seq_along(times)) {
  current_time <- times[t]
  
  # Background component
  bg <- grid_data |>
    mutate(
      bg = mu * trend.fun(current_time) *
        weekly.fun(current_time) *
        monthly.fun(current_time) * basevalue
    )
  
  # Event history within 90 days
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

for (t in seq_along(results_list)) {
  bg <- results_list[[t]]
  valid <- bg[bg$marks == 1, ]
  
  p <- ggplot() +
    geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +
    geom_tile(data = valid, aes(x = x, y = y, fill = lambda)) +
    scale_fill_gradientn(
      colours = c("darkblue", "white", "yellow", "orange", "red"),
      limits = c((min_lambda - 0.1*min_lambda),(max_lambda)),  # Fijamos la escala
      name = expression("Intensity " * lambda)
    ) +
    coord_sf() +
    labs(
      title = paste0("Intensity ", bg$year[1]),
      x = "X",
      y = "Y"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.title = element_blank(),        # Quita los títulos "X" y "Y"
      axis.text = element_blank(),         # Quita los números del eje
      axis.ticks = element_blank(),        # Quita las líneas de los ticks
      panel.grid = element_blank()         # Quita la grilla
    )
  
  ggsave(
    filename = paste0("~/EcuadorOpenData/City/Model6-SA2/PlotsLambda/intensity_", 2013 + t, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 600  # más calidad
  )

  }

##############################################################################################################

subset_indices <- which((2013 + seq_along(results_list)) %in% 2014:2024)

all_lambda_subset <- unlist(
  lapply(results_list[subset_indices], function(bg) bg$lambda[bg$marks == 1])
)

min_lambda_sset <- min(all_lambda_subset, na.rm = TRUE)
max_lambda_sset <- max(all_lambda_subset, na.rm = TRUE)

plots <- list()

for (i in subset_indices) {
  bg <- results_list[[i]]
  valid <- bg[bg$marks == 1, ]
  
  p <- ggplot() +
    geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +
    geom_tile(data = valid, aes(x = x, y = y, fill = lambda)) +
    scale_fill_gradientn(
      colours = c("darkblue", "white", "yellow", "orange", "red"),
      limits = c((min_lambda_sset - 0.1 * min_lambda_sset), max_lambda_sset),
      name = expression("Intensity " * lambda)
    ) +
    coord_sf() +
    labs(title = paste0("Intensity ", bg$year[1])) +
    theme_minimal() +
    theme(
      legend.position = "none",  # Para evitar múltiples leyendas
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  plots[[length(plots) + 1]] <- p
}

library(patchwork)

combined_plot <- wrap_plots(plots, ncol = 4) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)

# Guarda el resultado
ggsave(
  filename = "~/EcuadorOpenData/City/Model6-SA2/intensity.png",
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 600
)

#############################################################################################
# PLot Prob being a triggered event

subset_indices <- which((2013 + seq_along(results_list)) %in% 2014:2024)

all_lambda_subset <- unlist(
  lapply(results_list[subset_indices], function(bg) bg$ex[bg$marks == 1])
)

min_lambda_sset <- min(all_lambda_subset, na.rm = TRUE)
max_lambda_sset <- max(all_lambda_subset, na.rm = TRUE)

plots <- list()

for (i in subset_indices) {
  bg <- results_list[[i]]
  valid <- bg[bg$marks == 1, ]
  
  p <- ggplot() +
    geom_sf(data = city, fill = NA, colour = "black", linewidth = 0.7) +
    geom_tile(data = valid, aes(x = x, y = y, fill = ex)) +
    scale_fill_gradientn(
      colours = c("darkblue", "white", "yellow", "orange", "red"),
      limits = c((min_lambda_sset - 0.1 * min_lambda_sset), max_lambda_sset),
      name = expression("Triggering effect ")
    ) +
    coord_sf() +
    labs(title = paste0("Triggering ", bg$year[1])) +
    theme_minimal() +
    theme(
      legend.position = "none",  # Para evitar múltiples leyendas
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  plots[[length(plots) + 1]] <- p
}

library(patchwork)

combined_plot <- wrap_plots(plots, ncol = 4) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(combined_plot)

# Guarda el resultado
ggsave(
  filename = "~/EcuadorOpenData/City/Model6-SA2/triggers.png",
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 600
)



































