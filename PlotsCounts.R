library(ggplot2)
library(sf)
library(dplyr)
library(patchwork) # Para combinar
library(lubridate)

# Convertir días en fechas
events <- events %>%
  mutate(date = as.Date("2014-01-01") + days,
         year = year(date))

# Crear una paleta de colores por año
years <- sort(unique(events$year))
palette <- scales::hue_pal()(length(years))
names(palette) <- years

# Lista para guardar los plots
plot_list <- list()

# Crear los plots por año
for (y in years) {
  year_data <- events %>% filter(year == y)
  count <- nrow(year_data)
  plot_list[[as.character(y)]] <- ggplot() +
    geom_sf(data = city, fill = NA, colour = "grey70") +
    geom_point(data = year_data, aes(x = coorx, y = coory),
               colour = palette[as.character(y)], size = 0.6) +
    annotate("text", x = Inf, y = Inf, label = paste0("N = ", count),
             hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "bold") +
    labs(title = paste("Year:", y), x = "X", y = "Y") +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 🖼️ Seleccionar 4 años para mostrar juntos
years_to_plot <- c(2014, 2015,2017,2018,2020,2021, 2023,2024)  # Puedes cambiar estos años
selected_plots <- plot_list[as.character(years_to_plot)]

# 📦 Combinar usando patchwork
combined_plot <- wrap_plots(selected_plots, ncol = 4)

# 💾 Guardar en alta calidad
ggsave("Model6/homicides_by_year.png", combined_plot, width = 12, height = 8, dpi = 600)
