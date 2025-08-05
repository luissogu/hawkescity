# Definir rango alrededor de (0,0)
xlim <- c(-5000, 5000)
ylim <- c(-5000, 5000)

# Filtrar índices para x e y dentro del rango
ix <- which(excite.spatial.base.x >= xlim[1] & excite.spatial.base.x <= xlim[2])
iy <- which(excite.spatial.base.y >= ylim[1] & excite.spatial.base.y <= ylim[2])

# Subconjunto de los valores para x, y y la matriz de valores (suponiendo que la matriz está en la forma [x,y])
x_sub <- excite.spatial.base.x[ix]
y_sub <- excite.spatial.base.y[iy]
z_sub <- excite.spatial.basevalue[ix, iy]

# Graficar sólo la zona cercana a (0,0)
filled.contour(x_sub/1000, y_sub/1000, z_sub,
               main = paste0("Spatial triggering Loop: ", kk))

# Crear data.frame para ggplot
df <- expand.grid(x = x_sub, y = y_sub)
df$z <- as.vector(z_sub)

# Plot
ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("darkblue", "yellow", "orange", "red"),
    na.value = "white",
    name = expression(h)  # 👈 mu con subíndice b
  )+
  labs(title = paste0("Spatial triggering Loop: ", kk), x = "X", y = "Y") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14,
      family = "sans"
    )
  )
