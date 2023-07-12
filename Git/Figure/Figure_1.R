# #--------------------
# # DATA


# Make a map with mollweide projection
crs.proj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"

# Map
worldMap.proj <- spTransform(worldMap, CRS(crs.proj))

# Source/sink
source.sink.proj <- source.sink.xy.distrib[,c("Lon","Lat")]
coordinates(source.sink.proj) <- c("Lon", "Lat")
proj4string(source.sink.proj) <- crs(worldMap)
source.sink.proj <- spTransform(source.sink.proj, CRS(crs.proj))

# Sampling
sampling_coord <- data.frame(Lon = c(connectivity.final$Lon_from,connectivity.final$Lon_to),
                             Lat = c(connectivity.final$Lat_from,connectivity.final$Lat_to))
sampling_coord <- unique(sampling_coord)

coordinates(sampling_coord) <- c("Lon", "Lat")
proj4string(sampling_coord) <- crs(worldMap)
sampling_coord.proj <- spTransform(sampling_coord, CRS(crs.proj))

# Make background map
bb.west <- data.frame(Lon = -180,
                       Lat = seq(90,-90,-1))
bb.east <- data.frame(Lon = 180,
                      Lat = seq(90,-90,-1))
equator <- data.frame(Lon = seq(-180,180,1),
                      Lat = 0)
trop.capricorne <- data.frame(Lon = seq(-180,180,1),
                      Lat = -23.43624)
trop.cancer <- data.frame(Lon = seq(-180,180,1),
                          Lat = 23.43624)
trop.capricorne <- data.frame(Lon = seq(-180,180,1),
                              Lat = -23.43624)
polar.artic <- data.frame(Lon = seq(-180,180,1),
                          Lat = 66)
polar.antartic <- data.frame(Lon = seq(-180,180,1),
                          Lat = -66)
meridien.0 <- data.frame(Lon = 0,
                         Lat = seq(90,-90,-1))
meridien.west <- data.frame(Lon = -90,
                         Lat = seq(90,-90,-1))
meridien.east <- data.frame(Lon = 90,
                         Lat = seq(90,-90,-1))


coordinates(bb.west) <- c("Lon", "Lat")
proj4string(bb.west) <- crs(worldMap)
bb.west <- spTransform(bb.west, CRS(crs.proj))

coordinates(bb.east) <- c("Lon", "Lat")
proj4string(bb.east) <- crs(worldMap)
bb.east <- spTransform(bb.east, CRS(crs.proj))

coordinates(equator) <- c("Lon", "Lat")
proj4string(equator) <- crs(worldMap)
equator <- spTransform(equator, CRS(crs.proj))

coordinates(trop.capricorne) <- c("Lon", "Lat")
proj4string(trop.capricorne) <- crs(worldMap)
trop.capricorne <- spTransform(trop.capricorne, CRS(crs.proj))

coordinates(trop.cancer) <- c("Lon", "Lat")
proj4string(trop.cancer) <- crs(worldMap)
trop.cancer <- spTransform(trop.cancer, CRS(crs.proj))

coordinates(polar.artic) <- c("Lon", "Lat")
proj4string(polar.artic) <- crs(worldMap)
polar.artic <- spTransform(polar.artic, CRS(crs.proj))

coordinates(polar.antartic) <- c("Lon", "Lat")
proj4string(polar.antartic) <- crs(worldMap)
polar.antartic <- spTransform(polar.antartic, CRS(crs.proj))

coordinates(meridien.0) <- c("Lon", "Lat")
proj4string(meridien.0) <- crs(worldMap)
meridien.0 <- spTransform(meridien.0, CRS(crs.proj))

coordinates(meridien.west) <- c("Lon", "Lat")
proj4string(meridien.west) <- crs(worldMap)
meridien.west <- spTransform(meridien.west, CRS(crs.proj))

coordinates(meridien.east) <- c("Lon", "Lat")
proj4string(meridien.east) <- crs(worldMap)
meridien.east <- spTransform(meridien.east, CRS(crs.proj))


mapRegion <- ggplot() +
  geom_path(data = as.data.frame(bb.west) ,  aes(x = Lon, y = Lat)) +
  geom_path(data = as.data.frame(bb.east) ,  aes(x = Lon, y = Lat)) +
  geom_path(data = as.data.frame(equator) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(trop.capricorne) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(trop.cancer) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(polar.artic) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(polar.antartic) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(meridien.0) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(meridien.west) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_path(data = as.data.frame(meridien.east) ,  aes(x = Lon, y = Lat), linewidth = 0.1, alpha = 0.2) +
  geom_polygon(data = worldMap.proj , fill = "#CDCDCD", colour = "#CDCDCD" , linewidth=0.1 ,  aes(long, lat, group = group))
  
# Add sampled population and marine forest distribution

figure_1 <- mapRegion + # mapRegionNet
  geom_point(data = as.data.frame(source.sink.proj),  aes(x = Lon, y = Lat), size = 0.1, alpha = 0.6, colour = "#333333") +
  geom_point(data = as.data.frame(sampling_coord.proj) ,  aes(x = Lon, y = Lat), size = 0.1, alpha = 1, colour = "#ffd700") +
  annotate("text", alpha = 0.65, x = as.data.frame(equator)[1,1] , y = -as.data.frame(meridien.0)[1,2], hjust=0,vjust=0.5 , size = 2,
           label = paste0("Sampled populations\n",
                          "Marine forest distribution")) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Save
ggsave(filename = paste0("/Figure_1",".pdf"),
       plot = figure_1,
       device = cairo_pdf,
       path = saving_dir.paper.main,
       width = 16.7,
       height = 8.4,
       units = "cm"
)
