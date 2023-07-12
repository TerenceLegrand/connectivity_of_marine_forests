# For LDD
pld.period = 0
saving_dir.PD.LDD <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/correlation_final",".csv"), check.names=FALSE)
centrality.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/centrality_distrib",".csv"), check.names=FALSE)

# For mean in-house dispersal
pld.period = 7
saving_dir.PD.mean <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.mean <- read.csv(paste0(saving_dir.PD.mean,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.mean <- read.csv(paste0(saving_dir.PD.mean,"/correlation_final",".csv"), check.names=FALSE)
centrality.mean <- read.csv(paste0(saving_dir.PD.mean,"/centrality_distrib",".csv"), check.names=FALSE)


# -------------------------------
# Data processed

connectivity.final.LDD$nbr_step_diff <- connectivity.final.LDD$nbr_step - connectivity.final.mean$nbr_step

connected.ldd <- is.infinite(connectivity.final.LDD$nbr_step) != is.infinite(connectivity.final.mean$nbr_step)

sum(connected.ldd)/length(connected.ldd)

connectivity.final.LDD$nbr_step_ratio <- connectivity.final.mean$nbr_step/connectivity.final.LDD$nbr_step

connectivity.final.LDD$nbr_step[connectivity.final.LDD$nbr_step_ratio<1 & !is.infinite(connectivity.final.LDD$nbr_step_ratio)]
connectivity.final.mean$nbr_step[connectivity.final.LDD$nbr_step_ratio<1 & !is.infinite(connectivity.final.LDD$nbr_step_ratio)]

connectivity.final.LDD$nbr_step_ratio[is.infinite(connectivity.final.LDD$nbr_step_ratio)] <- NA
connectivity.final.mean$nbr_step[is.infinite(connectivity.final.mean$nbr_step)] <- NA
connectivity.final.LDD$nbr_step[is.infinite(connectivity.final.LDD$nbr_step)] <- NA

to_test <- connectivity.final.LDD$nbr_step
mean(to_test, na.rm = TRUE)
median(to_test, na.rm = TRUE)
ic_n <- sum(!is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

to_test <- connectivity.final.mean$nbr_step
mean(to_test, na.rm = TRUE)
median(to_test, na.rm = TRUE)
ic_n <- sum(!is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

to_test <- connectivity.final.LDD$nbr_step_ratio
mean(to_test, na.rm = TRUE)
median(to_test, na.rm = TRUE)
ic_n <- sum(!is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

# Boxplot


ggplot(connectivity.final.LDD, aes(PLD,nbr_step_ratio)) + 
  geom_boxplot() +
  coord_flip(clip = "off") +
  labs(x="PD", y = "nbr of steps") +
  stat_summary(fun=mean, geom="point", shape=23, size=1) +
  theme(axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6))

# Confidence interval
ic_n <- length(! is.na(connectivity.final.LDD$nbr_step_ratio))
ic_mean = mean(connectivity.final.LDD$nbr_step_ratio, na.rm = TRUE)
ic_sd = sd(connectivity.final.LDD$nbr_step_ratio, na.rm = TRUE)

margin <- qt(0.975,df=ic_n-1)*ic_sd/sqrt(ic_n)

# Get outliers

#outliers.mean <- connectivity.final.LDD$nbr_step_ratio > quantile(connectivity.final.LDD$nbr_step_ratio, 0.975, na.rm = TRUE)

lineConnections <- list()
lineStrenght <- numeric(0)

for( l.i in sort(connectivity.final.LDD$nbr_step_ratio, index.return=T, decreasing = TRUE, na.last = TRUE)$ix ){
  lineStrenght <- c(lineStrenght,(connectivity.final.LDD[l.i,"nbr_step_ratio"] ) )
  pointFrom <- c(as.numeric(as.character(connectivity.final.LDD[l.i,"Lon_from"])),as.numeric(as.character(connectivity.final.LDD[l.i,"Lat_from"])))
  pointTo <- c(as.numeric(as.character(connectivity.final.LDD[l.i,"Lon_to"])),as.numeric(as.character(connectivity.final.LDD[l.i,"Lat_to"])))
  routes_sl <- gcIntermediate(matrix(pointFrom,ncol=2),matrix(pointTo,ncol=2),n = 100, addStartEnd = TRUE, sp = TRUE, breakAtDateLine=TRUE)
  lineConnections = c(lineConnections,sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = l.i), match.ID = F))
}

lineConnectionsSp <- do.call(rbind, lineConnections)

# Extract important info

outliers.mean <- lineStrenght > quantile(lineStrenght, 0.975, na.rm = TRUE)
disconnected <- is.na(lineStrenght)
fastest.mean <- lineStrenght < 1

connected.ldd <- connected.ldd[sort(connectivity.final.LDD$nbr_step_ratio, index.return=T, decreasing = TRUE, na.last = TRUE)$ix]

# -------------------------------

# Panel A

MFC.PD <- read.xlsx(paste0(data.folder,"PD_litterature_review.xlsx"), check.names=FALSE, startRow = 1)

pd.in.house <- numeric(0)
pd.third.party <- numeric(0)
for (i in 1:length(MFC.PD$Species)) {
  if (MFC.PD$Dispersal.category[i] ==  "in-house") {
    pd.in.house <- c(pd.in.house,seq(as.numeric(MFC.PD$`Minimum.(d)`[i]),as.numeric(MFC.PD$`Maximum.(d)`[i]),1))
  } else {
    pd.third.party <- c(pd.third.party,seq(as.numeric(MFC.PD$`Minimum.(d)`[i]),as.numeric(MFC.PD$`Maximum.(d)`[i]),1))
  }
}

ldd.coeff <- 1/1000

x.pd <- seq(1,max.PD,0.05)
y.pd <- dnorm(x.pd, mean = mean(pd.in.house), sd = sd(pd.in.house))
y.pd[y.pd < (ldd.coeff/(sqrt(2*pi)*sd(pd.in.house)))] <- (ldd.coeff/(sqrt(2*pi)*sd(pd.in.house))) # min value is 1/1000 max value

ds <- data.frame(x = x.pd, y = y.pd)

pd.range <- data.frame(x.min = c(max(pd.in.house),min(pd.third.party),min(pd.in.house)),
                       x.max = c(min(c(max.PD,max(pd.third.party))),min(c(max.PD,max(pd.third.party))),max(pd.in.house)),
                       x.mean = c(mean(c(max(pd.in.house),(min(c(max.PD,max(pd.third.party)))))),mean(c(min(pd.third.party),(min(c(max.PD,max(pd.third.party)))))),mean(c(min(pd.in.house),max(pd.in.house)))),
                       y = c(0,0.5,1),
                       label = c("long distance dispersal","third-party dispersal","in-house dispersal"))



plot_density <- ggplot(ds, aes(x,y)) + 
  geom_line(colour = "black") +
  labs(x="Mean trajectory time (d)", y = "Normalised weighting factor") +
  scale_x_continuous(breaks = c(0,mean(pd.in.house),max(pd.in.house),max.PD), labels = as.character(c(0,round(mean(pd.in.house),2),max(pd.in.house),max.PD))) +
  scale_y_continuous(breaks = c(ldd.coeff/(sqrt(2*pi)*sd(pd.in.house)),1/(sqrt(2*pi)*sd(pd.in.house))), labels = c(0.001,1),expand = c(0.1,0)) +
  theme_classic() 

plot_range <- ggplot() +
  geom_linerange(data = pd.range, aes(xmin=x.min, xmax=x.max, y=y), colour = "#999999", linewidth = 0.5,inherit.aes = FALSE) +
  geom_point(data = pd.range, aes(x.min,y), colour = "#999999", size = 0.5) + 
  geom_point(data = pd.range, aes(x.max,y), colour = "#999999", size = 0.5) +
  geom_text(data = pd.range, aes(x.mean, y, label = label),  colour = "#999999",  hjust = 0.5, vjust = -0.5, size = 4) +
  theme_void() +
  coord_cartesian(ylim = c(0, 2))

# -------------------------------
# Panel B

step_nbr <- rbind(data.frame(step = connectivity.final.mean$nbr_step,
                               PD = "mean"),
                  data.frame(step = connectivity.final.LDD$nbr_step,
                             PD = "LDD"))

plot_boxplot <- ggplot(step_nbr, aes(PD, step)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 0.1, outlier.alpha = 0.1) +
  geom_boxplot(fill = "#CDCDCD") +
  coord_flip(clip = "off") +
  scale_y_log10(breaks = c(1,5,10,50,100,150)) +
  scale_x_discrete(labels=c('LDD','Fixed PD')) +
  labs(y = "Number of stepping stones") +
  theme_classic() + theme(axis.line.y = element_blank(), axis.title.y = element_blank()) 

# -------------------------------
# Panel C

crs.proj <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"

# Map
worldMap.proj <- spTransform(worldMap, CRS(crs.proj))

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


lineConnectionsSp.proj <- spTransform(lineConnectionsSp, CRS(crs.proj))

plot_connexions_LDD <- mapRegion + # mapRegionNet
  geom_sf(data = st_as_sf(lineConnectionsSp.proj [connected.ldd,]) , size= 1 , colour = color[3] ,alpha =0.1) +
  theme(legend.position = "none",plot.margin = margin(0,0,0,0,"cm")) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# -------------------------------
# Final figure

layout <- c(
  area(t = 0, l = 0, b = 3, r = 9),
  area(t = 4, l = 0, b = 10, r = 9),
  area(t = 4, l = 11, b = 10, r = 20),
  area(t = 11, l = 0, b = 30, r = 20)
)

figure_4 <- plot_range + plot_density + plot_boxplot + plot_connexions_LDD +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')

ggsave(filename = paste0("/Figure_4",".pdf"),
       plot = figure_4,
       device = cairo_pdf,
       path = saving_dir.paper.main,
       width = 21,
       height = 29.7,
       units = "cm"
)
