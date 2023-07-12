# Loop on all the datasets

correlation.final$Marker[correlation.final$Marker == "msat"] <- "microsatellites"

for (k in 1:length(correlation.final$Studies)) {
  
  study <- paste0(correlation.final$Studies[k],"_",correlation.final$Species[k],"_",correlation.final$Marker[k],correlation.final$Specificities[k],correlation.final$Diff.index[k])
  
  print(study)
  print(k)
  
  connectivity.study <- read.csv(paste0(saving_dir.studies.data,"/connectivity_",study,".csv"), check.names=FALSE)
  correlation.study <- read.csv(paste0(saving_dir.studies.data,"/correlation_",study,".csv"), check.names=FALSE)
  centrality <- read.csv(paste0(saving_dir.PD,"/centrality_distrib",".csv"), check.names=FALSE)
  
  # -------
  # Study info
  
  label_text <- paste0(correlation.final$species[study_k],'\n',correlation.final$Marker[study_k],", ",correlation.final$Diff.index[study_k],'\n',correlation.final$ref[study_k])
  
  info_study <- ggplot() +
    annotate("text", x = 1, y = 0,label = label_text, alpha = 1, hjust = 1, size = 3, angle = 90, colour = "#737476") +
    theme_void()
  
  # -------
  # Scatter predicted vs observed
  
  range.y <- range(connectivity.study$Differentiation, na.rm = TRUE)#range(connectivity.study$Differentiation)
  range.x <- range(connectivity.study$CCM, na.rm = TRUE) #range(connectivity.study$predict.mean.hc)
  
  
  scatter_plot_study <-  ggplot(connectivity.study, aes(x=CCM, y=Differentiation, group = 1)) +
    geom_point(size=0.5,color="#000000", alpha = 0.9) + xlim(range.x) + ylim(range.y) +
    geom_segment(aes(x = max(c(range.x[1],range.y[1])), xend = min(c(range.x[2],range.y[2])), y = max(c(range.x[1],range.y[1])), yend =min(c(range.x[2],range.y[2]))),color="grey",alpha = 0.5, linewidth=0.15) +
    invisible(geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=TRUE ,size=0.35, linetype = "longdash", formula = 'y ~ x')) +
    ylab(paste0("Observed genetic differentiation")) + xlab(paste0("Predicted genetic differentiation")) +
    annotate("text", alpha = 0.65, x = range.x[1], y = range.y[2], hjust=0,vjust=1 ,
             label = paste0("\nR²: ",format(round(correlation.final$r2.CCM[k], 3), nsmall = 3),
                            "\nDelta R²: ",format(round(correlation.final$delta_r2[k], 3), nsmall = 3),
                            "\np-value: ",format.pval(correlation.final$p.CCM[k], digits=3, eps = 0.001)), size = 2) +
    theme_plot +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          aspect.ratio=1)
  
  
  # -------
  # Map with populations and links
  
  source.sink.xy.i <- source.sink.xy[source.sink.xy$Pair %in% unique(c(connectivity.study$cell_from,connectivity.study$cell_to)),]
  
  x.min <- min(c(connectivity.study$Lon_from,connectivity.study$Lon_to)) 
  x.max <- max(c(connectivity.study$Lon_from,connectivity.study$Lon_to)) 
  y.min <- min(c(connectivity.study$Lat_from,connectivity.study$Lat_to)) 
  y.max <- max(c(connectivity.study$Lat_from,connectivity.study$Lat_to)) 
  
  x.min <- x.min - min(c(5, round( abs(diff(c(x.min , x.max))) + 1 )))
  x.max <- x.max + min(c(5, round( abs(diff(c(x.min , x.max))) + 1 )))
  y.min <- y.min - min(c(5, round( abs(diff(c(y.max , y.min))) + 1 )))
  y.max <- y.max + min(c(5, round( abs(diff(c(y.max , y.min))) + 1 )))
  
  if(abs(diff(c(x.max , x.min))) > abs(diff(c(y.max , y.min))) ) { 
    y.max <- ( y.min + (abs(diff(c(y.max , y.min))) / 2) ) + ( abs(diff(c(x.max , x.min))) / 2 )
    y.min <- ( y.min + (abs(diff(c(y.max , y.min))) / 2) ) - ( abs(diff(c(x.max , x.min))) / 2 )
  }
  
  if(abs(diff(c(x.max , x.min))) < abs(diff(c(y.max , y.min))) ) { 
    
    x.max <- ( x.min + (abs(diff(c(x.max , x.min))) / 2) ) + ( abs(diff(c(y.max , y.min))) / 2 )
    x.min <- ( x.min + (abs(diff(c(x.max , x.min))) / 2) ) - ( abs(diff(c(y.max , y.min))) / 2 )
    
  }
  
  worldMap <- ne_countries(scale = 10, returnclass = "sp")
  regionMap <- crop(worldMap,extent(x.min,x.max,y.min,y.max))
  
  lineConnections <- list()
  lineStrenght <- numeric(0)
  
  pathCoordinates_lon = numeric(0)
  pathCoordinates_lat = numeric(0)
  pathText <- data.frame()
  
  # Remove non-connected studies
  connectivity.study$Connectivity.distance[is.infinite(connectivity.study$nbr_step)] <- NA
  
  
  for( l.i in sort(connectivity.study$Connectivity.distance, index.return=T, decreasing = TRUE, na.last = FALSE)$ix ){
    lineStrenght <- c(lineStrenght,(connectivity.study[l.i,"Connectivity.distance"] ) )
    pointFrom <- c(as.numeric(as.character(connectivity.study[l.i,"Lon_from"])),as.numeric(as.character(connectivity.study[l.i,"Lat_from"]))) 
    pointTo <- c(as.numeric(as.character(connectivity.study[l.i,"Lon_to"])),as.numeric(as.character(connectivity.study[l.i,"Lat_to"]))) 
    routes_sl <- gcIntermediate(matrix(pointFrom,ncol=2),matrix(pointTo,ncol=2),n = 100, addStartEnd = TRUE, sp = TRUE, breakAtDateLine=TRUE)
    lineConnections = c(lineConnections,sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = l.i), match.ID = F))
    # Put the number of step on links
    pathCoordinates_lon  = rbind(pathCoordinates_lon,gcIntermediate(matrix(pointFrom,ncol=2),matrix(pointTo,ncol=2),n = 1, addStartEnd = FALSE, sp = FALSE, breakAtDateLine=FALSE)[1])
    pathCoordinates_lat  = rbind(pathCoordinates_lat,gcIntermediate(matrix(pointFrom,ncol=2),matrix(pointTo,ncol=2),n = 1, addStartEnd = FALSE, sp = FALSE, breakAtDateLine=FALSE)[2])
    pathText <- rbind(pathText,connectivity.study[l.i,"nbr_step"])
  }
  
  
  # remove no-connected pairs of site
  pathCoordinates_lon <- pathCoordinates_lon[!is.na(pathText)]
  pathCoordinates_lat <- pathCoordinates_lat[!is.na(pathText)]
  pathText <- pathText[!is.na(pathText),]
  
  lineConnectionsSp <- do.call(rbind, lineConnections)
  lineStrenght <- (lineStrenght - min(lineStrenght, na.rm = TRUE)) / max( lineStrenght - min(lineStrenght, na.rm = TRUE) , na.rm = TRUE)
  
  for( l.s in 1:length(lineStrenght) ) {
    if (! is.na(lineStrenght[l.s])) {
      if( lineStrenght[l.s] >= 0.8 ) { lineStrenght[l.s] <- 1 }
      if( lineStrenght[l.s] >= 0.6 & lineStrenght[l.s] < 0.8 ) { lineStrenght[l.s] <- 2 }
      if( lineStrenght[l.s] >= 0.4 & lineStrenght[l.s] < 0.6 ) { lineStrenght[l.s] <- 3 }
      if( lineStrenght[l.s] >= 0.2 & lineStrenght[l.s] < 0.4 ) { lineStrenght[l.s] <- 4 }
      if( lineStrenght[l.s] < 0.2 ) { lineStrenght[l.s] <- 5 }
    }
  }
  
  cent.i <- centrality$hc[centrality$cell %in% unique(c(connectivity.study$cell_from,connectivity.study$cell_to))]
  
  cent.i <- (max(cent.i)-cent.i)/(max(cent.i)-min(cent.i))
  
  myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414") # blue, green, yellow, orange, red
  
  # Connected
  map_link_study <- ggplot() + 
    geom_point(data = source.sink.xy.i, aes(x = Lon, y = Lat, color = 1),size=cent.i/2,alpha =1) +
    geom_polygon(data = regionMap, aes(x = long, y = lat, group = group), fill="#CDCDCD", colour = "#9E9E9E" , size=0.25 ) +
    geom_sf(data = st_as_sf(lineConnectionsSp[is.na(lineStrenght),]) , size= 1 , colour = "#E9E9E9",alpha =1) +
    geom_sf(data = st_as_sf(lineConnectionsSp) , size= 1 , colour = myColors[lineStrenght],alpha =1) +
    geom_point(data = source.sink.xy.i, aes(x = Lon, y = Lat,),size=cent.i/2,alpha =1) +
    scale_x_continuous("Longitude", limits = c(x.min , x.max), expand = c(0, 0)) +
    scale_y_continuous("Latitude", limits = c(y.min , y.max), expand = c(0, 0)) +
    binned_scale(aesthetics = "color", scale_name = "stepsn", name = "Oceanographic connectivity probability",
                 palette = function(x) myColors, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1), show.limits = TRUE, guide = "colorsteps") + 
    theme_map + theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      plot.background = element_blank(),
                      panel.background = element_blank())
  
  
  layout <- "
    AB
    "
  
  plot_annotation(title = paste(gsub("_"," ",correlation.study$Species),"--",correlation.study$Marker,"-",correlation.study$Diff.ind),
                  subtitle = correlation.final$ref[k])
  
  figure <- map_link_study + scatter_plot_study +
    plot_annotation(tag_levels = 'A', 
                    title = paste(gsub("_"," ",correlation.study$Species),"--",correlation.study$Marker,"-",correlation.study$Diff.ind),
                    subtitle = correlation.final$ref[k]) +
    plot_layout(design = layout, guides = "collect") & theme(legend.position = 'bottom')
  
  study_name <- paste0("ID_",correlation.final$data_id[k],"_",correlation.final$Species[k])
  
  
  ggsave(filename = paste0("/",study_name,".pdf"),
         plot = figure,
         device = cairo_pdf,
         path = saving_dir.studies.figure)
  
}

