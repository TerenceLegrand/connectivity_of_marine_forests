# -------------------------
# Load Connectivity results

#Connectivity.desc <- paste0("../Results [temporary right TL]/","/particlePairedConnectivityAveragedDistProba",n.season,"Table.desc")
Connectivity.desc <- paste0(connectivity.folder,"/particlePairedConnectivityAveragedDist",n.seasons,"Table.desc")
Connectivity <- attach.big.matrix(Connectivity.desc)
Connectivity <- data.table(Connectivity[,])

# Select only pair in distribution
if (useDistrib ==  TRUE) {
  Connectivity <- Connectivity[Connectivity$Pair.from %in% hex.distrib.range & Connectivity$Pair.to %in% hex.distrib.range,
                               c("Pair.from","Pair.to","Mean.Time","Mean.events","Distance")]
}

# -------------------------
# PLD distribution

if (pld.period == 0) {
  
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
    
  
  layout <- c(
    area(t = 0, l = 0, b = 4.5, r = 30),
    area(t = 5, l = 0, b = 20, r = 30)
  )
  # final plot arrangement
  plot_pd_weighting <- plot_range + plot_density + plot_layout(design = layout)
  
  
  ggsave(filename = "PD_weighting.pdf",
         plot = plot_pd_weighting,
         device = cairo_pdf,
         path = saving_dir.PD,
  )
  

  Connectivity$Mean.events.PD <- Connectivity$Mean.events
  
  for (pd.i in 1:max.PD) {
    Connectivity$Mean.events.PD[Connectivity$Mean.Time <= pd.i & Connectivity$Mean.Time > pd.i -1] <- Connectivity$Mean.events[Connectivity$Mean.Time <= pd.i & Connectivity$Mean.Time > pd.i -1]  * 
      mean(ds$y[ds$x <= pd.i & ds$x > pd.i -1])
  }
  
  
  Connectivity <- Connectivity[Connectivity$Mean.Time <= max.PD,]
  
  #plot(Connectivity$Mean.Time,Connectivity$Mean.events.PD)
  
  matrix.ev <- sparseMatrix(Connectivity$Pair.from,Connectivity$Pair.to, ,Connectivity$Mean.events.PD)
  
} else {
  
  # TL: filter on the pld.period given in the Pipeliner or at the beginning of this code
  Connectivity <- Connectivity[Connectivity$Mean.Time < pld.period + 1,]
  
  matrix.ev <- sparseMatrix(Connectivity$Pair.from,Connectivity$Pair.to, ,Connectivity$Mean.events)
  
}
  

# --------------------------
# For mean events 

sum.send <- rowSums(matrix.ev,na.rm = TRUE) # TL: sum horizontally, what cells send

sum.receive <- colSums(matrix.ev,na.rm = TRUE) # TL: sum vertically, what cells receive

# Backward
connectivity.matrix.back <- sweep(matrix.ev,2,sum.receive,FUN="/")

rownames(connectivity.matrix.back) <- colnames(connectivity.matrix.back) <- 1:length(connectivity.matrix.back[,1])

connectivity.matrix.back[is.nan(connectivity.matrix.back)] <- 0 # Remove NaN because 0/0 give NaN (many of colsum is 0)
connectivity.matrix.back[is.infinite(connectivity.matrix.back)] <- 0 # Remove Inf
connectivity.matrix.back <- as(connectivity.matrix.back, "CsparseMatrix") # recreate a sparse matrix

#Test to verify if normalisation went well
#if(sum(connectivity.matrix.back[,sample(source.sink.xy$Pair,1)]) != 1 ) { stop("Error :: 816")}

connectivity.matrix.back.inv <- t(connectivity.matrix.back)

#Test to verify if normalisation went well
#if(sum(connectivity.matrix.back.inv[sample(species_range,1),]) != 1 ) { stop("Error :: 816")}

# ---------------------------
# Compute centrality

connectivity.matrix.back.inv.list <- as.data.frame(summary(connectivity.matrix.back.inv)) # Transform matrix into list

colnames(connectivity.matrix.back.inv.list) <- c("from","to","weight")

# Create graph object from list
graph.back <- graph_from_data_frame(connectivity.matrix.back.inv.list,
                                    directed = TRUE,
                                    vertices = NULL)

graph.back.save <- graph.back

# Save proba weights 
graph.back.weight.proba <- E(graph.back)$weight

# Transform probabilities into distance using negative log transformation
E(graph.back)$weight <- -log(E(graph.back)$weight)

# Look for edges weight of 0 and remove them
graph.back <- delete_edges(graph.back, which(E(graph.back)$weight == 0))
if(sum(E(graph.back)$weight == 0) != 0 ) { stop("Error :: 361")}

# Run closeness centrality

cc =closeness(graph.back,
              vids = V(graph.back),
              mode = "out",
              weights = NULL,
              normalized = TRUE,
              cutoff = -1)

# Run harmonic centrality

hc = harmonic_centrality(graph.back,
                         vids = V(graph.back),
                         mode = "out",
                         weights = NULL,
                         normalized = TRUE,
                         cutoff = -1)

# Run betwennes centrality

bc = betweenness(graph.back,
                 v = V(graph.back),
                 directed = TRUE,
                 weights = NULL,
                 cutoff = -1)

if(sum(as.numeric(names(bc)) != as.numeric(names(hc))) != 0 ) { stop("Error :: 627")} # Test on centrality indices

centrality <- data.frame(cell =as.numeric(names(bc)),
                         cc = cc,
                         hc = hc,
                         bc = bc)

file_name <- paste0(saving_dir.PD,'/centrality_distrib','.csv')
write.table(centrality, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)

# Delete big matrices

rm(connectivity.matrix.back, connectivity.matrix.back.inv, matrix.ev,hc,bc)