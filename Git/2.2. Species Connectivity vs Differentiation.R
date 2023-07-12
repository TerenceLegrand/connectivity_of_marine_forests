# -----------------
# Compute shortest distance

dist.mean <- matrix(NA, nrow = length(pop.cells.dataset$Pair), ncol = length(pop.cells.dataset$Pair))
path_nbr.mean <- matrix(NA, nrow = length(pop.cells.dataset$Pair), ncol = length(pop.cells.dataset$Pair))

for (a in 1:length(pop.cells.dataset$Pair)) {
  for (b in 1:length(pop.cells.dataset$Pair)) {
    
    if (a != b & (sum(connectivity.matrix.back.inv.list$from == pop.cells.dataset[a,"Pair"]) * sum(connectivity.matrix.back.inv.list$to == pop.cells.dataset[b,"Pair"])) > 0 ) { # add an if condition on presence of the vertice in the graph, if not, prob = 0 and dist = inf

      path <- shortest_paths(
        graph.back,
        from = as.character(pop.cells.dataset[a,"Pair"]),
        to = as.character(pop.cells.dataset[b,"Pair"]),
        mode = "out",
        weights = NULL,
        algorithm = c("dijkstra"))
      
      path <- as.numeric(path$vpath[[1]])
      
      if (length(path) > 0) { # it means there is a path
        
        dist.mean[a,b] <- exp(-distances(
          graph.back,
          v = as.character(pop.cells.dataset[a,"Pair"]),
          to = as.character(pop.cells.dataset[b,"Pair"]),
          mode = "out",
          weights = NULL,
          algorithm = "dijkstra"
        ))
        
        path_nbr.mean[a,b] <- length(path)-1
        
        # EP = rep(path, each=2)[-1]
        # EP = EP[-length(EP)]
        # 
        # E(graph.back.save)[get.edge.ids(graph.back.save,(EP))] # In grpah.back.save is the probabilities
        # dist.mean[a,b] <- prod(E(graph.back.save)$weight[get.edge.ids(graph.back.save,(EP))],na.rm = TRUE)
      } else {
        dist.mean[a,b] <- 0
        path_nbr.mean[a,b] <- Inf
      } 
    } else {         
      dist.mean[a,b] <- 0
      path_nbr.mean[a,b] <- Inf
      }
  } 
}


# ----

# Direct probabilities to indirect probabilities (Pab == Pba)
# Formula from Legrand et al., 2022

# Pab = 1
# Pba = 1
# 1 - (1-Pab)*(1-Pba)
# Pab + Pba - Pab*Pba

dist.mean.ab <- dist.mean.ba <- dist.mean
dist.mean.ab[upper.tri(dist.mean, diag = TRUE)] <- 0
dist.mean.ba[lower.tri(dist.mean, diag = FALSE)] <- 0
dist.mean <- dist.mean.ab + t(dist.mean.ba) - dist.mean.ab*t(dist.mean.ba) # to deal with 1-low proba issues
dist.mean[upper.tri(dist.mean, diag = FALSE)] <- NA

# Find the shortest number of path

path_nbr.mean.min <- path_nbr.mean
for (a in 1:length(pop.cells.dataset$Pair)) {
  for (b in 1:length(pop.cells.dataset$Pair)) {
    if (a !=b) {
      path_nbr.mean.min[a,b] <- min(c(path_nbr.mean[a,b],path_nbr.mean[b,a]),na.rm = TRUE)
    }
  }
}
path_nbr.mean.min[upper.tri(path_nbr.mean.min, diag = TRUE)] <- NA
path_nbr.mean <- path_nbr.mean.min

# ----------------------------------
# Compute centrality differences

# 1. Using mean back probabilities

for (a in 1:length(pop.cells.dataset$Pair)) {
  if (sum(connectivity.matrix.back.inv.list$from == pop.cells.dataset[a,"Pair"]) * sum(connectivity.matrix.back.inv.list$to == pop.cells.dataset[a,"Pair"]) > 0) { # add an if condition on presence of the vertice in the graph, if not, centrality = 0
    pop.cells.dataset$hc[a] <- centrality$hc[centrality$cell == pop.cells.dataset$Pair[a]]
  } else  {
    pop.cells.dataset$hc[a] <- 0
  }
}

hc.mean <- matrix(NA, nrow = length(pop.cells.dataset$Pair), ncol = length(pop.cells.dataset$Pair))

for (a in 1:length(pop.cells.dataset$Pair)) {
  for (b in 1:a) {
    hc.mean[a,b] <- abs(diff(c(pop.cells.dataset$hc[a],pop.cells.dataset$hc[b])))
  }
}

diag(hc.mean) <- NA
hc.mean <- (max(hc.mean,na.rm = TRUE)-hc.mean)/abs(diff(range(hc.mean,na.rm = TRUE)))

# -----------------

# Find duplicated 

dup <- which(! duplicated(pop.cells.dataset$Pair))
pop.cells.dataset.dup <- pop.cells.dataset[dup,]

hc.mean.dup <- dist.max.dup <- dist.mean.dup <- popDifferentiation.dup <- path_nbr.mean.dup <- matrix(NA, nrow = length(dup), ncol = length(dup))
rownames(popDifferentiation.dup) <- colnames(popDifferentiation.dup) <- pop.cells.dataset.dup$Pair
rownames(dist.mean.dup) <- colnames(dist.mean.dup) <- pop.cells.dataset.dup$Pair
rownames(hc.mean.dup) <- colnames(hc.mean.dup) <- pop.cells.dataset.dup$Pair
rownames(path_nbr.mean.dup) <- colnames(path_nbr.mean.dup) <- pop.cells.dataset.dup$Pair

for (da in 1:length(dup)) {
  for (db in 1:length(dup)) {
    if (da > db ) {
      
      popDifferentiation.dup[da,db] <- mean(c(as.numeric(as.matrix(popDifferentiation[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]]])), # ab
                                              as.numeric(as.matrix(popDifferentiation[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]]]))), na.rm = TRUE) # ba

      dist.mean.dup[da,db] <- mean(c(as.numeric(as.matrix(dist.mean[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]]])), # ab
                                     as.numeric(as.matrix(dist.mean[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]]]))), na.rm = TRUE) # ba

      path_nbr.mean.dup[da,db] <- mean(c(as.numeric(as.matrix(path_nbr.mean[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]]])), # ab
                                     as.numeric(as.matrix(path_nbr.mean[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]]]))), na.rm = TRUE) # ba
      
      hc.mean.dup[da,db] <- mean(c(as.numeric(as.matrix(hc.mean[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]]])), # ab
                                   as.numeric(as.matrix(hc.mean[pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[db]],pop.cells.dataset$Pair %in% pop.cells.dataset$Pair[dup[da]]]))), na.rm = TRUE) # ba
      
    }
  }
}

# -----------------
# Computing Connectivity vs Differentiation

# Constructing dataframe with cell from, cell to, differentiation, connectivity 

size_matrix <- (ncol(popDifferentiation.dup)^2-ncol(popDifferentiation.dup))/2
mat_glmm <- matrix(nrow=size_matrix,ncol=17)
loop <- 0

for (a in 2:ncol(popDifferentiation.dup)) { 
  for (b in 1:(a-1)) {
    
    loop <- loop +1
    
    # Pair from cells
    mat_glmm[loop,1] <- pop.cells.dataset.dup$Pair[a]
    mat_glmm[loop,2] <- pop.cells.dataset.dup$Lon[a]
    mat_glmm[loop,3] <- pop.cells.dataset.dup$Lat[a]
    # Pair to cells
    mat_glmm[loop,4] <- pop.cells.dataset.dup$Pair[b]
    mat_glmm[loop,5] <- pop.cells.dataset.dup$Lon[b]
    mat_glmm[loop,6] <- pop.cells.dataset.dup$Lat[b]
    # Pair from sample
    mat_glmm[loop,7] <- pop.cells.dataset.dup$Sample.Lon[a]
    mat_glmm[loop,8] <- pop.cells.dataset.dup$Sample.Lat[a]
    # Pair to sample
    mat_glmm[loop,9] <- pop.cells.dataset.dup$Sample.Lon[b]
    mat_glmm[loop,10] <- pop.cells.dataset.dup$Sample.Lat[b]
    # Differentiation
    mat_glmm[loop,11] <- popDifferentiation.dup[a,b]
    # Connectivity mean proba
    mat_glmm[loop,12] <- dist.mean.dup[a,b]
    # Connectivity mean distance
    mat_glmm[loop,13] <- -log(dist.mean.dup[a,b] + 9e-324) #+ 9e-324 to deal with zero
    # HC diff
    mat_glmm[loop,14] <- hc.mean.dup[a,b]
    # HC from
    mat_glmm[loop,15] <- pop.cells.dataset.dup$hc[a]
    # HC to
    mat_glmm[loop,16] <- pop.cells.dataset.dup$hc[b]
    # Number of path
    mat_glmm[loop,17] <- path_nbr.mean.dup[a,b]
  }
}

colnames(mat_glmm) <- c('pop_A','pop_A_Lon','pop_A_Lat','pop_B','pop_B_Lon','pop_B_Lat','sample_A_Lon','sample_A_Lat','sample_B_Lon','sample_B_Lat','differentiation','connectivity.mean.proba','connectivity.mean.distance','HC_diff','HC_pop_A','HC_pop_B','nbr_step')
mat_glmm <- as.data.frame(mat_glmm)

connectivity.study <- data.frame(cell_from = mat_glmm$pop_A,
                                 Lon_from = mat_glmm$pop_A_Lon,
                                 Lat_from = mat_glmm$pop_A_Lat,
                                 cell_to = mat_glmm$pop_B,
                                 Lon_to = mat_glmm$pop_B_Lon,
                                 Lat_to = mat_glmm$pop_B_Lat,
                                 sample_Lon_from = mat_glmm$sample_A_Lon,
                                 sample_Lat_from = mat_glmm$sample_A_Lat,
                                 sample_Lon_to = mat_glmm$sample_B_Lon,
                                 sample_Lat_to = mat_glmm$sample_B_Lat,
                                 Differentiation = mat_glmm$differentiation,
                                 Connectivity.mean = mat_glmm$connectivity.mean.proba,
                                 Connectivity.distance = mat_glmm$connectivity.mean.distance,
                                 Connectivity.distance.inverse = 1/mat_glmm$connectivity.mean.distance,
                                 HC_diff = mat_glmm$HC_diff,
                                 HC_pop_A = mat_glmm$HC_pop_A,
                                 HC_pop_B = mat_glmm$HC_pop_B,
                                 nbr_step = mat_glmm$nbr_step,
                                 population = as.factor(paste0(species.studies[s],"_",mat_glmm$pop_A)),
                                 PLD = pld.period,
                                 Species = species,
                                 Study = species.studies[s],
                                 Marker = studies.marker[i],
                                 Specificity = studies.spec[e],
                                 Diff.index = studies.diffi[j],
                                 stringsAsFactors = FALSE )

# Prune the dataset
if(pruneData) { 
  connectivity.study <- connectivity.study[connectivity.study$Differentiation > 0 & connectivity.study$Differentiation < 1,] # removing Diff < 0 and Diff = 1 
}

# Transform Fst into Fst/(1-Fst)
if( transformFST & studies.diffi[j] == "Fst" ) { 
  connectivity.study$Differentiation <- connectivity.study$Differentiation/(1-connectivity.study$Differentiation)
}

# Filter uncanny data
connectivity.study$Connectivity.distance[is.infinite(connectivity.study$Connectivity.distance)] <- NA
connectivity.study <- connectivity.study[!is.na(connectivity.study$Differentiation),]

filter <- length(connectivity.study$Differentiation)

if (filter > 8) { # To have enough point to process the correlation
  
  # Processing simple pearson correlation and linear model (should go for linear mixed model)
  
  # GLMM
  
  # null model
  model.null <- lmer(Differentiation ~ 1 + (1|population), connectivity.study, REML=F)
  obs.pred.null <- data.frame(observed=(connectivity.study$Differentiation),predicted=predict(model.null),population=connectivity.study$population)
  model <- lm(observed ~ predicted, data = obs.pred.null)
  aic.SM <- AIC(model)
  r2.SM = summary(model)$adj.r.squared
  p.SM = anova(model)$`Pr(>F)`[1]
  
  # Mean connectivity
  model.mean <- lmer(Differentiation ~ Connectivity.distance + (1|population), connectivity.study, REML=F,na.action = na.omit)
  obs.pred <- data.frame(observed=(connectivity.study$Differentiation),predicted=predict(model.mean),population=connectivity.study$population)
  model <- lm(observed ~ predicted, data = obs.pred)
  aic.CM <- AIC(model)
  r2.CM = summary(model)$adj.r.squared
  p.CM = anova(model)$`Pr(>F)`[1]

  # Distance + Centrality
  model.mean.cent <- lmer(Differentiation ~ Connectivity.distance + HC_diff + (1|population), connectivity.study, REML=F)
  obs.pred <- data.frame(observed=(connectivity.study$Differentiation),predicted=predict(model.mean.cent),population=connectivity.study$population)
  model <- lm(observed ~ predicted, data = obs.pred)
  aic.CCM <- AIC(model)
  r2.CCM <- summary(model)$adj.r.squared
  p.CCM <- anova(model)$`Pr(>F)`[1]
  
  correlation.study <- data.frame(Species = species,
                                   Studies = species.studies[s],
                                   Marker = studies.marker[i],
                                   Specificities = studies.spec[e],
                                   Diff.index = studies.diffi[j],
                                   PLD = pld.period,
                                   nbr.pop = length(pop.cells.dataset$Pair),
                                   nbr.Pairs = nrow(connectivity.study),
                                   nbr.Pairs.Connected = sum(!is.infinite(connectivity.study$nbr_step)),
                                   aic.SM = aic.SM, # null model
                                   r2.SM = r2.SM,
                                   p.SM = p.SM,
                                   aic.CM = aic.CM, # mean connectivity
                                   r2.CM = r2.CM,
                                   p.CM= p.CM,
                                   aic.CCM = aic.CCM, # mean connectivity + hc 
                                   r2.CCM = r2.CCM,
                                   p.CCM= p.CCM,
                                   stringsAsFactors = FALSE)
  
  connectivity.study <- cbind(connectivity.study,data.frame(null.obs.pred = predict(model.null),
                                                            CM = predict(model.mean),
                                                            CCM = predict(model.mean.cent)))
  
  
} else { # If not, is NA value
  
  obs.pred <- data.frame(observed=NA,predicted=NA,population=NA)
  
  correlation.study <- data.frame(Species = species,
                                   Studies = species.studies[s],
                                   Marker = studies.marker[i],
                                   Specificities = studies.spec[e],
                                   Diff.index = studies.diffi[j],
                                   PLD = pld.period,
                                   nbr.pop = length(pop.cells.dataset$Pair),
                                   nbr.Pairs = nrow(connectivity.study),
                                   nbr.Pairs.Connected = sum(!is.infinite(connectivity.study$nbr_step)),
                                   aic.SM = NA, # null model
                                   r2.SM = NA,
                                   p.SM = NA,
                                   aic.CM = NA, # mean connectivity
                                   r2.CM = NA,
                                   p.CM= NA,
                                   aic.CCM = NA, # mean connectivity + hc 
                                   r2.CCM = NA,
                                   p.CCM= NA,
                                   stringsAsFactors = FALSE)
  
  
  if (length(connectivity.study$Differentiation) == 0) {connectivity.study[1,] <- NA}
  
  connectivity.study <- cbind(connectivity.study,data.frame(null.obs.pred = NA,
                                                            CM = NA,
                                                            CCM = NA))
  
  
  
}


# -----------------
# Global information 


# Ref
correlation.study$ref <- paste0(MFC.studies$First_author[MFC.studies$studies %in% correlation.study$Studies],
                                   " et al., ",as.character(MFC.studies$Year[MFC.studies$studies %in% correlation.study$Studies]))
# Marine realm

popCoordinates <- unique(cbind(c(connectivity.study$sample_Lon_from,connectivity.study$sample_Lon_to),
                               c(connectivity.study$sample_Lat_from,connectivity.study$sample_Lat_to)))
popCoordinates <- as.data.frame(popCoordinates)
colnames(popCoordinates) <- c("Lon","Lat")

# Latitudinal extent
correlation.study$lat_extent <- round(abs(diff(range(popCoordinates$Lat))),1)

# Sampling extent
popCoordinates.sf <- st_as_sf(popCoordinates, coords = c("Lon", "Lat"), crs = dt.projection)
popDistance <- st_distance(popCoordinates.sf, which = "Great Circle")

correlation.study$sampling_max_extent <- max(popDistance)
correlation.study$sampling_mean_distance <- as.numeric(mean(popDistance[lower.tri(popDistance, diag = FALSE)]))
correlation.study$sampling_sd_distance <- sd(popDistance[lower.tri(popDistance, diag = FALSE)])
correlation.study$sampling_cv_distance <- sd(popDistance[lower.tri(popDistance, diag = FALSE)])/as.numeric(mean(popDistance[lower.tri(popDistance, diag = FALSE)]))

# Marine biogeographical region 
dat <- data.frame(Longitude = popCoordinates$Lon,
                  Latitude = popCoordinates$Lat)
# Assignment modified according
coordinates(dat) <- ~ Longitude + Latitude
# Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
proj4string(dat) <- proj4string(realm.shp)

diff_biogeo <- over(dat, realm.shp)

correlation.study$ecoregions <- paste(unique(diff_biogeo$ECOREGION), collapse = ',')
correlation.study$ecoregions_nbr <- length(unique(diff_biogeo$ECOREGION))
correlation.study$province <- paste(unique(diff_biogeo$PROVINCE), collapse = ',')
correlation.study$province_nbr <- length(unique(diff_biogeo$PROVINCE))
correlation.study$REALM <- paste(unique(diff_biogeo$REALM), collapse = ',')
correlation.study$REALM_nbr <- length(unique(diff_biogeo$REALM))
correlation.study$ecoregions.cl <- paste(gsub("[^A-Z]","", unique(diff_biogeo$ECOREGION)), collapse = ',')
correlation.study$province.cl <- paste(gsub("[^A-Z]","", unique(diff_biogeo$PROVINCE)), collapse = ',')
correlation.study$REALM.cl <- paste(gsub("[^A-Z]","", unique(diff_biogeo$REALM)), collapse = ',')

# Phylogeny
correlation.study$genus <- MFC.taxo$taxGenus[MFC.taxo$name %in% gsub("_"," ",correlation.study$Species)]
correlation.study$family <- MFC.taxo$taxFamily[MFC.taxo$name %in% gsub("_"," ",correlation.study$Species)]
correlation.study$order <- MFC.taxo$taxOrder[MFC.taxo$name %in% gsub("_"," ",correlation.study$Species)]

# Obj study
correlation.study$categories <- MFC.studies$Categories[MFC.studies$studies %in% correlation.study$Studies]


# -----

correlation.final <- rbind(correlation.final,
                           correlation.study)

connectivity.final <- rbind(connectivity.final,
                            connectivity.study)

# Save data
#print("--> Save data")
file_name <- paste0(saving_dir.studies.data,"/","connectivity_",diff.name,'.csv')
#write.xlsx(connectivity.study, file = file_name, rowNames = TRUE, colNames = TRUE)
write.table(connectivity.study, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)

file_name <- paste0(saving_dir.studies.data,"/","correlation_",diff.name,'.csv')
#write.xlsx(correlation.study, file = file_name, rowNames = TRUE, colNames = TRUE)
write.table(correlation.study, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)

