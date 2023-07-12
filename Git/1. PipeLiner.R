rm( list=(ls()[ls()!="v"]) )
gc(reset=TRUE)
closeAllConnections()

source("0. Config.R")
source("Dependences/mainFunctions.R")


# ----

pipeLiner <- TRUE

ecologicalGroup <- "MarineForest"
project.name.results <- "GEB_results"

# load results from the Biophysical model (https://github.com/jorgeassis/biophysicalModelling)
load(paste0(connectivity.folder,"modelSummaryResults.RData"))
n.days.max <- global.simulation.parameters.res$particle.max.duration

max.PD <- n.days.max

spawn.p <- 1:12  # spawn.p <- c(6,7,8,9)
n.pld.period <- c(0,7) # O is considering PD distribution, 7 is considering mean in-house dispersal
n.seasons <- "" # c("YearRound","Spring","Summer","Autumn","Winter")
combinations <- expand.grid(season=n.seasons,pld.period=n.pld.period,stringsAsFactors = F)

# If yes, it extrapols stacked distribution information to source/sink hexagons data. (it takes a will, should be avoid if already did)
extractDistrib <- FALSE

# If yes, it use the stacked distribution of the ecological group to compute shortest path
useDistrib <- TRUE

# If yes, Fst is transformed into Fst/(1-Fst)
transformFST <- FALSE

# If yes, it removes pair of sites with differentiation == 0 and below
pruneData <- FALSE

# ------------------------
# Folder index

saving_dir <- file.path(results.folder,project.name.results)
dir.create(saving_dir, showWarnings = FALSE)

# Save general pipeline info into project folder

pipe.info <- data.frame(
  spawn.p = paste(spawn.p, collapse = "-"),
  n.pld.period = paste(n.pld.period, collapse = "-"),
  n.seasons = n.seasons,
  extractDistrib = extractDistrib,
  useDistrib = useDistrib,
  transformFST = transformFST,
  pruneData = pruneData
)
write.table(pipe.info, file = file.path(saving_dir,"PipeLiner_info.txt"), append = FALSE, sep = ",", dec = ".",
            row.names =FALSE, col.names = TRUE)

# ------------------------
# Load data

MFC.species <- read.xlsx(file.path(data.folder,"MFC_species_table_save.xlsx"), check.names=FALSE)
MFC.studies <- read.xlsx(file.path(data.folder,"MFC_studies_table_save.xlsx"), check.names=FALSE)
MFC.taxo <- read.xlsx(file.path(data.folder,"MFC_species_taxonomy.xlsx"), check.names=FALSE)

# Download GIS data

distrib.range=raster(file.path(data.folder,"SDM/Global/distributionRange.tif"))
distrib.richness=raster(file.path(data.folder,"SDM/Global/richnessDistribution.tif"))

# Dowload Spalding realm data
realm.shp <- shapefile(paste0(data.folder,"marine_realm/MEOW/meow_ecos.shp"))

# ------------------------
# Load cells

if( is.null(alternativeLandmass) ) { worldMap <- getMap(resolution = "high") } # TL: Modify in 0. Config
if( ! is.null(alternativeLandmass) ) { worldMap <- shapefile(alternativeLandmass) }

load(paste0(connectivity.folder,"sourceSinkSites.RData")) # TL: file resuming coordinates of cells center
colnames(source.sink.xy) <- c("Pair" , "Lon" , "Lat" , "Source" )
source.sink.xy <- source.sink.xy[source.sink.xy$Source == 1,]

hexagons.sourcesink.shp <- shapefile(paste0(connectivity.folder,"/sourceSinkSites.shp")) # TL: file resuming coordinates of cells in spatial polygons format (spatial data)
if( "SOURCE" %in% names(hexagons.sourcesink.shp) ) { hexagons.sourcesink.shp <- hexagons.sourcesink.shp[hexagons.sourcesink.shp$SOURCE==1,"ID"] } 

## ----
# Extrapole stacked Marine forest SDM to source/sink level

if (extractDistrib) {
  # We set a buffer of 15 km since the average edge length is ~ 10 km for hexagon = 5. So it will get info from the 6 neighbor hexagons
  # We need to do that to solve the missmatch between hex implementation and raster distribution 
  hex.distrib.richness <- matrix(0,nrow = length(hexagons.sourcesink.shp), ncol = 1)
  for (k in 1:length(hexagons.sourcesink.shp)) {
    print(paste0(" Hexagon loop ---> ",as.character(round((k/length(hexagons.sourcesink.shp)*100),3)),"%"))
    hex_i <- hexagons.sourcesink.shp@polygons[[k]]@Polygons[[1]]@coords
    hex_i <- hex_i[-1,]
    hex.distrib.richness[k,] <- mean(extract(distrib.richness, hex_i, buffer = 20000, weights = TRUE, fun = mean, na.rm = TRUE))
  }
  hex.distrib.range <- hexagons.sourcesink.shp$ID[which(hex.distrib.richness > 0)]
  hex.distrib.richness <- data.frame(id = hex.distrib.range, richness = hex.distrib.richness[hex.distrib.range])
  save(hex.distrib.range, file = file.path(data.folder,"SDM/Global/distributionRange.Rdata"))
  save(hex.distrib.richness, file = file.path(data.folder,"SDM/Global/richnessDistribution.Rdata"))
} else if (!extractDistrib) {
  load(file.path(data.folder,"SDM/Global/richnessDistribution.Rdata"))
  load(file.path(data.folder,"SDM/Global/distributionRange.Rdata"))
}

# Get distrib source.sink
if (useDistrib) {
source.sink.xy.distrib <- source.sink.xy[source.sink.xy$Pair %in% hex.distrib.range,]
} else if (!useDistrib) {
  source.sink.xy.distrib <- source.sink.xy
}

# Map
robinson <- CRS("+proj=robin +over")

if( is.null(alternativeLandmass) ) { worldMap <- getMap(resolution = "high") }
if( ! is.null(alternativeLandmass) ) { worldMap <- shapefile(alternativeLandmass) }

# ------------------------
# Loop over [1] PLD [2] Species [3] studies [4] Marker [5] Differentiation index [6] Specificity 

for (pld.i in 1:length(n.pld.period)) {
  
  pld.period <- n.pld.period[pld.i]
  
  saving_dir.PD <- paste0(saving_dir,"/PD_",as.character(pld.period))
  dir.create(saving_dir.PD, showWarnings = FALSE)
  
  saving_dir.studies <- file.path(saving_dir.PD,"studies")
  saving_dir.studies.data <- file.path(saving_dir.studies,"data")
  saving_dir.studies.figure <- file.path(saving_dir.studies,"figure")
  
  dir.create(saving_dir.studies, showWarnings = FALSE)
  dir.create(saving_dir.studies.data, showWarnings = FALSE)
  dir.create(saving_dir.studies.figure, showWarnings = FALSE)
  
  saving_dir.paper <- file.path(saving_dir.PD,"paper")
  saving_dir.paper.main <- file.path(saving_dir.paper,"main")
  saving_dir.paper.SI <- file.path(saving_dir.paper,"SI")
  
  dir.create(saving_dir.paper, showWarnings = FALSE)
  dir.create(saving_dir.paper.main, showWarnings = FALSE)
  dir.create(saving_dir.paper.SI, showWarnings = FALSE)
  
  # Implement dataframe
  connectivity.final <- data.frame()
  correlation.final <- data.frame()
  
  # ------------------------
  # Compute backward indirected graphs
  
  source("6.2. Species Connectivity Computation.R")
  
  # ------------------------
  #k = s = i = j = e = 1 for test
  for (k in 1:length(MFC.species$species)) { # loop on different species
    
    print(paste0(" SPECIES ---> ",as.character(round((k/length(MFC.species$species)*100),3)),"%",
                 " --> ",as.character(k),"/",as.character(length(MFC.species$species))," -- ",MFC.species$species[k],
          " -- PD = ",as.character(pld.period)))
    
    species <- gsub(" ","_",MFC.species$species[k])
    
    species.studies <- unlist(strsplit(MFC.species$studies_info[k],","))
    
    # Create export result folders 
    
    # saving_dir <- paste0(results.folder,"Correlation/","Species","/",species)
    # dir.create(saving_dir)
    
    
    for (s in 1:length(species.studies)) { # loop on different studies per species
      
      print(paste0("              STUDIES ---> ",as.character(round((s/length(species.studies)*100),3)),"%",
                   " --> ",as.character(s),"/",as.character(length(species.studies))," -- ",species.studies[s]))
      
      # Create export result folders
      # saving_dir <- paste0(results.folder,"Correlation/","Species","/",species,"/",species.studies[s])
      # dir.create(saving_dir)
      
      studies.marker <- unlist(strsplit(MFC.studies$diff_marker_info[MFC.studies$studies %in% species.studies[s]],","))
      
      studies.diffi <- unlist(strsplit(MFC.studies$diff_index_info[MFC.studies$studies %in% species.studies[s]],","))
      
      studies.spec <- unlist(strsplit(MFC.studies$diff_study_spec[MFC.studies$studies %in% species.studies[s]],","))
      studies.spec[!is.na(studies.spec)] <- paste0(studies.spec[!is.na(studies.spec)],"_")
      studies.spec[is.na(studies.spec)] <- ""
      studies.spec <- paste0("_",studies.spec)
      
      for (i in 1:length(studies.marker)) { # loop on different marker per studies per species
        
        print(paste0("                            MARKER ---> ",as.character(round((i/length(studies.marker)*100),3)),"%",
                     " --> ",as.character(i),"/",as.character(length(studies.marker))," -- ",studies.marker[i]))
        
        global.name <- paste0(species.studies[s],"_",species,"_",studies.marker[i])
        popCoordinates <- read.xlsx(paste0(data.folder,"/Global/",global.name,".xlsx"), check.names=FALSE)

        
        # Find duplicated sample sites and remove them (i.e. it happened when different marker spec is used for a same sampling scheme)
        popCoordinates <- popCoordinates[! duplicated(popCoordinates$Sample.Code),]
        
        # FIND SAMPLING CELLS
        
        ## source sink based on distrib range
        
        pop.distance.cells <- spDists(as.matrix(source.sink.xy.distrib[,2:3]),as.matrix(popCoordinates[,c("Lon","Lat")]),longlat = TRUE)
        #pop.distance.cells <- spDists(as.matrix(source.sink.xy[,2:3]),as.matrix(popCoordinates[,c("Lon","Lat")]),longlat = TRUE)
        
        pop.cells <- source.sink.xy.distrib[apply(pop.distance.cells,2,which.min),]
        #pop.cells <- source.sink.xy[apply(pop.distance.cells,2,which.min),]
        
        pop.cells$Sample.Code <- popCoordinates$Sample.Code
        pop.cells$Sample.Lon <- popCoordinates$Lon
        pop.cells$Sample.Lat <- popCoordinates$Lat
        
        for (j in 1:length(studies.diffi)) { # Loop on different diff index per marker per studies per species
          
          print(paste0("                                       DIFF INDEX ---> ",as.character(round((j/length(studies.diffi)*100),3)),"%",
                       " --> ",as.character(j),"/",as.character(length(studies.diffi))," -- ",studies.marker[j]))
          
          for (e in 1:length(studies.spec)) { # Loop on different study specificity per diff index per marker per studies per species
            
            
            diff.name <- paste0(species.studies[s],"_",species,"_",studies.marker[i],studies.spec[e],studies.diffi[j])
            popDifferentiation <- read.xlsx(paste0(data.folder,"Differentiation/",diff.name,".xlsx"),
                                            rowNames = TRUE,
                                            check.names= FALSE)
            
            
            pop.cells.dataset <- pop.cells[pop.cells$Sample.Code %in% colnames(popDifferentiation),]
            
            # Test to see if the data is OK
            if (sum(! rownames(popDifferentiation) %in% colnames(popDifferentiation)) !=0) {print(paste0(diff.name," --> Error 601"))} # Test if rowname and colname of differentiation matrix are differents
            if (sum(! colnames(popDifferentiation) %in% rownames(popDifferentiation) ) !=0) {print(paste0(diff.name," --> Error 601"))} # Test if rowname and colname of differentiation matrix are differents
            if (sum(! rownames(popDifferentiation) %in% popCoordinates$Sample.Code) !=0) {print(paste0(diff.name," --> Error 671"))} # Test to find if differentiation site names are contained in global files
            if (sum( duplicated(popCoordinates$Sample.Code)) !=0) {print(paste0(diff.name," --> Error 731"))} # Test to find if there are duplicate site names
            if (sum( is.na(popCoordinates$Sample.Code)) !=0) {print(paste0(diff.name," --> Error 231"))} # Test to find if there are NAs site names
            
            # -------------------------------- NOW, inside the loop with selected - species - study - marker - diff index - study specification
        
            source("6.2. Species Connectivity vs Differentiation.R")
            
            # ----------------------------
            
          }  # Closing loop on different study specificity per diff index per marker per studies per species
          
        } # Closing loop on different diff index per marker per studies per species
        
      } # Closing loop on different marker per studies per species
      
    } # Closing loop on different studies per species
    
  } # Closing loop on different species
  
  file_name <- paste0(saving_dir.PD,'/connectivity_final','.csv')
  write.table(connectivity.final, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)
  
  file_name <- paste0(saving_dir.PD,'/correlation_final','.csv')
  write.table(correlation.final, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)
  
}
