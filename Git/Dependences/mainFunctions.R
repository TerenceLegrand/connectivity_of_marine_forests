## ------------------------------------------------------------------------------------------------------------------
## PlankTonic
## Assis et al., 2018
## ------------------------------------------------------------------------------------------------------------------

# Add h3 package below
#packageurl <- "http://cran.r-project.org/src/contrib/Archive/h3jsr/h3jsr_1.2.3.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

#install.packages("http://cran.r-project.org/src/contrib/Archive/h3jsr/h3jsr_1.2.3.tar.gz", repos=NULL, type="source")


 
packages.to.use <- c("rnaturalearth",
                     "factoextra",
                     "rnaturalearth",
                     "geosphere",
                     "openxlsx",
                     "rgeos",
                     "av",
                     "gdata",
                     "dplyr",
                     "sf",
                     "countrycode", 
                     "spatialEco", 
                     "geosphere",
                     "httr",
                     "gstat",
                     "fasterize",
                     "spdep",
                     "rworldxtra",
                     "rworldmap",
                     "cleangeo",
                     "compiler",
                     "data.table",
                     "raster",
                     "rgdal",
                     "ncdf4",
                     "parallel",
                     "doParallel",
                     "rgeos",
                     "FNN",
                     "sqldf",
                     "igraph",
                     "reshape2",
                     "gdistance",
                     "ggplot2",
                     "bigmemory",
                     "devtools",
                     "spatstat",
                     "dismo",
                     "digest",
                     "h3js",
                     #"h3r",
                     "h3", #"h3-r",
                     "reshape2",
                     "stringr",
                     "h3jsr",
                     "lme4", # For linear mixed-effects model
                     "MuMIn",
                     "patchwork",
                     "sn",
                     "ggspatial") # For plotting panels

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  
  sink("/dev/null") 
  
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }

  if( ! package %in% rownames(installed.packages()) & package == "h3js" ) { devtools::install_github("saurfang/h3js") }
  if( ! package %in% rownames(installed.packages()) & package == "h3" ) { devtools::install_github("crazycapivara/h3-r") }
  if( ! package %in% rownames(installed.packages()) & package == "h3r" ) { devtools::install_github("scottmmjackson/h3r") }
  if( ! package %in% rownames(installed.packages()) & package == "h3r" ) { devtools::install_github("harryprince/h3r", ref="bug-fix/Makefile") }
  if( ! package %in% rownames(installed.packages()) & package == "h3jsr" ) { remotes::install_github("obrl-soil/h3jsr") }
  if( ! package %in% rownames(installed.packages()) & package == "rnaturalearthhires" ) { devtools::install_github("ropensci/rnaturalearthhires")  }
  if( ! package %in% rownames(installed.packages()) & package == "rnaturalearth" ) { devtools::install_github("ropensci/rnaturalearth"); install.packages("rnaturalearthdata")  }
  if( ! package %in% rownames(installed.packages()) & package == "spatialEco" ) { remotes::install_github("jeffreyevans/spatialEco")  }
  if( ! package %in% rownames(installed.packages()) & ! package %in% rownames(installed.packages()) ) { sink() ; stop("Error on package instalation") }
  if( ! package %in% rownames(installed.packages()) & package == "h3jsr" ) { remotes::install_version("h3jsr", version = "1.2.3", repos = "http://cran.us.r-project.org") } # Install previous version

  library(package, character.only = TRUE)
  
  sink()
}


## ---------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------
## Functions

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  while (format(date, format="%m") == m) { date <- date + 1 }
  return(as.integer(format(date - 1, format="%d")))
}

## -------------------

distinctColors <- function(n) {
  
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector_long <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  if(n >= 433) { col_vector <- sample(rep(col_vector_long,round(n/433)+1), n,replace=FALSE) }
  if(n >= 74 & n < 433) { col_vector <- sample(col_vector_long, n,replace=FALSE) }
  if(n <= 74) { col_vector <- sample(col_vector, n,replace=FALSE) }
  
  return(col_vector)
  
}



reclassVals <- function(vals,valsMin,valsMax) {
  
  vals <- vals - valsMin
  vals <- vals / valsMax
  
  vals[ vals >= 0.5 ] <- 2.25
  vals[ vals < 0.5 & vals >= 0.1 ] <- 0.5
  vals[ vals < 0.1 & vals >= 0.01 ] <- 0.25
  vals[ vals < 0.01 & vals >= 0.001 ] <- 0.1
  vals[ vals < 0.001 & vals >= 0.0001 ] <- 0.05
  vals[ vals < 0.0001 ] <- 0.01
  return(vals)

}

## -------------------

assignIDs <- function(coords) {
  
  coordinates(coords) <- ~Lon+Lat
  crs(coords) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  
  assignments <- data.frame(Name=rep("",length(coords)))
  loopNonInt <- TRUE
  colors.plot <- character(length(assignments$Name))
  colors.plot[which( assignments$Name == "")] <- "black"
  
  while(loopNonInt) {
    
    if(sum(assignments$Name != "") > 0) { 
      
      colors.names <- unique(assignments$Name)
      colors.names <- colors.names[colors.names != ""]
      colors.n <- length(colors.names)
      colors.n.val <- distinctColors(colors.n)
      
      for( i in 1:colors.n) {
        
        colors.plot[which( assignments$Name == colors.names[i])] <- colors.n.val[i]
        
    } }
    
    plot(coords,pch=20,col=colors.plot)
    poly <- spatstat::clickpoly(add=TRUE)
    p = Polygon(cbind(poly$bdry[[1]]$x,poly$bdry[[1]]$y))
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    crs(sps) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
    sps <- as(sps, "SpatialPolygons")
    
    assigned <- which(!is.na(over(coords,sps)))
    name <- readline("Name? > ")
    assignments[assigned,1] <- name
    break.i <- readline("Assign more? [y/n] > ")
    if( break.i != "y" ) { loopNonInt <- FALSE }
    
  }
    
  return( assignments )
  
}

## -------------------

distinctColors <- function(n) {
  library(RColorBrewer)
  
  col_vector <- 0
  errorManagment <- 0
  while(length(unique(col_vector)) != n) {
    
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector = sample(col_vector, n)
    errorManagment <- errorManagment + 1
    if(errorManagment > 1000) { stop("Error :: Cannot define distinctColors")}
  }

  
  return(col_vector)
}

## -------------------

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}

## -------------------

bezierCurve <- function(x, y, n=10)
{
  outx <- NULL
  outy <- NULL
  
  i <- 1
  for (t in seq(0, 1, length.out=n))
  {
    b <- bez(x, y, t)
    outx[i] <- b$x
    outy[i] <- b$y
    
    i <- i+1
  }
  
  return (list(x=outx, y=outy))
}

bez <- function(x, y, t)
{
  outx <- 0
  outy <- 0
  n <- length(x)-1
  for (i in 0:n)
  {
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
  }
  
  return (list(x=outx, y=outy))
}

## ------------------------------

distinctColors <- function(n) {
  
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector_long <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  if(n > 74) { col_vector <- sample(col_vector_long, n,replace=FALSE) }
  if(n <= 74) { col_vector <- sample(col_vector, n,replace=FALSE) }

  return(col_vector)
  
}


## ----------------------------
## ----------------------------

getLocation <- function(coordLon,coordLat) {
  
  construct.geocode.url <- function(coordLon,coordLat, return.call = "json", sensor = "false") {
    
    root <- "http://www.mapquestapi.com/geocoding/v1/reverse?key=mpIx2AWq4Lj9R0mDbW1hNWrPe1Jju4X9&"
    u <- paste(root, "location=", coordLat,",",coordLon,"", sep = "") #&includeRoadMetadata=true&includeNearestIntersection=true
    return(URLencode(u))
  }
  
  u <- construct.geocode.url(coordLon,coordLat)
  doc <- getURL(u)
  x <- fromJSON(doc,simplify = FALSE)
  
  
  return(x$results[[1]]$locations[[1]]$adminArea1)
  
}


## Clean Dump

clean.dump.files <- function(clean.dump.files,files,dump.folder) {
  
  if( clean.dump.files ) { 
    file.remove( list.files(dump.folder, full.names = TRUE, pattern = files) ) 
  }
  
}

## ---------------------------------------------------------------------------------------------------------------------

list.memory.used <- function(type) {
  
  if(type==1) { return(sum(gc()[,7]) ) }
  if(type==2) { return(sum(list.memory()$Size)) }
  
}

## ---------------------------------------------------------------------------------------------------------------------

list.memory <- cmpfun( function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 100) {
  
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
  
} )

## ---------------------------------------------------------------------------------------------------------------------

errors <- cmpfun( function(t.step.test) {
  

} )

## ---------------------------------------------------------------------------------------------------------------------

monitor.processes <- cmpfun( function (process.name) {  
  if( !exists("ptm") ) { ptm <- proc.time() }
  t.step.time <- proc.time() - ptm
  ptm <- proc.time()
  
  if( !exists("day") ) { day <- 0 }

  cat("P:" ,process.name, " | Day: ",day, " | Time Taken: " , round(t.step.time[3]) , " | Memory: ", round(sum(list.memory()$Size)) , "\n") 
  flush.console()
  ptm <<- proc.time()
} )           


## ---------------------------------------------------------------------------------------------------------------------

trim.by.distance <- function(xyDF,source.sink.dist,parallel) {
  
  coastline.pts.t <- xyDF

  if(parallel){
    
    seqListing <- round(seq(min(coastline.pts.t$y),max(coastline.pts.t$y),length.out =number.cores))
    parallelChunks <- data.frame(from = c(min(coastline.pts.t$y),seqListing[-c(1,length(seqListing))]), to = c(seqListing[-c(1,length(seqListing))] , max(coastline.pts.t$y) ) )
    
    cl.2 <- makeCluster(number.cores)
    registerDoParallel(cl.2)
    
    source.sink.xy.t <- foreach(section=1:nrow(parallelChunks), .combine=rbind, .verbose=FALSE, .packages=c("dismo","gstat","gdata","raster","data.table","bigmemory","FNN")) %dopar% { 
    
      coastline.pts.t.sec <- coastline.pts.t[coastline.pts.t$y >= parallelChunks[section,1] & coastline.pts.t$y <= parallelChunks[section,2] ,]
    
      source.sink.xy.t <- data.frame(matrix(NA,ncol=2,nrow=nrow(coastline.pts.t.sec)))
      
      iteraction <- 0
      
      while( nrow(coastline.pts.t.sec) > 0 ){
        
        iteraction <- iteraction + 1
        
        pt.i = coastline.pts.t.sec[1,,drop=FALSE]
        
        source.sink.xy.t[iteraction,] <- as.data.frame(pt.i)
        
        coastline.pts.t.i <- coastline.pts.t.sec
        colnames(coastline.pts.t.i) <- c("x","y")
        coordinates(coastline.pts.t.i) <- c("x","y")
        crs(coastline.pts.t.i) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        circle <- circles(pt.i, lonlat=TRUE, d=source.sink.dist*1000, dissolve=FALSE)
        circle <- geometry(circle)
        crs(circle) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        
        to.extract <- which(!is.na(over(coastline.pts.t.i,circle)))
        coastline.pts.t.sec <- coastline.pts.t.sec[-to.extract,]
        
      }
      
      source.sink.xy.t <- source.sink.xy.t[complete.cases(source.sink.xy.t),]
      return(source.sink.xy.t)

    }
      
    stopCluster(cl.2) ; rm(cl.2)
    
    to.remove <- which(duplicated(source.sink.xy.t))
    if(length(to.remove) > 0) { source.sink.xy.t <- source.sink.xy.t[-to.remove,]}
    
  }
  
  if( ! parallel){
    
  source.sink.xy.t <- data.frame(matrix(NA,ncol=2,nrow=nrow(coastline.pts.t)))
    
  iteractions <- nrow(coastline.pts.t)
  iteraction <- 0
  
  while( nrow(coastline.pts.t) > 0 ){
    
    iteraction <- iteraction + 1
    progress.percent <- 100 - round((nrow(coastline.pts.t) / iteractions) * 100)
    
    cat('\014')
    cat('\n')
    
    cat('\n',paste0(rep("-",100),collapse = ""))
    cat('\n',paste0(rep("-",progress.percent),collapse = ""),"||",progress.percent,"% (",iteraction,")")
    cat('\n',paste0(rep("-",100),collapse = ""))
    
    pt.i = coastline.pts.t[1,,drop=FALSE]
    
    source.sink.xy.t[iteraction,] <- as.data.frame(pt.i)
    
    coastline.pts.t.i <- coastline.pts.t
    colnames(coastline.pts.t.i) <- c("x","y")
    coordinates(coastline.pts.t.i) <- c("x","y")
    crs(coastline.pts.t.i) <- dt.projection
    circle <- circles(pt.i, lonlat=TRUE, d=source.sink.dist*1000, dissolve=FALSE)
    circle <- geometry(circle)
    crs(circle) <- dt.projection
    
    to.extract <- which(!is.na(over(coastline.pts.t.i,circle)))
    coastline.pts.t <- coastline.pts.t[-to.extract,]
    
  }
  
  source.sink.xy.t <- source.sink.xy.t[complete.cases(source.sink.xy.t),]
  
  }
  
  return(source.sink.xy.t)
  
}

## ---------------------------------------------------------------------------------------------------------------------

trim.by.distance.poly <- function(poly,source.sink.dist) {
    
    polyT <- poly
    polyT <- as(polyT,"SpatialLines")
    crds <- coordinates(as(polyT, 'SpatialPoints'))
    source.sink.xy.t <- matrix(NA,nrow=nrow(crds),ncol=2)
    colnames(source.sink.xy.t) <- c("Lon","Lat")
    source.sink.xy.t <- data.frame(source.sink.xy.t)
    source.sink.xy.t[1,] <- crds[which.max(crds[,2]), ]
    
    i <- 0
    
    while( !is.null(polyT) ){
    
      i <- i + 1
      cat("\r","Processing point:",i)  
      
      circletCentroid <- which(!is.na(source.sink.xy.t[,1]))
      circle <- circles(source.sink.xy.t[circletCentroid[length(circletCentroid)],], lonlat=TRUE, d=source.sink.dist*1000, dissolve=FALSE)
      circle <- geometry(circle)
      crs(circle) <- dt.projection
      
      intersectionPolys <- gIntersection(polyT, circle)
      
      if( !is.null(intersectionPolys) ) {
        
        pt.i <- coordinates(as(intersectionPolys, 'SpatialPoints'))
        pt.i <- data.frame(Lon=pt.i[2,1],lat=pt.i[2,2])
        
      }
      
      if( is.null(intersectionPolys) ) {
        
        crds <- coordinates(as(polyT, 'SpatialPoints'))
        pt.i <- data.frame(Lon=crds[which.max(crds[,2]), 1],lat=crds[which.max(crds[,2]), 2])

      }
      
      source.sink.xy.t[i,] <- pt.i
      polyT = gDifference(polyT,circle,byid=TRUE)
      gc(reset=TRUE)
      
    }

    source.sink.xy.t <- source.sink.xy.t[complete.cases(source.sink.xy.t),]

  return(source.sink.xy.t)
  
}

## ---------------------------------------------------------------------------------------------------------------------

# network.type = "Prob"
# comb = distance.probability
# crop.network = FALSE
# buffer = 10
# cells = source.sink.xy
# new.extent

produce.network <- function(network.type,comb,n.days,crop.network,buffer,cells,new.extent) {
  
  if(crop.network) {  final.cells <- which(   cells[,2] >= (new.extent[1] - buffer) & 
                                                cells[,2] <= (new.extent[2] + buffer) & 
                                                cells[,3] >= (new.extent[3] - buffer) & 
                                                cells[,3] <= (new.extent[4] + buffer) )   
  
  final.cells <- cells[final.cells,1]
  plot(cells[final.cells,2:3])
  
  }
  
  if( ! crop.network ) {  final.cells <- cells[,1]  }
  
  comb[ which(comb[,"Time.max"] > n.days) , "Probability" ] <- 0
  comb <- comb[,c("Pair.from","Pair.to","Probability")]
  comb <- comb[comb$Pair.from %in% as.vector(unlist(final.cells)) & comb$Pair.to %in% as.vector(unlist(final.cells)) ,]
  # comb <- comb[comb$Pair.from != comb$Pair.to,]
  comb <- as.data.frame( comb[ sort( as.vector(unlist(comb[,"Probability"])) , decreasing = TRUE, index.return =TRUE)$ix , ] )
  
  if( network.type == "Prob" ) {
    
    net.function <<- prod
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    
    # E(graph.obj)$weight = 1 - comb[,3] # The wheight has a negative impact on finding the closest path
    E(graph.obj)$weight = ifelse(-log(comb[,3]) == Inf,0,-log(comb[,3])) # Hock, Karlo Mumby, Peter J 2015
    graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
    graph.obj <- simplify(graph.obj)
    
  }
  
  if( network.type == "Time" ) {
    
    net.function <<- sum
    graph.obj <- graph.edgelist( cbind( as.character( comb[,1]) , as.character(comb[,2]) ) , directed = TRUE )
    E(graph.obj)$weight = comb[,3]
    graph.obj <- delete.edges(graph.obj, which(E(graph.obj)$weight ==0))
    graph.obj <- simplify(graph.obj)
    
  }
  
  return(list(comb,graph.obj))
  
}

## ---------------------------------------------------------------------------------------------------------------------

padlock <- function(dir,fun,id) {
  
  if( missing(id) ) { id <- NULL } 
  if( missing(dir) ) { dir <- NULL } 
  if( missing(fun) ) { fun <- NULL } 
  
  options(warn=-1)
  
  file.t <- paste0("Padlocker#",id,".Lk")
  
  ## ---------------------------
  
  if( fun == "LockAndHaltOn" ) {
  
    int <- 0
    
    repeat {
      
      if( ! padlock(dir,"isLocked") ) { 
        
        int <- int + 1
        
        padlock(dir,"Lock",id=id)
        
        Sys.sleep(sample(seq(1,2,length.out = 10),1))

        if( ! padlock(dir,"uniqueLocker",id=id) ) { padlock(dir,"Unlock",id=id) ; next }
        
        if( padlock(dir,"uniqueLocker",id=id) ) { break }
        
      }
      
      if( int > 9999 ) { stop("More than 9999 tries") }
      
  }
    
  }
  
  ## ---------------------------
    
  if( fun == "Lock" ) {
    
    write("Locked",file=paste0(dir,"/",file.t),append=FALSE)
    
  }
  
  ## ---------------------------
  
  if( fun == "Unlock" ) {
    
    file.remove( paste0(dir,"/",file.t) )
    
    if( id == -1 ) {
      
      file.remove( list.files(dir,"Padlocker#",full.names = TRUE) )
      
    }
  }
  
  ## ---------------------------
  
  if( fun == "isLocked" ) {
    
    Locked <- TRUE
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    
    if( length(fs) != 0 ) {
      Locked <- TRUE
    } 
    if( length(fs) == 0 ) {
      Locked <- FALSE
    } 
    
    return( Locked  )
    
  }
  
  ## ---------------------------
  
  if( fun == "uniqueLocker" ) {
    
    u.locker <- FALSE
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    
    if( length(fs) == 0) { u.locker <- FALSE  }
    if( length(fs) > 1) { u.locker <- FALSE  }
    if( length(fs) == 1 ) { 
      
      if( grepl(file.t,fs) ) { u.locker <- TRUE  }  
      
    }
    
    return( u.locker  )
    
  }
  
  ## ---------------------------
  
  if( fun == "countLockersAge" ) {
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    age <- as.numeric( difftime(Sys.time(), file.info(fs )$atime, units ="mins") )
    return( data.frame(fs,age)  )
    
  }
  
  ## ---------------------------
  
  if( fun == "Age" ) {
    
    fs <- list.files(dir,"Padlocker#",full.names = TRUE)
    age.i <- numeric(0)
    
    if( length(fs) > 0 ) {
      
      for( i in 1:length(fs)) {
        
        age <- as.numeric( difftime(Sys.time(), file.info(fs[i] )$atime, units ="mins") )
        age.i <- c(age.i,age)
        
      }
    }
    if( length(fs) == 0 ) {
      age.i <- 0
    }
    
    return(min(age.i))
    
  }
  
  ## ---------------------------
  
  options(warn=0)
  
}

## ---------------------------------------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------------------------------------
