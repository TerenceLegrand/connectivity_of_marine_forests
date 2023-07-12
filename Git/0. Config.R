closeAllConnections()
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)

## ------------------------------------
## Files and folders

#project.name <- "oceanographicConnectivity" 
project.folder <- "../" 
connectivity.folder <- "../../The Biogeography of Oceanographic Connectivity/Results/"
data.folder <- paste0(project.folder,"Data/")
results.folder <- paste0(project.folder,"Results/")
temp.folder <- paste0(project.folder,"tempFolder/")

# -----------------------------------
# Main region

sim.resolution <- 5 # >= 5 for HR # https://h3geo.org/docs/core-library/restable/
min.lon <- -180
max.lon <- 180
min.lat <- -90
max.lat <- 90
dt.projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# -----------------------------------
# Source Sink sites

sourceSinkLocationType <- "centroid" # centroid or peripheral along H3 polygons
alternativeLandmass <- NULL # NULL or .shp file for specific landmass
removeLandmassSourceSinkSites <- FALSE ## landmass regions are unwanted Source Sink sites

maskSourceSinkSites <- "../Data/Spatial/innerSeas_save.shp" ## NULL or .shp file to mask Source Sink sites
maskSourceSinkSitesType <- "exclude" ## include or exclude Source Sink sites from sim

addSourceSinkRegions <- NULL # "../connectivityNEPacific/Data/sourceSinkSites.shp" ## NULL or .shp file for additional Points or Polygon regions
addSourceSinkRegionsForceCoast <- NULL ## added regions are strictly coastal distributed

addSourceSinkRegionsBathymetry <- NULL ## NULL or range
bathymetryRasterFile <- NULL ## NULL or .tif raster layer

# -----------------------------------
# Traits

months.all <- 1:12 # 
from.day <- 1 ; to.day <- 31
from.year <- 2008 ; to.year <- 2009 # 

allow.back.to.origin <- FALSE                     # at t == t.start

n.hours.per.day <- 12                             # Needs recoding for diferent than 12 # how many tracks for each particle during a day
n.new.particles.per.day <- 1
remove.new.particles.last.days <- FALSE           # If last days (particle.max.duration) are not to deliver new particles  
remove.new.particles.last.days.n.days <- 30

longevity <- TRUE
particle.max.duration <- 180                      # Days allowed to travel


# -----------------------------------
# -----------------------------------
# Theme of ggplot

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=13) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)) ,
                   legend.title = element_blank() ,
                   legend.margin=margin(c(0.3,1,0.3,1), unit='lines') ,
                   legend.background = element_rect(fill="white", linewidth=0.2, linetype="solid",  colour ="#979797"))

theme_map <- theme(plot.background = element_rect(fill = "#FFFFFF", color = NA),
                   panel.background = element_rect(fill = "#F8F8F8", color = NA),
                   panel.border = element_rect(colour = "black", fill  =NA, linewidth = 0.1),
                   text = element_text(size=8),
                   axis.text.x =  element_text(size=6),
                   axis.text.y =  element_text(size=6))

theme_plot <- theme(plot.background = element_rect(fill = "#FFFFFF", color = NA),
                   panel.background = element_rect(fill = "#F8F8F8", color = NA),
                   panel.border = element_rect(colour = "black", fill  =NA, linewidth = 0.1),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   text = element_text(size=8),
                   axis.text.x =  element_text(size=6),
                   axis.text.y =  element_text(size=6))

theme_boxplot <- theme(plot.background = element_rect(fill = "#FFFFFF", color = NA),
                    panel.background = element_rect(fill = "#F8F8F8", color = NA),
                    panel.border = element_rect(colour = "black", fill  =NA, linewidth = 0.1),
                    text = element_text(size=8),
                    axis.text.x =  element_text(size=6),
                    axis.text.y =  element_text(size=6))