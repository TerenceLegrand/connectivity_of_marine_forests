# Choose PLD
pld.period = 0

# Saving dir

saving_dir.PD <- paste0(saving_dir,"/PD_",as.character(pld.period))
saving_dir.studies <- file.path(saving_dir.PD,"studies")
saving_dir.studies.data <- file.path(saving_dir.studies,"data")
saving_dir.studies.figure <- file.path(saving_dir.studies,"figure")
saving_dir.paper <- file.path(saving_dir.PD,"paper")
saving_dir.paper.main <- file.path(saving_dir.paper,"main")
saving_dir.paper.SI <- file.path(saving_dir.paper,"SI")

# Load results

connectivity.final <- read.csv(paste0(saving_dir.PD,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final <- read.csv(paste0(saving_dir.PD,"/correlation_final",".csv"), check.names=FALSE)
centrality <- read.csv(paste0(saving_dir.PD,"/centrality_distrib",".csv"), check.names=FALSE)

# For fixed pd
pld.period = 7

saving_dir.PD.fixed <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.pd.fixed  <- read.csv(paste0(saving_dir.PD.fixed,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.pd.fixed  <- read.csv(paste0(saving_dir.PD.fixed,"/correlation_final",".csv"), check.names=FALSE)

## -------------------------
# Process results

# Info

number.species <- length(unique(correlation.final$Species))

sampling_coord <- data.frame(Lon = c(connectivity.final$sample_Lon_from,connectivity.final$sample_Lon_to),
                             Lat = c(connectivity.final$sample_Lat_from,connectivity.final$sample_Lat_to))
sampling_coord <- unique(sampling_coord)

number.sampled.population <- length(sampling_coord$Lon)

correlation.final$Marker %>%  as.factor() %>% table() %>% `/` (length(correlation.final$Marker))

correlation.final$Diff.index %>%  as.factor() %>% table() %>% `/` (length(correlation.final$Diff.index))

# Nbr of pop

to_test <- correlation.final$nbr.pop
mean(to_test,na.rm = TRUE)
median(to_test,na.rm = TRUE)
ic_n <- length(! is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

# Sampling distance

to_test <- correlation.final$sampling_mean_distance/1000
mean(to_test,na.rm = TRUE)
median(to_test,na.rm = TRUE)
ic_n <- length(! is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

# Realm

(sum(correlation.final$REALM_nbr == 1)*100)/length(correlation.final$REALM_nbr)

(sum(correlation.final$REALM == 1)*100)/length(correlation.final$REALM_nbr)

(length(grep("Temperate Northern Atlantic",correlation.final$REALM))*100)/length(correlation.final$REALM_nbr)
(length(grep("Temperate Northern Pacific",correlation.final$REALM))*100)/length(correlation.final$REALM_nbr)
(length(grep("Temperate Australasia",correlation.final$REALM))*100)/length(correlation.final$REALM_nbr)


# % of connection allowed by LDD

correlation.final$nbr.Pairs.Connected.LDD <- round((((correlation.final$nbr.Pairs.Connected-correlation.final.pd.fixed$nbr.Pairs.Connected)*100)
                                                    /correlation.final$nbr.Pairs.Connected),1)


# First filter, remove the one that are NA
correlation.final.raw <- correlation.final
correlation.final.pd.fixed.raw <- correlation.final.pd.fixed
sum(! is.na(correlation.final$aic.SM))/length(correlation.final$aic.SM)
correlation.final <- correlation.final[! is.na(correlation.final$aic.SM),]
correlation.final.pd.fixed <- correlation.final.pd.fixed [! is.na(correlation.final.pd.fixed$aic.SM),]


# Transform NA to p-value = 1 and r² = 0
# Transform negative r² to r² = 0

correlation.final$p.SM[is.na(correlation.final$p.SM)] <- 1
correlation.final$p.CM[is.na(correlation.final$p.CM)] <- 1
correlation.final$p.CCM[is.na(correlation.final$p.CM)] <- 1

correlation.final$r2.CM[is.na(correlation.final$r2.CM)] <- 0
correlation.final$r2.CCM[is.na(correlation.final$r2.CCM)] <- 0

correlation.final$r2.CM[correlation.final$r2.CM < 0] <- 0
correlation.final$r2.CCM[correlation.final$r2.CCM < 0] <- 0

# For fixed

correlation.final.pd.fixed$r2.CM[is.na(correlation.final.pd.fixed$r2.CM)] <- 0
correlation.final.pd.fixed$r2.CCM[is.na(correlation.final.pd.fixed$r2.CCM)] <- 0

correlation.final.pd.fixed$r2.SM[correlation.final$r2.diff.null < 0] <- 0
correlation.final.pd.fixed$r2.CM[correlation.final.pd.fixed$r2.CM < 0] <- 0
correlation.final.pd.fixed$r2.CCM[correlation.final.pd.fixed$r2.CCM < 0] <- 0

correlation.final.pd.fixed$r2.SM[correlation.final.pd.fixed$r2.diff.null < 0] <- 0
correlation.final.pd.fixed$r2.CM[correlation.final.pd.fixed$r2.CM < 0] <- 0
correlation.final.pd.fixed$r2.CCM[correlation.final.pd.fixed$r2.CCM < 0] <- 0

# Which one are no-significative
correlation.final[!(!is.na(correlation.final$p.SM) & correlation.final$p.SM < 0.05),]
nbr_signif <- sum((!is.na(correlation.final$p.SM) & correlation.final$p.SM < 0.05))
(nbr_signif*100)/length(correlation.final$Species)

correlation.final[!(!is.na(correlation.final$p.CM) & correlation.final$p.CM < 0.05),]
nbr_signif <- sum((!is.na(correlation.final$p.CM) & correlation.final$p.CM < 0.05))
(nbr_signif*100)/length(correlation.final$Species)

correlation.final[!(!is.na(correlation.final$p.CCM) & correlation.final$p.CCM < 0.05),]
nbr_signif <- sum((!is.na(correlation.final$p.CCM) & correlation.final$p.CCM < 0.05))
(nbr_signif*100)/length(correlation.final$Species)

# Which one is the best model with AIC
to_test <- exp((apply(correlation.final[,c("aic.SM","aic.CM","aic.CCM")], 1, FUN = min)-correlation.final$aic.SM)/2)
mean(to_test)
ic_n <- length(to_test)
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

nbr_lowest_AIC <- sum(apply(correlation.final[,c("aic.SM","aic.CM","aic.CCM")], 1, FUN = min) == correlation.final$aic.SM)
(nbr_lowest_AIC*100)/length(correlation.final$aic.SM)

to_test <- exp((apply(correlation.final[,c("aic.SM","aic.CM","aic.CCM")], 1, FUN = min)-correlation.final$aic.CM)/2)
mean(to_test)
ic_n <- length(to_test)
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

nbr_lowest_AIC <- sum(apply(correlation.final[,c("aic.SM","aic.CM","aic.CCM")], 1, FUN = min) == correlation.final$aic.CM)
(nbr_lowest_AIC*100)/length(correlation.final$aic.CM)

to_test <- exp((apply(correlation.final[,c("aic.SM","aic.CM","aic.CCM")], 1, FUN = min)-correlation.final$aic.CCM)/2)
mean(to_test)
ic_n <- length(to_test)
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

nbr_lowest_AIC <- sum(apply(correlation.final[,c("aic.SM","aic.CM","aic.CCM")], 1, FUN = min) == correlation.final$aic.CCM)
(nbr_lowest_AIC*100)/length(correlation.final$aic.CCM)

nbr_lowest_AIC <- sum(apply(correlation.final[,c("aic.SM","aic.CCM")], 1, FUN = min) == correlation.final$aic.CCM)
(nbr_lowest_AIC*100)/length(correlation.final$aic.CCM)


# Global results

# Null
# r² 

to_test <- correlation.final$r2.SM

mean(to_test,na.rm = TRUE)
median(to_test,na.rm = TRUE)

ic_n <- length(! is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

# Connectivity
# r² 

to_test <- correlation.final$r2.CM
mean(to_test,na.rm = TRUE)
median(to_test,na.rm = TRUE)

ic_n <- length(! is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

# Centrality
# r² 

to_test <- correlation.final

mean(to_test$r2.CCM,na.rm = TRUE)
median(to_test$r2.CCM,na.rm = TRUE)

ic_n <- length(! is.na(to_test$r2.CCM))
ic_sd = sd(to_test$r2.CCM, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)


to_test <- correlation.final[(!is.na(correlation.final$p.CCM) & correlation.final$p.CCM < 0.05),]
min(to_test$r2.CCM)
to_test[which.min(to_test$r2.CCM),]

max(to_test$r2.CCM)
to_test[which.max(to_test$r2.CCM),]

# Diff with r²

correlation.final$delta_r2 <- correlation.final$r2.CCM-correlation.final$r2.SM

(sum(correlation.final$delta_r2 > 0)*100)/length(correlation.final$delta_r2)

to_test <- correlation.final
mean(to_test$delta_r2,na.rm = TRUE)
median(to_test$delta_r2,na.rm = TRUE)

quantile(to_test$delta_r2)
quantile(to_test$delta_r2,0.30)
quantile(to_test$delta_r2,0.9)

ic_n <- length(! is.na(to_test$delta_r2))
ic_sd = sd(to_test$delta_r2, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

to_test <- correlation.final[(!is.na(correlation.final$p.CCM) & correlation.final$p.CCM < 0.05),]
min(to_test$delta_r2)
to_test[which.min(to_test$delta_r2),]
max(to_test$delta_r2)
to_test[which.max(to_test$delta_r2),]

# Long distance dispersal

correlation.final$r2.improve.LDD <- correlation.final$r2.CCM-correlation.final.pd.fixed$r2.CCM

(sum(correlation.final$r2.improve.LDD> 0)*100)/length(correlation.final$r2.improve.LDD)

to_test <- correlation.final$r2.improve.LDD
mean(to_test,na.rm = TRUE)
median(to_test,na.rm = TRUE)
ic_n <- length(! is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)


AIC_disp <- data.frame(fixed = correlation.final.pd.fixed$aic.CCM,
                       LDD = correlation.final$aic.CCM)

# Fixed
to_test <- exp((apply(AIC_disp, 1, FUN = min)-AIC_disp$fixed)/2)
mean(to_test, na.rm = TRUE)
ic_n <- sum(!is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

nbr_lowest_AIC <- sum(apply(AIC_disp, 1, FUN = min) == AIC_disp$fixed, na.rm = TRUE)
(nbr_lowest_AIC*100)/sum(!is.na(to_test))

# LDD
to_test <- exp((apply(AIC_disp, 1, FUN = min)-AIC_disp$LDD)/2)
mean(to_test, na.rm = TRUE)
ic_n <- sum(!is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

nbr_lowest_AIC <- sum(apply(AIC_disp, 1, FUN = min) == AIC_disp$LDD, na.rm = TRUE)
(nbr_lowest_AIC*100)/sum(!is.na(to_test))

correlation.final$nbr.Pairs.Connected

(sum(!is.infinite(connectivity.final$nbr_step))*100)/length(connectivity.final$nbr_step)
(sum(!is.infinite(connectivity.final.pd.fixed$nbr_step))*100)/length(connectivity.final.pd.fixed$nbr_step)

((sum(!is.infinite(connectivity.final$nbr_step)) - sum(!is.infinite(connectivity.final.pd.fixed$nbr_step)))*100)/length(connectivity.final$nbr_step)

# Mean Fst

to_test <- connectivity.final$Differentiation[connectivity.final$Diff.index == "Fst"]
mean(to_test, na.rm = TRUE)
ic_n <- sum(!is.na(to_test))
ic_sd = sd(to_test, na.rm = TRUE)
qt(0.975,df=ic_n-1)*ic_sd /sqrt(ic_n)

# ---------------------------------
# Make Figure

# Figure 1
source("Figure/Figure_1.R")

# Figure 2
source("Figure/Figure_2.R")

# Figure 3
source("Figure/Figure_3.R")

# Figure 4
source("Figure/Figure_4.R")

# ---------------------------------
# Make SI

correlation.final.raw$p.SM[is.na(correlation.final.raw$p.SM) & ! is.na(correlation.final.raw$aic.SM)] <- 1
correlation.final.raw$p.CM[is.na(correlation.final.raw$p.CM) & ! is.na(correlation.final.raw$aic.SM)] <- 1
correlation.final.raw$p.CCM[is.na(correlation.final.raw$p.CM) & ! is.na(correlation.final.raw$aic.SM)] <- 1

correlation.final.raw$r2.SM[is.na(correlation.final.raw$r2.diff.null) & ! is.na(correlation.final.raw$aic.SM)] <- 0
correlation.final.raw$r2.CM[is.na(correlation.final.raw$r2.CM) & ! is.na(correlation.final.raw$aic.SM)] <- 0
correlation.final.raw$r2.CCM[is.na(correlation.final.raw$r2.CCM) & ! is.na(correlation.final.raw$aic.SM)] <- 0

correlation.final.raw$r2.SM[correlation.final.raw$r2.diff.null < 0 & ! is.na(correlation.final.raw$aic.SM)] <- 0
correlation.final.raw$r2.CM[correlation.final.raw$r2.CM < 0 & ! is.na(correlation.final.raw$aic.SM)] <- 0
correlation.final.raw$r2.CCM[correlation.final.raw$r2.CCM < 0 & ! is.na(correlation.final.raw$aic.SM)] <- 0

correlation.final.raw$species_short <- sub("_", " ", correlation.final.raw$Species)

correlation.final.raw[,c("p.SM","p.CM","p.CCM","r2.SM","r2.CM","r2.CCM",
                         "aic.SM","aic.CM","aic.CCM")] <- round(correlation.final.raw[,c("p.SM","p.CM","p.CCM","r2.SM","r2.CM","r2.CCM",
                                                                                         "aic.SM","aic.CM","aic.CCM")],2)

correlation.final.raw[,c("sampling_max_extent","sampling_mean_distance")] <- round(correlation.final.raw[,c("sampling_max_extent","sampling_mean_distance")]/1000,2)

correlation.final.raw$Marker[correlation.final.raw$Marker == "microsatellites"] <- "msat"

correlation.final.raw <- correlation.final.raw[order(correlation.final.raw$Species),]

# Table S2

Table_S2 <- data.frame(Species = correlation.final.raw$Species,
                       Marker = correlation.final.raw$Marker,
                       Diff.Index = correlation.final.raw$Diff.index,
                       Spec = correlation.final.raw$Specificities,
                       n = correlation.final.raw$nbr.pop,
                       Samp.Ext = correlation.final.raw$sampling_max_extent,
                       Samp.Mean.Dist = correlation.final.raw$sampling_mean_distance,
                       Lat.Ext = correlation.final.raw$lat_extent,
                       Realm = correlation.final.raw$REALM.cl,
                       Reference = correlation.final.raw$ref)

file_name <- paste0(saving_dir.paper.SI,'/Table_S2','.csv')
write.table(Table_S2, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)

# Table S3

Table_S3 <- data.frame(Species = correlation.final.raw$Species,
                                    Marker = correlation.final.raw$Marker,
                                    Diff.Index = correlation.final.raw$Diff.index,
                                    Spec = correlation.final.raw$Specificities,
                                    Pairs = correlation.final.raw$nbr.Pairs,
                                    Pairs.Co = correlation.final.raw$nbr.Pairs.Connected,
                                    AIC.SM = correlation.final.raw$aic.SM,
                                    AIC.CM = correlation.final.raw$aic.CM,
                                    AIC.CCM = correlation.final.raw$aic.CCM,
                                    pval.SM = correlation.final.raw$p.SM,
                                    pval.CM = correlation.final.raw$p.CM,
                                    pval.CCM = correlation.final.raw$p.CCM,
                                    rsquared.SM = correlation.final.raw$r2.SM,
                                    rsquared.CM = correlation.final.raw$r2.CM,
                                    rsquared.CCM = correlation.final.raw$r2.CCM,
                                    Reference = correlation.final.raw$ref)

file_name <- paste0(saving_dir.paper.SI,'/Table_S3','.csv')
write.table(Table_S3, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)

# Table S4

correlation.final.pd.fixed.raw$r2.CM[is.na(correlation.final.pd.fixed.raw$r2.CM) & ! is.na(correlation.final.pd.fixed.raw$aic.SM)] <- 0
correlation.final.pd.fixed.raw$r2.CCM[is.na(correlation.final.pd.fixed.raw$r2.CCM) & ! is.na(correlation.final.pd.fixed.raw$aic.SM)] <- 0

correlation.final.pd.fixed.raw$r2.CM[correlation.final.pd.fixed.raw$r2.CM < 0 & ! is.na(correlation.final.pd.fixed.raw$aic.SM)] <- 0
correlation.final.pd.fixed.raw$r2.CCM[correlation.final.pd.fixed.raw$r2.CCM < 0 & ! is.na(correlation.final.pd.fixed.raw$aic.SM)] <- 0

correlation.final.pd.fixed.raw$r2.CM[correlation.final.pd.fixed.raw$r2.CM < 0 & ! is.na(correlation.final.pd.fixed.raw$aic.SM)] <- 0
correlation.final.pd.fixed.raw$r2.CCM[correlation.final.pd.fixed.raw$r2.CCM < 0 & ! is.na(correlation.final.pd.fixed.raw$aic.SM)] <- 0

correlation.final.pd.fixed.raw[,c("p.SM","p.CM","p.CCM","r2.SM","r2.CM","r2.CCM",
                                  "aic.SM","aic.CM","aic.CCM")] <- round(correlation.final.pd.fixed.raw[,c("p.SM","p.CM","p.CCM","r2.SM","r2.CM","r2.CCM",
                                                                                                           "aic.SM","aic.CM","aic.CCM")],2)
correlation.final.pd.fixed.raw$Marker[correlation.final.pd.fixed.raw$Marker == "microsatellites"] <- "msat"

correlation.final.pd.fixed.raw <- correlation.final.pd.fixed.raw[order(correlation.final.pd.fixed.raw$Species),]

Table_S4 <- data.frame(Species = correlation.final.pd.fixed.raw$Species,
                       Marker = correlation.final.pd.fixed.raw$Marker,
                       Diff.Index = correlation.final.pd.fixed.raw$Diff.index,
                       Spec = correlation.final.pd.fixed.raw$Specificities,
                       Pairs = correlation.final.pd.fixed.raw$nbr.Pairs,
                       Pairs.Co = correlation.final.pd.fixed.raw$nbr.Pairs.Connected,
                       AIC.SM = correlation.final.pd.fixed.raw$aic.SM,
                       AIC.CM = correlation.final.pd.fixed.raw$aic.CM,
                       AIC.CCM = correlation.final.pd.fixed.raw$aic.CCM,
                       pval.SM = correlation.final.pd.fixed.raw$p.SM,
                       pval.CM = correlation.final.pd.fixed.raw$p.CM,
                       pval.CCM = correlation.final.pd.fixed.raw$p.CCM,
                       rsquared.SM = correlation.final.pd.fixed.raw$r2.SM,
                       rsquared.CM = correlation.final.pd.fixed.raw$r2.CM,
                       rsquared.CCM = correlation.final.pd.fixed.raw$r2.CCM,
                       Reference = correlation.final.pd.fixed.raw$ref)

file_name <- paste0(saving_dir.paper.SI,'/Table_S4','.csv')
write.table(Table_S4, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)

# Table S5

source("SI/SI_Table_S5.R")

# Table S6

# get hexagon name
for (k in 1:length(connectivity.final$cell_from)) {
  connectivity.final$cell_from_hex[k] <- hexagons.sourcesink.shp$HEX[hexagons.sourcesink.shp$ID == connectivity.final$cell_from[k]]
  connectivity.final$cell_to_hex[k] <- hexagons.sourcesink.shp$HEX[hexagons.sourcesink.shp$ID == connectivity.final$cell_to[k]]
}

Table_S6 <- data.frame(Species = connectivity.final$Species,
                              Marker = connectivity.final$Marker,
                              Diff.Index = connectivity.final$Diff.index,
                              Spec = connectivity.final$Specificity,
                              Site.From = connectivity.final$cell_from_hex,
                              Site.From.Lon = connectivity.final$Lon_from,
                              Site.From.Lat = connectivity.final$Lat_from,
                              Site.To = connectivity.final$cell_to_hex,
                              Site.To.Lon = connectivity.final$Lon_to,
                              Site.To.Lat = connectivity.final$Lat_to,
                              Sample.From.Lon = connectivity.final$sample_Lon_from,
                              Sample.From.Lat = connectivity.final$sample_Lat_from,
                              Sample.To.Lon = connectivity.final$sample_Lon_to,
                              Sample.To.Lat = connectivity.final$sample_Lat_to,
                              Gen.Diff = connectivity.final$Differentiation,
                              Con.Prob = connectivity.final$Connectivity.mean,
                              Con.Dist = connectivity.final$Connectivity.distance,
                              Con.Prob.fixed = connectivity.final.pd.fixed$Connectivity.mean,
                              Con.Dist.fixed = connectivity.final.pd.fixed$Connectivity.distance,
                              HC.From = connectivity.final$HC_pop_A,
                              HC.To = connectivity.final$HC_pop_B,
                              Step.Stones= connectivity.final$nbr_step,
                              Step.Stones.Fixed= connectivity.final.pd.fixed$nbr_step)

Table_S6$Species[Table_S6$Species == "Halidrys_dioica"] <- "Stephanocystis_dioica"

file_name <- paste0(saving_dir.paper.SI,'/Table_S6','.csv')
write.table(Table_S6, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)


# Figure S7

source("SI/SI_Figure_S7.R")
