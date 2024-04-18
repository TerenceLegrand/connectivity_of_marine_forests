# Load data

# For LDD
# 10^-3
pld.period = 0
saving_dir.PD.LDD <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/correlation_final",".csv"), check.names=FALSE)
correlation.final.LDD$LDD_coef <- "10^-3"
connectivity.final.LDD$LDD_coef <- "10^-3"

correlation.final.LDD.all <- correlation.final.LDD
connectivity.final.LDD.all <- connectivity.final.LDD

# 10^-4
pld.period = 1
saving_dir.PD.LDD <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/correlation_final",".csv"), check.names=FALSE)
correlation.final.LDD$LDD_coef <- "10^-4"
connectivity.final.LDD$LDD_coef <- "10^-4"

correlation.final.LDD.all <- rbind(correlation.final.LDD.all,correlation.final.LDD)
connectivity.final.LDD.all <- rbind(connectivity.final.LDD.all,connectivity.final.LDD)

# 10^-5
pld.period = 2
saving_dir.PD.LDD <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/correlation_final",".csv"), check.names=FALSE)
correlation.final.LDD$LDD_coef <- "10^-5"
connectivity.final.LDD$LDD_coef <- "10^-5"

correlation.final.LDD.all <- rbind(correlation.final.LDD.all,correlation.final.LDD)
connectivity.final.LDD.all <- rbind(connectivity.final.LDD.all,connectivity.final.LDD)

# 10^-2
pld.period = 3
saving_dir.PD.LDD <- paste0(saving_dir,"/PD_",as.character(pld.period))

connectivity.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/connectivity_final",".csv"), check.names=FALSE)
correlation.final.LDD <- read.csv(paste0(saving_dir.PD.LDD,"/correlation_final",".csv"), check.names=FALSE)
correlation.final.LDD$LDD_coef <- "10^-2"
connectivity.final.LDD$LDD_coef <- "10^-2"

correlation.final.LDD.all <- rbind(correlation.final.LDD.all,correlation.final.LDD)
connectivity.final.LDD.all <- rbind(connectivity.final.LDD.all,connectivity.final.LDD)

# Plot

correlation.final.LDD.all$LDD_coef <- as.factor(correlation.final.LDD.all$LDD_coef)
connectivity.final.LDD.all$LDD_coef <- as.factor(connectivity.final.LDD.all$LDD_coef)



give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

boxplot_AIC <-ggplot(correlation.final.LDD.all, aes(LDD_coef,aic.CCM)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 0.1, outlier.alpha = 0.1) +
  geom_boxplot(fill = "#CDCDCD") +
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
  #              position = position_dodge(width = 0.75)) +
  #coord_flip(clip = "off") +
  #scale_y_log10(breaks = c(1,5,10,50,100,150)) +
  #scale_x_discrete(labels=c('Long Distance Dispersal','Fixed Dispersal Time')) +
  labs(x ="LDD weighting factor",y = "AIC") +
  theme_classic() 


boxplot_r <-ggplot(correlation.final.LDD.all, aes(LDD_coef,r2.CCM)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 0.1, outlier.alpha = 0.1) +
  geom_boxplot(fill = "#CDCDCD") +
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
  #              position = position_dodge(width = 0.75)) +
  #coord_flip(clip = "off") +
  #scale_y_log10(breaks = c(1,5,10,50,100,150)) +
  #scale_x_discrete(labels=c('Long Distance Dispersal','Fixed Dispersal Time')) +
  labs(x ="LDD weighting factor",y = "R²") +
  theme_classic()

boxplot_step <- ggplot(connectivity.final.LDD.all, aes(LDD_coef,nbr_step)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 0.1, outlier.alpha = 0.1) +
  geom_boxplot(fill = "#CDCDCD") +
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
  #              position = position_dodge(width = 0.75)) +
  #coord_flip(clip = "off") +
  scale_y_log10(breaks = c(1,5,10,50,100,150)) +
  #scale_x_discrete(labels=c('Long Distance Dispersal','Fixed Dispersal Time')) +
  labs(x ="LDD weighting factor",y = "Number of stepping stones") +
  theme_classic() 



layout <- "
    AB
    CD
    "

figure_LDD <- boxplot_AIC + boxplot_r + boxplot_step +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A')

ggsave(filename = paste0("/Figure_LDD_coeff",".pdf"),
       plot = figure_LDD,
       device = cairo_pdf,
       path = saving_dir.paper.SI
)


# Stat

connectivity.final.LDD.all$nbr_step[is.infinite(connectivity.final.LDD.all$nbr_step)] <- NA

correlation.final.LDD.all%>%group_by(LDD_coef)%>%summarise(Average=mean(aic.CCM,na.rm= TRUE ))
stats::kruskal.test(aic.CCM ~ LDD_coef, data = correlation.final.LDD.all)

correlation.final.LDD.all%>%group_by(LDD_coef)%>%summarise(Average=mean(r2.CCM,na.rm= TRUE ))
stats::kruskal.test(r2.CCM ~ LDD_coef, data = correlation.final.LDD.all)

connectivity.final.LDD.all%>%group_by(LDD_coef)%>%summarise(Average=mean(nbr_step,na.rm= TRUE ))
stats::kruskal.test(nbr_step ~ LDD_coef, data = connectivity.final.LDD.all)
