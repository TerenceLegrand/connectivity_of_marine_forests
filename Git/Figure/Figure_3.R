color <- c("#ffd700","#fa8775","#cd34b5")

# Forest plot

correlation.final$count <- seq(1,length(correlation.final$Species),1)
correlation.final$Species_i <- paste(correlation.final$Species_accepted,correlation.final$count)

# Genus order

correlation.final$Species_i <- factor(correlation.final$Species_i, levels = correlation.final$Species_i[order(correlation.final$genus,correlation.final$r2.CCM, decreasing = TRUE)])
#results$species <- factor(results$species, levels = results$species[order(results$r2.cent, decreasing = FALSE)])
correlation.final$species_plot <- paste0(substr(correlation.final$genus,1,1),". ",sub(".* ","", correlation.final$Species_accepted))
correlation.final$species_plot <- correlation.final$species_plot[order(correlation.final$genus,correlation.final$r2.CCM, decreasing = TRUE)]

correlation.final$signif_SM_plot <- "*"
correlation.final$signif_SM_plot[!correlation.final$p.SM <= 0.05] <- ""

correlation.final$signif_CM_plot <- "*"
correlation.final$signif_CM_plot[!correlation.final$p.CM <= 0.05] <- ""

correlation.final$signif_CCM_plot <- "*"
correlation.final$signif_CCM_plot[!correlation.final$p.CCM <= 0.05] <- ""

shadded <- correlation.final$genus[order(correlation.final$genus,correlation.final$r2.CCM, decreasing = TRUE)] 
shadded <- which(! duplicated(shadded))
shadded <- cbind(shadded,c(shadded[-1],length(correlation.final$genus)+1))
shadded <- shadded[c(TRUE, FALSE),]
shadded <- shadded -0.5

forest_plot <- correlation.final %>%
  ggplot(aes(y = factor(Species_i))) + 
  geom_point(aes(r2.SM,Species_i), colour = color[1]) + 
  geom_point(aes(r2.CM,Species_i), colour = color[2]) + 
  geom_point(aes(r2.CCM,Species_i), colour = color[3]) +
  annotate("rect", xmin = -.1, xmax = 1.1, ymin = shadded[,1], ymax = shadded[,2], alpha = .2) +
  geom_hline(yintercept = seq(1,length(correlation.final$Species_i),1), linetype="dotted", linewidth = 0.1) +
  geom_linerange(aes(xmin=r2.SM, xmax=r2.CM), colour = "#666666", linewidth = 1) +
  geom_linerange(aes(xmin=r2.SM, xmax=r2.CCM), colour = "#666666", linewidth = 1) +
  geom_point(aes(r2.SM,Species_i), colour = color[1]) + 
  geom_point(aes(r2.CM,Species_i), colour = color[2]) + 
  geom_point(aes(r2.CCM,Species_i), colour = color[3]) +
  theme_classic() +
  scale_x_continuous("R²",expand = c(0,0)) +
  scale_y_discrete(labels=correlation.final$species_plot) +
  theme(axis.title.y = element_blank(),axis.text.y = element_text(size=6, face = "italic"))


## ------------
# Boxplot

plot_results <- rbind(data.frame(model = "Spatial Model",
                                 species = correlation.final$Species,
                                 studies = correlation.final$Studies,
                                 r2 = correlation.final$r2.SM),
                      data.frame(model = "Connectivity Model",
                                 species = correlation.final$Species,
                                 studies = correlation.final$Studies,
                                 r2 = correlation.final$r2.CM),
                      data.frame(model = "Connectivity and Correlation Model",
                                 species = correlation.final$Species,
                                 studies = correlation.final$Studies,
                                 r2 = correlation.final$r2.CCM))

plot_results$model <- factor(plot_results$model, levels = c("Connectivity and Correlation Model","Connectivity Model","Spatial Model"))

boxplot <- ggplot(plot_results) +
  geom_boxplot(aes(r2,model),colour = rev(color), fill = rev(color), alpha = 0.4) +
  geom_jitter(aes(r2,model),width = 0, size = 0.2, alpha = 0.1) +
  annotate("rect", xmin = -.1, xmax = 1.1, ymin = 1, ymax = 1, alpha = 0) +
  scale_x_continuous(expand = c(0,0), position = "top") +
  scale_y_discrete(labels=c("CCM","CM","SM")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_text(hjust=1, size = 8, angle=0, vjust = 0.5, face = "bold"))

## ------------
# Traits

correlation.final$Marker[correlation.final$Marker == "microsatellites"] <- "msat"
correlation.final$nbr.PairsRegression.connected <- round((correlation.final$nbr.Pairs.Connected/correlation.final$nbr.Pairs)*100,0)
correlation.final$sampling_mean_distance <- round(correlation.final$sampling_mean_distance/1000,0)

good_order <- matrix(NA,ncol = 1,nrow = length(correlation.final$Species))
order_forest <- order(correlation.final$genus,correlation.final$r2.CCM, decreasing = FALSE)
for (i in 1:length(correlation.final$Species)) {
  good_order[order_forest[i]] <- i
}

correlation.final$data_id <- as.character(good_order )

traits <- correlation.final %>%
  ggplot() +
  geom_text(aes(x = 0, y = Species_i, label = signif_SM_plot),  hjust = 0, size = 4,  colour = color[1]) +
  geom_text(aes(x = 0.01, y = Species_i, label = signif_CM_plot),  hjust = 0, size = 4,  colour = color[2]) +
  geom_text(aes(x = 0.02, y = Species_i, label = signif_CCM_plot),  hjust = 0, size = 4,  colour = color[3]) +
  geom_text(aes(x = 0.035, y = Species_i, label = data_id),  hjust = 0, size = 2) +
  geom_text(aes(x = 0.06, y = Species_i, label = Marker),  hjust = 0, size = 2) +
  geom_text(aes(x = 0.15, y = Species_i, label = Diff.index),  hjust = 0, size = 2) +
  geom_text(aes(x = 0.21, y = Species_i, label = nbr.pop),  hjust = 0, size = 2) +
  geom_text(aes(x = 0.26, y = Species_i, label = nbr.Pairs.Connected.LDD),  hjust = 0, size = 2) +
  geom_text(aes(x = 0.31, y = Species_i, label = sampling_mean_distance),  hjust = 0, size = 2) +
  #geom_text(aes(x = 0.36, y = Species_i, label = lat_extent),  hjust = 0, size = 2) +
  geom_text(aes(x = 0.36, y = Species_i, label = ref),  hjust = 0, size = 2) +
  #geom_text(aes(x = 1, y = 65, label = "ref"),  hjust = 0, size = 2, fontface = "bold") +
  annotate("rect", xmin = 0, xmax = 0.55, ymin = shadded[,1], ymax = shadded[,2], alpha = .2) +
  theme_void()

title <- 
  ggplot() +
  geom_text(aes(x = 0.03, y = 1, label = "ID"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  geom_text(aes(x = 0.06, y = 1, label = "Marker"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  geom_text(aes(x = 0.15, y = 1, label = "Index"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  geom_text(aes(x = 0.21, y = 1, label = "Nbr pop"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  geom_text(aes(x = 0.26, y = 1, label = "LDD connections (%)"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  geom_text(aes(x = 0.31, y = 1, label = "Mean sampling distance (km)"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  geom_text(aes(x = 0.36, y = 1, label = "Reference"),  hjust = 0, size = 2, fontface = "bold",angle = 50) +
  annotate("rect", xmin = 0, xmax = 0.55, ymin = 1, ymax = 1, alpha = 0) +
  theme_void()

## ------------
# Genus

genus.seg <- correlation.final$genus[order(correlation.final$genus,correlation.final$r2.CCM, decreasing = TRUE)] 
genus.seg <- which(! duplicated(genus.seg))
genus.seg <- cbind(genus.seg-0.2,c(genus.seg[-1]-0.8,length(correlation.final$genus)+0.2))

genus <- correlation.final %>%
  ggplot() +
  geom_text(aes(x = 1, y = Species_i, label = ''),  hjust = 1, size = 4) +
  annotate("text", x = 1, y = rowMeans(genus.seg),label = rev(sort(unique(correlation.final$genus))), alpha = 1, hjust = 1, size = 3) +
  annotate("segment", x = 1.1, xend = 1.1, y = genus.seg[,1], yend = genus.seg[,2], alpha = 1) +
  theme_void() +
  coord_cartesian(xlim = c(0, 2)) 

  
## ------------
# Combine

height = 30
width = 20

layout <- c(
  patchwork::area(t = 0, l = 5, b = 4.5, r = 12),
  patchwork::area(t = 0, l = 13, b = 8, r = 20),
  patchwork::area(t = 5, l = 0, b = 30, r = 5), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  patchwork::area(t = 5, l = 5, b = 30, r = 12),
  patchwork::area(t = 5, l = 13, b = 30, r = 20)
)


# final plot arrangement
figure_3 <- boxplot + title + genus +forest_plot + traits +  plot_layout(design = layout)

ggsave(filename = paste0("/Figure_3",".pdf"),
       plot = figure_3,
       device = cairo_pdf,
       path = saving_dir.paper.main,
       width = 21,
       height = 25.7,
       units = "cm"
)
