# SI Sentivity tests

sensitivity.test <- correlation.final[correlation.final$p.CCM <= 0.05,]

sensitivity <- data.frame()

model.results <- c("r2.CCM","r2.SM","delta_r2")

for (k in 1:length(model.results)) {
  
  weight <- model.results[k]
  
  characteristics <- c("Marker","Diff.index","genus","family","order")
  
  for (i in 1:length(characteristics)) {
    
    group <- characteristics[i]
    
    test <- sensitivity.test[,c(group,weight)]
    colnames(test) <- c("group","weight")
    
    c_summary <- test %>% 
      group_by(group) %>% 
      tally()
    
    Kruskal.Wallis.results <- kruskal.test(weight ~ group, data = test)
    
    sensitivity <- rbind(sensitivity,data.frame(Model = weight,
                                                Feature = group,
                                                Test = "Kruskal-Wallis",
                                                n = length(test$group),
                                                pvalue = Kruskal.Wallis.results$p.value,
                                                adj.rsquared = NA))
  }
  
  # Realm
  
  realm.i <- unique(str_trim(tolower(unlist(strsplit(sensitivity.test$REALM,","))), side = c("both")))
  
  test <- data.frame()
  for (i in 1:length(realm.i)) {
    realm.results <- data.frame(group = rep(realm.i[i],length(grep(realm.i[i],tolower(sensitivity.test$REALM)))),
                                weight = sensitivity.test[grep(realm.i[i],tolower(sensitivity.test$REALM)),weight])
    test <- rbind(test,
                  realm.results)
    
  }
  
  
  c_summary <- test %>% 
    group_by(group) %>% 
    tally()
  
  Kruskal.Wallis.results <- kruskal.test(weight ~ group, data = test)
  
  
  sensitivity <- rbind(sensitivity,data.frame(Model = weight,
                                              Feature = "Realm",
                                              Test = "Kruskal-Wallis",
                                              n = length(test$group),
                                              pvalue = Kruskal.Wallis.results$p.value,
                                              adj.rsquared = NA))
  
  # Linear models 
  
  characteristics <- c("nbr.pop","sampling_mean_distance","lat_extent","sampling_max_extent")
  
  for (i in 1:length(characteristics)) {
    
    group <- characteristics[i]
    test <- sensitivity.test[,c(group,weight)]
    colnames(test) <- c("group","weight")
    
    
    model <- lm(weight ~ group, data = test)
    
    sensitivity <- rbind(sensitivity,data.frame(Model = weight,
                                                Feature = group,
                                                Test = "Linear model",
                                                n = length(test$group),
                                                pvalue = anova(model)$`Pr(>F)`[1],
                                                adj.rsquared = summary(model)$adj.r.squared))
  }
}


sensitivity[,c("pvalue","adj.rsquared")] <- round(sensitivity[,c("pvalue","adj.rsquared")],2)

file_name <- paste0(saving_dir.paper.SI,'/Table_S5','.csv')
write.table(sensitivity, file = file_name, sep = ",", eol = "\n", row.names = FALSE, col.names = TRUE)
