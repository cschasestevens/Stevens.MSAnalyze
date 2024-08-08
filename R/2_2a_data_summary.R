#### Data transformation ####


# Data Distribution Function

dst <- function(d1) {
  
  ## density plot (data frame, mdata + 1, plot.title)
  
  dst <- function(df,trf1,md) {
    dsty <- as.data.frame(colMeans(df[(md + 1):ncol(df)]))
    
    names(dsty) <- c("Value")
    
    p <- ggplot(dsty, 
                aes(x=`Value`)) +
      
      # Density Plot
      
      geom_density(color = 'darkslategrey', 
                   fill = col1b[[3]]) +
      
      # Plot Theme
      
      labs(title = paste(trf1,"Data Distribution"), 
           y= "Density") +
      
      thm.univ +
      
      theme(plot.margin = unit(c(rep(0.2,4)),"cm"))
    
    return(p)
    
  }
  
  
  # log2-transform
  
  lg2.trf <- function(df,md,grp){
    
    # Log2-transform dataset
    
    d.log2 <- cbind(df[,c(1:md)],
                    as.data.frame(lapply(df[,c((md + 1):ncol(df))], 
                                         function(x) log2(x))))
    
    d.log2[[name.grp]] <- factor(d.log2[[name.grp]],
                                 levels = unique(d.log2[[name.grp]]))
    
    return(d.log2)
    
  }
  
  d1[["Log2 Datasheet"]] <- lg2.trf(d1[["Datasheet"]],
                                    data.md,
                                    name.grp)
  
  
  # Pareto scaling (for PCA/PLS-DA)
  
  d1[["Pareto-scaled Datasheet"]] <- cbind(d1[["Log2 Datasheet"]][1:data.md],
                                           pareto_scale(d1[["Log2 Datasheet"]][(data.md + 1):
                                                                                ncol(d1[["Log2 Datasheet"]])]))
  
  
  p.dist <- ggarrange(dst(d1[["Datasheet"]],
                          "Normal",
                          data.md),
                      dst(d1[["Log2 Datasheet"]],
                          "Log2 Transformed",
                          data.md),
                      dst(d1[["Pareto-scaled Datasheet"]],
                          "Pareto Scaled",
                          data.md),
                      nrow = 1,
                      ncol = 3)
  
  
  # Visualize data distribution
  
  ggsave(paste(save.file.in,"data_dist.png",sep = ""),
         p.dist,
         width = 15,
         height = 5,
         dpi = 600)
  
  
  return(d1)
  
  
  
}




