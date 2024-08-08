#### Data transformation ####


# Data Distribution Function

fun.dst <- function(d1,md1,grp1,ext1) {
  
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
    
    d.log2[[grp1]] <- factor(d.log2[[grp1]],
                                 levels = unique(d.log2[[grp1]]))
    
    return(d.log2)
    
  }
  
  d1[["Log2 Datasheet"]] <- lg2.trf(d1[["Data Import"]][["Datasheet"]],
                                    md1,
                                    grp1)
  
  
  # Pareto scaling (for PCA/PLS-DA)
  
  d1[["Pareto-scaled Datasheet"]] <- cbind(d1[["Log2 Datasheet"]][1:md1],
                                           pareto_scale(d1[["Log2 Datasheet"]][(md1 + 1):
                                                                                ncol(d1[["Log2 Datasheet"]])]))
  
  
  p.dist <- ggarrange(dst(d1[["Data Import"]][["Datasheet"]],
                          "Normal",
                          md1),
                      dst(d1[["Log2 Datasheet"]],
                          "Log2 Transformed",
                          md1),
                      dst(d1[["Pareto-scaled Datasheet"]],
                          "Pareto Scaled",
                          md1),
                      nrow = 1,
                      ncol = 3)
  
  d1[["Distribution Plot"]] <- p.dist
  
  
  # Visualize data distribution
  
  ggsave(paste(ext1,"data_dist.png",sep = ""),
         p.dist,
         width = 15,
         height = 5,
         dpi = 600)
  
  
  return(d1)
  
  
  
}




