#### Correlation Analysis ####

# names for each list element

names(list.data$`Multivariate - Correlation`$`Input - Filter`) <- cor.comp


# Compute correlation matrices in parallel

## CLUSTER SET-UP: creates a cluster for parallelization of script execution using 75% of available cores
##                -all qc functions are nested in cluster evaluation function to load in each instance
##                -include specific libraries that are needed to run functions in clusterEvalQ()

## Runs in parallel using 75% of available cores or number of cores for each list element

clus1 <- makeCluster(ifelse(core.num > length(cor.comp),
                            length(cor.comp),
                            core.num*0.75))

### Plot Save function

fun.cor.plot.output <- function(ext.name,c1,p1,d.w,d.h,res1) {
  
  ggsave(paste(save.file.in,"Cor/",
               ext.name,".png",sep = ""),
         ggcorrplot::ggcorrplot(c1,
                                method = "circle",
                                type = "lower",
                                show.diag = T,
                                hc.order = T,
                                colors = c(col2a[[1]],col2a[[6]],col2a[[12]]),
                                outline.color = "white",
                                p.mat = p1,
                                insig = "blank"),
         width = d.w,
         height = d.h,
         dpi = res1)
  
}

### Pass libraries, functions, and variables to cluster

clusterEvalQ(clus1,{
  library(corrplot)
  library(ggplot2)
  library(ggcorrplot)
  
  # Insert functions and variables as needed
  
  fun.cor.plot.output <- function(ext.name,c1,p1,d.w,d.h,res1) {
    
    ggsave(paste(save.file.in,"Cor/",
                 ext.name,".png",sep = ""),
           ggcorrplot::ggcorrplot(c1,
                                  method = "square",
                                  show.diag = T,
                                  hc.order = T,
                                  colors = c(col2a[[1]],col2a[[6]],col2a[[12]]),
                                  outline.color = "white",
                                  p.mat = p1,
                                  insig = "blank",
                                  ggtheme = (thm.univ +
                                    thm.leg.title.y +
                                    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm")))),
           width = d.w,
           height = d.h,
           dpi = res1)
    
  }
  
})

clusterExport(clus1,varlist = c("cor.comp","save.file.in",
                                "list.data","col2a","thm.univ",
                                "thm.leg.title.y"))


## Generate each matrix

list.data$`Multivariate - Correlation`[["Matrix All"]] <- list("Correlation" = parLapply(clus1,
                                                                                         list.data$`Multivariate - Correlation`$`Input - Filter`,
                                                                                         function(x) round(cor(x,
                                                                                                               method = "spearman"),digits = 2)
),
"P-Value" = parLapply(clus1,
                      list.data$`Multivariate - Correlation`$`Input - Filter`,
                      function(x) cor.mtest(x,
                                            conf.level = 0.95)$p)
)


# Save plots

## Replace NA with 0

list.data$`Multivariate - Correlation`$`Matrix All`$Correlation <- lapply(list.data$`Multivariate - Correlation`$`Matrix All`$Correlation,
                                                                          function(x) replace(x,is.na(x),0)
)

list.data$`Multivariate - Correlation`$`Matrix All`$`P-Value` <- lapply(list.data$`Multivariate - Correlation`$`Matrix All`$`P-Value`,
                                                                        function(x) replace(x,is.na(x),0)
)




