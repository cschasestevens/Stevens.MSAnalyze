#### ChemRICH (adapted from Barupal et. al. 2017) ####

### Annotations

cr.anno <- list.data$`Data Import`$Annotations[,c("Index","Metabolite name")]

names(cr.anno) <- c("Index","Name")

cr.anno <- left_join(cr.anno,
                     list.data$`Annotation Summary`$Merged,
                     by = "Name",
                     multiple = "first")


### P-values and Fold Changes

cr.stat <- list.data$Plots$Input


### df list for generating all ChemRICH; Required columns are Index,Name,Subclass,Saturation,Chain Type

cr.comp <- unique(cr.stat$Comp.new)

fun.t.cr <- function(c,x) {
  
  cr1 <- dplyr::filter(cr.stat,
                       Comp.new == cr.comp[[x]] &
                         Label == c) %>%
    t() %>%
    as.data.frame()
  
  names(cr1) <- cr1[1,]
  
  return(cr1)
  
}

list.data[["Multivariate - ChemRICH"]] <- list("Input" = vector("list",
                                                                length = length(cr.comp)))

for(i in 1:length(list.data$`Multivariate - ChemRICH`$Input)) {
  
  list.data$`Multivariate - ChemRICH`$Input[[i]] <- cbind(fun.t.cr("pv",i),
                                                          fun.t.cr("fc",i))
  
  list.data$`Multivariate - ChemRICH`$Input[[i]] <- data.frame(cr.anno[,c("Index","Name",
                                                                          "Subclass","Saturation","Chain Type")],
                                                               as.data.frame(lapply(list.data$`Multivariate - ChemRICH`$Input[[i]][-c(1:3),],
                                                                                    function(x) as.numeric(x)
                                                               )
                                                               )
  )
  
  list.data$`Multivariate - ChemRICH`$Input[[i]] <- list.data$`Multivariate - ChemRICH`$Input[[i]][!is.na(list.data$`Multivariate - ChemRICH`$Input[[i]]$Subclass),]
  
}

names(list.data$`Multivariate - ChemRICH`$Input) <- cr.comp


### Conduct ks-tests for each by class, saturation, and chain length

fun.cr.input <- function(name.var) {
  
  cr.input <- lapply(cr.comp,
                     function(j) {
                       
                       
                       # Identify Groups for KS-test
                       
                       ks <- unique(list.data$`Multivariate - ChemRICH`$Input[[j]][[name.var]])
                       
                       ks.df <- list.data$`Multivariate - ChemRICH`$Input[[j]]
                       
                       
                       # Perform KS-test
                       
                       ks2 <- lapply(ks,
                                     function(x) {
                                       
                                       ks3 <- ks.test(dplyr::select(filter(ks.df,
                                                                           .data[[name.var]] == x),
                                                                    6),
                                                      "punif",
                                                      alternative = "greater")
                                       
                                       ks3 <- data.frame(`Cluster` = x,
                                                         `p-value` = ks3[["p.value"]])
                                       
                                     })
                       
                       
                       # Convert to ChemRICH Input
                       
                       cr.input.fun <- function(df,
                                                df2,
                                                v1,
                                                v2) {
                         
                         cr.plot <- data.frame(`Index` = seq(1:length(df)))
                         
                         ### KS-result
                         
                         cr.ks.result <- data.frame(`Cluster` = as.character(t(as.data.frame(lapply(df,
                                                                                                    "[",
                                                                                                    paste(v1))))),
                                                    `p-value` = t(as.data.frame(lapply(df,
                                                                                       "[",
                                                                                       "p.value"))))
                         
                         cr.ks.result[["FDR"]] <- p.adjust(cr.ks.result$p.value,
                                                           method = "fdr")
                         
                         ### Cluster Size
                         
                         cr.clust <- df2 %>%
                           dplyr::count(as.character(df2[[paste(v2)]]))
                         
                         names(cr.clust) <- c("Cluster",
                                              "Cluster Size")
                         
                         ### Inc/Dec number
                         
                         cr.rat <- data.frame(Cluster = as.character(df2[[paste(v2)]]),
                                              
                                              Ratio = ifelse(df2[[7]] > 0,
                                                             "Inc",
                                                             "Dec")
                         )
                         
                         cr.rat.inc <- dplyr::count(dplyr::filter(cr.rat,
                                                                  Ratio == "Inc"),
                                                    Cluster)
                         
                         cr.rat.dec <- dplyr::count(dplyr::filter(cr.rat,
                                                                  Ratio == "Dec"),
                                                    Cluster)
                         
                         names(cr.rat.inc) <- c("Cluster",
                                                "Increased")
                         
                         names(cr.rat.dec) <- c("Cluster",
                                                "Decreased")
                         
                         ## Combine for final plot input and ratio calculation
                         
                         cr.plot <- cbind(cr.plot,
                                          left_join(cr.clust,
                                                    left_join(cr.ks.result,
                                                              full_join(cr.rat.inc,
                                                                        cr.rat.dec,
                                                                        by = "Cluster"),
                                                              by = "Cluster"),
                                                    by = "Cluster"))
                         
                         cr.plot[is.na(cr.plot)] <- 0
                         
                         
                         cr.plot[["Ratio"]] <- ifelse(cr.plot$Decreased == 0,
                                                      1,
                                                      (cr.plot$Increased)/(cr.plot$Increased +
                                                                             cr.plot$Decreased)
                         ) 
                         
                         return(cr.plot)
                         
                       }
                       
                       
                       cr.in <- cr.input.fun(ks2,
                                             ks.df,
                                             "Cluster",
                                             name.var)
                       
                       return(cr.in)
                       
                     })
  
  names(cr.input) <- cr.comp
  
  return(cr.input)
  
  
}


