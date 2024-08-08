#### Output Function for ChemRICH ####

# ChemRICH plot

fun.cr.plot <- function(df,
                        xlab1,
                        p.title) {
  
  # Main Plot
  
  cr.plot <- ggplot(df,
                    aes(x = Cluster,
                        y = -log10(FDR),
                        fill = Ratio)) +
    ## Plot points
    
    geom_point(shape = 21,
               col = "black",
               aes(size = df[["Cluster Size"]]),
               alpha = 0.7) +
    ## Significance line
    
    geom_hline(yintercept = 1.3,
               linetype = "dashed") +
    ## Data scale (only if needed)
    
    # scale_y_continuous(limits = c(0,
    #                               10),
    #                    breaks = c(0,2,4,
    #                               6,8,10)) +
    
    ## Theme
    thm.mult +
    ## Axis labels
    xlab(xlab1) +
    ylab("-log10(p)") +
    ## Gradient and size scaling
    
    scale_fill_gradientn(name = "Increased Ratio",
                         colors = c('midnightblue','royalblue4','royalblue2',
                                    'lightskyblue',"azure","white",
                                    "navajowhite",'goldenrod1',"orange",
                                    "firebrick1",'firebrick4'),
                         breaks = c(c(0,.1,.2,.3,.4,
                                      0.5,.6,.7,.8,.9,
                                      1)),
                         labels = c("0","","0.2","","0.4",
                                    "","0.6","","0.8","","1.0")
    ) +
    scale_size_area(name = "Cluster Size",
                    max_size = 16)
  
  # Save
  
  ggsave(paste(save.chem.in,
               p.title,
               ".png",
               sep = ""),
         width = 12,
         height = 10,
         dpi = 700)
  
}


fun.cr.plot.output <- function(list.cr,var1) {
  
  lapply(cr.comp,
         function(x) {
           
          tryCatch({fun.cr.plot(list.data$`Multivariate - ChemRICH`[[list.cr]][[x]],
                                var1,
                                paste(var1,"/",
                                      gsub("\\.","v",x)))},
                   error = function(e) {message("Check Plot Input"); print(e)},
                   finally = {ggplot() +
                       theme_void()})
           
         })
  
}




