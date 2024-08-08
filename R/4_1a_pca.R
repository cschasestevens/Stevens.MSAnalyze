#### Principal Components Analysis (PCA) ####

# Prepare input for PCA (data frame, grouping variable, point overlay, qc1, qc2, mdata + 1)

fun.pca <- function(df, 
                md1,
                grp1) {
  
  # Separate metadata 
  
  pca.in <- prcomp(df[,(md1 + 1):ncol(df)], 
                   scale. = T) 
  
  vrc <- summary(pca.in)$importance[2,]
  
  PC1 <- pca.in$x[, 1]
  
  PC2 <- pca.in$x[, 2]
  
  p9 <- cbind(df, 
              PC1, 
              PC2)
  
  
  # Assign PC1 and 2 labels
  
  pc1_lab <- function(val, digits = 2, format = 'f') {
    paste0('PC1 (', formatC(100*val, 
                            format = format, 
                            digits = digits),
           '%)')
  }
  
  pc2_lab <- function(val, digits = 2, format = 'f') {
    paste0('PC2 (', formatC(100*val, 
                            format = format, 
                            digits = digits),
           '%)')
  }
  
  
  # Generate plot
  
  pcp <- ggplot(p9, 
                aes(x=PC1, 
                    y=PC2, 
                    fill = .data[[grp1]])) +
    
    geom_point(shape=21, 
               col="black", 
               size = 3) +
    
    labs(y = pc2_lab(vrc[2]), 
         x = pc1_lab(vrc[1])) +
    
    thm.mult +
    
    scale_fill_manual(values = col1a)
  
  return(list("Input" = pca.in,
              "Plot - PCA" = pcp))
  
}




