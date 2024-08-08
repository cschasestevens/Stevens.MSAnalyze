#### PLS-DA with CV plot ####

fun.pls <- function (df,md,grp) {
  
  # Prepare data
  
  pls.in <- df
  
  g1 <- pls.in[[grp]]
  
  data.pls <- pls.in[,-c(1:md)]
  
  pls_in <- plsda(data.pls, 
                  g1, 
                  ncomp = 15)
  
  set.seed(234)
  
  
  # cross-validation
  
  pls_cv <- perf(pls_in, 
                 folds = 15, 
                 progressBar = T,
                 auc = T, 
                 nrepeat = 15)
  
  
  # error rate
  
  cv1 <- pls_cv$error.rate$overall %>%
    as.data.frame() %>% 
    mutate(val = 'Overall')
  
  cv1$comp <- row.names(cv1)
  
  cv2 <- as.data.frame(pls_cv$error.rate$BER) %>%
    as.data.frame() %>% 
    mutate(val = 'BER')
  
  cv2$comp <- row.names(cv2)
  
  cv_r1 <- rbind(cv1,cv2)
  
  cv_r1$comp <- sub('.*comp', '', cv_r1$comp) %>% 
    as.numeric()
  
  cv_p1 <- cv_r1 %>% 
    pivot_longer(max.dist:mahalanobis.dist,
                 names_to = 'out',
                 values_to = 'Error Rate')
  
  
  # standard deviation
  
  cv3 <- as.data.frame(pls_cv$error.rate.sd$overall) %>%
    as.data.frame() %>% 
    mutate(val = 'Overall')
  
  colnames(cv3) <- paste(colnames(cv3), 
                         'sd', 
                         sep = '.')
  
  cv3$comp <- row.names(cv3)
  
  
  cv4 <- as.data.frame(pls_cv$error.rate.sd$BER) %>%
    as.data.frame() %>% 
    mutate(val = 'BER')
  
  colnames(cv4) <- paste(colnames(cv4), 
                         'sd', 
                         sep = '.')
  
  cv4$comp <- row.names(cv4)
  
  cv_r2 <- rbind(cv3,cv4)
  
  cv_r2$comp <- sub('.*comp', '', cv_r1$comp) %>% 
    as.numeric()
  
  cv_p3 <- cv_r2 %>% 
    pivot_longer(max.dist.sd:mahalanobis.dist.sd,
                 names_to = 'out',
                 values_to = 'Error Rate')
  
  
  
  # create cv plot (scatter)
  
  cv_p2 <- ggplot(cv_p1, aes(x = comp, 
                             y = `Error Rate`, 
                             color = out)) + 
    geom_line(aes(linetype = val)) +
    geom_point() +
    labs(x = 'Component') +
    thm.mult +
    scale_x_continuous(breaks = c(1,2,3,
                                  4,5,6,
                                  7,8,9,
                                  10)) +
    geom_errorbar(aes(ymin = `Error Rate` - cv_p3$`Error Rate`,
                      ymax = `Error Rate` - cv_p3$`Error Rate`),
                  width = .2,
                  position = position_dodge(0.05)) +
    
    scale_color_manual(values = col1a)
  
  
  # extract vars
  
  axes <- pls_in$variates[['X']]
  g <- g1
  
  pl_in <-  as.data.frame(cbind(axes, g))
  pl_in$g <- as.factor(pl_in$g)
  pl_in$comp1 <- as.numeric(pl_in$comp1)
  pl_in$comp2 <- as.numeric(pl_in$comp2)
  
  
  # format labels for axes
  ax <- pls_in$prop_expl_var[['X']]
  
  
  
  pl1_lab <- function(val, digits = 2, format = 'f') {
    paste0('X-Variate 1 (', formatC(100*val, 
                                    format = format, 
                                    digits = digits),
           '%)')
  }
  
  pl2_lab <- function(val, digits = 2, format = 'f') {
    paste0('X-Variate 2 (', formatC(100*val, 
                                    format = format, 
                                    digits = digits),
           '%)')
  }
  
  
  # Create PLS-DA plot
  
  plda <- ggplot(pl_in, aes(x=comp1, 
                            y=comp2, 
                            fill = g)) +
    
    stat_ellipse(geom = "polygon", 
                 col= "black", 
                 alpha =0.25, 
                 show.legend = F) +
    
    geom_point(shape=21, 
               col="black", 
               size = 3) +
    
    labs(y = pl2_lab(ax[2]), 
         x = pl1_lab(ax[1]),
         fill = grp) +
    
    thm.mult +
    
    scale_fill_manual(values = col1a)
  
  
  return(list("Input - PLSDA" = pls_in,
              "Plot - PLSCV" = cv_p2,
              "Plot - PLSDA" = plda))
  
}


# Create VIP Plot

fun.vip <- function(df,com1,com.sel) {
  
  # Plot input
  
  p.vip <- dplyr::arrange(dplyr::mutate(dplyr::slice_max(df,
                                                        .data[[com.sel]],
                                                        n = 25),
                                       Name = row.names(dplyr::slice_max(df,
                                                                         .data[[com.sel]],
                                                                         n = 25))),
                         desc(.data[[com.sel]]))
  
  p.vip$Name <- factor(p.vip$Name,
                      levels = rev(p.vip$Name))
  
  # Create plot
  
  p.vip2 <- ggplot(p.vip,
                  aes(x = .data[[com.sel]],
                      y = Name,
                      fill = .data[[com.sel]])) +
    geom_col() +
    
    scale_fill_gradientn(colors = col1b) +
    
    labs(x = com1,
         fill = "VIP Score")
  
  return(p.vip2)
  
}

















