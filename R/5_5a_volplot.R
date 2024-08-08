#### Single Volcano Plot ####
  
  
  # Plot function
  
fun.vol <- function (df,
                       group_fc, 
                       group_p, 
                       title,
                       p.name) {
  
  
  v_1 <- EnhancedVolcano(df, 
                         lab = df[["Label"]],
                         title = element_blank(),
                         subtitle = element_blank(), 
                         caption = element_blank(), 
                         x= group_fc, 
                         y= group_p,
                         pCutoff = 0.005, 
                         FCcutoff = 0.25,
                         cutoffLineType = 'twodash',
                         legendLabels = c('NS','Fold Change',
                                          'p-value','FC+P'), 
                         legendLabSize = 12,
                         labFace = 'bold', 
                         col = col1b,
                         colAlpha = 0.7,
                         legendIconSize = 4,
                         pointSize = 2,
                         border = 'full',
                         borderWidth = 1.5,
                         legendPosition = 'right',
                         labSize = 3,
                         drawConnectors = T,
                         typeConnectors = "open",
                         min.segment.length = unit(1,
                                                   "mm")
  ) +
    
    thm.univ +
    labs(color='Key') +
    ggtitle(title) +
    
    theme(plot.margin = unit(c(.2,.2,.2,.2),"cm"))
  
  # v1 <- v_1 +
  #   ggplot2::coord_cartesian(xlim=c(-5,5)) + 
  #   # ylim = c(0,100)) +
  #   ggplot2::scale_x_continuous(breaks=seq(-5,5,1))
  
  # v1 <- move_layers(v1,
  #                   "GeomPolygon",
  #                   position = "bottom")
  
  ggsave(paste(save.file.in,"Vol/",
               p.name,".png",sep = ""),
         v_1,
         width = 8,
         height = 8,
         dpi = 700)
  
}
  
  
  




