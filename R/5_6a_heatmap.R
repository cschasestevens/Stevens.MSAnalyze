#### Heatmap ####

# Perform hierarchical clustering on DF

## Split into p-values and fold change

hc.data <- list("P" = dplyr::select(dplyr::filter(list.data$Plots$Input_filt,
                                                  Label == "pv"),
                                    -c("Comp.orig","Label")
                                    ),
             "FC" = dplyr::select(dplyr::filter(list.data$Plots$Input_filt,
                                                Label == "fc"),
                                  -c("Comp.orig","Label")
                                  )
             )

hc.data.grp <- as.matrix(hc.data$FC[2:ncol(hc.data$FC)])

row.names(hc.data.grp) <- hc.data$FC$Comp.new

hc.data.met <- t(hc.data.grp)


## Compute distance and clustering

set.seed(1234)

fun.hclust <- function(mat) {
  
  d.dist <- dist(scale(mat),
                 method = "euclidean")
  
  d.hcls <- hclust(d.dist,
                   method = "ward.D2")
  
  d.dend <- as.dendrogram(d.hcls)
  
  d.dend2 <- dendro_data(d.dend,type = "rectangle")
  
  d.dend.seg <- d.dend2$segments
  
  d.dend.lab <- d.dend2$labels
  
  return(list(d.dist,
              d.hcls,
              d.dend,
              d.dend.seg,
              d.dend.lab)
         )
  
  }

### Treatment groups

hc.grp <- fun.hclust(hc.data.grp)

### Metabolites/lipids

hc.met <- fun.hclust(hc.data.met)


## View dendrogram

fun.dend.row <- function(hc,
                     num.k,
                     exp.by) {
  
  p.dend <- hc[[3]] %>%
    set("branches_k_color",
        k = num.k,
        value = col1a)
  
  p.dend2 <- ggplot(as.ggdend(p.dend),
                    labels = F) +
    thm.mult +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = unit(c(0,0,0,0),
                             "cm"),
          axis.ticks.length.x = unit(0,"cm"),
          axis.ticks.length.y = unit(0,"cm")
    ) +
    coord_flip() +
    scale_x_reverse(expand = exp.by) +
    scale_y_reverse(expand = c(0.0,0.0))
  
}

fun.dend.col <- function(hc,
                         num.k,
                         exp.by) {
  
  p.dend <- hc[[3]] %>%
    set("branches_k_color",
        k = num.k,
        value = col1a)
  
  p.dend2 <- ggplot(as.ggdend(p.dend),
                    labels = F) +
    thm.mult +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = unit(c(0,0,0,0),
                             "cm"),
          axis.ticks.length.x = unit(0,"cm"),
          axis.ticks.length.y = unit(0,"cm")
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = exp.by)
  
}

### Save for incorporating into final heatmap

fun.dend.sav <- function(title1,p1,d.w,d.h) {
  
  ggsave(paste(save.file.in,
               title1,".png",
               sep = ""),
         p1,
         width = d.w,
         height = d.h,
         dpi = 700)
  }


## validation of clustering/cluster number

### Function (scale nboot accordingly to limit runtime {start at n=50})

fun.hclust.val <- function(hc,b.num) {
  
  hc.pval <- pvclust(hc,
                     method.dist = "cor",
                     method.hclust = "average",
                     nboot = b.num,
                     parallel = T)
  
  return(hc.pval)
  
}


## Final Heatmap

### order rows and columns by hclust results
### determine optimal cluster number by pvclust


## Input heatmap dataframe

fun.heat.in <- function(k.num.grp,k.num.met) {
  
  p.heat.in <- as.data.frame(hc.data.grp[hc.grp[[5]]$label,])
  
  p.heat.var <- data.frame(variable = names(p.heat.in),
                           new = list.data$`Data Import`$Annotations[[name.metab]],
                           Cluster2 = cutree(hc.met[[3]],k = k.num.met,
                                             order_clusters_as_data = T))
  
  
  hc.cut <- cutree(hc.grp[[3]],k = k.num.grp,
                   order_clusters_as_data = F)
  
  p.heat.in <- data.frame(Group = factor(row.names(p.heat.in),
                                         levels = unique(row.names(p.heat.in))),
                          Cluster = factor(hc.cut,
                                           levels = unique(hc.cut)),
                          p.heat.in)
  
  p.heat.in <- p.heat.in %>%
    dplyr::select(Group,Cluster,hc.met[[5]]$label)
  
  
  p.heat.in <- reshape2::melt(p.heat.in,
                              id.vars = c("Group","Cluster"))
  
  p.heat.in <- left_join(p.heat.in,
                         p.heat.var,
                         by = "variable")
  
  p.heat.in$new <- factor(p.heat.in$new,
                          levels = unique(p.heat.in$new))
  
  p.heat.in$Cluster2 <- factor(p.heat.in$Cluster2,
                               levels = unique(p.heat.in$Cluster2))
  
  return(p.heat.in)
  
}

## Add metabolite annotation clusters





fun.heat.out <- function(df,anno1,anno2) {
  
  p.dot2 <- ggplot(df,
                   aes(x = .data[[anno1]],
                       y = Group,
                       fill = value)) +
    geom_tile() +
    thm.mult +
    labs(fill = "Log2 Fold Change") +
    theme(plot.margin = unit(c(0,0,
                               0,0),
                             "cm"),
          # Axes
          axis.ticks.y = element_blank(),
          # axis.text.x = element_text(face = 'bold',
          #                            size = 8,
          #                            angle = 45,
          #                            hjust = 1,
          #                            vjust = 1),
          axis.text.x = element_text(face = "italic",
                                     size = 10),
          axis.text.y = element_text(face = "bold",
                                     size = 6),
          axis.title.x = element_text(face = "bold",
                                      size = 14),
          axis.title.y = element_text(face='bold',
                                      size = 14),
          
          # Strip
          strip.background = element_rect(fill = 'slategray2'),
          strip.text = element_text(face = 'bold',
                                    size = 8)) +
    
    theme(legend.text = element_text(size = 12))+
    scale_fill_gradientn(colors = col2a) +
    
    facet_grid2(
      # select data (text before ~ facets y, text after ~ facets x)
      Cluster ~ Cluster2,
                scales = "free",
                space = "free",
      switch = "y",
                strip = strip_themed(background_y = elem_list_rect(fill = col1a[1:length(unique(df[["Cluster"]])
                                                                                         )]
                                                                   ),
                                     background_x = elem_list_rect(fill = col1a[1:length(unique(df[["Cluster2"]]))])
                                     )) +

    theme(panel.spacing.y = unit(0,"cm"),
          panel.spacing.x = unit(0,"cm"),
          axis.ticks.length = unit(0,"cm"),
          axis.ticks = element_blank()) +

    scale_y_discrete(expand = c(0.0,0)) +
    xlab(anno2)
  
  return(p.dot2)
  
}


## Annotation Legends

fun.heat.legs <-function(df,var1,var2) {
  
  anno.row <- as_ggplot(get_legend(ggplot(df,
                     aes(x = Group,
                         y = .data[[var1]],
                         fill = .data[[var1]])) +
    geom_col(show.legend = T,
             width = 1) +
    thm.univ +
    labs(x = "",
         y = var1,
         fill = var1) +
    theme(axis.title.y = element_text(angle = 360,
                                      vjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,0,0,0),
                             "cm")) +
    scale_fill_manual(values = col1a) +
    thm.leg.title.y))
  
  
  anno.col <- as_ggplot(get_legend(ggplot(df,
                                          aes(x = Group,
                                              y = .data[[var2]],
                                              fill = .data[[var2]])) +
                                     geom_col(show.legend = T,
                                              width = 1) +
                                     thm.univ +
                                     labs(x = "",
                                          y = var2,
                                          fill = var2) +
                                     theme(axis.title.y = element_text(angle = 360,
                                                                       vjust = 0.5),
                                           axis.text.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           plot.margin = unit(c(0,0,0,0),
                                                              "cm")) +
                                     scale_fill_manual(values = col1a) +
                                     thm.leg.title.y))
  
  return(list("anno.row" = anno.row,
              "anno.col" = anno.col)
         )
  
}










