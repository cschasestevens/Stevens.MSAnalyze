#### Annotation Summary ####

# Count samples per group (if sample types are present in df)

fun.count.samp <- function(df) {
  
  grps <- df %>% 
    dplyr::select(sampleType) %>%
    dplyr::count(sampleType)
  
  grps$n <- as.character(grps$n)
  
  grps$Samp <- paste(grps$n,
                          grps$sampleType,
                          sep = " ")
  
  return(grps$Samp)
  
  }





# Add Class information to annotation list

fun.anno.merge <- function(anno.lab) {
  
  list.data <- c(list.data,
                 list("Annotation Summary" = 
                        list("Annotation List" = 
                               as.data.frame(list.data[["Data Import"]][["Annotations"]][[anno.lab]])
                             )
                      )
                 )
  
  
  names(list.data[[2]][[1]]) <- c("Name")
  
  list.data[[2]][["Merged"]] <- left_join(list.data[[2]][[1]],
                                          list.ref,
                                          by = "Name")
  
  list.data[[2]][["Missing"]] <- list.data[[2]][["Merged"]][is.na(list.data[[2]][["Merged"]]$`Major Class`),]
  
  nrow(list.data[[2]][["Missing"]])
  
  list.data[[2]][["Filter"]] <- list.data[[2]][["Merged"]][!is.na(list.data[[2]][["Merged"]]$`Major Class`),]
  
  list.data[[2]][["Summary"]] <- list.data[[2]][["Filter"]] %>%
    dplyr::select(`Major Class`,Class,Subclass,`Subclass 2`)
  
  return(list.data)
  
}

list.data <- fun.anno.merge(anno.lab1)





# Create Network

## Nodes

### Input df

anno.node.in <- list.data[[2]][["Summary"]]

anno.node.in <- anno.node.in %>%
  dplyr::arrange(`Major Class`,Class,Subclass,`Subclass 2`)

anno.node.in <- as.data.frame(lapply(anno.node.in,
                function(x) factor(x,levels = unique(x))))


### Annotation count and node-level assignment

fun.anno.cnt <- function(df) {
  
  fun.anno.cnt2 <- function(df1,var1,lev1) {
    
    cnt1 <- df1 %>%
      dplyr::count(.data[[var1]])
    
    cnt1[["Level"]] <- lev1
    
    names(cnt1) <- c("Node","N","Level")
    
    return(cnt1)
    
  }
  
  rbind(fun.anno.cnt2(df,"Major.Class",1),
        fun.anno.cnt2(df,"Class",2),
        fun.anno.cnt2(df,"Subclass",3),
        fun.anno.cnt2(df,"Subclass.2",4)
        )
  
}


### Annotation class proportions

fun.anno.prp <- function(df,l1,gvar,ovar) {
  
  cnt2 <- df %>%
    dplyr::count(.[,l1]) %>%
    dplyr::group_by(.data[[gvar]]) %>%
    dplyr::mutate(Prop = round((n/sum(n))*100,
                               digits = 2)) %>%
    dplyr::ungroup(.data[[gvar]]) %>%
    dplyr::select(all_of(c("Major.Class",ovar,"Prop")))
  
  names(cnt2) <- c("Major.Class","Node","Prop")
  
  return(cnt2)
  
}


### Run functions and append to node input

cnt.maj <- anno.node.in %>%
  dplyr::count(`Major.Class`) %>%
  dplyr::mutate(Prop = round((n/sum(n)*100),
                             digits = 2),
                Node = `Major.Class`) %>%
  dplyr::select(`Major.Class`,Node,Prop)

  
anno.node.weight <- rbind(cnt.maj,
                          fun.anno.prp(anno.node.in,c("Major.Class","Class"),
                                       "Major.Class","Class"),
                          fun.anno.prp(anno.node.in,c("Major.Class","Class","Subclass"),
                                       "Class","Subclass"),
                          fun.anno.prp(anno.node.in,c("Major.Class","Class","Subclass","Subclass.2"),
                                       "Subclass","Subclass.2"))


anno.node <- fun.anno.cnt(anno.node.in)


### Assign node sizes and label positions for plot generation

anno.node <- anno.node %>%
  dplyr::group_by(Level) %>%
  dplyr::mutate(End = 2*pi*cumsum(N)/sum(N),
                Start = dplyr::lag(End,
                                   default = 0),
                Mid = 0.5*(Start + End))

anno.node <- left_join(anno.node,
                       anno.node.weight,
                       by = "Node")

anno.node[["Label"]] <- ifelse(anno.node$Level == 4,
                               paste(gsub(".*-","",anno.node$Node),", ",
                                     anno.node$N,", ",anno.node$Prop,"%",sep = ""),
                               ifelse(anno.node$Level == 3,
                                      paste(gsub("d.*","d",anno.node$Node),", ",
                                            anno.node$N,", ",anno.node$Prop,"%",sep = ""),
                                      paste(anno.node$Node,", ",
                                            anno.node$N,", ",anno.node$Prop,"%",sep = "")))


anno.node[["Alt"]] <- ifelse(anno.node$Level == 4,
                             paste(gsub(".*-","",anno.node$Node),sep = ""),
                             ifelse(anno.node$Level == 3,
                                    paste(gsub("d.*","d",anno.node$Node),sep = ""),
                                    paste(anno.node$Node,sep = "")))

anno.node[["Label"]] <- ifelse((anno.node$Prop < (1/nrow(unique(anno.node[anno.node$Level == 4,"Alt"])))*100 &
                                 anno.node$Level == 4) |
                                 (anno.node$Prop < (1/nrow(unique(anno.node[anno.node$Level == 3,"Alt"])))*100 &
                                    anno.node$Level == 3) |
                                 (anno.node$Prop < (1/nrow(unique(anno.node[anno.node$Level == 2,"Alt"])))*100 &
                                    anno.node$Level == 2) |
                                 (anno.node$Prop < (1/nrow(unique(anno.node[anno.node$Level == 1,"Alt"])))*100 &
                                    anno.node$Level == 1),
                               "",
                               ifelse(anno.node$N < 10,
                                      "",
                                      anno.node$Label))





## Edges

### Generate edges based on node levels

fun.anno.edge <- function(df,
                        var1,
                        var2) {
  
  edge1 <- data.frame(df[var1],
                      df[var2])
  
  names(edge1) <- c("from","to")
  
  return(edge1)
  
}

anno.edge.in <- vector("list",
                length = ncol(anno.node.in))

for(i in 1:(length(anno.edge.in) - 1)) {
  
  anno.edge.in[[i]] <- fun.anno.edge(anno.node.in,
                                  i,
                                  i + 1)
  
}

### Bind list to obtain final edge list and remove duplicate edges

anno.edge.in <- dplyr::bind_rows(anno.edge.in)

anno.edge.in <- anno.edge.in[!duplicated(anno.edge.in),]


## Create Network

p.anno <- graph_from_data_frame(anno.edge.in,
                                  vertices = anno.node)





# Plot network as sunburst diagram

p.anno.sun <- function(net1,var1,var.lab) {
  
  # Create graph
  
  GO.graph.out <- ggraph(net1,
                         layout = "partition",
                         circular = T) +
    
    # Sunburst specific parameters
    
    geom_node_arc_bar(aes(fill = .data[[var1]],
                          r0 = ((Level - 1)/4)^1.25,
                          r = ((Level)/4)^1.25,
                          start = Start,
                          end = End,
                          alpha = Prop
    ),
    show.legend = T
    ) +
    
    # Color scale for distinguishing major classes
    
    scale_fill_manual(values = col1a) +
    
    # Plot labels
    
    geom_shadowtext(aes(x = (((Level)/4)^1.25+((Level - 1)/4)^1.25)*
                    sin(Mid)/2, 
                  y = (((Level)/4)^1.25+((Level - 1)/4)^1.25)*
                    cos(Mid)/2, 
                  label = Label),
                  bg.color = "black",
                  color = "white",
                  alpha = 1,
                  bg.r = 0.2) +
    
    # Sunburst theme parameters
    
    theme(panel.background = element_rect(fill = "white"),
          legend.key.size = unit(1,
                                 "cm"),
          legend.title = element_text(size = 18,
                                      face = "bold",
                                      vjust = 0.2),
          plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm"),
          legend.position = c("right"),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 20,
                                    face = "bold",
                                    hjust = 0.5)
          
    ) +
    
    # Plot labels
    
    labs(fill = var.lab)
  
  # Save plot in chosen directory
  
  ggsave(paste(save.file.in,"anno_sum.png",sep = ""),
         GO.graph.out,
         width = 14,
         height = 12,
         dpi = 600)
  
  }





# Add node information and network to list data

list.data[["Annotation Summary"]][["Nodes"]] <- anno.node

list.data[["Annotation Summary"]][["Network"]] <- p.anno


# Remove redundant objects from global environment

remove(anno.edge.in,anno.edge.weight,list.edge,
       p.dist,anno.node.in,anno.node.weight,cnt.maj)










