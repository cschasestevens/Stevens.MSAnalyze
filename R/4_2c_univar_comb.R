#### Combine Fold Change & Univariate Statistics ####

list.data.p <- list.data$`Univariate - OW-ANOVA`$`Raw P`

list.data.fc <- data.frame(Name = names(list.data$`Descriptive Statistics`$`Log2 Fold Change`[2:ncol(list.data$`Descriptive Statistics`$`Log2 Fold Change`)]),
                           lapply(fun.t(list.data$`Descriptive Statistics`$`Log2 Fold Change`),
                                  function(x) round(as.numeric(x),
                                                    digits = 2)))

names(list.data.p) <- c("Name",
                        paste(names(list.data.p[2:ncol(list.data.p)]),
                              ".pv",sep = ""))

names(list.data.fc) <- c("Name",
                         paste(names(list.data.fc[2:ncol(list.data.fc)]),
                               ".fc",sep = ""))



list.data[["Plots"]] <- list("Input" = dplyr::left_join(list.data.p,
                                                        list.data.fc,
                                                        by = "Name"))

## Transpose Comparison names to match between p-value and fold change

fun.t2 <- function(df,lab1,lab2) {
  
  pfc <- as.data.frame(names(df))
  
  pfc2 <- as.data.frame(gsub(lab1,"",pfc[[1]]))
  
  names(pfc2) <- c("Comp")
  
  pfc2 <- dplyr::mutate(pfc2,
                        g1 = gsub("\\..*",
                                  "",
                                  pfc2[[1]]),
                        g2 = gsub("^.*\\.",
                                  "",
                                  pfc2[[1]]))
  
  pfc2 <- fun.t(pfc2)
  
  pfc3 <- as.data.frame(lapply(pfc2,
                               function(x) sort(x)))
  
  pfc3 <- as.data.frame(t(pfc3))
  
  pfc3 <- data.frame(Comp.orig = as.data.frame(names(df)),
                     Comp.new = paste(pfc3$V1,
                                      pfc3$V2,
                                      sep = "."),
                     Label = lab2)
  
  names(pfc3) <- c("Comp.orig","Comp.new","Label")
  
  return(pfc3)
  
}

list.data.comps <- base::rbind(fun.t2(list.data.p[-1],
                                      ".pv",
                                      "pv"),
                               fun.t2(list.data.fc[-1],
                                      ".fc",
                                      "fc"))

list.data.plots <- fun.t(list.data$Plots$Input)

list.data.plots <- data.frame(Comp.orig = row.names(list.data.plots),
                              as.data.frame(lapply(list.data.plots,
                                                   function(x) as.numeric(x)
                              )
                              )
)






