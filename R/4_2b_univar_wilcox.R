#### Wilcox Rank-sum Test ####

# Function

fun.wilcox <- function(df,g1,g2) {
  
  d.wil <- dplyr::filter(df,
                        .data[[name.grp]] == g1 |
                          .data[[name.grp]] == g2)
  
  d.wil2 <- lapply(d.wil[(data.md + 1):ncol(d.wil)],
                  function(x) wilcox.test(x ~ d.wil[[name.grp]]))
  
  d.wil3 <- as.data.frame(lapply(d.wil2,
                                function(x) x[["p.value"]]))
  
  d.wil3 <- data.frame(Name = names(d.wil3),
                      `p-value` = as.numeric(t(d.wil3)))
  
  return(d.wil3)
  
  }






