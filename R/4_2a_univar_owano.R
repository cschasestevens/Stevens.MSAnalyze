#### One-way ANOVA with Tukey's Post-Hoc analysis and Benjamini-Hochberg correction

## One-way ANOVA Function

### Function

fun.ow.ano <- function(df, 
                   md) {
  # Input data
  
  d.in <- df
  
  m.name <- d.in[,-c(1:md)] %>%
    names()
  
  
  # Run OW-ANOVA with Tukey's post-hoc analysis
  
  d.ow <- lapply(d.in[(md + 1):ncol(d.in)],
                 function(x) TukeyHSD(aov(x ~ d.in[[name.grp]]), 
                                      conf.level = 0.95))
  
  
  # extract post-hoc result
  
  d.ow.post <- as.data.frame(lapply(lapply(d.ow, 
                                           function(x) lapply(x, 
                                                              '[',
                                                              missing_arg(),
                                                              "p adj")),
                                    "["))
  
  names(d.ow.post) <- m.name
  
  
  # FDR correction
  
  d.ow.fdr <- as.data.frame(apply(d.ow.post,
                                  1,
                                  function(x) as.data.frame(p.adjust(x, 
                                                                     method= "fdr"))))
  
  names(d.ow.fdr) <- row.names(d.ow.post)
  
  
  # Return list of raw and FDR-adjusted p-values
  
  return(list("Raw P" = data.frame("Name" = row.names(t(d.ow.post)),
                                   t(d.ow.post)),
              "FDR P" = data.frame("Name" = row.names(d.ow.fdr),
                                   d.ow.fdr),
              "Sig Raw" = data.frame("Significant" = colSums(as.data.frame(t(d.ow.post)) < 0.05)),
              "Sig FDR" = data.frame("Significant" = colSums(d.ow.fdr < 0.05)
                                        )
              )
         )
  
}




