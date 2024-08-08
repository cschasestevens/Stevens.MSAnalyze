#### Group-wise Fold Change Function ####

fun.fold <- function(df,md,grp) {
  
  # Group Combinations  
  
  fold.comb <- as.data.frame(combn(levels(as.factor(df[[grp]])),2))
  
  # Metabolite Group Means
  
  fold.mean <- aggregate(df[,(md + 1):ncol(df)],
                         by = list(df[[grp]]),
                         FUN = mean)
  
  # Group Fold Changes
  
  fun.comb <- function(x,y) {x / y}
  
  fold.all <- lapply(fold.comb,
                     function(x) {
                       
                       fun.comb(fold.mean[fold.mean[[1]] == x[[1]],2:ncol(fold.mean)],
                                fold.mean[fold.mean[[1]] == x[[2]],2:ncol(fold.mean)])
                       
                     })
  
  # Return Combined List
  
  fold.all <- data.table::rbindlist(fold.all)
  
  fold.all[["Comparison"]] <- paste(fold.comb[1,],
                                    fold.comb[2,],
                                    sep = "-")
  
  fold.all <- dplyr::select(fold.all,.data[["Comparison"]],
                            everything())
  
  return(fold.all)
  
}






fun.fold2 <- function(l1.num,grp) {
  
  fc <- fun.fold(list.data[[l1.num]][[1]],
                 data.md,
                 grp)
  
  data.fold <- fc
  
  fc.log2 <- cbind(data.fold[[1]],
                   as.data.frame(apply(data.fold[,2:ncol(data.fold)],
                                       2,
                                       function(x) log2(x))
                   )
  )
  
  names(fc.log2) <- c("Comparison",
                      names(fc.log2[2:ncol(fc.log2)]))
  
  
  ## Check for missing or infinite values
  
  m1 <- fc.log2[is.na(fc.log2),]
  
  inf1 <- fc.log2[apply(fc.log2,
                        1,
                        function(x) is.infinite(x)),]
  
  fc.list <- list("Input" = list.data[[l1.num]][[1]],
                  "Fold Change" = fc,
                  "Log2 Fold Change" = fc.log2,
                  "Missing" = m1,
                  "Infinite" = inf1)
  
  return(fc.list)
  
}





