#### Group-wise and Class-wise Fold Change Function ####

# Transpose

fun.t <- function(df) {
  
  df.t <- df %>%
    t() %>%
    as.data.frame()
  
  names(df.t) <- as.character(df.t[1,])
  
  df.t <- df.t[-1,]
  
  return(df.t)
  
}





# Class-wise fold changes

fun.fold3 <- function(df,md,grp,list.cls) {
  
  # Metabolite Group Means
  
  fold.mean <- aggregate(df[,(md + 1):ncol(df)],
                         by = list(df[[grp]]),
                         FUN = mean)
  
  # Transpose and add class column
  
  fold.mean <- fun.t(fold.mean)
  
  fold.mean <- as.data.frame(lapply(fold.mean,
                                    function(x) as.numeric(x))) %>%
    dplyr::mutate(Class = rownames(fold.mean)) %>%
    dplyr::select(Class, everything())
  
  # Filter classes prior to fold change calculation
  
  fold.mean <- fold.mean[grepl(list.cls,fold.mean[["Class"]]),]
  
  # Group Combinations  
  
  fold.comb <- as.data.frame(combn(levels(as.factor(fold.mean[["Class"]])),2))
  
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





# Log2-transform and missing/infinite check

fun.fold4 <- function(l1.num,grp,list.cls2) {
  
  fc <- fun.fold3(list.data[[l1.num]][[1]],
                  data.md,
                  grp,
                  list.cls2)
  
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


