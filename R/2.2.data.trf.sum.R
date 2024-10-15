#' Data Transformation and Quality Check
#'
#' Performs data imputation, log2-transformation, and
#' quality checks on an untargeted dataset.
#'
#' @param ld Data list generated from ms_input().
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_data_check(d)
#'
#' @export
ms_data_check <- function(ld) {
  # input data from ms_input()
  ld1 <- ld
  # data imputation
  ## replace zeroes or NA with 1/10th of
  ## the lowest non-zero value present
  ld1[["imputed"]] <- ld1[["data"]]
  ld1[["imputed"]] <- setNames(
    as.data.frame(
      lapply(
        seq.int(1, ncol(ld1[["imputed"]]), 1),
        function(x) {
          ld1[["imputed"]][[x]][is.na(ld1[["imputed"]][[x]])] <- 0
          d1 <- ld1[["imputed"]][[x]]
          d1 <- ifelse(
            d1 == 0,
            round(0.1 * min(d1[d1 > 0]), digits = 0),
            d1
          )
          return(d1)
        }
      )
    ),
    c(names(ld1[["imputed"]]))
  )
  # Log2-transformation
  ld1[["data.log2"]] <- setNames(
    as.data.frame(
      lapply(
        seq.int(
          1,
          ncol(ld1[["imputed"]]),
          1
        ),
        function(x) {
          log2(ld1[["imputed"]][[x]])
        }
      )
    ),
    c(names(ld1[["imputed"]]))
  )
  # Pareto scaling
  ld1[["data.pareto"]] <- setNames(
    as.data.frame(
      apply(
        apply(
          as.matrix(ld1[["data.log2"]]),
          2,
          function(x) x - mean(x)
        ),
        2,
        function(y) round(y / sqrt(sd(y)), digits = 2)
      )
    ),
    c(names(ld1[["data.log2"]]))
  )
  # Data distributions
  dst <- function(df, n1) {
    dsty <- as.data.frame(colMeans(df))
    names(dsty) <- c("Value")
    p <- ggplot2::ggplot(
      dsty,
      ggplot2::aes(x = Value) # nolint 
    ) +
      # Density Plot
      ggplot2::geom_density(
        color = "darkslategrey",
        fill = col_univ()[[1]] # nolint
      ) +
      # Plot Theme
      ggplot2::labs(
        title = paste(n1, "Data Distribution"),
        y = "Density"
      ) +
      ms_theme1() +
      ggplot2::theme(
        plot.margin = ggplot2::unit(c(rep(0.5, 4)), "cm")
      )
    return(p)
  }
  lp1 <- list(
    "plot.dist" = ggpubr::ggarrange(
      dst(ld1[["data"]], "Input"),
      dst(ld1[["data.log2"]], "Log2"),
      dst(ld1[["data.pareto"]], "Pareto Scaled"),
      nrow = 1,
      ncol = 3
    )
  )

  return(
    list(
      "data" = ld1,
      "plots" = lp1
    )
  )
}
