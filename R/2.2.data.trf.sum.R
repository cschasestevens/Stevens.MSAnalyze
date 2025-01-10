#' Data Transformation and Quality Check
#'
#' Performs data imputation, log2-transformation, and
#' quality checks on an untargeted dataset.
#'
#' @param ld Data list generated from ms_input().
#' @param cl_var Variable for annotating sample correlation heatmap.
#' @param samp_id Variable containing sample IDs.
#' @return List containing formatted input data and metadata
#' for downstream analysis.
#' @examples
#'
#' # ms_data_check(d)
#'
#' @export
ms_data_check <- function(ld, cl_var, samp_id) {
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
  # Data normalization
  head(d[["data"]][1:10])
  head(d_prev)
?scale
## reverse median-centering log2 normalization
### test for first sample
# median centering log2 for raw data
t(d_prev[1, 6:ncol(d_prev)])
median(t(d_prev[1, 6:ncol(d_prev)]))
sort(t(d_prev[1, 6:ncol(d_prev)]))
sort(log2(t(d_prev[1, 6:ncol(d_prev)]) / median(t(d_prev[1, 6:ncol(d_prev)]))))

sort(log2(2 ^ t(d[["data"]][1, ])))
median(2 ^ t(d[["data"]][1, ]))
### mTIC normalization
d1.tic <- dplyr::select(
  dplyr::left_join(
    d1.tic,
    setNames(
      aggregate(
        d1.tic[["mTIC"]],
        list(
          d1.tic[["Group"]]
        ),
        function(x) 
          sum(x)
      ),
      c(
        "Group","mTIC.sum"
      )
    ),
    by = "Group"
  ),
  c(
    "mTIC.sum",
    everything()
  )
)

# Group TCC
d1.tic[,6:ncol(d1.tic)] <- as.data.frame(
  lapply(
    d1.tic[,6:ncol(
      d1.tic
    )],
    function(x) 
      (x/d1.tic[["TCC.million"]])*
      median(
        d1.tic[["TCC.million"]]
      )
  )
)

# Group mTIC
d1.tic[,7:ncol(d1.tic)] <- as.data.frame(
  lapply(
    d1.tic[,7:ncol(
      d1.tic
      )],
    function(x) 
      (x/d1.tic[["mTIC.sum"]])*
      median(
        d1.tic[["mTIC.sum"]]
        )
    )
  )

d1.tic <- data.frame(
  "mTIC.norm" = apply(
    d1.tic[,6:ncol(
      d1.tic
    )],
    1,
    function(x) 
      sum(x)
  ),
  d1.tic
)

p.tic2 <- ggplot(d1.tic, 
                 aes(x = .data[["Index"]], 
                     y = .data[["TCC.million"]])) +
  geom_line(alpha = 0.5) +
  geom_point(aes(color = as.factor(.data[["Group"]])),
             shape=16,
             size = 1,
             alpha = 0.5) +
  geom_smooth(color = "firebrick1",alpha = 0.5) +
  labs(y = "TCC (million cells)", 
       x = "Sample Index") +
  thm.mult +
  scale_color_manual(values = col1)


ggsave(
  "analysis/plot.compare.cell.counts.png",
  ggarrange(p.tic1,p.tic2,ncol = 2,common.legend = T),
  width = 14,
  height = 6,
  dpi = 600
)


### RSD before and after norm
## calculate group RSDs for each compound
sd(d1.tic[d1.tic$Group == "0.5.No BE","CE.16.1."])/
mean(d1.tic[d1.tic$Group == "0.5.No BE","CE.16.1."])



q.rsd <- lapply(
  d1.tic[,8:ncol(d1.tic)],
  function(x) 
    aggregate(
      x,
      list(d1.tic[["Group"]]),
      function(y) 
        (sd(y)/
           mean(y)*
           100
        )
    )
)

q.rsd <- q.rsd %>%
  purrr::reduce(
    dplyr::left_join,
    by = "Group.1"
  )

names(q.rsd) <- c("Group.RSD",names(d1.tic[,8:ncol(d1.tic)]))

p.rsd <- reshape2::melt(
  q.rsd,
  id.vars = "Group.RSD"
)

p.rsd[["value"]] <- round(
  as.numeric(
    p.rsd[["value"]]
  ),
  digits = 2
)

p.rsd <- p.rsd[order(p.rsd[["Group.RSD"]]),]
p.rsd[["ID"]] <- seq.int(1,nrow(p.rsd),1)
aggregate(
  p.rsd[["value"]],
  list(
    p.rsd[["Group.RSD"]]
  ),
  function(x) 
    median(x)
)


p.rsd3 <- ggplot(
  p.rsd, 
  aes(
    x = ID, 
    y = value
  )
) +
  geom_point(
    aes(
      color = .data[["Group.RSD"]]),
    shape=16, 
    size = 1,
    alpha = 0.5
  ) +
  # geom_text_repel(data = subset(p.rsd, 
  #                               value > 50),
  #                 aes(label = p.rsd[p.rsd$value > 50,"Name"]),
  #                 size = 3,
  #                 segment.size = 0.1,
  #                 segment.color = "grey",
  #                 show.legend = F,
  #                 bg.color = "black",
  #                 color = "white") +
  geom_hline(
    yintercept = max(aggregate(
      p.rsd[["value"]],
      list(
        p.rsd[["Group.RSD"]]
      ),
      function(x) 
        median(x)
    )[["x"]]),
    linetype = "dashed",
    color = "firebrick1") +
  labs(
    y = "Group % RSD", 
    x = "Index"
  ) +
  thm.mult +
  scale_color_manual(values = col1)


ggsave(
  "analysis/plot.qrsd.mtic.png",
  ggarrange(p.rsd2,p.rsd3,ncol = 2,common.legend = T),
  width = 14,
  height = 6,
  dpi = 400
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

  ## Sample correlation heatmap
  # Input
  h1 <- ld1[["data.pareto"]]
  h1md <- data.frame(
    "Group" = as.factor(ld1[["meta"]][[cl_var]]),
    "ID" = ld1[["meta"]][[samp_id]]
  )
  h2 <- t(
    magrittr::set_rownames(
      as.matrix(h1),
      h1md[["ID"]]
    )
  )
  ld1[["data.samp.corr"]] <- round(
    cor(
      h2,
      method = "spearman"
    ),
    digits = 2
  )
  # Heatmap colors
  fun_hm_col <- circlize::colorRamp2(
    c(-1, 0, 0.5, 1),
    colors = col_grad()[c(1, 3, 6, 12)]
  )
  set.seed(1234)
  fun_hm_bar <- list(
    cl_var = setNames(
      col_univ()[1:length(levels(h1md[["Group"]]))], # nolint
      as.character(levels(h1md[["Group"]]))
    )
  )
  # Annotations
  hm_anno_list <- list(
    hm_col <- ComplexHeatmap::HeatmapAnnotation( # nolint
      `Group` = h1md[["Group"]],
      col = fun_hm_bar,
      show_annotation_name = FALSE,
      show_legend = TRUE
    )
  )
  # Plot
  h_out <- ComplexHeatmap::Heatmap(
    ld1[["data.samp.corr"]],
    col = fun_hm_col,
    name = "Correlation",
    top_annotation = hm_anno_list[[1]],
    show_column_names = TRUE,
    show_row_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    heatmap_width = ggplot2::unit(16, "cm"),
    heatmap_height = ggplot2::unit(16, "cm"),
    column_title = "Sample Correlation"
  )

  lp1 <- list(
    "plot.dist" = ggpubr::ggarrange(
      dst(ld1[["data"]], "Input"),
      dst(ld1[["data.log2"]], "Log2"),
      dst(ld1[["data.pareto"]], "Pareto Scaled"),
      nrow = 1,
      ncol = 3
    ),
    "plot.cor" = h_out
  )

  return(
    list(
      "data" = ld1,
      "plots" = lp1
    )
  )
}
