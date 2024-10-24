#' Multivariate ANOVA
#'
#' Conducts an ANOVA for each metadata variable with Tukey's post-hoc
#' analysis for each compound and combines with fold changes
#' for each group comparison.
#'
#' @param matpv Data matrix for conducting ANOVA.
#' @param matfc Data matrix for calculating fold change.
#' @param md Metadata from ms_input().
#' @param md_var Metadata variables to test by ANOVA.
#' @param an Reference annotation data frame.
#' @param fc_class1 Logical indicating if compound classes rather
#' than individual compounds should be tested for significant differences.
#' @param grp_class1 If fc_class1 is TRUE, which variable should be used
#' to group individual compounds into separate compound classes?
#' @return A data frame containing the combined ANOVA and fold change results
#' for the specified variable(s).
#' @examples
#'
#' # ms_stat_anova(
#' #   matpv = d[["data"]][["data.pareto"]],
#' #   matfc = d[["data"]][["imputed"]],
#' #   md = d[["data"]][["meta"]],
#' #   md_var = c("Organ", "Sex", "Dose"),
#' #   an = d[["data"]][["anno"]],
#' #   fc_class1 = FALSE,
#' #   grp_class1 = "Subclass"
#' # )
#'
#' @export
ms_stat_anova <- function(matpv, matfc, md, md_var, an, fc_class1, grp_class1) { # nolint
  if(Sys.info()[["sysname"]] != "Windows") { # nolint
    if(fc_class1 == FALSE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      # Check for normality
      d1_check <- as.data.frame(
        apply(
          as.matrix(d1),
          2,
          function(x) mvnormtest::mshapiro.test(t(as.matrix(x)))[[2]]
        )
      )
      print(
        paste(
          length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ]),
          "of",
          nrow(d1_check),
          "compounds have normally distributed intensities."
        )
      )
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * 0.75),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[["Name"]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[["Name"]]
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * 0.75),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[["Name"]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[["Name"]]
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * 0.75),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[["Name"]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[["Name"]]
          )
        )
      }

      # Add fold change
      fold_all <- ms_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      head(d_out)
      write.table(
        d_out,
        "analysis/table.stats.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }

    if(fc_class1 == TRUE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      d1 <- aggregate(
        t(as.matrix(d1)),
        list(an1[[grp_class1]]),
        function(x) mean(x)
      )
      d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
      # Check for normality
      d1_check <- as.data.frame(
        apply(
          as.matrix(d1),
          2,
          function(x) mvnormtest::mshapiro.test(t(as.matrix(x)))[[2]]
        )
      )
      print(
        paste(
          length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ]),
          "of",
          nrow(d1_check),
          "classes have normally distributed intensities."
        )
      )
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * 0.75),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * 0.75),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            parallel::mclapply(
              mc.cores = ceiling(parallel::detectCores() * 0.75),
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }

      # Add fold change
      fold_all <- ms_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      write.table(
        d_out,
        "analysis/table.stats.class.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }
  }
  if(Sys.info()[["sysname"]] == "Windows") { # nolint
    if(fc_class1 == FALSE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      # Check for normality
      d1_check <- as.data.frame(
        apply(
          as.matrix(d1),
          2,
          function(x) mvnormtest::mshapiro.test(t(as.matrix(x)))[[2]]
        )
      )
      print(
        paste(
          length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ]),
          "of",
          nrow(d1_check),
          "compounds have normally distributed intensities."
        )
      )
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[["Name"]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[["Name"]]
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[["Name"]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[["Name"]]
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = an1[["Name"]][[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            an1[["Name"]]
          )
        )
      }

      # Add fold change
      fold_all <- ms_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      head(d_out)
      write.table(
        d_out,
        "analysis/table.stats.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }

    if(fc_class1 == TRUE){ # nolint
      # Determine interactions between selected variables
      # Load data
      d1 <- matpv
      md1 <- md
      an1 <- an
      d1 <- aggregate(
        t(as.matrix(d1)),
        list(an1[[grp_class1]]),
        function(x) mean(x)
      )
      d1 <- setNames(as.data.frame(t(d1[, 2:ncol(d1)])), d1[[1]])
      # Check for normality
      d1_check <- as.data.frame(
        apply(
          as.matrix(d1),
          2,
          function(x) mvnormtest::mshapiro.test(t(as.matrix(x)))[[2]]
        )
      )
      print(
        paste(
          length(d1_check[d1_check[[1]] > (0.05 / nrow(d1_check)), ]),
          "of",
          nrow(d1_check),
          "classes have normally distributed intensities."
        )
      )
      ## ANOVA
      if(length(md_var) == 3) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 7, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]]) *
                                  as.factor(md1[[md_var[[3]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 2) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 3, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]]) *
                                  as.factor(md1[[md_var[[2]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }
      if(length(md_var) == 1) { # nolint
        d_mult <- dplyr::bind_rows(
          setNames(
            lapply(
              seq.int(1, ncol(d1), 1),
              function(y) {
                d2 <- dplyr::bind_rows(
                  lapply(
                    seq.int(1, 1, 1),
                    function(x) {
                      as.data.frame(
                        as.matrix(
                          TukeyHSD(
                            aov(
                              d1[[y]] ~
                                as.factor(md1[[md_var[[1]]]])
                            )
                          )[[x]]
                        )
                      )
                    }
                  )
                )
                d2 <- dplyr::mutate(
                  d2,
                  "Name" = names(d1)[[y]],
                  "Comparison" = rownames(d2)
                )
                d2 <- d2[, c("Name", "Comparison", "p adj")]
                d2[["FDR"]] <- p.adjust(d2[["p adj"]], method = "fdr")
                return(d2)
              }
            ),
            names(d1)
          )
        )
      }

      # Add fold change
      fold_all <- ms_stat_fc( # nolint
        mat1 = matfc,
        md = md1,
        md_var = md_var,
        an = an1,
        fc_class = fc_class1,
        grp_class = grp_class1
      )

      # Fix anova comparisons to match fold change comparisons
      dm_fix <- data.frame(
        "Comparison" = unique(d_mult[["Comparison"]]),
        "Comparison.fc" = unlist(
          lapply(
            seq.int(1, length(unique(d_mult[["Comparison"]])), 1),
            function(x) {
              fx1 <- unique(d_mult[["Comparison"]])[[x]]
              fx2 <- ifelse(
                fx1 %in% unique(fold_all[["Comparison.fc"]]) == FALSE,
                paste(
                  unlist(strsplit(fx1, "-"))[[2]],
                  unlist(strsplit(fx1, "-"))[[1]],
                  sep = "-"
                ),
                fx1
              )
              return(fx2)
            }
          )
        )
      )

      d_out <- dplyr::left_join(
        dplyr::left_join(
          d_mult,
          dm_fix,
          by = "Comparison"
        ),
        fold_all,
        by = c("Comparison.fc", "Name")
      )
      write.table(
        d_out,
        "analysis/table.stats.class.txt",
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }
  }
  return(d_out)
}
