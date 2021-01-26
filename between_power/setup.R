library("funfact") ## devtools::install_github("dalejbarr/funfact")
library("tidyverse")
library("parallel")

make_data <- function(ns, ni, sv,
                      eff_A = 0, eff_B = 0,
                      r_sub = .6, r_itm = .6) {
  my_design <- list(ivs = list(A = 2L,
                               B = 2L),
                    between_subj = c("B"),
                    between_item = c("B"),
                    n_item = ni)

  params <- funfact::gen_pop(my_design, ns)
  params$fixed[] <- c(2000, eff_A, eff_B, 0)
  params$subj_rfx[, ] <- c(100^2, r_sub * 100 * sv,
                           r_sub * 100 * sv, sv^2)
  params$item_rfx[, ] <- c(100^2, r_itm * 100 * sv,
                           r_itm * 100 * sv, sv^2)
  params$err_var <- 300^2

  funfact::sim_norm(my_design, ns, params) %>%
    funfact::with_dev_pred(iv_names = c("A", "B"))
}

args <- if (!interactive()) {
          commandArgs(trailingOnly = TRUE)
        } else {
          c(24L, 24L, 10L)
        }

design_tbl <- crossing(
  tibble(nsubj = as.integer(args[1]),
         nitems = as.integer(args[2])),
  tibble(eff_B = 50),
  tibble(svar = seq(0, 120, 20)))
  
nmc <- as.integer(args[3]) # number of monte carlo simulations
sample_info <- paste0(args[1], "_", args[2], "_", args[3])
## print(sample_info)

source("start_cluster.R")
