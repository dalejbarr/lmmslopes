options(tidyverse.quiet = TRUE)

args <- if (!interactive()) {
          commandArgs(trailingOnly = TRUE)
        } else {
          as.character(c(24L, 24L, 10L, "A"))
        }

if (length(args) < 4L) {
  stop("must supply at least 4 arguments to command line:\n",
       "  'nsubj', 'nitems', 'nmc', 'A/B', (uncorr)")
}

source("setup.R")
source("start_cluster.R")
  
nmc <- as.integer(args[3]) # number of monte carlo simulations
sample_pfx <- paste0(args[1], "_", args[2], "_", args[3], "_", args[4])
## print(sample_info)

corr <- TRUE
if (length(args) == 5L) {
  if (args[5] == "uncorr") {
    corr <- FALSE
  }
}
if (corr) {
  sample_info <- paste0(sample_pfx, "_corr")
  todo <- if (args[4] == "A") {
            design_tbl_corr_A
          } else {
            design_tbl_corr_B
          }
} else {
  sample_info <- paste0(sample_pfx, "_uncorr")
  todo <- if (args[4] == "A") {
            design_tbl_uncorr_A
          } else {
            design_tbl_uncorr_B
          }
}

.junk <- clusterCall(cl, function() {
  library(dplyr)
  library(tidyr)
})
clusterExport(cl, c("fit5", "make_data_A", "make_data_AB", "compare_mods",
                    "allfit", "tryFit", "tryUpdate", "get_chisq"))

res <- pmap(todo, function(nsubj, nitems, eff_A, eff_B,
                           svar_subj, svar_item,
                           func, alpha, test) {
  message("running ", nmc, " runs of ",
          nsubj, "S, ", nitems, "I; ",
          "eff = c(A=", eff_A, ", B=", eff_B, "); ",
          "slope variance of (", svar_subj, ", ", svar_item, "); ",
          "'", func, "'...")
  parSapply(cl, seq_len(nmc), function(.x) {
    repeat {
      dat <- do.call(func,
                     list(ns = nsubj, ni = nitems,
                          sv_subj = svar_subj, sv_item = svar_item,
                          eff_A = eff_A, eff_B = eff_B))
      mods <- allfit(dat, test)
      conv <- sapply(mods, function(x) x[["converged"]])
      if (all(conv)) break;
    }
    fit5(mods, alpha = alpha, test = test) 
  })
})

dtbl <- todo
dtbl[["results"]] <- res
fname <- paste0("results_", sample_info, "_lmer.rds")

saveRDS(dtbl, fname)
message("results saved to ", fname)

stopCluster(cl)
