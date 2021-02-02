options(tidyverse.quiet = TRUE)

args <- if (!interactive()) {
          commandArgs(trailingOnly = TRUE)
        } else {
          c(24L, 24L, 10L)
        }

if (length(args) < 3L) {
  stop("must supply three arguments to command line: 'nsubj', 'nitems', 'nmc'")
}

source("setup.R")
source("start_cluster.R")
  
nmc <- as.integer(args[3]) # number of monte carlo simulations
sample_pfx <- paste0(args[1], "_", args[2], "_", args[3])
## print(sample_info)

corr <- TRUE
if (length(args) == 4L) {
  if (args[4] == "uncorr") {
    corr <- FALSE
  }
}
if (corr) {
  sample_info <- paste0(sample_pfx, "_corr")
  todo <- design_tbl_corr_A
} else {
  sample_info <- paste0(sample_pfx, "_uncorr")
  todo <- design_tbl_uncorr_A
}

.junk <- clusterCall(cl, function() {
  library(dplyr)
  library(tidyr)
})
clusterExport(cl, c("fit5", "make_data_A", "compare_mods",
                    "tryFit", "tryUpdate", "get_chisq"))

res <- pmap(todo, function(nsubj, nitems, eff_A, eff_B, svar_subj, svar_item,
                           func, alpha, test) {
  message("running ", nmc, " runs of ",
          nsubj, "S, ", nitems, "I; ",
          "eff = c(A=", eff_A, ", B=", eff_B, "); ",
          "slope variance of (", svar_subj, ", ", svar_item, "); ",
          "'", func, "'...")
  parSapply(cl, seq_len(nmc), function(.x) {
    dat <- do.call(func,
                   list(ns = nsubj, ni = nitems,
                        sv_subj = svar_subj, sv_item = svar_item,
                        eff_A = eff_A, eff_B = eff_B))
    fit5(dat, alpha = alpha, test = test)
  })
})

dtbl <- todo
dtbl[["results"]] <- res
fname <- paste0("results_", sample_info, "_lmer.rds")

saveRDS(dtbl, fname)
message("results saved to ", fname)

stopCluster(cl)
