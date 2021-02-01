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
  todo <- design_tbl_corr
} else {
  sample_info <- paste0(sample_pfx, "_uncorr")
  todo <- design_tbl_uncorr
}

.junk <- clusterCall(cl, function() {library(dplyr)})
clusterExport(cl, c("fit5", "make_data", "tryFit",
                    "tryUpdate", "get_chisq"))

res <- pmap(todo, function(nsubj, nitems, eff_B, svar_subj, svar_item) {
  message("running ", nmc, " runs of ",
          nsubj, " subjects with ", nitems, " items and ",
          "slope variance of (", svar_subj, ", ", svar_item, ")...")
  parSapply(cl, seq_len(nmc), function(.x) {
    dat <- make_data(nsubj, nitems, eff_B = eff_B,
                     sv_subj = svar_subj, sv_item = svar_item)
    fit5(dat)
  })
})

dtbl <- todo
dtbl[["results"]] <- res
fname <- paste0("results_", sample_info, "_lmer.rds")

saveRDS(dtbl, fname)
message("results saved to ", fname)

stopCluster(cl)
