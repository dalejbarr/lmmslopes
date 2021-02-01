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

design_tbl <- crossing(
  tibble(nsubj = as.integer(args[1]),
         nitems = as.integer(args[2])),
  tibble(eff_B = 120),
  tibble(svar = seq(0, 120, 20)))
  
nmc <- as.integer(args[3]) # number of monte carlo simulations
sample_info <- paste0(args[1], "_", args[2], "_", args[3])
## print(sample_info)

.junk <- clusterCall(cl, function() {library(dplyr)})
clusterExport(cl, c("fit5", "make_data", "tryFit",
                    "tryUpdate", "get_chisq", "better_model_LRT"))

res <- pmap(design_tbl, function(nsubj, nitems, eff_B, svar) {
  message("running ", nmc, " runs of ",
          nsubj, " subjects with ", nitems, " items and ",
          "slope variance of ", svar, "...")
  parSapply(cl, seq_len(nmc), function(.x) {
    dat <- make_data(nsubj, nitems, eff_B = eff_B, sv = svar)
    fit5(dat)
  })
})

dtbl <- design_tbl
dtbl[["results"]] <- res
fname <- paste0("results_", sample_info, "_lmer.rds")

saveRDS(dtbl, fname)
message("results saved to ", fname)

stopCluster(cl)
