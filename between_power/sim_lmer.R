options(tidyverse.quiet = TRUE)
source("setup.R")

tryFit <- function(tf.formula, tf.data, ...) {
  converged <- TRUE
  w.handler <- function(w) {
    converged <- FALSE
    invokeRestart("muffleWarning")
  }
  arg.list <- c(list(formula=tf.formula, data=tf.data), list(...))
  list(value=withCallingHandlers(tryCatch(
         do.call(lme4::lmer, arg.list),
         error=function(e) e),
         warning=w.handler),
       converged=converged)
}

tryUpdate <- function(tf.formula, tf.model, ...) {
  converged <- TRUE
  w.handler <- function(w) {
    converged <- FALSE
    invokeRestart("muffleWarning")
  }
  arg.list <- c(list(object = tf.model, formula. = tf.formula), 
                list(...))
  list(value = withCallingHandlers(tryCatch(
         do.call("update", arg.list), 
         error = function(e) e), warning = w.handler),
       converged = converged)
}

get_chisq <- function(bigger, smaller) {
  deviance(smaller) - deviance(bigger)
}

better_model_LRT <- function(bigger, smaller, alpha = .2) {
  chi <- get_chisq(bigger, smaller)
  df_diff <- attr(logLik(bigger), "df") - attr(logLik(smaller), "df")
  pval <- pchisq(chi, df_diff, lower.tail = FALSE)
  if (pval <= alpha) {
    return(1L) ## smaller significantly worse, stay
  } else {
    return(2L) ## no sig diff, shift
  }
}

fit5 <- function(mcr.data) {
  mods <- c(max = Y ~ BB2 + (AA2 | subj_id) + (AA2 | item_id),
            nrc = Y ~ BB2 + (AA2 || subj_id) + (AA2 || item_id),
            zis = Y ~ BB2 + (AA2 || subj_id) + (1 | item_id),
            zss = Y ~ BB2 + (1 | subj_id) + (AA2 || item_id),
            rio = Y ~ BB2 + (1 | subj_id) + (1 | item_id))
  res <- sapply(mods, tryFit, mcr.data, REML = FALSE, simplify = FALSE)
  conv <- sapply(res, function(x) x[["converged"]])
  mod_res <- sapply(res, function(x) x[["value"]])
  result <- c(Max = NA_real_,
              AIC = NA_real_, AIC.ix = NA_real_,
              LRT = NA_real_, LRT.ix = NA_real_)

  ## maximal
  if (conv[["max"]]) {
    mod2 <- tryUpdate(. ~ . -BB2, mod_res[["max"]])
    result["Max"] <- get_chisq(mod_res[["max"]], mod2[["value"]])
  } else {}

  ## AIC competition
  aic_val <- sapply(res[conv], function(x) AIC(x[["value"]]))
  aic_win <- which.min(aic_val)
  mod_aic <- res[conv][[aic_win]]
  mod_aic2 <- tryUpdate(. ~ . -BB2, mod_aic[["value"]])
  
  if (mod_aic2[["converged"]]) {
    result["AIC"] <- get_chisq(mod_aic[["value"]], mod_aic2[["value"]])
    result["AIC.ix"] <- aic_win
  } else {}

  ## LRT competition
  paths <- list(c("max", "nrc", "zis", "rio"),
                c("max", "nrc", "zss", "rio"))
  wins <- sapply(paths, function(path) {
    mconv <- path[conv[path]]
    winner <- mconv[1]
    while(length(mconv) > 1) {
      mod_bet <- better_model_LRT(mod_res[[mconv[1]]], mod_res[[mconv[2]]])
      if (mod_bet == 1L) {
        winner <- mconv[1L]

        mconv <- winner ## stop
      } else {
        winner <- mconv[2L]
        mconv <- mconv[-1L] ## continue
      }
    }
    return(winner)
  })
  lrt_win_name <-
    wins[which.min(sapply(wins, function(x) which(names(mods) == x)))]
  lrt_win <- match(lrt_win_name, names(mods))
  mod2 <- tryUpdate(. ~ . -BB2, mod_res[[lrt_win]])
  if (mod2[["converged"]]) {
    result["LRT"] <- get_chisq(mod_res[[lrt_win]], mod2[["value"]])
    result["LRT.ix"] <- lrt_win
  } else {}
  
  return(result)
}

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
