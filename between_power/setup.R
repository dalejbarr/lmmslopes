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

fit5 <- function(mcr.data, alpha = .2) {
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
  comps <- list(nrc = list(m = c("max", "nrc"), df = 2L),
                zis = list(m = c("nrc", "zis"), df = 1L),
                zss = list(m = c("nrc", "zss"), df = 1L),
                rio = list(m = c("zss", "rio"), df = 1L))

  lrt_ps <- sapply(comps, function(.x) {
    ms <- .x[["m"]]
    this_df <- .x[["df"]]
    m1 <- res[[ms[1]]]$value
    m2 <- res[[ms[2]]]$value
    devdiff <- deviance(m2) - deviance(m1)
    c(chisq = devdiff,
      p = pchisq(devdiff, df = this_df, lower.tail = FALSE))
  })

  back_selection <- which(lrt_ps["p", ] < alpha)

  winner <- if (length(back_selection) == 0L) {
              ## random intercept wins
              "rio"
            } else {
              ## calculate winner
              names(mods)[min(back_selection)]
            }

  mod2 <- tryUpdate(. ~ . -BB2, res[[winner]]$value)
  if (mod2[["converged"]]) {
    result["LRT"] <- get_chisq(res[[winner]]$value, mod2[["value"]])
    result["LRT.ix"] <- which(names(mods) == winner)
  } else {}
  
  result
}

