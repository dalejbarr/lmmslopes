library("funfact") ## devtools::install_github("dalejbarr/funfact")
library("tidyverse")
library("parallel")

## for the design with one ws/wi factor (A)
make_data_A <- function(ns, ni, sv_subj, sv_item,
                        eff_A = 0, eff_B = 0,
                        r_sub = .6, r_itm = .6) {

  my_design <- list(ivs = list(A = 2L),
                    n_item = ni)

  params <- funfact::gen_pop(my_design, ns)
  params$fixed[] <- c(2000, eff_A)
  params$subj_rfx[, ] <- c(100^2, r_sub * 100 * sv_subj,
                           r_sub * 100 * sv_subj, sv_subj^2)
  params$item_rfx[, ] <- c(100^2, r_itm * 100 * sv_item,
                           r_itm * 100 * sv_item, sv_item^2)
  params$err_var <- 300^2

  funfact::sim_norm(my_design, ns, params) %>%
    funfact::with_dev_pred(iv_names = c("A"))
}

## design with A (ws/wi) and B (bs/bi)
make_data_AB <- function(ns, ni, sv_subj, sv_item,
                      eff_A = 0, eff_B = 0,
                      r_sub = .6, r_itm = .6) {
  
  my_design <- list(ivs = list(A = 2L,
                               B = 2L),
                    between_subj = c("B"),
                    between_item = c("B"),
                    n_item = ni)

  params <- funfact::gen_pop(my_design, ns)
  params$fixed[] <- c(2000, eff_A, eff_B, 0)
  params$subj_rfx[, ] <- c(100^2, r_sub * 100 * sv_subj,
                           r_sub * 100 * sv_subj, sv_subj^2)
  params$item_rfx[, ] <- c(100^2, r_itm * 100 * sv_item,
                           r_itm * 100 * sv_item, sv_item^2)
  params$err_var <- 300^2

  funfact::sim_norm(my_design, ns, params) %>%
    funfact::with_dev_pred(iv_names = c("A", "B"))
}

## compare the models named in 'comps'; return matrix
compare_mods <- function(comps, modfits) {
  sapply(comps, function(.x) {
    ms <- .x[["m"]]
    this_df <- .x[["df"]]
    m1 <- modfits[[ms[1]]]$value
    m2 <- modfits[[ms[2]]]$value
    devdiff <- deviance(m2) - deviance(m1)
    c(chisq = devdiff,
      p = pchisq(devdiff, df = this_df, lower.tail = FALSE))
  })
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

fit5 <- function(mcr.data, alpha = .2, test = "A") {
  if (test == "B") {
    mods <- c(max = Y ~ AA2 * BB2 + (AA2 | subj_id) + (AA2 | item_id),
              nrc = Y ~ AA2 * BB2 + (AA2 || subj_id) + (AA2 || item_id),
              zis = Y ~ AA2 * BB2 + (AA2 || subj_id) + (1 | item_id),
              zss = Y ~ AA2 * BB2 + (1 | subj_id) + (AA2 || item_id),
              rio = Y ~ AA2 * BB2 + (1 | subj_id) + (1 | item_id))
  } else {
    mods <- c(max = Y ~ AA2 + (AA2 | subj_id) + (AA2 | item_id),
              nrc = Y ~ AA2 + (AA2 || subj_id) + (AA2 || item_id),
              zis = Y ~ AA2 + (AA2 || subj_id) + (1 | item_id),
              zss = Y ~ AA2 + (1 | subj_id) + (AA2 || item_id),
              rio = Y ~ AA2 + (1 | subj_id) + (1 | item_id))
  }

  res <- sapply(mods, tryFit, mcr.data, REML = FALSE, simplify = FALSE)
  conv <- sapply(res, function(x) x[["converged"]])
  mod_res <- sapply(res, function(x) x[["value"]])
  result <- c(Max = NA_real_,
              AIC = NA_real_, AIC.ix = NA_real_,
              LRT = NA_real_, LRT.ix = NA_real_)

  ## maximal
  if (test == "A") {
    mod2 <- tryUpdate(. ~ . -AA2, mod_res[["max"]])
  } else if (test == "B") {
    mod2 <- tryUpdate(. ~ . -BB2, mod_res[["max"]])
  } else {
    stop("'test' must be 'A' or 'B'")
  }

  result["Max"] <- get_chisq(mod_res[["max"]], mod2[["value"]])

  ## AIC competition
  aic_val <- sapply(res, function(x) AIC(x[["value"]]))
  aic_win <- which.min(aic_val)
  mod_aic <- res[[aic_win]]
  if (test == "A") {
    mod_aic2 <- tryUpdate(. ~ . -AA2, mod_aic[["value"]])
  } else if (test == "B") {
    mod_aic2 <- tryUpdate(. ~ . -BB2, mod_aic[["value"]])
  } else {
    stop("'test' must be 'A' or 'B'")
  }
  
  result["AIC"] <- get_chisq(mod_aic[["value"]], mod_aic2[["value"]])
  result["AIC.ix"] <- aic_win

  ## LRT competition
  ## compare max to nrc
  lrt_comp <- compare_mods(list(list(m = c("max", "nrc"), df = 2L)), res)

  if (lrt_comp["p", 1] > alpha) { ## reject max, test nrc
    comps <- list(zis = list(m = c("nrc", "zis"), df = 1L),
                  zss = list(m = c("nrc", "zss"), df = 1L),
                  rio_zis = list(m = c("zis", "rio"), df = 1L),
                  rio_zss = list(m = c("zss", "rio"), df = 1L))
    result2 <- compare_mods(comps, res)
    if (result2["p", "zis"] > alpha) { ## reject nrc; test zis against rio
      if (result2["p", "rio_zis"] > alpha) {
        winner <- "rio"
      } else {
        winner <- "zis"
      }
    } else { ## retain nrc; test nrc against zss/rio
      if (result2["p", "zss"] > alpha) { ## reject nrc
        if (result2["p", "rio_zss"] > alpha) {
          winner <- "rio"
        } else {
          winner <- "zss"
        }
      } else {
        winner <- "nrc"
      }
    }
  } else { ## accept max
    winner <- "max"
  }
  
  if (test == "A") {
    mod2 <- tryUpdate(. ~ . -AA2, res[[winner]]$value)
  } else if (test == "B") {
    mod2 <- tryUpdate(. ~ . -BB2, res[[winner]]$value)
  } else {
    stop("'test' must be 'A' or 'B'")
  }
  
  result["LRT"] <- get_chisq(res[[winner]]$value, mod2[["value"]])
  result["LRT.ix"] <- which(names(mods) == winner)
  
  result
}

design_tbl_corr_A <- crossing(
  tibble(nsubj = as.integer(args[1]),
         nitems = as.integer(args[2])),
  tibble(eff_A = c(0, 25),
         eff_B = 0),
  tibble(svar_subj = seq(0, 120, 20),
         svar_item = seq(0, 120, 20))) %>%
  mutate(func = "make_data_A",
         alpha = .2,
         test = "A")

design_tbl_uncorr_A <- crossing(
  tibble(nsubj = as.integer(args[1]),
         nitems = as.integer(args[2])),
  tibble(eff_A = c(0, 25),
         eff_B = 0),
  tibble(svar_subj = seq(0, 120, 20)),
  tibble(svar_item = seq(0, 120, 20))) %>%
  mutate(func = "make_data_A",
         alpha = .2,
         test = "A")

## testing B: eff_B = c(0, 120)
