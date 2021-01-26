source("setup.R")

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
    mods <- c(max = Resp ~ Cond + (Cond | SubjID) + (Cond | ItemID),
              nrc = Resp ~ Cond + (Cond || SubjID) + (Cond || ItemID),
              zis = Resp ~ Cond + (Cond || SubjID) + (1 | ItemID),
              zss = Resp ~ Cond + (1 | SubjID) + (Cond || ItemID),
              rio = Resp ~ Cond + (1 | SubjID) + (1 | ItemID))
    res <- sapply(mods, tryFit, mcr.data, REML = FALSE, simplify = FALSE)
    conv <- sapply(res, function(x) x[["converged"]])
    mod_res <- sapply(res, function(x) x[["value"]])
    result <- c(Max = NA_real_,
                AIC = NA_real_, AIC.ix = NA_real_,
                LRT = NA_real_, LRT.ix = NA_real_)

    ## maximal
    if (conv["max"]) {
        mod2 <- tryUpdate(. ~ . -Cond, mod_res[["max"]])
        result["Max"] <- get_chisq(mod_res[["max"]], mod2[["value"]])
    } else {}
    
    ## AIC competition
    aic_val <- sapply(res[conv], function(x) AIC(x[["value"]]))
    aic_win <- which.min(aic_val)
    mod_aic <- res[conv][[aic_win]]
    mod_aic2 <- tryUpdate(. ~ . -Cond, mod_aic[["value"]])
    
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
    mod2 <- tryUpdate(. ~ . -Cond, mod_res[[lrt_win]])
    if (mod2[["converged"]]) {
        result["LRT"] <- get_chisq(mod_res[[lrt_win]], mod2[["value"]])
        result["LRT.ix"] <- lrt_win
    } else {}
    
    return(result)
}

clusterExport(cl, c("fit5", "tryUpdate", "get_chisq", "better_model_LRT"))

parms_fixed[["eff"]] = 0
type_I <- mcRun("fit5", mcr.cluster = cl,
                mcr.datFn = "simdat2",
                mcr.constant = parms_fixed,
                mcr.varying = full_set)

saveRDS(type_I, paste0("type_I_", sample_info, "_lmer.rds"))

parms_fixed[["eff"]] = 25
pow <- mcRun("fit5", mcr.cluster = cl,
             mcr.datFn = "simdat2",
             mcr.constant = parms_fixed,
             mcr.varying = full_set)

saveRDS(pow, paste0("power_", sample_info, "_lmer.rds"))

stopCluster(cl)
