library("simgen") ## devtools::install_github("dalejbarr/simgen")

simdat2 <- function(mcr.params) {
  mkSigma <- function(v1, v2, r) {
    cv <- r * sqrt(v1) * sqrt(v2)
    return(matrix(c(v1, cv, cv, v2), ncol = 2))
  }
  nsubj <- mcr.params[["nsubj"]]
  nitem <- mcr.params[["nitem"]]
  exp.list <- expand.grid(list(ItemID = rep(factor(1:nitem), each = 2),
                               SubjID = factor(1:nsubj)))
  exp.list$Cond <- c(-.5, .5)
  subj <- cbind(SubjID = factor(1:nsubj),
                MASS::mvrnorm(nsubj, 
                              mu = c(0, 0),
                              Sigma = mkSigma(mcr.params[["t00"]],
                                              mcr.params[["t11"]], 
                                              mcr.params[["rsub"]])))
  colnames(subj) <- c("SubjID", "ri.subj", "rs.subj")

  item <- cbind(ItemID = factor(1:nitem),
                MASS::mvrnorm(nitem, 
                              mu = c(0, 0),
                              Sigma = mkSigma(mcr.params[["w00"]],
                                              mcr.params[["w11"]], 
                                              mcr.params[["ritm"]])))
  colnames(item) <- c("ItemID", "ri.item", "rs.item")
  
  x <- merge(subj,
             merge(exp.list, item))[, c("SubjID", "ItemID", "Cond",
                                        colnames(subj)[-1], colnames(item)[-1])]
  x <- x[order(x$SubjID, x$ItemID), ]
  rownames(x) <- 1:nrow(x)
  x$err <- rnorm(nrow(x), mean = 0, sd = sqrt(mcr.params[["err"]]))
  x$Resp <- mcr.params[["int"]] + x$ri.subj + x$ri.item + mcr.params[["eff"]] * 
    x$Cond + x$err + x$rs.subj * x$Cond + x$rs.item * x$Cond
  x <- x[, c("SubjID", "ItemID", "Cond", "Resp")]
  x$SubjID <- factor(x$SubjID)
  return(x)
}

args <- commandArgs(trailingOnly = TRUE)

## parameters that are fixed across all simulations
parms_fixed <- list(
    nsubj = as.integer(args[1]),
    nitem = as.integer(args[2]), # sample size
    int = 2000, eff = 0, # fixed effect
    err = 300^2, # residual variance
    miss = 0, pMin = 0, pMax = 0, # missing data
    t00 = 100^2, # by-subj random intercept
    w00 = 100^2) # by-item random intercept

n_increments = 7L ## 
parms_vary <-
  cbind(t11 = seq(0, 120, length.out = n_increments)^2, # by-subj slope
        rsub = .6, 
        w11 = seq(0, 120, length.out = n_increments)^2, # by-item slope
        ritm = .6) 

nmc <- as.integer(args[3]) # number of monte carlo simulations
sample_info <- paste0(args[1], "_", args[2], "_", args[3])
print(sample_info)

saveRDS(parms_vary, "parms_varying.rds")

## each parameter setting is run 10000 times, a total of 100 million sims
full_set <- parms_vary[rep(seq_len(nrow(parms_vary)), each = nmc), ]

source("start_cluster.R")
initializeCluster(cl)
