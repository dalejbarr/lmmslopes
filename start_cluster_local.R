## TODO: customize the line below for your own computing cluster
ncores <- if ((parallel::detectCores() - 2L) < 1L) {
            1L
          } else {
            parallel::detectCores() - 2L
          }

if (Sys.getenv("USER") == "daleb" &&
    grepl("daleb-pc$", Sys.getenv("STY"))) {
  ncores <- rep(c("localhost", "chatter", "gossip", "yap"), c(1L, 6L, 6L, 6L))
  message("initializing ", length(ncores), " cores on ",
          paste(unique(ncores), collapse = ", "))
}

cl <- parallel::makePSOCKcluster(ncores,
                                 homogeneous = FALSE)
