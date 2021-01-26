## TODO: customize the line below for your own computing cluster
ncores <- if ((parallel::detectCores() - 2L) < 1L) {
            1L
          } else {
            parallel::detectCores() - 2L
          }

cl <- parallel::makeCluster(ncores)
