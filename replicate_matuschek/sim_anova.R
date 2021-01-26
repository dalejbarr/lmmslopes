source("setup.R")

parms_fixed[["eff"]] = 0
type_I <- mcRun("fitanova",
                mcr.fnArgs = list(wsbi = FALSE),
                mcr.cluster = cl,
                mcr.datFn = "simdat2",
                mcr.constant = parms_fixed,
                mcr.varying = full_set)

saveRDS(type_I, paste0("type_I_", sample_info, "_anova.rds"))

parms_fixed[["eff"]] = 25
pow <- mcRun("fitanova",
             mcr.fnArgs = list(wsbi = FALSE),             
             mcr.cluster = cl,
             mcr.datFn = "simdat2",
             mcr.constant = parms_fixed,
             mcr.varying = full_set)

saveRDS(pow, paste0("power_", sample_info, "_anova.rds"))

stopCluster(cl)

