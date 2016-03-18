library("dplyr")
library("tidyr")
library("simgen")
library("ggplot2")

path <- "orig_10000"

performance_chi <- function(x, cutoff = 3.85) {
    sum(x >= cutoff, na.rm = TRUE) / sum(!is.na(x))
}

performance_p <- function(x, cutoff = .05) {
    sum(x <= cutoff, na.rm = TRUE) / sum(!is.na(x))
}

## figure out which model was selected
get_model <- function(x, path) {
    which_model <- function(y) {
        mx <- dat[y[["from"]]:y[["to"]], c("AIC.ix", "LRT.ix")]
        mods <- c("Maximal", rep("Reduced", 3), "Intercept")
        mod_mx <- apply(mx, 2, function(mm) mods[mm])
        ff <- data_frame(selection = rep(c("AIC", "LRT"), each = nrow(mod_mx)),
                         model = c(mod_mx)) %>% count(selection, model)
        all <- expand.grid(selection = c("AIC", "LRT"),
                           model = c("Intercept", "Maximal", "Reduced"),
                           stringsAsFactors = FALSE)
        left_join(all, ff, c("selection", "model")) %>%
        mutate(n = ifelse(is.na(n), 0, n)) %>%
        group_by(selection) %>%
        mutate(p = n / sum(n)) %>% ungroup() %>% arrange(selection, model)
    }
    dat <- readRDS(paste0(path, "/", x[["fname"]]))
    nmc <- x[["nmc"]]
    slopes <- data_frame(slope = as.integer(sqrt(pvary[, "t11"])),
                         from = (seq_len(nrow(pvary)) - 1) * nmc + 1,
                         to = seq_len(nrow(pvary)) * nmc)
    slopes %>% group_by(slope) %>% do(which_model(.)) %>% ungroup()
}

## type I error and power function
get_results <- function(x, pvary, type_I, cutoffs = NULL, path) {
    calc_prop <- function(y, dat) {
        result <- NA_real_
        if (x[["type"]] == "anova") {
            cols <- c(`min_F'` = "pmf", F1xF2 = "pmax")
            fn <- performance_p
            c_max <- .05
            q <- .05
            c_fn <- ">"
        } else { # lmer
            cols <- c(Maximal = "Max", AIC = "AIC", LRT = "LRT")
            fn <- performance_chi
            c_max <- 3.85
            q <- .95
            c_fn <- "<"
        }
        r0 <- y[["from"]][1]
        r1 <- y[["to"]][1]
        result <- apply(dat[r0:r1, cols], 2, fn)
        if (type_I) {
            cutoff <- apply(dat[y[["from"]]:y[["to"]], cols], 2, quantile,
                            probs = q)
            bind_rows(
                data_frame(method = names(cols), var = "type_I",
                           value = result),
                data_frame(method = names(cols), var = "cutoff",
                           value = ifelse(do.call(c_fn, list(cutoff, c_max)),
                                          cutoff, cutoff)) )
                           ## value = ifelse(do.call(c_fn, list(cutoff, c_max)),
                           ##    c_max, cutoff)) )
        } else { ## power
            c_vals <- spread(y, method, cutoff)[names(cols)] %>% unlist()
            res2 <- mapply(function(cix, cutf) {
                apply(dat[r0:r1, cix, drop = FALSE], 2, fn, cutoff = cutf)
            }, cols, c_vals, USE.NAMES = FALSE)
            bind_rows(
                data_frame(method = names(cols), var = "uncorrected",
                           value = result),
                data_frame(method = names(cols), var = "corrected",
                           value = res2))
        }
    }
    nmc <- x[["nmc"]]
    dat <- readRDS(paste0(path, "/", x[["fname"]]))
    s1 <- data_frame(slope = as.integer(sqrt(pvary[, "t11"])),
                     from = (seq_len(nrow(pvary)) - 1) * nmc + 1,
                     to = seq_len(nrow(pvary)) * nmc)
    if (!type_I) {
        slopes <- merge(x %>% select(-fname), s1) %>%
        inner_join(cutoffs, c("nsubj", "nitem", "type", "slope")) %>%
        select(slope, from, to, method, cutoff)
    } else {
        slopes <- s1
    }
    slopes %>% group_by(slope) %>% do(calc_prop(., dat))
}

plot_results <- function(pdat, fname, nmc, ptitle, ylim = NULL,
                         cbPalette = c("#F0E442", "#009E73", "#56B4E9",
                                       "#E69F00", "#000000"),
                         hgt = 4.5) {
    all_results <- pdat %>%
    rename(analysis = var) %>%
    mutate(method = factor(method,
                           levels = rev(c("Maximal", "AIC", "LRT", "F1xF2", "min_F'"))),
           nitem = factor(nitem, rev(sort(unique(nitem)))),
           nsubj = factor(nsubj, rev(sort(unique(nsubj)))))

    line_ann <- all_results %>%
    filter(analysis == "type_I") %>% group_by(nsubj, nitem) %>% slice(1) %>%
    ungroup() %>%
    select(nsubj, nitem, analysis) %>%
    mutate(yintercept = .05)
    
    g <- ggplot(all_results,
                aes(slope, value, color = method, shape = method)) +
    geom_hline(data = line_ann, aes(yintercept = yintercept),
               linetype = 2, color = "#999999") +
    geom_point() + geom_line() +
    facet_grid(analysis ~ nsubj + nitem) +
    scale_color_manual(values = cbPalette) +
    scale_shape_manual(values = rev(c(19, 24, 25, 5, 0))) +
    scale_x_discrete(limits = seq(0, 120, 20)) +
    coord_cartesian(xlim = c(-10, 130), ylim = ylim) +
    theme_bw() +
    theme(legend.key = element_blank(),
          strip.background = element_rect(colour = NA)) +
    labs(x = "SD of random slopes", y = "proportion significant",
         title = paste0(ptitle, " (", nmc, " for each setting)"))
    
    ggsave(fname, g, width = 8, height = hgt)
    return(g)
}

pvary <- readRDS(paste0(path, "/", "parms_varying.rds"))

#########################################
## type I error and power performance

todo <-
  data_frame(fname = list.files(path,
                                pattern = "^(type_I)|(power)_.+\\.rds$"),
             analysis = sub("(type_I|power).+" ,"\\1", fname),
             nsubj = paste0(sub("(type_I|power)_([0-9]+)_.+", "\\2", fname),
                            " subjects"),
             nitem = paste0(sub("(type_I|power)_([0-9]+)_([0-9]+)_.+", "\\3",
                                fname), " items"),
             nmc = as.integer(sub("(type_I|power)_([0-9]+)_([0-9]+)_([0-9].+)_(.+)\\.rds$", "\\4", fname)),
             type = sub("(type_I|power)_([0-9]+)_([0-9]+)_([0-9]+)_(.+)\\.rds$", "\\5", fname))

nmc <- todo[["nmc"]][1]
stopifnot(length(unique(todo[["nmc"]])) == 1)

all_type_I <- todo %>%
  filter(analysis == "type_I") %>%
  group_by(nsubj, nitem, type) %>%
  do(get_results(., pvary, TRUE, path = path)) %>% ungroup() 

cutoffs_full <- filter(all_type_I, var == "cutoff") %>% select(-var) %>%
  rename(cutoff = value)

cutoffs_anti <- cutoffs_full %>%
  inner_join(data_frame(method = c("min_F'", "F1xF2",
                                   "Maximal", "AIC", "LRT"),
                        c_max = rep(c(.05, 3.85), c(2, 3))), "method") %>%
  mutate(cutoff = ifelse(method %in% c("min_F'", "F1xF2"),
                                ifelse(cutoff > c_max, c_max, cutoff),
                                ifelse(cutoff < c_max, c_max, cutoff))) %>%
  select(-c_max)

t1 <- filter(all_type_I, var == "type_I")

g_typeI <- plot_results(t1, paste0(path, "/sim_type_I.pdf"),
                        nmc, "Type I Error", c(0, .10))

power_c2n <- todo %>%
    filter(analysis == "power") %>%
    group_by(nsubj, nitem, type) %>%
    do(get_results(., pvary, FALSE, cutoffs_full, path = path)) %>%
       ungroup() %>%
    mutate(var = ifelse(var == "corrected", "corrected_to_nominal", var))

power_cfa <- todo %>%
    filter(analysis == "power") %>%
    group_by(nsubj, nitem, type) %>%
    do(get_results(., pvary, FALSE, cutoffs_anti, path = path)) %>%
       ungroup() %>%
    filter(var == "corrected") %>%
    mutate(var = "corrected_for_anticonservativity")

power_cfa %>% select(-type, -var) %>% spread(method, value) %>%
   mutate(maxdiff = ifelse(LRT > AIC, LRT, AIC) - Maximal) %>%
   arrange(desc(maxdiff))

all_power <- bind_rows(power_c2n, power_cfa) %>%
mutate(var = factor(var, c("uncorrected", "corrected_for_anticonservativity",
                           "corrected_to_nominal")))

g_power <- plot_results(all_power,
                        paste0(path, "/sim_power.pdf"), nmc, "Power",
                        hgt = 8)

########################
## which model was selected
todo_wh <- data_frame(fname = list.files(path, pattern = "lmer\\.rds$"),
           analysis = sub("(type_I|power).+" ,"\\1", fname),
           nsubj = paste0(sub("(type_I|power)_([0-9]+)_.+", "\\2", fname),
               " subjects"),
           nitem = paste0(sub("(type_I|power)_([0-9]+)_([0-9]+)_.+", "\\3",
               fname), " items"),
           nmc = as.integer(sub("(type_I|power)_([0-9]+)_([0-9]+)_([0-9]+)_.+",
               "\\4", fname)))

mods <- todo_wh %>%
    group_by(nsubj, nitem, analysis) %>%
    do(get_model(., path)) %>% ungroup() %>%
    mutate(analysis = factor(analysis, levels = c("type_I", "power")),
        nsubj = factor(nsubj, levels = rev(sort(unique(nsubj)))),
        nitem = factor(nitem, levels = rev(sort(unique(nitem)))),
        model = factor(model, levels = c("Maximal", "Reduced", "Intercept"))) %>%
    group_by(nsubj, nitem, model, slope, selection, model) %>%
    summarize(n = sum(n), p = mean(p)) %>% ungroup()

mods %>% filter(slope == 0) %>% mutate(nonint = model != "Intercept") %>%
   group_by(nsubj, nitem, nonint, selection) %>% summarize(p = sum(p))

mods %>% filter(slope == max(slope)) %>% mutate(nonint = model != "Intercept") %>%
   group_by(nsubj, nitem, nonint, selection) %>% summarize(p = sum(p))

gg <- ggplot(mods,
       aes(slope, p, color = model, shape = model)) +
    geom_point() + geom_line() +
    facet_grid(selection ~ nsubj + nitem) +
    scale_color_manual(values = c("#F0E442", "#009E73", "#56B4E9", "#E69F00", "#000000")) +
    scale_shape_manual(values = c(19, 24, 25)) +
    scale_x_discrete(limits = seq(0, 120, 20)) +
    coord_cartesian(xlim = c(-10, 130)) +
    theme_bw() +
    theme(legend.key = element_blank(),
          strip.background = element_rect(colour = NA)) +
    labs(x = "SD of random slopes", y = "proportion selected",
         title = paste0("simulation results (", nmc, " for each setting)"))
ggsave(paste0(path, "/sim_model.pdf"), gg, width = 8, height = 4.5)
