library("dplyr")
library("purrr")
library("tidyr")
library("ggplot2")

mod_stats <- function(m, prefix = "") {
  ff <- factor(as.integer(m), levels = 1:5)
  ff2 <- as.vector(table(ff))
  if (prefix != "") {
    names(ff2) <- paste(prefix, paste0("m", 1:5), sep = "_")
  } else {
    names(ff2) <- paste0("m", 1:5)
  }
  ff2
}

performance_chi <- function(x, cutoff = 3.85) {
  sum(x >= cutoff, na.rm = TRUE) / sum(!is.na(x))
}

## type I error and power function
get_results <- function(x) {
  nmc <- ncol(x)
  pvals <- apply(x[c("Max", "AIC", "LRT"), ], 1, performance_chi)
  mods_AIC <- mod_stats(x["AIC.ix", ], "AIC")
  mods_LRT <- mod_stats(x["LRT.ix", ], "LRT")

  bind_cols(tibble(nmc = nmc),
            as_tibble(as.list(pvals)),
            as_tibble(as.list(mods_AIC)),
            as_tibble(as.list(mods_LRT)))
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


if (interactive()) {
  rfile <- 
} else {
  rfile <- commandArgs(TRUE)[1]
}

res <- bind_rows(readRDS("results_30_10_10000_A_uncorr_lmer.rds"),
                 readRDS("results_50_20_10000_A_uncorr_lmer.rds"))

pow_results <- res %>%
  mutate(map_df(results, get_results)) %>%
  select(-results) %>%
  mutate(simulation = factor(sprintf("%d subjects, %d items",
                                     nsubj, nitems),
                             levels = c("50 subjects, 20 items",
                                        "30 subjects, 10 items")))

pbys <- pow_results %>%
  select(nsubj:LRT, simulation) %>%
  pivot_longer(cols = c(Max, AIC, LRT),
               names_to = "model",
               values_to = "power") %>%
  select(-nsubj, -nitems, -eff_B, -test)

modsel <- pow_results %>%
  select(-nmc, -Max, -AIC, -LRT) %>%
  pivot_longer(cols = AIC_m1:LRT_m5,
               names_to = c("criterion", "model"),
               names_sep = "_",
               values_to = "count")

corr <- pbys %>%
  filter(near(svar_subj, svar_item))

corr %>%
  rename(slopes_var = svar_subj, p = power) %>%
  mutate(technique = case_when(model == "Max" ~ "maximal",
                               model == "AIC" ~ "model selection (AIC)",
                               TRUE ~ "model selection (LRT)"),
         analysis = if_else(near(eff_A, 0), "type I error", "power")) %>%
  select(analysis, simulation, slopes_var, technique, p) %>%
  arrange(simulation) %>%
  saveRDS(file = "../manuscript/data/results_lmm.rds")

ggplot(corr %>% filter(near(eff_A, 0)),
       aes(svar_subj, power, color = model)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 120, 20)) +
  coord_cartesian(ylim = c(0, .1)) +
  facet_wrap(~simulation)

ggplot(corr %>% filter(near(eff_A, 25)),
       aes(svar_subj, power, color = model)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 120, 20)) +
  facet_wrap(~simulation)

##ggplot(modsel %>% filter(near(eff_A, 0)),
ggplot(modsel %>% filter(near(svar_subj, svar_item), !near(eff_A, 0)),
       aes(svar_item, count, color = model)) +
  geom_point(aes(shape = model), alpha = .4, size = 2) +
  geom_line(alpha = .4) +
  scale_x_continuous(breaks = seq(0, 120, 20)) +
  facet_wrap(simulation~criterion)
