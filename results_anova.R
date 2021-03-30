library("dplyr")
library("purrr")
library("tidyr")
library("ggplot2")

get_results_anova <- function(x) {
  tibble(minF = x["pmf", ],
         f1f2 = x["pmax", ])
}

rfile1 <- "results_30_10_10000_A_corr_anova.rds"
rfile2 <- "results_50_20_10000_A_corr_anova.rds"

res <- bind_rows(readRDS(rfile1),
                 readRDS(rfile2))

pow_results <- res %>%
  mutate(d = map(results, get_results_anova)) %>%
  select(-results) %>%
  unnest(c(d)) %>%
  mutate(minF_sig = minF < .05,
         f1f2_sig = f1f2 < .05) %>%
  group_by(nsubj, nitems, eff_A, svar_subj, svar_item) %>%
  summarize(p_minf = mean(minF_sig),
            p_f1f2 = mean(f1f2_sig),
            .groups = "drop") %>%
  mutate(simulation = factor(sprintf("%d subjects, %d items", nsubj, nitems),
                             levels = c("50 subjects, 20 items",
                                        "30 subjects, 10 items"))) %>%
  pivot_longer(c(p_minf, p_f1f2),
               names_to = "technique", values_to = "p") %>%
  mutate(technique = recode(technique,
                            p_minf = "min-F'",
                            p_f1f2 = "F1xF2"),
         analysis = factor(if_else(near(eff_A, 0), "type I error", "power"),
                           levels = c("type I error", "power")))

pow_results %>%
  select(analysis, simulation, slopes_var = svar_subj, technique, p) %>%
  arrange(analysis, simulation, slopes_var) %>%
  saveRDS("results_all_anova.rds")


ggplot(pow_results %>% filter(analysis == "power"),
       aes(svar_subj, p, color = technique)) +
  geom_hline(yintercept = .05, linetype = 2) +
  geom_point() +
  geom_line() +
  labs(x = "random slope variances", y = "Pr(reject H0)") +
  facet_wrap(~simulation)
