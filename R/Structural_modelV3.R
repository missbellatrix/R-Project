# 03_structural_model_v3.R -----------------------------------------------
# Phase 2 rerun with the revised measurement structure (v3).
#   PE   performance expectancy (PE1-PE20)
#   DIFF difficulties using AI   (CU1-CU3)   [CU4 dropped: comfort item]
#   ATU  attitude                (ATU1-ATU5)
#   RISK concerns about AI       (AUP1-AUP3) [exploratory: low reliability]
#   MIUA motivation to keep using AI
# Demographics are covariates (female, private, S1).
#
# Expected signs keep the ORIGINAL theory-based direction; they are not
# changed after looking at the data.
source("R/Setup.R")

dat <- readRDS("data/clean_indonesia.rds")

scales_v3 <- list(
  PE   = paste0("PE", 1:20),
  DIFF = c("CU1", "CU2", "CU3"),
  ATU  = paste0("ATU", 1:5),
  RISK = c("AUP1", "AUP2", "AUP3"),
  MIUA = paste0("MIUA", 1:5)
)
items_v3 <- unlist(scales_v3, use.names = FALSE)
meas <- make_measurement(scales_v3)

# 1. Models -------------------------------------------------------------------
struct_full <- '
  ATU  ~ a1*PE + a2*DIFF + a3*RISK + d1*female + d2*private + d3*s1
  MIUA ~ c1*PE + c2*DIFF + c3*RISK + b*ATU + g1*female + g2*private + g3*s1

  ind_PE   := a1*b
  ind_DIFF := a2*b
  ind_RISK := a3*b
  tot_PE   := c1 + a1*b
  tot_DIFF := c2 + a2*b
  tot_RISK := c3 + a3*b
'
struct_med <- '
  ATU  ~ a1*PE + a2*DIFF + a3*RISK + d1*female + d2*private + d3*s1
  MIUA ~ b*ATU + g1*female + g2*private + g3*s1

  ind_PE   := a1*b
  ind_DIFF := a2*b
  ind_RISK := a3*b
'
model_full <- paste(meas, struct_full, sep = "\n")
model_med  <- paste(meas, struct_med,  sep = "\n")

# 2. Primary estimation (WLSMV, ordinal items) --------------------------------
fit_full <- sem(model_full, data = dat, ordered = items_v3, estimator = "WLSMV")
fit_med  <- sem(model_med,  data = dat, ordered = items_v3, estimator = "WLSMV")

fit_tbl <- bind_rows(fit_tab(fit_full, "v3 partial mediation"),
                     fit_tab(fit_med,  "v3 full mediation"))
print(as.data.frame(fit_tbl))
print(lavTestLRT(fit_full, fit_med))
print(lavInspect(fit_full, "rsquare")[c("ATU", "MIUA")])

lat <- lavInspect(fit_full, "cor.lv")
print(round(lat, 3))
if (any(abs(lat[lower.tri(lat)]) > 0.85)) warning("Latent correlation > .85")

# 3. Hypothesis table ----------------------------------------------------------
hyp <- tribble(
  ~label,     ~hypothesis,                 ~expected,
  "a1",       "H1: PE -> ATU",             "+",
  "a2",       "H2: DIFF -> ATU",           "-",
  "a3",       "H3: RISK -> ATU",           "?",
  "b",        "H4: ATU -> MIUA",           "+",
  "c1",       "H5: PE -> MIUA (direct)",   "+",
  "c2",       "H6: DIFF -> MIUA (direct)", "-",
  "c3",       "H7: RISK -> MIUA (direct)", "?",
  "ind_PE",   "H8a: PE -> ATU -> MIUA",    "+",
  "ind_DIFF", "H8b: DIFF -> ATU -> MIUA",  "-",
  "ind_RISK", "H8c: RISK -> ATU -> MIUA",  "?",
  "g1",       "H9a: female -> MIUA",       "?",
  "g2",       "H9b: private -> MIUA",      "?",
  "g3",       "H9c: S1 -> MIUA",           "?"
)

pe <- parameterEstimates(fit_full, standardized = TRUE, ci = TRUE) %>% as_tibble()
std <- standardizedSolution(fit_full, type = "std.all") %>%
  as_tibble() %>%
  select(lhs, op, rhs, std_all = est.std)

hyp_tbl <- pe %>%
  filter(label %in% hyp$label) %>%
  select(label, lhs, op, rhs, est, se, z, pvalue, ci.lower, ci.upper) %>%
  left_join(hyp, by = "label") %>%
  left_join(std, by = c("lhs", "op", "rhs")) %>%
  mutate(supported = case_when(
    expected == "+" ~ est > 0 & pvalue < .05,
    expected == "-" ~ est < 0 & pvalue < .05,
    TRUE            ~ pvalue < .05
  )) %>%
  arrange(match(label, hyp$label)) %>%
  select(hypothesis, expected, est, se, z, pvalue, ci.lower, ci.upper,
         std_all, supported)

print(hyp_tbl %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 20)

print(pe %>% filter(label %in% c("tot_PE", "tot_DIFF", "tot_RISK")) %>%
        select(label, est, se, pvalue, ci.lower, ci.upper))

# 4. Robustness A: MLR + FIML ---------------------------------------------------
fit_mlr <- sem(model_full, data = dat, estimator = "MLR", missing = "fiml")

path_compare <- standardizedSolution(fit_full, type = "std.all") %>%
  as_tibble() %>%
  filter(op == "~") %>%
  select(lhs, rhs, wlsmv = est.std) %>%
  left_join(standardizedSolution(fit_mlr, type = "std.all") %>%
              as_tibble() %>%
              filter(op == "~") %>% select(lhs, rhs, mlr = est.std),
            by = c("lhs", "rhs")) %>%
  mutate(diff = round(wlsmv - mlr, 3))
print(path_compare)

# 5. Robustness B: bootstrapped indirect effects (slow, off by default) -------
run_boot <- FALSE
if (run_boot) {
  set.seed(2026)
  fit_boot <- sem(model_full, data = dat, estimator = "ML", missing = "fiml",
                  se = "bootstrap", bootstrap = 2000,
                  parallel = if (.Platform$OS.type == "windows") "snow" else "multicore",
                  ncpus = max(1, parallel::detectCores() - 1))
  boot_tbl <- parameterEstimates(fit_boot, boot.ci.type = "perc", ci = TRUE) %>%
    as_tibble() %>%
    filter(op == ":=") %>%
    select(label, est, se, ci.lower, ci.upper)
  print(boot_tbl)
  write_csv(boot_tbl, "output/indirect_effects_bootstrap_v3.csv")
}

# 6. Save ------------------------------------------------------------------------
write_csv(hyp_tbl,      "output/hypothesis_table_v3.csv")
write_csv(fit_tbl,      "output/structural_fit_v3.csv")
write_csv(path_compare, "output/wlsmv_vs_mlr_paths_v3.csv")
saveRDS(fit_full,   "output/fit_structural_full_v3.rds")
saveRDS(model_full, "output/model_full_syntax_v3.rds")

r2 <- lavInspect(fit_full, "rsquare")[c("ATU", "MIUA")]
write_csv(tibble(construct = names(r2), r2 = round(as.numeric(r2), 3)),
          "output/r2_v3.csv")

library(semPlot)
png("output/fig_structural_v3.png", width = 2000, height = 1200, res = 250)
semPaths(fit_full, what = "std", whatLabels = "std", structural = TRUE,
         layout = "tree2", residuals = FALSE, intercepts = FALSE,
         thresholds = FALSE, nCharNodes = 0, sizeLat = 9,
         edge.label.cex = 1.1, edge.color = "black")
dev.off()
