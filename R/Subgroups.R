# 04_subgroups.R ---------------------------------------------------------
# Phase 3: do the findings differ across student groups?
#
# Design choices (driven by group sizes: male 166, D3 91, private 61, ...):
#  A. Construct-score differences by group (ANOVA, eta-squared).
#  B. Measurement invariance for gender, education level and field group,
#     using PE parcels and MLR. Items are treated as continuous here as an
#     approximation, because ordinal WLSMV models are unstable in small groups.
#     Public/private is NOT tested (n = 61 is too small for an 18-indicator CFA).
#  C. Multi-group PATH analysis on composite (mean) scores for all four
#     groupings. This is the only workable option for n = 61, and it ignores
#     measurement error. Label it EXPLORATORY in the paper.
#  Missing data: students with a missing item are dropped from these models
#  (about 10 of 535).
source("R/Setup.R")

dat <- readRDS("data/clean_indonesia.rds")

sc <- list(
  PE   = paste0("PE", 1:20),
  DIFF = c("CU1", "CU2", "CU3"),
  ATU  = paste0("ATU", 1:5),
  RISK = c("AUP1", "AUP2", "AUP3"),
  MIUA = paste0("MIUA", 1:5)
)
cons <- paste0(names(sc), "_s")   # *_s = composite score (avoids clashing with latent names)

# Composite scores, PE parcels and grouping factors ---------------------------
dat <- dat %>%
  mutate(
    PE_s   = rowMeans(across(all_of(sc$PE)),   na.rm = TRUE),
    DIFF_s = rowMeans(across(all_of(sc$DIFF)), na.rm = TRUE),
    ATU_s  = rowMeans(across(all_of(sc$ATU)),  na.rm = TRUE),
    RISK_s = rowMeans(across(all_of(sc$RISK)), na.rm = TRUE),
    MIUA_s = rowMeans(across(all_of(sc$MIUA)), na.rm = TRUE),
    PEp1 = rowMeans(across(all_of(paste0("PE", 1:4))),   na.rm = TRUE),
    PEp2 = rowMeans(across(all_of(paste0("PE", 5:8))),   na.rm = TRUE),
    PEp3 = rowMeans(across(all_of(paste0("PE", 9:12))),  na.rm = TRUE),
    PEp4 = rowMeans(across(all_of(paste0("PE", 13:16))), na.rm = TRUE),
    PEp5 = rowMeans(across(all_of(paste0("PE", 17:20))), na.rm = TRUE),
    g_gender = factor(if_else(female == 1, "Female", "Male"),
                      levels = c("Male", "Female")),
    g_uni    = factor(if_else(private == 1, "Private", "Public"),
                      levels = c("Public", "Private")),
    g_edu    = factor(if_else(s1 == 1, "S1", "D3"), levels = c("D3", "S1")),
    g_field  = factor(field_grp)
  )

group_vars <- c(gender = "g_gender", university = "g_uni",
                education = "g_edu", field = "g_field")
print(map(group_vars, ~ table(dat[[.x]])))

# =============================================================================
# A. Construct means by group
# =============================================================================
means_tbl <- imap_dfr(group_vars, function(gv, nm) {
  dat %>%
    group_by(group = as.character(.data[[gv]])) %>%
    summarise(n = n(), across(all_of(cons), ~ round(mean(.x, na.rm = TRUE), 3)),
              .groups = "drop") %>%
    mutate(grouping = nm, .before = 1)
})
print(means_tbl)

# eta-squared: .01 small, .06 medium, .14 large. For two groups,
# Cohen's d is roughly 2 * sqrt(eta2 / (1 - eta2)).
anova_tbl <- imap_dfr(group_vars, function(gv, nm) {
  map_dfr(cons, function(v) {
    s  <- summary(aov(as.formula(paste(v, "~", gv)), data = dat))[[1]]
    ss <- s[["Sum Sq"]]
    tibble(grouping = nm, construct = v,
           F = s[["F value"]][1], p = s[["Pr(>F)"]][1], eta2 = ss[1] / sum(ss))
  })
}) %>%
  group_by(grouping) %>%
  mutate(p_holm = p.adjust(p, method = "holm")) %>%
  ungroup() %>%
  mutate(across(c(F, p, eta2, p_holm), ~ round(.x, 4)))
print(anova_tbl, n = 40)
# Mean comparisons are only meaningful if the measurement is invariant (part B).

# =============================================================================
# B. Measurement invariance (gender, education, field)
# =============================================================================
meas_inv <- '
  PE   =~ PEp1 + PEp2 + PEp3 + PEp4 + PEp5
  DIFF =~ CU1 + CU2 + CU3
  ATU  =~ ATU1 + ATU2 + ATU3 + ATU4 + ATU5
  RISK =~ AUP1 + AUP2 + AUP3
  MIUA =~ MIUA1 + MIUA2 + MIUA3 + MIUA4 + MIUA5
'

run_invariance <- function(gv) {
  d <- dat %>% filter(!is.na(.data[[gv]]))
  conf   <- cfa(meas_inv, data = d, group = gv, estimator = "MLR")
  metric <- cfa(meas_inv, data = d, group = gv, estimator = "MLR",
                group.equal = "loadings")
  scalar <- cfa(meas_inv, data = d, group = gv, estimator = "MLR",
                group.equal = c("loadings", "intercepts"))
  get <- function(f, lab) {
    x <- fitMeasures(f, c("cfi.robust", "rmsea.robust", "srmr"))
    tibble(step = lab, cfi = x[["cfi.robust"]],
           rmsea = x[["rmsea.robust"]], srmr = x[["srmr"]])
  }
  out <- bind_rows(get(conf, "configural"), get(metric, "metric (loadings)"),
                   get(scalar, "scalar (+ intercepts)")) %>%
    mutate(d_cfi = cfi - lag(cfi), d_rmsea = rmsea - lag(rmsea),
           ok = case_when(is.na(d_cfi) ~ NA,
                          TRUE ~ d_cfi >= -0.01 & d_rmsea <= 0.015)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  list(table = out, fits = list(conf = conf, metric = metric, scalar = scalar))
}

# Criteria: a step is acceptable if CFI falls by no more than .01 and RMSEA
# rises by no more than .015 (Chen 2007).
inv_res <- imap(group_vars[c("gender", "education", "field")], function(gv, nm) {
  cat("\n===== invariance:", nm, "=====\n")
  r <- tryCatch(run_invariance(gv), error = function(e) {
    message("Invariance failed for ", nm, ": ", e$message); NULL })
  if (!is.null(r)) print(r$table)
  r
})
inv_tbl <- imap_dfr(compact(inv_res), ~ mutate(.x$table, grouping = .y, .before = 1))
# If metric fails, loadings differ by group; if only scalar fails, group MEAN
# comparisons are not trustworthy. Report this next to the comparisons in A.

# =============================================================================
# C. Multi-group path analysis on composite scores
# =============================================================================
paths <- c("a1", "a2", "a3", "c1", "c2", "c3", "b")
path_desc <- c(a1 = "PE -> ATU", a2 = "DIFF -> ATU", a3 = "RISK -> ATU",
               c1 = "PE -> MIUA", c2 = "DIFF -> MIUA", c3 = "RISK -> MIUA",
               b  = "ATU -> MIUA")

labelled <- function(nm, G) paste0("c(", paste0(nm, "_", 1:G, collapse = ","), ")*")

mg_syntax <- function(G) {
  paste0(
    "ATU_s  ~ ", labelled("a1", G), "PE_s + ", labelled("a2", G), "DIFF_s + ",
    labelled("a3", G), "RISK_s\n",
    "MIUA_s ~ ", labelled("c1", G), "PE_s + ", labelled("c2", G), "DIFF_s + ",
    labelled("c3", G), "RISK_s + ", labelled("b", G), "ATU_s\n")
}
plain_syntax <- "
  ATU_s  ~ PE_s + DIFF_s + RISK_s
  MIUA_s ~ PE_s + DIFF_s + RISK_s + ATU_s
"

run_paths <- function(gv) {
  d  <- dat %>% filter(!is.na(.data[[gv]])) %>%
    select(all_of(c(gv, cons)))
  G  <- nlevels(d[[gv]])
  lv <- levels(d[[gv]])
  
  # group.label fixes the group ORDER to the factor levels. Without it lavaan
  # orders groups by first appearance in the data, which can mislabel groups.
  fit_free <- sem(mg_syntax(G), data = d, group = gv, group.label = lv,
                  estimator = "MLR")
  fit_eq   <- sem(plain_syntax, data = d, group = gv, group.label = lv,
                  estimator = "MLR", group.equal = "regressions")
  
  overall <- lavTestLRT(fit_free, fit_eq)   # all seven paths equal across groups?
  
  # Which paths differ? One Wald test per path, Holm-adjusted across the seven.
  wald <- map_dfr(paths, function(p) {
    con <- paste(paste0(p, "_", 1:(G - 1), " == ", p, "_", 2:G), collapse = "\n")
    w <- lavTestWald(fit_free, constraints = con)
    tibble(path = path_desc[[p]], chisq = w$stat, df = w$df, p = w$p.value)
  }) %>% mutate(p_holm = p.adjust(p, "holm"))
  
  std <- standardizedSolution(fit_free) %>%
    as_tibble() %>%
    filter(op == "~") %>%
    transmute(group = lv[group], path = paste(rhs, "->", lhs),
              beta = est.std, se, pvalue, ci.lower, ci.upper)
  
  list(overall = overall, wald = wald, std = std, fit = fit_free,
       n = table(d[[gv]]))
}

path_res <- imap(group_vars, function(gv, nm) {
  cat("\n===== paths:", nm, "=====\n")
  r <- tryCatch(run_paths(gv), error = function(e) {
    message("Path model failed for ", nm, ": ", e$message); NULL })
  if (!is.null(r)) {
    print(r$n)
    cat("\nAll paths equal across groups?\n"); print(r$overall)
    cat("\nPath by path (Holm-adjusted):\n")
    print(r$wald %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
    cat("\nStandardized paths by group:\n")
    print(r$std %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 40)
  }
  r
})

# =============================================================================
# D. Save
# =============================================================================
write_csv(means_tbl, "output/group_construct_means.csv")
write_csv(anova_tbl, "output/group_anova.csv")
write_csv(inv_tbl,   "output/invariance_tests.csv")
walls <- imap_dfr(compact(path_res), ~ mutate(.x$wald, grouping = .y, .before = 1))
stds  <- imap_dfr(compact(path_res), ~ mutate(.x$std,  grouping = .y, .before = 1))
write_csv(walls, "output/subgroup_path_tests.csv")
write_csv(stds,  "output/subgroup_paths_std.csv")
saveRDS(compact(path_res), "output/subgroup_path_fits.rds")