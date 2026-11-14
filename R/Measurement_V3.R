# 02b_measurement_v3.R ---------------------------------------------------
# Phase 1c: revised measurement models (based on the questionnaire wording).
#   v2: DIFF (CU1-3) / ATU + AUP4, AUP5 / RISK (AUP1-3) / PE / MIUA
#   v3: same, but ATU = ATU1-5 only (AUP4, AUP5 dropped: they behaved like
#       intention items and pushed the ATU-MIUA correlation to .85)
# CU4 is dropped in both: it is a positively worded comfort item.
source("R/Setup.R")

dat <- readRDS("data/clean_indonesia.rds")

scales_v2 <- list(
  PE   = paste0("PE", 1:20),
  DIFF = c("CU1", "CU2", "CU3"),
  ATU  = c(paste0("ATU", 1:5), "AUP4", "AUP5"),
  RISK = c("AUP1", "AUP2", "AUP3"),
  MIUA = paste0("MIUA", 1:5)
)
scales_v3 <- scales_v2
scales_v3$ATU <- paste0("ATU", 1:5)

# AVE, composite reliability and Fornell-Larcker matrix for any fitted CFA
validity_tab <- function(fit) {
  l <- standardizedSolution(fit) %>% as_tibble() %>% filter(op == "=~")
  ave <- l %>%
    group_by(lhs) %>%
    summarise(k   = n(),
              AVE = mean(est.std^2),
              CR  = sum(est.std)^2 / (sum(est.std)^2 + sum(1 - est.std^2)),
              .groups = "drop")
  lat <- lavInspect(fit, "cor.lv")
  fl  <- lat
  diag(fl) <- sqrt(ave$AVE[match(rownames(lat), ave$lhs)])
  list(ave_cr = ave %>% mutate(across(c(AVE, CR), ~ round(.x, 3))),
       fornell_larcker = round(fl, 3))
}

# 1. Fit both versions -------------------------------------------------------
fits <- imap(list(v2 = scales_v2, v3 = scales_v3), ~ cfa(
  make_measurement(.x), data = dat,
  ordered = unlist(.x, use.names = FALSE), estimator = "WLSMV"))

fit_summary <- imap_dfr(fits, ~ fit_tab(.x, .y))
print(as.data.frame(fit_summary))          # includes SRMR

# 2. Validity checks ---------------------------------------------------------
# Fornell-Larcker: the diagonal (sqrt AVE) must exceed the numbers beside it.
for (nm in names(fits)) {
  cat("\n=====", nm, "=====\n")
  print(validity_tab(fits[[nm]]))
}

# HTMT (want < .85 strict, < .90 lenient)
htmt_v3 <- tryCatch(
  semTools::htmt(make_measurement(scales_v3), dat, missing = "listwise"),
  error = function(e) { message("htmt failed: ", e$message); NULL })
if (!is.null(htmt_v3)) print(round(htmt_v3, 3))

# 3. Can ATU and MIUA be told apart? ----------------------------------------
m3 <- make_measurement(scales_v3)
it3 <- unlist(scales_v3, use.names = FALSE)
fit_free <- cfa(m3, data = dat, ordered = it3, estimator = "WLSMV", std.lv = TRUE)
fit_fix  <- cfa(paste(m3, "ATU ~~ 1*MIUA", sep = "\n"), data = dat,
                ordered = it3, estimator = "WLSMV", std.lv = TRUE)
print(lavTestLRT(fit_free, fit_fix))
# Significant = the two constructs are statistically separable (still report
# the high correlation honestly).

# 4. Save ----------------------------------------------------------------------
write_csv(fit_summary, "output/cfa_fit_v2_v3.csv")
write_csv(validity_tab(fits$v3)$ave_cr, "output/reliability_validity_v3.csv")
saveRDS(fits$v3, "output/fit_cfa_v3.rds")
vt <- validity_tab(fits$v3)
write.csv(vt$fornell_larcker, "output/fornell_larcker_v3.csv")
if (!is.null(htmt_v3)) write.csv(round(htmt_v3, 3), "output/htmt_v3.csv")
write_csv(as_tibble(lavInspect(fits$v3, "cor.lv"), rownames = "construct"),
          "output/latent_correlations_v3.csv")
standardizedSolution(fits$v3) %>% as_tibble() %>% filter(op == "=~") %>%
  select(construct = lhs, item = rhs, loading = est.std, se, pvalue) %>%
  write_csv("output/cfa_loadings_v3.csv")