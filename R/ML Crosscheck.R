# 05_ml_crosscheck.R -----------------------------------------------------
# Phase 4: machine-learning cross-check of the SEM findings.
#
# What this DOES: asks whether the same predictors matter, in a similar order,
# when a flexible model (random forest, gradient boosting) predicts motivation,
# and whether non-linear / interaction structure is being missed by the SEM.
# What this does NOT do: validate causal paths or mediation.
#
# Outcome: MIUA composite (mean of MIUA1-5). Predictors: construct composite
# scores (v3 structure) + demographics. ATU is a MEDIATOR in the SEM, so models
# are run with and without it.
#
# One-time install:
# install.packages(c("tidymodels", "ranger", "xgboost", "vip", "shapviz",
#                    "mclust", "cluster"))
source("R/Setup.R")
suppressPackageStartupMessages({
  library(tidymodels)
  library(vip)
  library(shapviz)
})
conflicted::conflicts_prefer(vip::vi)
tidymodels_prefer()
# Optional speed-up:
# library(doParallel); registerDoParallel(max(1, parallel::detectCores() - 1))

set.seed(2026)

# 1. Data ---------------------------------------------------------------------
dat <- readRDS("data/clean_indonesia.rds")
sc <- list(
  PE   = paste0("PE", 1:20),
  DIFF = c("CU1", "CU2", "CU3"),
  ATU  = paste0("ATU", 1:5),
  RISK = c("AUP1", "AUP2", "AUP3"),
  MIUA = paste0("MIUA", 1:5)
)

ml <- dat %>%
  transmute(
    MIUA = rowMeans(across(all_of(sc$MIUA)), na.rm = TRUE),
    PE   = rowMeans(across(all_of(sc$PE)),   na.rm = TRUE),
    DIFF = rowMeans(across(all_of(sc$DIFF)), na.rm = TRUE),
    ATU  = rowMeans(across(all_of(sc$ATU)),  na.rm = TRUE),
    RISK = rowMeans(across(all_of(sc$RISK)), na.rm = TRUE),
    female, private, s1, chatgpt,
    field_grp = factor(field_grp),
    CU4, AUP4, AUP5          # items dropped from the v3 constructs (extended set)
  ) %>%
  drop_na()
cat("Rows used:", nrow(ml), "\n")

feature_sets <- list(
  all          = c("PE", "DIFF", "ATU", "RISK", "female", "private", "s1",
                   "chatgpt", "field_grp"),
  no_ATU       = c("PE", "DIFF", "RISK", "female", "private", "s1",
                   "chatgpt", "field_grp"),
  constructs   = c("PE", "DIFF", "ATU", "RISK"),
  demographics = c("female", "private", "s1", "chatgpt", "field_grp"),
  extended     = c("PE", "DIFF", "ATU", "RISK", "female", "private", "s1",
                   "chatgpt", "field_grp", "CU4", "AUP4", "AUP5")
)
feat_all <- feature_sets$all

# 2. Models and cross-validation ------------------------------------------------
# 10-fold CV repeated 3 times. Reported R-squared is OUT-OF-SAMPLE.
folds   <- vfold_cv(ml, v = 10, repeats = 3)
metrics <- metric_set(rsq, rmse)

make_rec <- function(feats, norm = TRUE) {
  r <- recipe(reformulate(feats, response = "MIUA"), data = ml) %>%
    step_dummy(all_nominal_predictors()) %>%
    step_zv(all_predictors())
  if (norm) r <- r %>% step_normalize(all_numeric_predictors())
  r
}

lin_spec <- linear_reg() %>% set_engine("lm")
rf_spec  <- rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression")
xgb_spec <- boost_tree(trees = tune(), tree_depth = tune(),
                       learn_rate = tune(), min_n = tune()) %>%
  set_engine("xgboost") %>%
  set_mode("regression")

xgb_grid <- crossing(trees = c(150, 300), tree_depth = c(2, 3),
                     learn_rate = c(0.03, 0.1), min_n = c(5, 15))

best_metrics <- function(res, label) {
  b <- select_best(res, metric = "rsq")
  collect_metrics(res) %>% filter(.config == b$.config) %>% mutate(model = label)
}

run_models <- function(feats) {
  rec_lin  <- make_rec(feats, norm = TRUE)
  rec_tree <- make_rec(feats, norm = FALSE)
  p <- ncol(bake(prep(rec_tree), new_data = NULL)) - 1
  rf_grid <- crossing(mtry = unique(pmin(c(2, 3, 5), p)), min_n = c(5, 10, 20))
  
  lin <- workflow(rec_lin, lin_spec) %>%
    fit_resamples(folds, metrics = metrics)
  rf  <- workflow(rec_tree, rf_spec) %>%
    tune_grid(folds, grid = rf_grid, metrics = metrics)
  xgb <- workflow(rec_tree, xgb_spec) %>%
    tune_grid(folds, grid = xgb_grid, metrics = metrics)
  
  list(
    metrics = bind_rows(collect_metrics(lin) %>% mutate(model = "linear"),
                        best_metrics(rf,  "random forest"),
                        best_metrics(xgb, "xgboost")),
    rf = rf, xgb = xgb
  )
}

res <- imap(feature_sets, function(feats, nm) {
  cat("\n--- feature set:", nm, "---\n")
  run_models(feats)
})

cmp <- imap_dfr(res, ~ mutate(.x$metrics, feature_set = .y)) %>%
  filter(.metric %in% c("rsq", "rmse")) %>%
  select(feature_set, model, .metric, mean, std_err) %>%
  pivot_wider(names_from = .metric, values_from = c(mean, std_err)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  arrange(feature_set, model)
print(cmp, n = 30)

null_rmse <- sd(ml$MIUA) * sqrt((nrow(ml) - 1) / nrow(ml))
cat("\nBaseline RMSE (predict the mean):", round(null_rmse, 3), "\n")
# Reading it: if random forest / xgboost are no better than linear, the
# relationships are essentially linear and additive, so the SEM is not missing
# important non-linearity. If they are clearly better, look at SHAP dependence.

# 3. What matters? ------------------------------------------------------------
# 3a. Random forest permutation importance (final fit on all rows)
rec_tree_all <- make_rec(feat_all, norm = FALSE)
rf_final <- workflow(rec_tree_all, rf_spec) %>%
  finalize_workflow(select_best(res$all$rf, metric = "rsq")) %>%
  fit(ml)
imp_rf <- vi(extract_fit_parsnip(rf_final)) %>% mutate(Importance = round(Importance, 4))
print(imp_rf)
# Plot built with ggplot2 directly (vip() can return a non-ggplot object here)
p_imp <- imp_rf %>%
  mutate(Variable = reorder(Variable, Importance)) %>%
  ggplot(aes(Importance, Variable)) +
  geom_col(fill = "blue") +
  labs(y = NULL, title = "Random forest permutation importance") +
  theme_minimal()
ggsave("output/ml_rf_importance.png", p_imp, width = 7, height = 4)

# 3b. Leave-one-feature-out: how much out-of-sample R-squared is lost?
# Correlated predictors (PE, ATU, RISK) share credit, so a small drop does NOT
# mean a variable is unimportant, only that others compensate.
rf_default <- rand_forest(trees = 500, min_n = 5) %>%
  set_engine("ranger") %>% set_mode("regression")
cv_rsq <- function(feats, spec, norm) {
  workflow(make_rec(feats, norm), spec) %>%
    fit_resamples(folds, metrics = metric_set(rsq)) %>%
    collect_metrics() %>% pull(mean)
}
drops <- list(PE = "PE", DIFF = "DIFF", ATU = "ATU", RISK = "RISK",
              demographics = c("female", "private", "s1", "chatgpt", "field_grp"))
loco <- bind_rows(
  tibble(model = "linear",
         dropped = names(drops),
         rsq_without = map_dbl(drops, ~ cv_rsq(setdiff(feat_all, .x), lin_spec, TRUE)),
         rsq_full = cv_rsq(feat_all, lin_spec, TRUE)),
  tibble(model = "random forest",
         dropped = names(drops),
         rsq_without = map_dbl(drops, ~ cv_rsq(setdiff(feat_all, .x), rf_default, FALSE)),
         rsq_full = cv_rsq(feat_all, rf_default, FALSE))
) %>% mutate(loss = round(rsq_full - rsq_without, 3),
             across(c(rsq_without, rsq_full), ~ round(.x, 3)))
print(loco)

# 3c. SHAP values from the tuned gradient-boosting model
xgb_final <- workflow(rec_tree_all, xgb_spec) %>%
  finalize_workflow(select_best(res$all$xgb, metric = "rsq")) %>%
  fit(ml)
X_shap <- bake(extract_recipe(xgb_final), new_data = ml, all_predictors())
sv <- shapviz(extract_fit_engine(xgb_final), X_pred = data.matrix(X_shap), X = X_shap)

shap_imp <- tibble(feature = colnames(get_shap_values(sv)),
                   mean_abs_shap = colMeans(abs(get_shap_values(sv)))) %>%
  arrange(desc(mean_abs_shap)) %>%
  mutate(mean_abs_shap = round(mean_abs_shap, 4))
print(shap_imp)
ggsave("output/ml_shap_beeswarm.png", sv_importance(sv, kind = "beeswarm"),
       width = 7, height = 5)
ggsave("output/ml_shap_dependence.png", sv_dependence(sv, v = c("ATU", "PE")),
       width = 9, height = 4)
# SHAP shows how each predictor moves the prediction. Compare the ORDER of
# importance with the standardized SEM paths (ATU > PE > RISK > DIFF, small
# demographics). Order agreeing is the cross-check; exact sizes will differ.

# 4. Respondent segments (clustering) ---------------------------------------------
cl_vars <- c("PE", "DIFF", "ATU", "RISK", "MIUA")
Z <- as.data.frame(scale(ml[cl_vars]))

sil <- map_dfr(2:6, function(k) {
  km <- kmeans(Z, centers = k, nstart = 50)
  s  <- cluster::silhouette(km$cluster, dist(Z))
  tibble(k = k, avg_silhouette = round(mean(s[, "sil_width"]), 3))
})
print(sil)
# Silhouette guide: above .50 clear clusters; .26 to .50 weak; below .25 none.
# A low value means respondents form a CONTINUUM, not distinct types.

k_use <- sil$k[which.max(sil$avg_silhouette)]   # override, e.g. k_use <- 3
km <- kmeans(Z, centers = k_use, nstart = 50)
ml$cluster <- factor(km$cluster)

profile <- ml %>%
  group_by(cluster) %>%
  summarise(n = n(),
            across(all_of(cl_vars), ~ round(mean(.x), 2)),
            female = round(mean(female), 2), private = round(mean(private), 2),
            s1 = round(mean(s1), 2), chatgpt = round(mean(chatgpt), 2),
            .groups = "drop")
print(profile)

# Model-based clustering (latent profile analysis); mclust:: avoids masking purrr
mc <- Mclust(Z, G = 1:6)
cat("\nmclust best model:", mc$modelName, "with", mc$G, "profiles\n")
print(mc$BIC)

# Are clusters related to demographics? (Cramer's V: .1 small, .3 medium)
cramer <- function(a, b) {
  t <- suppressWarnings(chisq.test(table(a, b)))
  tibble(chisq_p = round(t$p.value, 4),
         cramers_v = round(sqrt(unname(t$statistic) /
                                  (sum(t$observed) * (min(dim(t$observed)) - 1))), 3))
}
assoc <- bind_rows(
  cramer(ml$cluster, ml$s1)        %>% mutate(variable = "education (S1)"),
  cramer(ml$cluster, ml$private)   %>% mutate(variable = "university (private)"),
  cramer(ml$cluster, ml$female)    %>% mutate(variable = "gender (female)"),
  cramer(ml$cluster, ml$field_grp) %>% mutate(variable = "field group")
)
print(assoc)

pc <- prcomp(Z)
p_pca <- as_tibble(pc$x[, 1:2]) %>% mutate(cluster = ml$cluster) %>%
  ggplot(aes(PC1, PC2, colour = cluster)) + geom_point(alpha = 0.6) + theme_minimal()
ggsave("output/ml_clusters_pca.png", p_pca, width = 6, height = 5)

# 5. Save --------------------------------------------------------------------------
write_csv(cmp,      "output/ml_model_comparison.csv")
write_csv(imp_rf,   "output/ml_rf_importance.csv")
write_csv(loco,     "output/ml_leave_one_out.csv")
write_csv(shap_imp, "output/ml_shap_importance.csv")
write_csv(sil,      "output/ml_cluster_silhouette.csv")
write_csv(profile,  "output/ml_cluster_profile.csv")
write_csv(assoc,    "output/ml_cluster_demographics.csv")