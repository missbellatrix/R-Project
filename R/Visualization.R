source("R/Setup.R")
dat <- readRDS("data/clean_indonesia.rds")

# Mean score per construct (PE, CU, ATU, AUP, MIUA)
score_cols <- paste0(names(scales), "_score")
for (k in names(scales)) dat[[paste0(k, "_score")]] <- rowMeans(dat[scales[[k]]], na.rm = TRUE)

# 1. Which AI tools do students use?
p1 <- dat %>% count(ai_tool_lab) %>%
  ggplot(aes(reorder(ai_tool_lab, n), n)) +
  geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "AI tools used by students", x = NULL, y = "Respondents")

# 2. Respondents by field of study
p2 <- dat %>% count(field_lab) %>%
  ggplot(aes(reorder(field_lab, n), n)) +
  geom_col(fill = "darkorange") + coord_flip() +
  labs(title = "Respondents by field of study", x = NULL, y = "Respondents")

# 3. Construct scores by gender
p3 <- dat %>%
  mutate(gender_lab = if_else(female == 1, "Female", "Male")) %>%
  pivot_longer(all_of(score_cols), names_to = "construct", values_to = "score") %>%
  ggplot(aes(construct, score, fill = gender_lab)) +
  geom_boxplot() +
  labs(title = "Construct scores by gender", x = NULL, y = "Mean score (1-4)", fill = NULL)

# 4. Correlation heatmap of constructs
p4 <- cor(dat[score_cols], use = "pairwise.complete.obs") %>%
  as.data.frame() %>% rownames_to_column("a") %>%
  pivot_longer(-a, names_to = "b", values_to = "r") %>%
  ggplot(aes(a, b, fill = r)) +
  geom_tile() + geom_text(aes(label = round(r, 2))) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits = c(-1, 1)) +
  labs(title = "Correlation between constructs", x = NULL, y = NULL)

print(p1); print(p2); print(p3); print(p4)
ggsave("output/viz_ai_tools.png", p1, width = 7, height = 5)
ggsave("output/viz_field.png",    p2, width = 7, height = 5)
ggsave("output/viz_gender.png",   p3, width = 8, height = 5)
ggsave("output/viz_corr.png",     p4, width = 6, height = 5)