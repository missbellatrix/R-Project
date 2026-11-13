source("R/Setup.R")
dat <- readRDS("data/clean_indonesia.rds")

# ---------------------------------------------------------------
# Construct scores (mean of items per construct)
# ---------------------------------------------------------------
score_cols <- paste0(names(scales), "_score")
for (k in names(scales)) {
  dat[[paste0(k, "_score")]] <- rowMeans(dat[scales[[k]]], na.rm = TRUE)
}

# ---------------------------------------------------------------
# Helper: pie chart with % labels and custom colours
# ---------------------------------------------------------------
pie_chart <- function(x, title, colors) {
  d <- tibble(group = x) %>%
    count(group) %>%
    mutate(pct = n / sum(n),
           label = paste0(group, "\n", round(pct * 100), "%"))
  ggplot(d, aes(x = "", y = n, fill = group)) +
    geom_col(width = 1, colour = "white", linewidth = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5) +
    scale_fill_manual(values = colors) +
    labs(title = title) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# ---------------------------------------------------------------
# Demographic charts
# ---------------------------------------------------------------

# Fig 1: Gender distribution
p1 <- pie_chart(
  if_else(dat$female == 1, "Female", "Male"),
  "Gender distribution",
  c("Male" = "#55B5FF", "Female" = "#FF8070")
)

# Fig 2: Type of universities
p2 <- pie_chart(
  if_else(dat$private == 1, "Private University", "Public University"),
  "Type of universities",
  c("Private University" = "#FFC000", "Public University" = "#9966E6")
)

# Education level (D3 vs S1)
p3 <- pie_chart(
  if_else(dat$s1 == 1, "Undergraduate Degree (S1)", "Associate Degree (D3)"),
  "Education level",
  c("Associate Degree (D3)" = "#FF5C8A", "Undergraduate Degree (S1)" = "#80FFD4")
)

# Fig 4: Fields of study distribution
field_order <- c("Language", "Economics", "Natural sciences", "Social sciences",
                 "Education", "Mathematics", "Art", "Engineering")
field_cols  <- c("Language" = "#80FFD4", "Economics" = "#FF5500",
                 "Natural sciences" = "#2E7D32", "Social sciences" = "#1A9BD7",
                 "Education" = "#FFB3F5", "Mathematics" = "#55B5FF",
                 "Art" = "#FF5C8A", "Engineering" = "#8B3A0F")

p4 <- dat %>%
  count(field_lab) %>%
  mutate(field_lab = factor(field_lab, levels = field_order)) %>%
  ggplot(aes(field_lab, n, fill = field_lab)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = field_cols) +
  labs(title = "Fields of study distribution", x = NULL, y = "Respondents") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

# AI tools used by students
p5 <- dat %>%
  count(ai_tool_lab) %>%
  ggplot(aes(reorder(ai_tool_lab, n), n)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip() +
  labs(title = "AI tools used by students", x = NULL, y = "Respondents") +
  theme_minimal()

# ---------------------------------------------------------------
# Construct charts
# ---------------------------------------------------------------

# Construct scores by gender
p6 <- dat %>%
  mutate(gender_lab = if_else(female == 1, "Female", "Male")) %>%
  pivot_longer(all_of(score_cols), names_to = "construct", values_to = "score") %>%
  ggplot(aes(construct, score, fill = gender_lab)) +
  geom_boxplot() +
  labs(title = "Construct scores by gender", x = NULL,
       y = "Mean score (1-4)", fill = NULL) +
  theme_minimal()

# Correlation heatmap of constructs
p7 <- cor(dat[score_cols], use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  rownames_to_column("a") %>%
  pivot_longer(-a, names_to = "b", values_to = "r") %>%
  ggplot(aes(a, b, fill = r)) +
  geom_tile() +
  geom_text(aes(label = round(r, 2))) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", limits = c(-1, 1)) +
  labs(title = "Correlation between constructs", x = NULL, y = NULL) +
  theme_minimal()

# ---------------------------------------------------------------
# Show all plots
# ---------------------------------------------------------------
print(p1); print(p2); print(p3); print(p4)
print(p5); print(p6); print(p7)

# ---------------------------------------------------------------
# Save all plots to the output folder
# ---------------------------------------------------------------
ggsave("output/viz_gender_pie.png",     p1, width = 5, height = 5)
ggsave("output/viz_university_pie.png", p2, width = 5, height = 5)
ggsave("output/viz_education_pie.png",  p3, width = 5, height = 5)
ggsave("output/viz_fields_bar.png",     p4, width = 8, height = 5)
ggsave("output/viz_ai_tools.png",       p5, width = 7, height = 5)
ggsave("output/viz_gender_boxplot.png", p6, width = 8, height = 5)
ggsave("output/viz_corr.png",           p7, width = 6, height = 5)