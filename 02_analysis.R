library(lme4)
library(emmeans)

# Test substrate x allele effects
## Data formatting

act_long <- CYP3A4.act %>%
  pivot_longer(
    cols = starts_with("CYP3A4"),   # allele columns
    names_to = "allele",
    values_to = "activity"
  )

act_long <- act_long %>%
  filter(!is.na(activity))

act_long$allele <- factor(act_long$allele)

df_subset <- act_long %>%
  group_by(allele, substrate) %>%
  filter(n() >= 2) %>%
  ungroup()

model_subset <- lmer(activity ~ allele * substrate + (1 | PMID), data = df_subset, REML=FALSE)
model_subset_no_int <- lmer(activity ~ allele + substrate + (1 | PMID), data = df_subset, REML=FALSE)
anova(model_subset_no_int, model_subset)

## Diagnostics
plot(model_subset)
qqnorm(resid(model_subset))
qqline(resid(model_subset))

# Repeat analysis with restricted dataset and log-transformed activity values
df_sub_multi <- df_subset %>%
  group_by(allele) %>%
  filter(n_distinct(substrate) >= 2) %>%
  ungroup() %>%
  droplevels()

model_sub_multi <- lmer(
  log(activity) ~ allele * substrate + (1 | PMID),
  data = df_sub_multi,
  REML = FALSE
)

model_sub_multi_no_int <- lmer(
  log(activity) ~ allele + substrate + (1 | PMID),
  data = df_sub_multi,
  REML = FALSE
)

anova(model_sub_multi_no_int, model_sub_multi)

## Diagnostics
isSingular(model_sub_multi)
df.residual(model_sub_multi)
plot(model_sub_multi)
qqnorm(resid(model_sub_multi))
qqline(resid(model_sub_multi))

# Estimated marginal means on the log scale
emm_sub <- emmeans(
  model_sub_multi,
  ~ substrate | allele
)

# Pairwise substrate comparisons within each allele
sub_contrasts <- pairs(
  emm_sub,
  adjust = "holm"
)

# Obtain effect estimates and 95% confidence intervals
results <- as.data.frame(
  summary(
    sub_contrasts,
    infer = c(TRUE, TRUE),
    level = 0.95
  )
)

# Convert log-scale differences to activity ratios
results <- results %>%
  mutate(
    activity_ratio = exp(estimate),
    ratio_lower = exp(lower.CL),
    ratio_upper = exp(upper.CL)
  )

# Forest plot of results

## Data formatting
plot_data <- results %>%
  filter(
    is.finite(activity_ratio),
    is.finite(ratio_lower),
    is.finite(ratio_upper),
    activity_ratio > 0,
    ratio_lower > 0,
    ratio_upper > 0
  ) %>%
  mutate(
    allele = as.character(allele),
    comparison = as.character(contrast)
  )

plot_data_grouped <- plot_data %>%
  mutate(
    allele = as.character(allele),
    
    panel = case_when(
      allele %in% c(
        "CYP3A4.3",
        "CYP3A4.9",
        "CYP3A4.10",
        "CYP3A4.11",
        "CYP3A4.12"
      ) ~ "CYP3A4.3, .9, .10, .11, .12",
      
      allele %in% c(
        "CYP3A4.2",
        "CYP3A4.7"
      ) ~ "CYP3A4.2 and .7",
      
      allele == "CYP3A4.16" ~ "CYP3A4.16",
      allele == "CYP3A4.18" ~ "CYP3A4.18"
    ),
    
    # Include allele in row labels so comparisons remain identifiable
    comparison_label = paste0(
      allele, ": ", comparison
    ),
    
    panel = factor(
      panel,
      levels = c(
        "CYP3A4.3, .9, .10, .11, .12",
        "CYP3A4.2 and .7",
        "CYP3A4.16",
        "CYP3A4.18"
      )
    )
  )

plot_data_grouped <- plot_data_grouped %>%
  mutate(
    p_label = case_when(
      is.na(p.value) ~ "NA",
      p.value < 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", sprintf("%.3f", p.value))
    )
  )

# Plot
forest_plot_grouped <- ggplot(
  plot_data_grouped,
  aes(
    x = activity_ratio,
    y = comparison_label
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "gray50"
  ) +
  
  geom_errorbar(
    aes(
      xmin = ratio_lower,
      xmax = ratio_upper
    ),
    orientation = "y",
    width = 0.2,
    color = "#5B83AD",
    linewidth = 0.7
  ) +
  
  geom_point(
    color = "#5B83AD",
    size = 2.5
  ) +
  
  # Add p-values
  geom_text(
    aes(
      x = p_x,
      label = p_label
    ),
    hjust = 0,
    size = 3.2,
    color = "gray25"
  ) +
  
  facet_grid(
    panel ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_x_log10(
    expand = expansion(mult = c(0.05, 0.4))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Estimated activity ratio (log scale)",
    y = NULL
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    strip.background = element_rect(
      fill = "gray95",
      color = "gray80"
    ),
    strip.text.y = element_text(
      face = "bold",
      angle = 0
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    plot.margin = margin(10, 80, 10, 10)
  )

forest_plot_grouped
