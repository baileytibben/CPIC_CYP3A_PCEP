library(lme4)
#test model system and substrate effects
##data formatting
act_long <- CYP3A4.act %>%
  pivot_longer(
    cols = starts_with("CYP3A4"),   # allele columns
    names_to = "allele",
    values_to = "activity"
  )

act_long <- act_long %>%
  filter(!is.na(activity))

act_long$allele <- factor(act_long$allele)

##substrate
df_subset <- act_long %>%
  group_by(allele, substrate) %>%
  filter(n() >= 2) %>%
  ungroup()

model_subset <- lmer(activity ~ allele * substrate + (1 | PMID), data = df_subset)
model_subset_no_int <- lmer(activity ~ allele + substrate + (1 | PMID), data = df_subset)
anova(model_subset_no_int, model_subset)

##system
df_subset2 <- act_long %>%
    group_by(allele, system) %>%
    filter(n() >=2) %>%
    ungroup()

model_subset2 <- lmer(activity ~ allele * system + (1 | PMID), data = df_subset2)
model_subset2_no_int <- lmer(activity ~ allele + system + (1 | PMID), data = df_subset2)
anova(model_subset2_no_int, model_subset2)
