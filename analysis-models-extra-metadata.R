library(tidyverse)
library(phyloseq)
library(vegan)
library(mgcv)
library(gratia)
library(emmeans)
library(reshape2)

phage_metadata <- readRDS("rds/metadata.combined.subset.rds") %>% 
  mutate(age_months = floor(age_days/30)) %>% 
  mutate(diet = na_if(diet, "")) %>% 
  mutate(diet = fct_drop(diet))

# number of samples with delivery mode, sex, and diet metadata
phage_metadata %>% filter(!is.na(delivery_mode)) %>% dim() # 500 samples with delivery_mode metadata
phage_metadata %>% filter(!is.na(sex)) %>% dim() # 774 samples with sex metadata
phage_metadata %>% filter(!is.na(diet)) %>% dim() # 645 samples with diet metadata

phage_metadata %>% filter(!is.na(delivery_mode)) %>% filter(!is.na(sex)) %>% filter(!is.na(diet)) %>% dim()
# 280 samples have all 3 metadata variables

### ALPHA DIVERSITY VOTU - ADDITIONAL METADATA ###
ps.phage <- readRDS("rds/ps.phage.rds")

alpha.diversity.votu <- estimate_richness(ps.phage, measures=c("Shannon", "Observed")) %>% 
  mutate(Pielou = Shannon / log(Observed))
alpha.diversity.votu$run <- rownames(alpha.diversity.votu)
alpha.diversity.votu <- alpha.diversity.votu %>% 
  left_join(phage_metadata %>% select(run, depth, age_days, study, delivery_mode, subjectID, sex, diet)) %>% 
  mutate(age_months = floor(age_days/30))

alpha.diversity.votu.subset.sex <- alpha.diversity.votu %>% filter(!is.na(sex))
alpha.diversity.votu.subset.delivery_mode <- alpha.diversity.votu %>% filter(!is.na(delivery_mode))
alpha.diversity.votu.subset.diet <- alpha.diversity.votu %>% filter(!is.na(diet))


# pielou
n.sex <- nrow(alpha.diversity.votu.subset.sex)
n.delivery_mode <- nrow(alpha.diversity.votu.subset.delivery_mode)
n.diet <- nrow(alpha.diversity.votu.subset.diet)
alpha.diversity.votu.subset.sex.squeezed <- alpha.diversity.votu.subset.sex %>% mutate(Pielou_squeezed = (Pielou * (n.sex - 1) + 0.5) / n.sex)
alpha.diversity.votu.subset.delivery_mode.squeezed <- alpha.diversity.votu.subset.delivery_mode %>% mutate(Pielou_squeezed = (Pielou * (n.delivery_mode - 1) + 0.5) / n.delivery_mode)
alpha.diversity.votu.subset.diet.squeezed <- alpha.diversity.votu.subset.diet %>% mutate(Pielou_squeezed = (Pielou * (n.diet - 1) + 0.5) / n.diet)

model.pielou.sex.gam <- gam(Pielou_squeezed ~ s(age_months, by = sex) + log10(depth) + sex + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                        family = betar(),
                        method = "REML",
                        data = alpha.diversity.votu.subset.sex.squeezed)

model.pielou.delivery_mode.gam <- gam(Pielou_squeezed ~ s(age_months, by = delivery_mode) + log10(depth) + delivery_mode + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                            family = betar(),
                            method = "REML",
                            data = alpha.diversity.votu.subset.delivery_mode.squeezed)

model.pielou.diet.gam <- gam(Pielou_squeezed ~ s(age_months) + log10(depth) + diet + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                            family = betar(),
                            method = "REML",
                            data = alpha.diversity.votu.subset.diet.squeezed)

# shannon
model.shannon.sex.gam <- gam(Shannon ~ s(age_months, by = sex) + log10(depth) + sex + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                         data = alpha.diversity.votu.subset.sex,
                         family = gaussian(),
                         method = "REML")

model.shannon.delivery_mode.gam <- gam(Shannon ~ s(age_months, by = delivery_mode) + log10(depth) + delivery_mode + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                             data = alpha.diversity.votu.subset.delivery_mode,
                             family = gaussian(),
                             method = "REML")

model.shannon.diet.gam <- gam(Shannon ~ s(age_months) + log10(depth) + diet + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                             data = alpha.diversity.votu.subset.diet,
                             family = gaussian(),
                             method = "REML")

#richness
model.richness.sex.gam <- gam(Observed ~ s(age_months, by = sex) + log10(depth) + sex + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                          data = alpha.diversity.votu.subset.sex,
                          family = nb(),
                          method = "REML")

model.richness.delivery_mode.gam <- gam(Observed ~ s(age_months, by = delivery_mode) + log10(depth) + delivery_mode + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                              data = alpha.diversity.votu.subset.delivery_mode,
                              family = nb(),
                              method = "REML")

model.richness.diet.gam <- gam(Observed ~ s(age_months) + log10(depth) + diet + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                              data = alpha.diversity.votu.subset.diet,
                              family = nb(),
                              method = "REML")

library(plyr)
alpha.diversity.plot.data <- model.shannon.sex.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  rbind.fill(model.shannon.diet.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(.by = "diet")
  ) %>% 
  rbind.fill(model.shannon.delivery_mode.gam %>% 
               smooth_estimates() %>% 
               add_confint()
  ) %>% 
  mutate(metric = "shannon") %>% 
  rbind.fill(model.richness.sex.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "richness")
  ) %>% 
  rbind.fill(model.richness.diet.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "richness") %>% 
               mutate(.by = "diet")
  ) %>% 
  rbind.fill(model.richness.delivery_mode.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "richness")
  ) %>% 
  rbind.fill(model.pielou.sex.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "pielou")
  ) %>% 
  rbind.fill(model.pielou.diet.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "pielou") %>% 
               mutate(.by = "diet")
  ) %>% 
  rbind.fill(model.pielou.delivery_mode.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "pielou")
  ) %>% 
  filter(.type == "TPRS")

### TEST IF SEX/DELIVERY MODE SIGNIFICANTLY AFFECT ANY OF THE MODELS ###
emm.pielou.sex <- emmeans(model.pielou.sex.gam, 
                          pairwise ~ sex | age_months, 
                          at = list(age_months = seq(0, 36, by = 3)),
                          type = "response", data = alpha.diversity.votu.subset.sex)

emm.pielou.delivery_mode <- emmeans(model.pielou.delivery_mode.gam, 
                                    pairwise ~ delivery_mode | age_months, 
                                    at = list(age_months = seq(0, 36, by = 3)),
                                    type = "response", data = alpha.diversity.votu.subset.delivery_mode)

emm.richness.sex <- emmeans(model.richness.sex.gam, 
                            pairwise ~ sex | age_months, 
                            at = list(age_months = seq(0, 36, by = 3)),
                            type = "response", data = alpha.diversity.votu.subset.sex)

emm.richness.delivery_mode <- emmeans(model.richness.delivery_mode.gam, 
                            pairwise ~ delivery_mode | age_months, 
                            at = list(age_months = seq(0, 36, by = 3)),
                            type = "response", data = alpha.diversity.votu.subset.delivery_mode)

emm.shannon.sex <- emmeans(model.shannon.sex.gam, 
                          pairwise ~ sex | age_months, 
                          at = list(age_months = seq(0, 36, by = 3)),
                          type = "response", data = alpha.diversity.votu.subset.sex)

emm.shannon.delivery_mode <- emmeans(model.shannon.delivery_mode.gam, 
                                    pairwise ~ delivery_mode | age_months, 
                                    at = list(age_months = seq(0, 36, by = 3)),
                                    type = "response", data = alpha.diversity.votu.subset.delivery_mode)

summary(emm.pielou.sex$contrasts) %>% filter(p.value <= 0.05) # no significant differences
summary(emm.pielou.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # no significant differences
summary(model.pielou.diet.gam) #dietmix p-value = 0.007

summary(emm.richness.sex$contrasts) %>% filter(p.value <= 0.05) # significant at 21, 24, 27 mo (F > M)
summary(emm.richness.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # significant at 15, 18, 21, 24 mo (CS < VG)
summary(model.richness.diet.gam) # dietformula p-value = 0.005, dietsolid_food p-value = 0.006

summary(emm.shannon.sex$contrasts) %>% filter(p.value <= 0.05) # no significant differences
summary(emm.shannon.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # significant at 18, 21, 24, 27 mo (CS < VG)
summary(model.shannon.diet.gam) # dietmix p-value = 0.031


### VIROME DEVELOPMENTAL VELOCITY (PHF) ###
dist.phf.wuni <- readRDS("rds/dist.phf.wuni.rds")

dist.phf.wuni.melt <- dist.phf.wuni %>% 
  as.matrix() %>% 
  melt() %>% 
  `colnames<-`(c("sample1", "sample2", "dist")) %>% 
  filter(!(sample1 == sample2)) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, subjectID, depth, delivery_mode, sex, diet), by = c("sample1" = "run")) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, subjectID, depth, delivery_mode, sex, diet), by = c("sample2" = "run"), suffix = c(".sample1", ".sample2")) %>% 
  filter(study.sample1 %in% c("Lim2015", "Zhao2017", "Liang2020a", "Beller2022", "Walters2023", "Garmaeva2024")) %>%
  filter(study.sample2 %in% c("Lim2015", "Zhao2017", "Liang2020a", "Beller2022", "Walters2023", "Garmaeva2024")) %>% 
  filter(subjectID.sample1 == subjectID.sample2) %>% 
  filter(age_months.sample1 < age_months.sample2) %>% 
  group_by(sample1) %>% 
  filter(age_months.sample2 == min(age_months.sample2)) %>%
  ungroup() %>% 
  mutate(dT = age_months.sample2 - age_months.sample1,
         vel = dist / dT,
         mean_log_depth = (log10(depth.sample1) + log10(depth.sample2))/2,
         study  = study.sample1) %>%
  select(dist, dT, vel, study, mean_log_depth, subjectID = subjectID.sample1, age_months.sample2, delivery_mode = delivery_mode.sample1, sex = sex.sample1, diet = diet.sample2)

dist.phf.wuni.melt.subset.sex <- dist.phf.wuni.melt %>% filter(!is.na(sex))
dist.phf.wuni.melt.subset.delivery_mode <- dist.phf.wuni.melt %>% filter(!is.na(delivery_mode))
dist.phf.wuni.melt.subset.diet <- dist.phf.wuni.melt %>% filter(!is.na(diet))

table(dist.phf.wuni.melt.subset.delivery_mode$delivery_mode)
table(dist.phf.wuni.melt.subset.sex$sex)
table(dist.phf.wuni.melt.subset.diet$diet)

model.velocity.sex.gam <- gam(vel ~ s(age_months.sample2, by = sex) + mean_log_depth + sex +
                            s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                          family = tw(),
                          method = "REML",
                          data = dist.phf.wuni.melt.subset.sex)

model.velocity.delivery_mode.gam <- gam(vel ~ s(age_months.sample2, by = delivery_mode) + mean_log_depth + delivery_mode +
                                s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                              family = tw(),
                              method = "REML",
                              data = dist.phf.wuni.melt.subset.delivery_mode)

model.velocity.diet.gam <- gam(vel ~ s(age_months.sample2) + mean_log_depth + diet +
                                s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                              family = tw(),
                              method = "REML",
                              data = dist.phf.wuni.melt.subset.diet)

velocity.plot.data <- model.velocity.sex.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  rbind.fill(model.velocity.delivery_mode.gam %>% 
               smooth_estimates() %>% 
               add_confint()) %>% 
  rbind.fill(model.velocity.diet.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(.by = "diet")) %>% 
  mutate(metric = "velocity") %>%
  #filter(.smooth == "s(age_months.sample2)") %>% 
  filter(.type == "TPRS")

emm.velocity.sex <- emmeans(model.velocity.sex.gam, 
                           pairwise ~ sex | age_months.sample2, 
                           at = list(age_months = seq(0, 24, by = 3)),
                           type = "response", data = dist.phf.wuni.melt.subset.sex)

emm.velocity.delivery_mode <- emmeans(model.velocity.delivery_mode.gam, 
                                     pairwise ~ delivery_mode | age_months.sample2, 
                                     at = list(age_months = seq(0, 24, by = 3)),
                                     type = "response", data = dist.phf.wuni.melt.subset.delivery_mode)

summary(emm.velocity.sex$contrasts) %>% filter(p.value <= 0.05) # significant at 5 mo
summary(emm.velocity.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # significant at 8 mo
summary(model.velocity.diet.gam) # dietformula p-value = 0.0002

### VIROME DEVELOPMENTAL VELOCITY (VOTU) ###
dist.bray <- readRDS("rds/dist.bray.rds")

dist.bray.melt <- dist.bray %>% 
  as.matrix() %>% 
  melt() %>% 
  `colnames<-`(c("sample1", "sample2", "dist")) %>% 
  filter(!(sample1 == sample2)) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, subjectID, depth, delivery_mode, sex, diet), by = c("sample1" = "run")) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, subjectID, depth, delivery_mode, sex, diet), by = c("sample2" = "run"), suffix = c(".sample1", ".sample2")) %>% 
  filter(study.sample1 %in% c("Lim2015", "Zhao2017", "Liang2020a", "Beller2022", "Walters2023", "Garmaeva2024")) %>%
  filter(study.sample2 %in% c("Lim2015", "Zhao2017", "Liang2020a", "Beller2022", "Walters2023", "Garmaeva2024")) %>% 
  filter(subjectID.sample1 == subjectID.sample2) %>% 
  filter(age_months.sample1 < age_months.sample2) %>% 
  group_by(sample1) %>% 
  filter(age_months.sample2 == min(age_months.sample2)) %>%
  ungroup() %>% 
  mutate(dT = age_months.sample2 - age_months.sample1,
         vel = dist / dT,
         mean_log_depth = (log10(depth.sample1) + log10(depth.sample2))/2,
         study  = study.sample1) %>%
  select(dist, dT, vel, study, mean_log_depth, subjectID = subjectID.sample1, age_months.sample2, delivery_mode = delivery_mode.sample1, sex = sex.sample1, diet = diet.sample2)

dist.bray.melt.subset.sex <- dist.bray.melt %>% filter(!is.na(sex))
dist.bray.melt.subset.delivery_mode <- dist.bray.melt %>% filter(!is.na(delivery_mode))
dist.bray.melt.subset.diet <- dist.bray.melt %>% filter(!is.na(diet))

model.velocity.votu.sex.gam <- gam(vel ~ s(age_months.sample2, by = sex) + mean_log_depth + sex +
                                s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                              family = tw(),
                              method = "REML",
                              data = dist.bray.melt.subset.sex)

model.velocity.votu.delivery_mode.gam <- gam(vel ~ s(age_months.sample2, by = delivery_mode) + mean_log_depth + delivery_mode +
                                          s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                                        family = tw(),
                                        method = "REML",
                                        data = dist.bray.melt.subset.delivery_mode)

model.velocity.votu.diet.gam <- gam(vel ~ s(age_months.sample2) + mean_log_depth + diet +
                                 s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                               family = tw(),
                               method = "REML",
                               data = dist.bray.melt.subset.diet)

emm.velocity.votu.sex <- emmeans(model.velocity.votu.sex.gam, 
                            pairwise ~ sex | age_months.sample2, 
                            at = list(age_months = seq(0, 24, by = 3)),
                            type = "response", data = dist.bray.melt.subset.sex)

emm.velocity.votu.delivery_mode <- emmeans(model.velocity.votu.delivery_mode.gam, 
                                      pairwise ~ delivery_mode | age_months.sample2, 
                                      at = list(age_months = seq(0, 24, by = 3)),
                                      type = "response", data = dist.bray.melt.subset.delivery_mode)

summary(emm.velocity.votu.sex$contrasts) %>% filter(p.value <= 0.05) # significant at 5 mo
summary(emm.velocity.votu.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # significant at 8 mo
summary(model.velocity.votu.diet.gam) # dietformula p-value = 0.001


velocity.votu.plot.data <- model.velocity.votu.sex.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  rbind.fill(model.velocity.votu.delivery_mode.gam %>% 
               smooth_estimates() %>% 
               add_confint()) %>% 
  rbind.fill(model.velocity.votu.diet.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(.by = "diet")) %>% 
  mutate(metric = "velocity.votu") %>%
  #filter(.smooth == "s(age_months.sample2)") %>% 
  filter(.type == "TPRS")


### BETA DISPERSION ###
dist.bray <- readRDS("rds/dist.bray.rds")
dist.bray.metadata <- phage_metadata %>% 
  select(run, study, depth, age_months) %>% 
  filter(run %in% (dist.bray %>% as.matrix() %>% colnames()))

dist.phf.wuni <- readRDS("rds/dist.phf.wuni.rds")
dist.phf.wuni.metadata <- phage_metadata %>% 
  select(run, study, depth, age_months) %>% 
  filter(run %in% (dist.phf.wuni %>% as.matrix() %>% colnames()))

dispersion_model.votu <- betadisper(dist.bray, group = dist.bray.metadata$age_months)
dispersion_model.phf <- betadisper(dist.phf.wuni, group = dist.phf.wuni.metadata$age_months)

dispersion_data.votu <- data.frame(distance_to_centroid = dispersion_model.votu$distances) %>% 
  rownames_to_column("run") %>% 
  left_join(phage_metadata, by = c("run" = "run"))

dispersion_data.phf <- data.frame(distance_to_centroid = dispersion_model.phf$distances) %>% 
  rownames_to_column("run") %>% 
  left_join(phage_metadata, by = c("run" = "run"))

dispersion_data.votu.subset.sex <- dispersion_data.votu %>% filter(!is.na(sex))
dispersion_data.votu.subset.delivery_mode <- dispersion_data.votu %>% filter(!is.na(delivery_mode))
dispersion_data.votu.subset.diet <- dispersion_data.votu %>% filter(!is.na(diet))

dispersion_data.phf.subset.sex <- dispersion_data.phf %>% filter(!is.na(sex))
dispersion_data.phf.subset.delivery_mode <- dispersion_data.phf %>% filter(!is.na(delivery_mode))
dispersion_data.phf.subset.diet <- dispersion_data.phf %>% filter(!is.na(diet))

model_dispersion_gamm.votu.sex <- gam(
  distance_to_centroid ~ s(age_months, by = sex) + sex + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.votu.subset.sex,
  method = "REML"
)

model_dispersion_gamm.votu.delivery_mode <- gam(
  distance_to_centroid ~ s(age_months, by = delivery_mode) + delivery_mode + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.votu.subset.delivery_mode,
  method = "REML"
)

model_dispersion_gamm.votu.diet <- gam(
  distance_to_centroid ~ s(age_months) + diet + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.votu.subset.diet,
  method = "REML"
)

model_dispersion_gamm.phf.sex <- gam(
  distance_to_centroid ~ s(age_months, by = sex) + sex + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.phf.subset.sex,
  method = "REML"
)

model_dispersion_gamm.phf.delivery_mode <- gam(
  distance_to_centroid ~ s(age_months, by = delivery_mode) + delivery_mode + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.phf.subset.delivery_mode,
  method = "REML"
)

model_dispersion_gamm.phf.diet <- gam(
  distance_to_centroid ~ s(age_months) + diet + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.phf.subset.diet,
  method = "REML"
)

dispersion.plot.data <- model_dispersion_gamm.phf.sex %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(metric = "dispersion.phf") %>%
  rbind.fill(model_dispersion_gamm.phf.delivery_mode %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "dispersion.phf")) %>% 
  rbind.fill(model_dispersion_gamm.phf.diet %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "dispersion.phf") %>%
               mutate(.by = "diet")) %>% 
  rbind.fill(model_dispersion_gamm.votu.sex %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "dispersion.votu")) %>% 
  rbind.fill(model_dispersion_gamm.votu.delivery_mode %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "dispersion.votu")) %>% 
  rbind.fill(model_dispersion_gamm.votu.diet %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(metric = "dispersion.votu") %>%
               mutate(.by = "diet")) %>%
  filter(.type == "TPRS")

emm.dispersion.votu.sex <- emmeans(model_dispersion_gamm.votu, 
                          pairwise ~ sex | age_months, 
                          at = list(age_months = seq(0, 24, by = 3)),
                          type = "response", data = dispersion_data.votu.subset)

emm.dispersion.votu.delivery_mode <- emmeans(model_dispersion_gamm.votu, 
                                    pairwise ~ delivery_mode | age_months, 
                                    at = list(age_months = seq(0, 24, by = 3)),
                                    type = "response", data = dispersion_data.votu.subset)

emm.dispersion.phf.sex <- emmeans(model_dispersion_gamm.phf, 
                            pairwise ~ sex | age_months, 
                            at = list(age_months = seq(0, 24, by = 3)),
                            type = "response", data = dispersion_data.phf.subset)

emm.dispersion.phf.delivery_mode <- emmeans(model_dispersion_gamm.phf, 
                                      pairwise ~ delivery_mode | age_months, 
                                      at = list(age_months = seq(0, 24, by = 3)),
                                      type = "response", data = dispersion_data.phf.subset)


summary(emm.dispersion.votu.sex$contrasts) %>% filter(p.value <= 0.05) # no significant differences
summary(emm.dispersion.votu.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # no significant differences

summary(emm.dispersion.phf.sex$contrasts) %>% filter(p.value <= 0.05) # significant at 6, 9, 12, 15, 18, 21 mo
summary(emm.dispersion.phf.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # significant at 9 mo


### TEMPERATE PHAGES ###
bacphlip_out <- read_tsv("bacphlip.tsv") %>% 
  mutate(lifecycle = ifelse(Temperate > 0.5, "temperate", "virulent")) %>% 
  select(virus, lifecycle)

ps.phage.relabund <- readRDS("rds/ps.phage.relabund.rds")

otu_table.melt <- otu_table(ps.phage.relabund) %>% melt() #%>% filter(value > 0)
colnames(otu_table.melt) <- c("virus", "run", "relabund")

otu_table.melt.lifecycle <- otu_table.melt %>% 
  left_join(phage_metadata %>% select(run, age_months, study, delivery_mode, sex, depth, subjectID, diet)) %>% 
  left_join(bacphlip_out, by = c("virus" = "virus"))

temperate.relabund.summary <- otu_table.melt.lifecycle %>% 
  filter(lifecycle == "temperate") %>% 
  group_by(run, lifecycle) %>% 
  dplyr::summarize(total_relabund = sum(relabund), age_months, study, delivery_mode, sex, depth, subjectID, diet) %>% 
  distinct() %>% 
  ungroup()

temperate.relabund.summary.filt.sex <- temperate.relabund.summary %>% filter(!is.na(sex))
temperate.relabund.summary.filt.delivery_mode <- temperate.relabund.summary %>% filter(!is.na(delivery_mode))
temperate.relabund.summary.filt.diet <- temperate.relabund.summary %>% filter(!is.na(diet))

n.sex <- nrow(temperate.relabund.summary.filt.sex)
n.delivery_mode <- nrow(temperate.relabund.summary.filt.delivery_mode)
n.diet <- nrow(temperate.relabund.summary.filt.diet)
temperate.relabund.summary.filt.sex.squeezed <- temperate.relabund.summary.filt.sex %>% mutate(relabund_squeezed = (total_relabund * (n.sex - 1) + 0.5) / n.sex)
temperate.relabund.summary.filt.delivery_mode.squeezed <- temperate.relabund.summary.filt.delivery_mode %>% mutate(relabund_squeezed = (total_relabund * (n.delivery_mode - 1) + 0.5) / n.delivery_mode)
temperate.relabund.summary.filt.diet.squeezed <- temperate.relabund.summary.filt.diet %>% mutate(relabund_squeezed = (total_relabund * (n.diet - 1) + 0.5) / n.diet)

model.temperate.sex.gam <- gam(relabund_squeezed ~ s(age_months, by = sex) + sex + log10(depth) +
                             s(age_months, study, bs = "fs", m = 1)  + s(subjectID, bs = "re"),
                           family = betar(),
                           method = "REML",
                           data = temperate.relabund.summary.filt.sex.squeezed)

model.temperate.delivery_mode.gam <- gam(relabund_squeezed ~ s(age_months, by = delivery_mode) + delivery_mode + log10(depth) +
                             s(age_months, study, bs = "fs", m = 1)  + s(subjectID, bs = "re"),
                           family = betar(),
                           method = "REML",
                           data = temperate.relabund.summary.filt.delivery_mode.squeezed)

model.temperate.diet.gam <- gam(relabund_squeezed ~ s(age_months) + diet + log10(depth) +
                             s(age_months, study, bs = "fs", m = 1)  + s(subjectID, bs = "re"),
                           family = betar(),
                           method = "REML",
                           data = temperate.relabund.summary.filt.diet.squeezed)

temperate.plot.data <- model.temperate.sex.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  rbind.fill(model.temperate.delivery_mode.gam %>% 
               smooth_estimates() %>% 
               add_confint()) %>% 
  rbind.fill(model.temperate.diet.gam %>% 
               smooth_estimates() %>% 
               add_confint() %>% 
               mutate(.by = "diet")) %>% 
  mutate(metric = "temperate.abund") %>% 
  filter(.type == "TPRS")

emm.temperate.sex <- emmeans(model.temperate.sex.gam, 
                                  pairwise ~ sex | age_months, 
                                  at = list(age_months = seq(0, 24, by = 3)),
                                  type = "response", data = temperate.relabund.summary.filt.sex.squeezed)

emm.temperate.delivery_mode <- emmeans(model.temperate.delivery_mode.gam, 
                                            pairwise ~ delivery_mode | age_months, 
                                            at = list(age_months = seq(0, 24, by = 3)),
                                            type = "response", data = temperate.relabund.summary.filt.delivery_mode.squeezed)

summary(emm.temperate.sex$contrasts) %>% filter(p.value <= 0.05) # no significant differences
summary(emm.temperate.delivery_mode$contrasts) %>% filter(p.value <= 0.05) # significant at 0 mo (CS > VG)
summary(model.temperate.diet.gam) # dietformula p-value = 0.003 (LOWER than milk)



# fig.s3.data <- readRDS("rds/fig.s3.data.rds")
fig.s3.data <- rbind.fill(alpha.diversity.plot.data, velocity.plot.data, velocity.votu.plot.data, temperate.plot.data, dispersion.plot.data) %>% 
  mutate(age_months = ifelse(metric == "velocity.phf", age_months.sample2, age_months)) %>%
  mutate(age_months = ifelse(metric == "velocity.votu", age_months.sample2, age_months)) %>% 
  mutate(.smooth = str_replace_all(.smooth, ".sample2", "")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):delivery_modeVG"), "Vaginal (delivery mode)")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):delivery_modeCS"), "Cesarean section (delivery mode)")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):sexF"), "Female (sex)")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):sexM"), "Male (sex)")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months)"), "Diet")) %>%
  mutate(.smooth = factor(.smooth, 
                          levels = c("Cesarean section (delivery mode)", "Vaginal (delivery mode)", "Male (sex)", "Female (sex)", "Diet"))) %>% 
  mutate(metric = factor(metric, 
                         levels = c("richness", "shannon", "pielou", "velocity.phf", "velocity.votu", "dispersion.phf", "dispersion.votu", "temperate.abund")))

# saveRDS(fig.s3.data, "rds/fig.s3.data.rds")
metric.labels <- c(sex = "Sex",
                   diet = "Diet",
                   delivery_mode = "Delivery mode",
                   richness = "Richness",
                   shannon = "Shannon",
                   pielou = "Pielou",
                   velocity.phf = "VDV (PHF)",
                   velocity.votu = "VDV (vOTU)",
                   dispersion.phf = "Beta disp.\n(PHF)",
                   dispersion.votu = "Beta disp.\n(vOTU)",
                   temperate.abund = "Temperate\nrel. abund.")

signif.months <- read_csv("fig-s3-gamms-significant.csv") %>% 
  mutate(.by = factor(.by)) %>% 
  mutate(metric = factor(metric))

# original model data
orig.model.richness.gam <- readRDS("rds/model.richness.gam.minimal.rds")
orig.model.shannon.gam <- readRDS("rds/model.shannon.gam.minimal.rds")
orig.model.pielou.gam <- readRDS("rds/model.pielou.gam.minimal.rds")

orig.model.velocity.gam.votu <- readRDS("rds/model.velocity.gam.votu.rds")
orig.model.velocity.gam.phf <-  readRDS("rds/model.velocity.gam.phf.rds")

orig.model.dispersion.gam.votu <- readRDS("rds/model_dispersion_gamm.bray.votu.rds")
orig.model.dispersion.gam.phf <- readRDS("rds/model_dispersion_gamm.wuni.phf.rds")

orig.model.temperate.gam <- readRDS("rds/model.temperate.gam.rds")


orig.model.smooths <- 
  orig.model.richness.gam %>% 
  smooth_estimates() %>% 
  #add_confint() %>% 
  mutate(metric = "richness") %>% 
  rbind.fill(orig.model.shannon.gam %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "shannon")) %>% 
  rbind.fill(orig.model.pielou.gam %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "pielou")) %>% 
  rbind.fill(orig.model.velocity.gam.phf %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "velocity.phf")) %>% 
  rbind.fill(orig.model.velocity.gam.votu %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "velocity.votu")) %>% 
  rbind.fill(orig.model.dispersion.gam.phf %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "dispersion.phf")) %>% 
  rbind.fill(orig.model.dispersion.gam.votu %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "dispersion.votu")) %>% 
  rbind.fill(orig.model.temperate.gam %>% 
               smooth_estimates() %>% 
               #add_confint() %>% 
               mutate(metric = "temperate.abund")) %>% 
  filter(.type == "TPRS") %>% 
  mutate(age_months = ifelse(metric == "velocity.phf", age_months.sample2, age_months)) %>%
  mutate(age_months = ifelse(metric == "velocity.votu", age_months.sample2, age_months)) %>% 
  mutate(model = "original.model") %>% 
  select(-.by) %>% 
  crossing(.by = c("delivery_mode", "diet", "sex")) %>% 
  mutate(metric = factor(metric, 
                         levels = c("richness", "shannon", "pielou", "velocity.phf", "velocity.votu", "dispersion.phf", "dispersion.votu", "temperate.abund")))


fig.s3.plot <- fig.s3.data %>% 
  ggplot(aes(x = age_months, y = .estimate, color = .smooth)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth), linewidth = 0.1, alpha = 0.1) +
  geom_rug(data = signif.months, 
           aes(x = age_months),
           linewidth = 0.75,
           sides = "b",
           color = "black", 
           length = unit(2, "mm"), 
           inherit.aes = F) +
  geom_line(data = orig.model.smooths, aes(x = age_months, y = .estimate), color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.8) +
  theme_bw() + 
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Age (months)",
       y = "Partial effect",
       fill = "Metadata variable",
       color = "Metadata variable") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), limits = c(0, 36)) +
  facet_grid(metric ~ .by, scales = "free_y",
             labeller = as_labeller(metric.labels)) +
  theme(legend.position = "bottom") +
  guides(col = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

fig.s3.plot

ggsave("~/Desktop/SuppFigure3A.pdf", plot = fig.s3.plot, width = 6.5, height = 9, limitsize = FALSE)


### DIET SIGNIFICANCE ###
emm.pielou.diet <- pairs(emmeans(model.pielou.diet.gam, 
              ~ diet, 
              type = "link", 
              data = alpha.diversity.votu.subset.diet)) %>% 
  summary() %>% mutate(metric = "Pielou")

emm.shannon.diet <- pairs(emmeans(model.shannon.diet.gam, 
                                 ~ diet, 
                                 type = "link", 
                                 data = alpha.diversity.votu.subset.diet)) %>% 
  summary() %>% mutate(metric = "Shannon")

emm.richness.diet <- pairs(emmeans(model.richness.diet.gam, 
                                 ~ diet, 
                                 type = "link", 
                                 data = alpha.diversity.votu.subset.diet)) %>% 
  summary() %>% mutate(metric = "Richness")

emm.velocity.diet <- pairs(emmeans(model.velocity.diet.gam, 
                                   ~ diet, 
                                   type = "link", 
                                   data = dist.phf.wuni.melt.subset.diet)) %>% 
  summary() %>% mutate(metric = "VDV (PHF)")


emm.velocity.votu.diet <- pairs(emmeans(model.velocity.votu.diet.gam, 
                                   ~ diet, 
                                   type = "link", 
                                   data = dist.bray.melt.subset.diet)) %>% 
  summary() %>% mutate(metric = "VDV (vOTU)")

emm.dispersion.phf.diet <- pairs(emmeans(model_dispersion_gamm.phf.diet, 
                                         ~ diet, 
                                         type = "link", 
                                         data = dispersion_data.phf.subset.diet)) %>% 
  summary() %>% mutate(metric = "Beta disp. (PHF)")


emm.dispersion.votu.diet <- pairs(emmeans(model_dispersion_gamm.votu.diet, 
                                         ~ diet, 
                                         type = "link", 
                                         data = dispersion_data.votu.subset.diet)) %>% 
  summary() %>% mutate(metric = "Beta disp. (vOTU)")


emm.temperate.diet <- pairs(emmeans(model.temperate.diet.gam, 
                                     ~ diet, 
                                     type = "link", 
                                     data = temperate.relabund.summary.filt.diet.squeezed)) %>% 
  summary() %>% mutate(metric = "Temperate rel. abund.")

emm.all_models.diet <- 
  rbind(emm.richness.diet, emm.shannon.diet, emm.pielou.diet,
      emm.velocity.diet, emm.velocity.votu.diet,
      emm.dispersion.phf.diet, emm.dispersion.votu.diet, emm.temperate.diet)

# saveRDS(emm.all_models.diet, "rds/emm.all_models.diet.rds")

emm.diet.estimate <- emm.all_models.diet %>%
  select(metric, contrast, estimate) %>%
  pivot_wider(names_from = contrast, values_from = estimate) %>%
  column_to_rownames("metric") %>%
  as.matrix()

emm.diet.pval <- emm.all_models.diet %>%
  select(metric, contrast, p.value) %>%
  pivot_wider(names_from = contrast, values_from = p.value) %>%
  column_to_rownames("metric") %>%
  as.matrix()

#col_fun <- colorRamp2(seq(-0.5, 0.5, length.out = 50), viridisLite::viridis(50, option = "viridis"))
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

library(ComplexHeatmap)
fig.s3.heatmap <- Heatmap(emm.diet.estimate, 
        col = col_fun,
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "white", lwd = 1),
        column_title = "Pairwise comparisons by diet",
        column_names_rot = 45,
        row_names_side = "left",
        row_title = "Metric",
        cell_fun = function(j, i, x, y, width, height, fill) {
          p = emm.diet.pval[i, j]
          if(p < 0.05) {
            star = ifelse(p < 0.01, "**", "*")
            grid.text(star, x, y, gp = gpar(fontsize = 12, fontface = "bold"))
          }
        },
        
        heatmap_legend_param = list(
          title = "Estimate",
          at = c(-0.5, 0, 0.5),
          labels = c("Decrease", "", "Increase")
        ))

fig.s3.heatmap

pdf("~/Desktop/SuppFigure3B.pdf", height = 4, width = 5)
draw(fig.s3.heatmap)
dev.off()

fig.s3a <- fig.s3.plot
fig.s3b <- ggdraw() + draw_image(image_read_pdf("~/Desktop/SuppFigure3B.pdf", density = 600))

fig.s3 <- 
  plot_grid(
    fig.s3a,
    plot_grid(
      fig.s3b,
      ncol = 1, 
      nrow = 2, 
      rel_heights = c(0.5, 1)
    ),
    labels = c("A", "B"),
    label_size = 22,
    ncol = 2,
    nrow = 1,
    rel_widths = c(1, 0.75)
  )

ggsave("~/Desktop/SuppFigure3.pdf", plot = fig.s3, width = 10, height = 9, limitsize = FALSE)




