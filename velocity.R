library(tidyverse)
library(phyloseq)
library(vegan)
library(gratia)
library(mgcv)
library(reshape2)

### VIROME DEVELOPMENTAL VELOCITY (PHF) ###
phage_metadata <- readRDS("rds/metadata.combined.subset.rds") %>% mutate(age_months = floor(age_days/30))
dist.phf.wuni <- readRDS("rds/dist.phf.wuni.rds")
# dist.phf.bray <- readRDS("rds/dist.phf.bray.rds")

dist.phf.wuni.melt <- dist.phf.wuni %>% 
  as.matrix() %>% 
  melt() %>% 
  `colnames<-`(c("sample1", "sample2", "dist")) %>% 
  filter(!(sample1 == sample2)) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, subjectID, depth, delivery_mode, sex), by = c("sample1" = "run")) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, subjectID, depth, delivery_mode, sex), by = c("sample2" = "run"), suffix = c(".sample1", ".sample2")) %>% 
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
  select(dist, dT, vel, study, mean_log_depth, subjectID = subjectID.sample1, age_months.sample2, delivery_mode = delivery_mode.sample1, sex = sex.sample1)

model.velocity.gam <- gam(vel ~ s(age_months.sample2) + mean_log_depth + 
                               s(age_months.sample2, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                             family = tw(),
                             method = "REML",
                             data = dist.phf.wuni.melt)

summary(model.velocity.gam)
gam.check(model.velocity.gam)

velocity.plot <- model.velocity.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(metric = "velocity") %>% 
  filter(.smooth == "s(age_months.sample2)") %>% 
  ggplot(aes(x = age_months.sample2, y = .estimate, colour = metric)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = metric), linewidth = 0.1, alpha = 0.1) +
  scale_colour_manual(values = c("#F89C74FF")) +
  scale_fill_manual(values = c("#F89C74FF")) +
  theme_bw() + 
  labs(title = "Virome developmental velocity",
       x = "Age (months)",
       y = "Partial effect") +
  theme(aspect.ratio = 0.5,
        legend.title = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 36, by = 6))

velocity.plot


### BETA DISPERSION ###
library(vegan)
dist.bray <- distance(ps.phage.relabund, method = "bray")
dist.bray.metadata <- phage_metadata %>% 
  select(run, study, depth, age_months) %>% 
  filter(run %in% (dist.bray %>% as.matrix() %>% colnames()))

# dist.phf.wuni <- readRDS("rds/dist.phf.wuni.rds")
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

model_dispersion_gamm.votu <- gam(
  distance_to_centroid ~ s(age_months) + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.votu,
  method = "REML"
)

model_dispersion_gamm.phf <- gam(
  distance_to_centroid ~ s(age_months) + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
  family = tw(),
  data = dispersion_data.phf,
  method = "REML"
)

# saveRDS(model_dispersion_gamm, "rds/model_dispersion_gamm.bray.votu.rds")
# model_dispersion_gamm.votu <- readRDS("rds/model_dispersion_gamm.bray.votu.rds")
# model_dispersion_gamm.phf <- readRDS("rds/model_dispersion_gamm.wuni.phf.rds")

summary(model_dispersion_gamm.votu)
gam.check(model_dispersion_gamm.votu)

summary(model_dispersion_gamm.phf)
gam.check(model_dispersion_gamm.phf)

dispersion.plot <- 
  model_dispersion_gamm.phf %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(dataset = "phf") %>% 
  rbind(model_dispersion_gamm.votu %>% 
          smooth_estimates() %>% add_confint() %>% 
          filter(.smooth == "s(age_months)") %>% 
          mutate(dataset = "votu")) %>% 
  filter(.smooth == "s(age_months)") %>% 
  ggplot(aes(x = age_months, y = .estimate, color = dataset)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = dataset), linewidth = 0.1, alpha = 0.1) +
  scale_colour_manual(values = c("#87C55FFF", "#66C5CCFF"), labels = c("PHF\n(p=0.003)","vOTU\n(p=0.23)")) +
  scale_fill_manual(values = c("#87C55FFF", "#66C5CCFF"), labels = c("PHF\n(p=0.003)","vOTU\n(p=0.23)")) +
  theme_bw() + 
  labs(title = "Beta dispersion",
       x = "Age (months)",
       y = "Partial effect") +
  theme(aspect.ratio = 0.5,
        legend.title = element_blank(),
        legend.key.spacing.y = unit(0.15, 'cm'),
        legend.position = "inside",
        legend.position.inside = c(0.22, 0.12) ,
        legend.box.background = element_rect(color = "black")
        ) +
  guides(colour = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = seq(0, 36, by = 6)) 
  
dispersion.plot
