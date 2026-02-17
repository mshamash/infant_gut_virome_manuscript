library(tidyverse)
library(phyloseq)

#### PREVALENCE OF PHFs ####
# ps.phage.relabund.phf <- readRDS("rds/ps.phage.relabund.phf.rds")

ps_melt.phf <- psmelt(ps.phage.relabund.phf) %>% 
  mutate(
    age_months = floor(age_days/30),
    
    age_bin = case_when(
      age_months == 0 ~ "0 mo.",
      age_months %in% 1:3 ~ "1-3 mo.",
      age_months %in% 4:6 ~ "4-6 mo.",
      age_months %in% 7:12 ~ "7-12 mo.",
      TRUE ~ ">1 yr."
    )
  )

prevalence_data.phf <- ps_melt.phf %>%
  group_by(age_bin, subjectID, OTU) %>%
  summarise(is_present_in_infant = any(Abundance > 0),
            .groups = "drop") %>% 
  group_by(age_bin, OTU) %>%
  summarise(n_infants_in_bin = n(),
            n_positives = sum(is_present_in_infant),
            prevalence = n_positives / n_infants_in_bin,
            .groups = "drop")

# iphop_out.split <- readRDS("rds/iphop_out.split.rds")

top_phages.phf <- prevalence_data.phf %>%
  group_by(OTU) %>%
  summarise(mean_prev = mean(prevalence)) %>%
  slice_max(mean_prev, n = 10) %>%
  left_join(iphop_out.split %>% select(Virus, family, `Host genome`), by = c("OTU" = "Host genome")) %>% 
  select(-Virus) %>% 
  unique()
  
top.prevalent.phf.plot <- prevalence_data.phf %>%
  filter(OTU %in% top_phages.phf$OTU) %>%
  left_join(top_phages.phf) %>%
  mutate(age_bin = factor(age_bin, levels = c("0 mo.", "1-3 mo.", "4-6 mo.", "7-12 mo.", ">1 yr."))) %>% 
  ggplot(aes(x = age_bin, y = prevalence, color = family, group = family)) +
    geom_line(size = 0.7) +
    geom_point(size = 1) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "10 most prevalent PHFs",
         y = "Prevalence\n(fraction of infants)",
         x = "Age bin",
         color = "PHF") +
    theme_bw() +
    theme(legend.position = "right",# axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 11),
          plot.background = element_rect(colour = "black", fill="white", linewidth=0.5),
          legend.text = element_text(face = "italic")) +
    paletteer::scale_colour_paletteer_d("ggthemes::Classic_10_Medium")

prevalence.plot.data.phf <- prevalence_data.phf %>%
  left_join(iphop_out.split %>% select(Virus, family, `Host genome`), by = c("OTU" = "Host genome")) %>% 
  select(-Virus) %>% 
  unique() %>% 
  group_by(age_bin) %>% 
  arrange(desc(prevalence), .by_group = TRUE) %>% 
  mutate(rank = row_number()) %>%
  ungroup() %>% 
  mutate(age_bin = factor(age_bin, levels = c("0 mo.", "1-3 mo.", "4-6 mo.", "7-12 mo.", ">1 yr.")))

prevalence.plot.cutoffs.phf <- 
  prevalence.plot.data.phf %>% 
  group_by(age_bin) %>% 
  summarise(p50_rank = max(rank[prevalence >= 0.5], -Inf)) %>% 
  filter(p50_rank > 0)

prevalence.rank.phf.plot <- 
  prevalence.plot.data.phf %>% 
  ggplot(aes(x = rank, y = prevalence, color = age_bin)) +
    geom_point(aes(size = prevalence > 0.5)) +
    scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 0.5), guide = "none") +
    geom_line() +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    paletteer::scale_colour_paletteer_d("LaCroixColoR::Apricot") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "PHF prevalence",
         y = "Prevalence\n(fraction of infants)",
         x = "Rank",
         color = "Age bin")

top.prevalent.phf.plot
prevalence.rank.phf.plot


#### PREVALENCE OF vOTUs ####
# ps.phage <- readRDS("rds/ps.phage.rds")

ps_melt <- psmelt(ps.phage) %>% 
  mutate(
    age_months = floor(age_days/30),
    
    age_bin = case_when(
      age_months == 0 ~ "0 mo.",
      age_months %in% 1:3 ~ "1-3 mo.",
      age_months %in% 4:6 ~ "4-6 mo.",
      age_months %in% 7:12 ~ "7-12 mo.",
      TRUE ~ ">1 yr."
    )
  )

prevalence_data <- ps_melt %>%
  group_by(age_bin, subjectID, OTU) %>%
  summarise(is_present_in_infant = any(Abundance > 0),
            .groups = "drop") %>% 
  group_by(age_bin, OTU) %>%
  summarise(n_infants_in_bin = n(),
            n_positives = sum(is_present_in_infant),
            prevalence = n_positives / n_infants_in_bin,
            .groups = "drop")

top_phages <- prevalence_data %>%
  group_by(OTU) %>%
  summarise(mean_prev = mean(prevalence)) %>%
  slice_max(mean_prev, n = 10)

top.prevalent.plot <- 
  prevalence_data %>%
  filter(OTU %in% top_phages$OTU) %>%
  left_join(top_phages) %>%
  mutate(age_bin = factor(age_bin, levels = c("0 mo.", "1-3 mo.", "4-6 mo.", "7-12 mo.", ">1 yr."))) %>% 
  ggplot(aes(x = age_bin, y = prevalence, color = OTU, group = OTU)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  labs(title = "10 most prevalent vOTUs",
       y = "Prevalence\n(fraction of infants)",
       x = "Age bin",
       color = "vOTU") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 11),
        plot.background = element_rect(colour = "black", fill="white", linewidth = 0.5),
        legend.text = element_text(face = "italic")) +
  paletteer::scale_colour_paletteer_d("ggthemes::Classic_10_Medium")

prevalence.plot.data <- 
  prevalence_data %>%
  group_by(age_bin) %>% 
  arrange(desc(prevalence), .by_group = TRUE) %>% 
  mutate(rank = row_number()) %>%
  ungroup() %>% 
  mutate(age_bin = factor(age_bin, levels = c("0 mo.", "1-3 mo.", "4-6 mo.", "7-12 mo.", ">1 yr.")))

prevalence.plot.cutoffs <- 
  prevalence.plot.data %>% 
  group_by(age_bin) %>% 
  summarise(p50_rank = max(rank[prevalence >= 0.5], -Inf)) %>% 
  filter(p50_rank > 0)

prevalence.rank.plot <- 
  prevalence.plot.data %>% 
  ggplot(aes(x = rank, y = prevalence, color = age_bin)) +
  geom_point(aes(size = prevalence > 0.5)) +
  scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 0.5), guide = "none") +
  geom_line() +
  geom_vline(data = prevalence.plot.cutoffs,
             aes(xintercept = p50_rank, color = age_bin),
             linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  paletteer::scale_colour_paletteer_d("LaCroixColoR::Apricot") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "vOTU prevalence",
       y = "Prevalence\n(fraction of infants)",
       x = "Rank",
       color = "Age bin")

top.prevalent.plot
prevalence.rank.plot
