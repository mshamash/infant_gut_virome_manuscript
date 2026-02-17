library(tidyverse)
library(reshape2)
library(phyloseq)
library(gratia)
library(vegan)

bacphlip_out <- read_tsv("bacphlip.tsv") %>% 
  mutate(lifecycle = ifelse(Temperate > 0.5, "temperate", "virulent")) %>% 
  select(virus, lifecycle)

phage_metadata <- read_tsv("metadata.tsv")
# ps.phage.relabund <- readRDS("rds/ps.phage.relabund.rds")

otu_table.melt <- otu_table(ps.phage.relabund) %>% melt() #%>% filter(value > 0)
colnames(otu_table.melt) <- c("virus", "run", "relabund")

otu_table.melt.lifecycle <- otu_table.melt %>% 
  left_join(phage_metadata %>% select(run, age_months, study, delivery_mode, sex, depth, subjectID)) %>% 
  left_join(bacphlip_out, by = c("virus" = "virus"))

temperate.relabund.summary <- otu_table.melt.lifecycle %>% 
  filter(lifecycle == "temperate") %>% 
  group_by(run, lifecycle) %>% 
  summarise(total_relabund = sum(relabund), age_months, study, delivery_mode, sex, depth, subjectID) %>% 
  distinct() %>% 
  ungroup()


n <- nrow(temperate.relabund.summary)
temperate.relabund.summary.squeezed <- temperate.relabund.summary %>% mutate(relabund_squeezed = (total_relabund * (n - 1) + 0.5) / n)

model.temperate.gam <- gam(relabund_squeezed ~ s(age_months) + log10(depth) +
                             s(age_months, study, bs = "fs", m = 1)  + s(subjectID, bs = "re"),
                          family = betar(),
                          method = "REML",
                          data = temperate.relabund.summary.squeezed)
# saveRDS(model.temperate.gam, "rds/model.temperate.gam.minimal.rds")
# model.temperate.gam <- readRDS("rds/model.temperate.gam.rds")

summary(model.temperate.gam)
gam.check(model.temperate.gam)

temperate.plot <- model.temperate.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(metric = "temperate.abund") %>% 
  filter(.smooth == "s(age_months)") %>% 
  ggplot(aes(x = age_months, y = .estimate, colour = metric)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = metric), linewidth = 0.1, alpha = 0.1) +
  scale_colour_manual(values = c("#FE88B1FF")) +
  scale_fill_manual(values = c("#FE88B1FF")) +
  theme_bw() + 
  labs(title = "Temperate phage relative abundance",
       x = "Age (months)",
       y = "Partial effect") +
  theme(aspect.ratio = 0.5,
        legend.title = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 36, by = 6))

temperate.plot


# phf overview of temperate phages
# iphop_out.split <- readRDS("rds/iphop_out.split.rds")

temperate.phf.overview <- otu_table.melt.lifecycle %>% 
  filter(relabund > 0) %>% 
  select(virus, age_months) %>% 
   mutate(
     age_bin = case_when(
       age_months == 0 ~ "0 mo.",
       age_months %in% 1:3 ~ "1-3 mo.",
       age_months %in% 4:6 ~ "4-6 mo.",
       age_months %in% 7:12 ~ "7-12 mo.",
       TRUE ~ ">1 yr."
     )
   ) %>% 
  distinct() %>% 
  left_join(iphop_out.split %>% select(Virus, family), by = c("virus" = "Virus")) %>% 
  filter(!is.na(family))

temperate.votu.overview <- temperate.phf.overview %>% 
  group_by(age_bin) %>% 
  summarize(n_temperate_votus = n())

temperate.phf.rel.props <- temperate.phf.overview %>% 
  group_by(age_bin, family) %>% 
  summarize(n_tempetate_vots_per_phf = n()) %>% 
  left_join(temperate.votu.overview) %>% 
  mutate(temperate.phf.rel.prop = n_tempetate_vots_per_phf / n_temperate_votus)


temperate.phf.plot <- temperate.phf.rel.props %>% 
  mutate(age_bin = factor(age_bin, levels = c("0 mo.", "1-3 mo.", "4-6 mo.", "7-12 mo.", ">1 yr."))) %>% 
  ungroup() %>% 
  mutate(family_lumped = fct_lump_n(family, n = 11, w = temperate.phf.rel.prop, other_level = "Others")) %>%
  group_by(age_bin, family_lumped) %>% 
  summarize(temperate.phf.rel.prop = sum(temperate.phf.rel.prop)) %>% 
  ggplot(aes(x = age_bin, y = temperate.phf.rel.prop, fill = family_lumped, group = family_lumped)) + 
    geom_area(colour = "black", linewidth = 0.1) +
    #geom_bar(stat = "identity") +
    theme_bw() +
    #theme(legend.position = "none") +
    #scale_fill_brewer(palette = "Paired", name = "PHF") +
    paletteer::scale_fill_paletteer_d("rcartocolor::Pastel", name = "PHF") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.001)) +
    ylab("Proportion of detected temperate vOTUs") +
    xlab("Age bin") +
    theme(legend.text = element_text(face = "italic"),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


temperate.phf.plot
