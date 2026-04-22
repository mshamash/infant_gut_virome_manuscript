library(tidyverse)

phage_metadata <- readRDS("rds/metadata.combined.subset.rds") %>% 
  mutate(age_months = floor(age_days/30))

genomad.taxonomy <- read_tsv("votus-3kb_taxonomy.tsv")

genomad.taxonomy.filt <- genomad.taxonomy %>% 
  select(seq_name, lineage) %>% 
  filter(!str_detect(lineage, "Caudoviricetes|Microviridae|Inoviridae")) %>% 
  mutate(family = str_extract(lineage, "[^;]+viridae")) %>% 
  filter(!is.na(family))

genomad.taxonomy.filt$family %>% unique()
other.virus.family.counts <- as.data.frame(genomad.taxonomy.filt$family %>% 
                                             table())
colnames(other.virus.family.counts) <- c("family", "n")



ps.phage.relabund <- readRDS("rds/ps.phage.relabund.rds")
otu_table.melt <- otu_table(ps.phage.relabund) %>% melt()
colnames(otu_table.melt) <- c("virus", "run", "relabund")

otu_table.melt.filt <- 
  otu_table.melt %>% 
  filter(virus %in% genomad.taxonomy.filt$seq_name)

# Presence / absence heatmap
other.virus.presence_absence <- otu_table.melt.filt %>% 
  left_join(genomad.taxonomy.filt, by = c("virus" = "seq_name")) %>% 
  group_by(run, family) %>% 
  dplyr::summarize(fam.relabund = sum(relabund)) %>% 
  mutate(fam.present = ifelse(fam.relabund > 0, 1, 0)) %>% 
  select(run, family, fam.present) %>% 
  pivot_wider(names_from = family, values_from = fam.present, values_fill = 0) %>% 
  column_to_rownames("run") %>% 
  as.matrix() %>% 
  t()


phage_metadata.arranged <- phage_metadata %>% filter(run %in% colnames(other.virus.presence_absence)) %>% arrange(age_months)
samples.ordered <- phage_metadata.arranged$run

col_fun <- colorRamp2(c(0, 1), c("white", "black"))

library(ComplexHeatmap)
library(paletteer)

study.cols.hex <- paletteer_d("ggthemes::Classic_Green_Orange_12", n = 12)

all.months <- unique(phage_metadata.arranged$age_months)
selective_labels <- ifelse(all.months <= 12, 
                           as.character(all.months), 
                           ifelse(all.months %% 3 == 0, all.months, ""))

age.strip <- HeatmapAnnotation(
  Age = anno_block(
    gp = gpar(fill = "transparent", col = NA), 
    labels = selective_labels,
    labels_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  Study = phage_metadata.arranged$study,
  col = list(Study = setNames(as.character(study.cols.hex), sort(unique(phage_metadata.arranged$study)))),
  annotation_name_side = "left"
)

other.virus.presence_absence.ordered <- other.virus.presence_absence[, samples.ordered]


host_groups <- c(
  # https://ictv.global/virus-properties
  
  # Animal
  "Picornaviridae" = "Animal", "Astroviridae" = "Animal", "Anelloviridae" = "Animal", 
  "Caliciviridae" = "Animal", "Adenoviridae" = "Animal", "Parvoviridae" = "Animal",
  "Papillomaviridae" = "Animal", "Polyomaviridae" = "Animal", "Poxviridae" = "Animal", 
  "Herpesviridae" = "Animal", "Circoviridae" = "Animal", "Retroviridae" = "Animal", 
  "Adintoviridae" = "Animal", "Iridoviridae" = "Animal", "Sedoreoviridae" = "Animal",
  
  # Plants
  "Solemoviridae" = "Plant", "Geminiviridae" = "Plant", "Nanoviridae" = "Plant", 
  "Tombusviridae" = "Plant", "Virgaviridae" = "Plant", "Closteroviridae" = "Plant",
  
  # Fungi/Protist
  "Totiviridae" = "Fungi/Protist", "Endornaviridae" = "Fungi/Protist", 
  "Mimiviridae" = "Fungi/Protist", "Phycodnaviridae" = "Fungi/Protist",
  "Genomoviridae" = "Fungi/Protist",
  
  # Archaea
  "Fuselloviridae" = "Archaea", "Smacoviridae" = "Archaea", 
  "Bicaudaviridae" = "Archaea", "Sphaerolipoviridae" = "Archaea",
  
  # Virophage
  "Lavidaviridae" = "Virophage"
)

host_cols <- c(
  "Animal"        = "#7D9DC6FF",
  "Plant"         = "#59A55DFF",
  "Fungi/Protist" = "#EFDB56FF",
  "Archaea"       = "#CA4D2AFF",
  "Virophage"     = "#ECA23FFF"
)

row_grouping <- host_groups[rownames(other.virus.presence_absence.ordered)]

host_anno <- rowAnnotation(
  Host = row_grouping,
  col = list(Host = host_cols),
  show_annotation_name = TRUE
)

counts_vector <- other.virus.family.counts$n[match(rownames(other.virus.presence_absence.ordered), other.virus.family.counts$family)]

family.labels <- paste0(rownames(other.virus.presence_absence.ordered), " (n = ", counts_vector, ")")

other.virus.heatmap <- 
  Heatmap(other.virus.presence_absence.ordered,
        col = col_fun,
        show_heatmap_legend = F,
        name = "virome",
        
        # COLUMNS
        cluster_columns = F,
        show_column_names = F,
        column_split = phage_metadata.arranged$age_months, 
        column_title = "Age (months)",
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_gap = unit(1, "mm"),
        top_annotation = age.strip,
        
        # ROWS
        cluster_rows = T,
        cluster_row_slices = F,
        show_row_dend = F,
        show_row_names = T,
        row_labels = family.labels,
        row_split = row_grouping,
        row_title_rot = 90,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        row_gap = unit(2, "mm"),
        row_title = "Virus family",
        row_names_side = "left",
        left_annotation = host_anno
        )

ComplexHeatmap::draw(other.virus.heatmap)

pdf("~/Desktop/SuppFigure6A.pdf", width = 18, height = 7)
ComplexHeatmap::draw(other.virus.heatmap)
for(i in seq_along(unique(phage_metadata.arranged$age_months))) {
  decorate_annotation("Age", slice = i, {
    grid.lines(x = c(0.01, 0.99), y = c(0.05, 0.05), 
               gp = gpar(col = "black", lwd = 1))
  })
}
dev.off()


# Relative abundance by host
other.virus.abund <- 
  otu_table.melt.filt %>% 
  group_by(run) %>% 
  dplyr::summarize(other.virus.relabund = sum(relabund)) %>% 
  left_join(phage_metadata, by = c("run" = "run"))

summary(other.virus.abund$other.virus.relabund)

other.virus.abund %>% 
  ggplot(aes(x = age_months, y = other.virus.relabund)) +
  geom_point() +
  facet_wrap(~study)

other.virus.abund.byhost <- otu_table.melt.filt %>% 
  left_join(genomad.taxonomy.filt %>% select(seq_name, family),
          by = c("virus" = "seq_name")) %>% 
  left_join(as.data.frame(host_groups) %>% rownames_to_column("virus"), 
            by = c("family" = "virus")) %>% 
  left_join(phage_metadata %>% select(run, age_months, study, depth, subjectID, diet, delivery_mode, sex)) %>% 
  group_by(run, host_groups, age_months, study, depth, subjectID, diet, delivery_mode, sex) %>% 
  dplyr::summarize(host_groups.relabund = sum(relabund))

n.rows <- nrow(other.virus.abund.byhost)
other.virus.abund.byhost.squeezed <- other.virus.abund.byhost %>% 
  mutate(host_groups.relabund.squeezed = (host_groups.relabund * (n.rows - 1) + 0.5) / n.rows) %>% 
  mutate(host_groups = as.factor(host_groups))

model.other.virus.gam <- gam(host_groups.relabund.squeezed ~ 
                               s(age_months, by = host_groups) + host_groups + log10(depth) + 
                               s(age_months, study, bs = "fs", m = 1) + 
                               s(subjectID, bs = "re"),
                            family = betar(),
                            method = "REML",
                            data = other.virus.abund.byhost.squeezed)
summary(model.other.virus.gam)
plot(model.other.virus.gam, pages=1, se=TRUE)
#saveRDS(model.other.virus.gam, "rds/model.other.virus.gam.rds")
#model.other.virus.gam <- readRDS("rds/model.other.virus.gam.rds")

other.virus.abundance.plot <- 
  model.other.virus.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  filter(.type == "TPRS") %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):host_groupsAnimal"), "Animal")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):host_groupsPlant"), "Plant")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):host_groupsArchaea"), "Archaea")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):host_groupsFungi/Protist"), "Fungi/Protist")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):host_groupsVirophage"), "Virophage")) %>% 
  mutate(.smooth = factor(.smooth, 
                          levels = c("Animal", "Archaea", "Fungi/Protist", "Plant", "Virophage"))) %>% 
  ggplot(aes(x = age_months, y = .estimate, colour = .smooth)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth), linewidth = 0.1, alpha = 0.1) +
  scale_colour_manual(values = c("#7D9DC6FF", "#59A55DFF",  "#EFDB56FF", "#CA4D2AFF", "#ECA23FFF")) +
  scale_fill_manual(values = c("#7D9DC6FF", "#59A55DFF",  "#EFDB56FF", "#CA4D2AFF", "#ECA23FFF")) +
  theme_bw() + 
  labs(title = "Virus relative abundance (excl. bacteriophages)",
       x = "Age (months)",
       y = "Partial effect",
       colour = "Host",
       fill = "Host") +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), limits = c(0,38))

### PLANT VIRUS ABUNDANCE AND DIET ###
other.virus.abund.byhost.filt.diet <- other.virus.abund.byhost %>% 
  filter(!is.na(diet)) %>% 
  filter(host_groups == "Plant")

n.rows.diet <- nrow(other.virus.abund.byhost.filt.diet)
other.virus.abund.byhost.filt.diet.squeezed <- other.virus.abund.byhost.filt.diet %>% 
  mutate(host_groups.relabund.squeezed = (host_groups.relabund * (n.rows.diet - 1) + 0.5) / n.rows.diet)

model.other.virus.gam.plant_virus.diet <- 
  gam(host_groups.relabund.squeezed ~ 
        s(age_months) + diet + log10(depth) + 
        s(age_months, study, bs = "fs", m = 1) + 
        s(subjectID, bs = "re"),
      family = betar(),
      method = "REML",
      data = other.virus.abund.byhost.filt.diet.squeezed)

emm.plant_virus.diet <- pairs(emmeans(model.other.virus.gam.plant_virus.diet, 
                                   ~ diet, 
                                   type = "link", 
                                   data = other.virus.abund.byhost.filt.diet.squeezed)) %>% 
  summary()
#no significant difference in plant (or animal) virus abundance by diet



### ANIMAL VIRUS ABUNDANCE AND DELIVERY MDOE ###

# Delivery mode
other.virus.abund.byhost.filt.delivery_mode <- other.virus.abund.byhost %>% 
  filter(!is.na(delivery_mode)) %>% 
  filter(host_groups == "Animal")
n.rows.delivery_mode <- nrow(other.virus.abund.byhost.filt.delivery_mode)
other.virus.abund.byhost.filt.delivery_mode.squeezed <- other.virus.abund.byhost.filt.delivery_mode %>% 
  mutate(host_groups.relabund.squeezed = (host_groups.relabund * (n.rows.delivery_mode - 1) + 0.5) / n.rows.delivery_mode)

model.other.virus.gam.animal_virus.delivery_mode <- 
  gam(host_groups.relabund.squeezed ~ 
        s(age_months, by = delivery_mode) + delivery_mode + log10(depth) + 
        s(age_months, study, bs = "fs", m = 1) + 
        s(subjectID, bs = "re"),
      family = betar(),
      method = "REML",
      data = other.virus.abund.byhost.filt.delivery_mode.squeezed)
emm.animal_virus.delivery_mode <- emmeans(model.other.virus.gam.animal_virus.delivery_mode, 
                                       pairwise ~ delivery_mode | age_months, 
                                       at = list(age_months = seq(0, 24, by = 3)),
                                       type = "response", data = other.virus.abund.byhost.filt.delivery_mode.squeezed)
summary(emm.animal_virus.delivery_mode$contrasts) %>% filter(p.value <= 0.05)
# cs > vg (p < 0.05) at 6, 9, 12 mo

# Sex
other.virus.abund.byhost.filt.sex <- other.virus.abund.byhost %>% 
  filter(!is.na(sex)) %>% 
  filter(host_groups == "Animal")
n.rows.sex <- nrow(other.virus.abund.byhost.filt.sex)
other.virus.abund.byhost.filt.sex.squeezed <- other.virus.abund.byhost.filt.sex %>% 
  mutate(host_groups.relabund.squeezed = (host_groups.relabund * (n.rows.sex - 1) + 0.5) / n.rows.sex)

model.other.virus.gam.animal_virus.sex <- 
  gam(host_groups.relabund.squeezed ~ 
        s(age_months, by = sex) + sex + log10(depth) + 
        s(age_months, study, bs = "fs", m = 1) + 
        s(subjectID, bs = "re"),
      family = betar(),
      method = "REML",
      data = other.virus.abund.byhost.filt.sex.squeezed)
emm.animal_virus.sex <- emmeans(model.other.virus.gam.animal_virus.sex, 
                                pairwise ~ sex | age_months, 
                                at = list(age_months = seq(0, 24, by = 3)),
                                type = "response", data = other.virus.abund.byhost.filt.sex.squeezed)
summary(emm.animal_virus.sex$contrasts) %>% filter(p.value <= 0.05)
# no significant differences by sex


# Diet
other.virus.abund.byhost.filt.diet <- other.virus.abund.byhost %>% 
  filter(!is.na(diet)) %>% 
  filter(host_groups == "Animal")
n.rows.diet <- nrow(other.virus.abund.byhost.filt.diet)
other.virus.abund.byhost.filt.diet.squeezed <- other.virus.abund.byhost.filt.diet %>% 
  mutate(host_groups.relabund.squeezed = (host_groups.relabund * (n.rows.diet - 1) + 0.5) / n.rows.diet)

model.other.virus.gam.animal_virus.diet <- 
  gam(host_groups.relabund.squeezed ~ 
        s(age_months) + diet + log10(depth) + 
        s(age_months, study, bs = "fs", m = 1) + 
        s(subjectID, bs = "re"),
      family = betar(),
      method = "REML",
      data = other.virus.abund.byhost.filt.diet.squeezed)

emm.animal_virus.diet <- pairs(emmeans(model.other.virus.gam.animal_virus.diet, 
                                      ~ diet, 
                                      type = "link", 
                                      data = other.virus.abund.byhost.filt.diet.squeezed)) %>% 
  summary()
# no signficiant differences by diet

### Plot animal virus gamm for delivery mode
animal_virus.delivery_mode.plot <- 
  model.other.virus.gam.animal_virus.delivery_mode %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  filter(.type == "TPRS") %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):delivery_modeCS"), "Cesarean\nsection")) %>% 
  mutate(.smooth = str_replace_all(.smooth, fixed("s(age_months):delivery_modeVG"), "Vaginal")) %>% 
  ggplot(aes(x = age_months, y = .estimate, color = .smooth)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .smooth), linewidth = 0.1, alpha = 0.1) +
  geom_rug(data = data.frame(age_months = c(6, 9, 12)), 
           aes(x = age_months),
           linewidth = 0.75,
           sides = "b",
           color = "black", 
           length = unit(2, "mm"), 
           inherit.aes = F) +
  theme_bw() + 
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Animal virus relative abundance by delivery mode",
       x = "Age (months)",
       y = "Partial effect",
       fill = "Delivery mode",
       color = "Delivery mode") +
  scale_x_continuous(breaks = seq(0, 36, by = 6), limits = c(0, 36)) +
  theme(legend.position = "right")

### NUMBER OF VIRUSES PER HOST CATEGORY ###
genomad.taxonomy.filt %>% 
  select(seq_name, family) %>% 
  left_join(as.data.frame(host_groups) %>% rownames_to_column("virus"), 
            by = c("family" = "virus")) %>% 
  group_by(host_groups) %>% 
  dplyr::summarize(n = n())



