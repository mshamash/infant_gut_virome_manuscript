library(tidyverse)
library(phyloseq)
library(reshape2)
library(vegan)
library(cowplot)
library(mgcv)
library(gratia)
library(flextable)

phage_cov <- read_tsv("phage-coverage-combined.filt.tsv") #1x depth coverage and 75% breadth coverage filter done in bash
# 1,869 out of 1,893 samples remaining (991 out of 999 infants)

phage_cov.summary <- phage_cov %>%
  group_by(sample, contig) %>%
  summarize(phage_cov = meandepth)

phage_cov.totals <- phage_cov.summary %>% 
  group_by(sample) %>% 
  summarize(sample_total_cov = sum(phage_cov)) # get total coverage per sample

phage_cov.joined <- inner_join(phage_cov.summary, phage_cov.totals, by = c("sample" = "sample"))

phage_cov.joined$relabund_cov <- (phage_cov.joined$phage_cov / phage_cov.joined$sample_total_cov)

phage_metadata <- read_tsv("metadata.tsv")
rownames(phage_metadata) <- phage_metadata$run

phage_cov.joined <- left_join(phage_cov.joined, phage_metadata, by = c("sample" = "run"))

# saveRDS(phage_cov.joined, "rds/phage_cov.joined.rds")
# phage_cov.joined <- readRDS("rds/phage_cov.joined.rds")

votu.table <- as.data.frame(pivot_wider(phage_cov.joined[, c("sample", "contig", "phage_cov")], id_cols = "contig", names_from = "sample", values_from = "phage_cov"))
rownames(votu.table) <- votu.table$contig
votu.table <- votu.table[, -1]
votu.table <- floor(votu.table)
votu.table[is.na(votu.table)] <- 0

PHAGE_OTU <- otu_table(votu.table, taxa_are_rows = T)

PHAGE_METADATA <- sample_data(phage_metadata)

ps.phage <- phyloseq(PHAGE_OTU, PHAGE_METADATA)
ps.phage.relabund <- transform_sample_counts(ps.phage, function(x) x / sum(x) )
# saveRDS(ps.phage, "rds/ps.phage.rds")
# saveRDS(ps.phage.relabund, "rds/ps.phage.relabund.rds")
# ps.phage <- readRDS("rds/ps.phage.rds")
# ps.phage.relabund <- readRDS("rds/ps.phage.relabund.rds")

# PHF / PHAGE-HOST FAMILY / IPHOP IMPORT

iphop_out <- read_csv("Host_prediction_to_genome_m90-combined.csv")

length(unique(iphop_out$Virus)) # 36,800 unique phage contigs with host assignment, out of 49,745 total detected contigs (74%)

iphop_out.dedup <- iphop_out %>% group_by(Virus) %>% arrange(desc(`Confidence score`)) %>% slice(1)

iphop_out.dedup <- iphop_out.dedup[, c("Virus", "Host taxonomy", "Host genome")]

iphop_out.split <- iphop_out.dedup %>%
  separate(col = `Host taxonomy`, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

iphop_out.split$domain <- gsub("d__", "", iphop_out.split$domain)
iphop_out.split$phylum <- gsub("p__", "", iphop_out.split$phylum)
iphop_out.split$class <- gsub("c__", "", iphop_out.split$class)
iphop_out.split$order <- gsub("o__", "", iphop_out.split$order)
iphop_out.split$family <- gsub("f__", "", iphop_out.split$family)
iphop_out.split$genus <- gsub("g__", "", iphop_out.split$genus)
iphop_out.split$species <- gsub("s__", "", iphop_out.split$species)

iphop_out.split$phylum <- sapply(strsplit(iphop_out.split$phylum, "_"), `[`, 1)
iphop_out.split$class <- sapply(strsplit(iphop_out.split$class, "_"), `[`, 1)
iphop_out.split$order <- sapply(strsplit(iphop_out.split$order, "_"), `[`, 1)
iphop_out.split$genus <- sapply(strsplit(iphop_out.split$genus, "_"), `[`, 1)
iphop_out.split$species <- sapply(strsplit(iphop_out.split$species, "_"), `[`, 1)

iphop_out.split$phylum <- gsub("Proteobacteria", "Pseudomonadota", iphop_out.split$phylum)
iphop_out.split$phylum <- gsub("Desulfobacterota", "Thermodesulfobacteriota", iphop_out.split$phylum)
iphop_out.split$phylum <- gsub("Actinobacteriota", "Actinomycetota", iphop_out.split$phylum)
iphop_out.split$phylum <- gsub("Firmicutes", "Bacillota", iphop_out.split$phylum)

# saveRDS(iphop_out.split, "rds/iphop_out.split.rds")
# iphop_out.split <- readRDS("rds/iphop_out.split.rds")

taxtable <- iphop_out.split %>% ungroup() %>% select(domain, phylum, class, order, family, genus, species, `Host genome`) %>% distinct() %>% as.data.frame()

#solve duplicate RS_GCF_001455345.1 issue
taxtable <- taxtable %>% group_by(`Host genome`) %>% arrange(desc(`Host genome`)) %>% slice(1) %>% as.data.frame()

rownames(taxtable) <- taxtable$`Host genome`
taxtable <- taxtable %>% select(domain, phylum, class, order, family, genus, species, `Host genome`)
# saveRDS(taxtable, "rds/taxtable.rds")
# taxtable <- readRDS("rds/taxtable.rds")

TAXONOMY <- tax_table(as.matrix(taxtable))

TREE.gtdb.r202 <- read_tree("bac120_r202.tree")

ps.phage.relabund.melted <- ps.phage.relabund %>% psmelt()

ps.phage.relabund.melted <- ps.phage.relabund.melted %>% left_join(iphop_out.split, by = c("OTU" = "Virus"))

ps.phage.relabund.melted.filt <- ps.phage.relabund.melted %>%
  ungroup() %>%
  group_by(Sample, `Host genome`) %>%
  summarize(relabund = sum(Abundance)) %>%
  filter(!is.na(`Host genome`)) %>%
  filter(relabund > 0) %>%
  select(c(Sample, relabund, `Host genome`)) %>%
  pivot_wider(id_cols = "Host genome", names_from = "Sample", values_from = "relabund") %>%
  as.data.frame()

rownames(ps.phage.relabund.melted.filt) <- ps.phage.relabund.melted.filt$`Host genome`
ps.phage.relabund.melted.filt <- ps.phage.relabund.melted.filt[, -1]
ps.phage.relabund.melted.filt[is.na(ps.phage.relabund.melted.filt)] <- 0
saveRDS(ps.phage.relabund.melted.filt, "rds/ps.phage.relabund.melted.filt.rds")
ps.phage.relabund.melted.filt <- readRDS("rds/ps.phage.relabund.melted.filt.rds")




PHAGE_OTU.genome <- otu_table(ps.phage.relabund.melted.filt, taxa_are_rows = T)
ps.phage.relabund.hosts <- phyloseq(PHAGE_OTU.genome, PHAGE_METADATA, TREE.gtdb.r202, TAXONOMY)
ps.phage.relabund.phf <- ps.phage.relabund.hosts %>% 
  tax_glom(taxrank = "family")
# saveRDS(ps.phage.relabund.phf, "rds/ps.phage.relabund.phf.rds")
# ps.phage.relabund.phf <- readRDS("rds/ps.phage.relabund.phf.rds")

### BETA DIVERSITY PHF ###

ord.phf.pcoa.wuni <- ordinate(ps.phage.relabund.phf, method="PCoA", distance="wunifrac")
# saveRDS(ord.phf.pcoa.wuni, "rds/ord.phf.pcoa.wuni.rds")
# ord.phf.pcoa.wuni <- readRDS("rds/ord.phf.pcoa.wuni.rds")

ord.phf.pcoa.axis1.variance <- round(ord.phf.pcoa.wuni$values$Relative_eig[1]*100, 1)
ord.phf.pcoa.axis2.variance <- round(ord.phf.pcoa.wuni$values$Relative_eig[2]*100, 1)

pcoa.phf.wuni.age <- 
  ord.phf.pcoa.wuni$vectors %>% 
  as.data.frame() %>% 
  rownames_to_column("run") %>% 
  left_join(phage_metadata) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = age_months)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low = "lightblue", high = "navyblue", breaks = seq(0, 38, by=6), limits = c(-1, 39)) +
  xlab(paste("Axis.1 [", ord.phf.pcoa.axis1.variance, "%]", sep = "")) +
  ylab(paste("Axis.2 [", ord.phf.pcoa.axis2.variance, "%]", sep = "")) +
  labs(color = "Age (months)") +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pcoa.phf.wuni.age

pcoa.phf.wuni.study <- 
  ord.phf.pcoa.wuni$vectors %>% 
  as.data.frame() %>% 
  rownames_to_column("run") %>% 
  left_join(phage_metadata) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = study)) +
  geom_point() +
  theme_bw() +
  paletteer::scale_colour_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  xlab(paste("Axis.1 [", ord.phf.pcoa.axis1.variance, "%]", sep = "")) +
  ylab(paste("Axis.2 [", ord.phf.pcoa.axis2.variance, "%]", sep = "")) +
  labs(color = "Study") +
  guides(colour = guide_legend(nrow = 3)) +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pcoa.phf.wuni.study


title.votu <- ggdraw() + draw_label("vOTU", fontface = 'bold') + theme(plot.margin = margin(0, 0, 0, 40))
title.phf <- ggdraw() + draw_label("PHF", fontface = 'bold') + theme(plot.margin = margin(0, 0, 0, 40))

legend_pcoa.phf.wuni.age <- get_legend(pcoa.phf.wuni.age)
legend_pcoa.phf.wuni.study <- get_legend(pcoa.phf.wuni.study)

rm_legend <- function(p){p + theme(legend.position = "none")}

plot_grid(
  plot_grid(
    title.votu,
    title.phf,
    rm_legend(pcoa.votu.bray.age), rm_legend(pcoa.phf.wuni.age),
    rm_legend(pcoa.votu.bray.study), rm_legend(pcoa.phf.wuni.study),
    nrow = 3, 
    ncol = 2,
    rel_heights = c(0.1, 1, 1),
    labels = c("", "", "A", "B", "C", "D")
    ),
  plot_grid(
    legend_pcoa.votu.bray.age, legend_pcoa.votu.bray.study, 
    ncol = 2, 
    nrow = 1, 
    rel_widths = c(0.5,1)
    ),
  ncol = 1,
  rel_heights = c(1, 0.22)
)

ordination.plots <- plot_grid(
  plot_grid(
    rm_legend(pcoa.phf.wuni.age),
    rm_legend(pcoa.phf.wuni.study),
    nrow = 1, 
    ncol = 2
  ),
  plot_grid(
    legend_pcoa.phf.wuni.age, legend_pcoa.phf.wuni.study, 
    ncol = 2, 
    nrow = 1, 
    rel_widths = c(0.5,1)
  ),
  ncol = 1,
  rel_heights = c(1, 0.22)
)

dist.phf.wuni <- distance(ps.phage.relabund.phf, method = "wunifrac")
# saveRDS(dist.phf.wuni, "rds/dist.phf.wuni.rds")
# dist.phf.wuni <- readRDS("rds/dist.phf.wuni.rds")
dist.phf.wuni.metadata <- phage_metadata %>% filter(run %in% (dist.phf.wuni %>% as.matrix() %>% colnames()))
set.seed(123)
permanova.phf.wuni <- adonis2(dist.phf.wuni ~ log10(depth) + study + age_months,
                          data = dist.phf.wuni.metadata,
                          by = "terms",
                          permutations = 999,
                          strata = dist.phf.wuni.metadata$subjectID)
# saveRDS(permanova.phf.wuni, "rds/permanova.phf.wuni.minimal.rds")
# permanova.phf.wuni <- readRDS("rds/permanova.phf.wuni.rds")
permanova.phf.wuni

set.seed(123)
permanova.phf.wuni.infant <- adonis2(dist.phf.wuni ~ log10(depth) + study + subjectID, 
                              data = dist.phf.wuni.metadata, 
                              by = "terms",
                              permutations = 999)
permanova.phf.wuni.infant
# saveRDS(permanova.phf.wuni.infant, "rds/permanova.phf.wuni.infant.rds")


permanova.table <- data.frame(
  `Effect` = c("log10(depth)", "study", "age_months", "subjectID"),
  `R2` = c("0.009", "0.136", "0.001", "0.488"),
  `P.value` = c("0.988", "0.011 *", "0.032 *", "-"))

permanova.flextable <- flextable(permanova.table) %>% 
  add_header_row(colwidths = c(3),
                 values = c("Weighted UniFrac distance ~ log10(depth) + study + age_months")) %>% 
  theme_vanilla() %>% 
  add_footer_lines("PERMANOVA (adonis2) with 999 permutations") %>% 
  autofit(add_w = 0.8) %>% 
  flextable::compose(i = 2, j = 2, part = "header", value = as_paragraph("R", as_sup("2"))) %>% 
  flextable::compose(i = 2, j = 3, part = "header", value = as_paragraph("P-value"))

permanova.flextable
permanova.flextable.grob <- gen_grob(permanova.flextable, fit = "fixed", just = "center")


### ALPHA DIVERSITY VOTU ###
alpha.diversity.votu <- estimate_richness(ps.phage, measures=c("Shannon", "Observed")) %>% 
  mutate(Pielou = Shannon / log(Observed))
alpha.diversity.votu$run <- rownames(alpha.diversity.votu)
alpha.diversity.votu <- alpha.diversity.votu %>% 
  left_join(phage_metadata %>% select(run, depth, age_days, study, delivery_mode, subjectID, sex)) %>% 
  mutate(age_months = floor(age_days/30))

# pielou
n <- nrow(alpha.diversity.votu)
alpha.diversity.votu.squeezed <- alpha.diversity.votu %>% mutate(Pielou_squeezed = (Pielou * (n - 1) + 0.5) / n)
model.pielou.gam <- gam(Pielou_squeezed ~ s(age_months) + log10(depth) +
                             s(age_months, study, bs = "fs", m = 1)  + s(subjectID, bs = "re"),
                          family = betar(),
                          method = "REML",
                          data = alpha.diversity.votu.squeezed)
# saveRDS(model.pielou.gam, "rds/model.pielou.gam.minimal.rds")
# model.pielou.gam <- readRDS("rds/model.pielou.gam.minimal.rds")

summary(model.pielou.gam)
plot(model.pielou.gam, pages=1, se=TRUE)
gam.check(model.pielou.gam)

# shannon
model.shannon.gam <- gam(Shannon ~ s(age_months) + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                          data = alpha.diversity.votu,
                          family = gaussian(),
                          method = "REML")
# saveRDS(model.shannon.gam, "rds/model.shannon.gam.minimal.rds")
# model.shannon.gam <- readRDS("rds/model.shannon.gam.minimal.rds")

summary(model.shannon.gam)
plot(model.shannon.gam, pages=1, se=TRUE)
gam.check(model.shannon.gam)

#richness

model.richness.gam <- gam(Observed ~ s(age_months) + log10(depth) + s(age_months, study, bs = "fs", m = 1) + s(subjectID, bs = "re"),
                          data = alpha.diversity.votu,
                          family = nb(),
                          method = "REML")
# saveRDS(model.richness.gam, "rds/model.richness.gam.minimal.rds")
# model.richness.gam <- readRDS("rds/model.richness.gam.minimal.rds")

summary(model.richness.gam)
plot(model.richness.gam, pages=1, se=TRUE)
gam.check(model.richness.gam)

alpha.diversity.plot <- 
  model.shannon.gam %>% 
  smooth_estimates() %>% 
  add_confint() %>% 
  mutate(metric = "shannon") %>% 
  rbind(model.richness.gam %>% 
          smooth_estimates() %>% 
          add_confint() %>% 
          mutate(metric = "richness")
          ) %>% 
  rbind(model.pielou.gam %>% 
          smooth_estimates() %>% 
          add_confint() %>% 
          mutate(metric = "pielou")
  ) %>% 
  filter(.smooth == "s(age_months)") %>% 
  ggplot(aes(x = age_months, y = .estimate, color = metric)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = metric), linewidth = 0.1, alpha = 0.1) +
  scale_colour_manual(values = c("#DCB0F2FF", "#9EB9F3FF", "#D3B484FF"), labels = c("Pielou\n(p=0.15)", "Richness\n(p=0.04)","Shannon\n(p=0.85)")) +
  scale_fill_manual(values = c("#DCB0F2FF", "#9EB9F3FF", "#D3B484FF"), labels = c("Pielou\n(p=0.15)", "Richness\n(p=0.04)","Shannon\n(p=0.85)")) +
  theme_bw() + 
  labs(title = "Alpha diversity",
       x = "Age (months)",
       y = "Partial effect") +
  theme(legend.title = element_blank(),
        legend.key.spacing.y = unit(0.15, 'cm'),
        legend.position = "bottom"
        ) +
  scale_x_continuous(breaks = seq(0, 36, by = 6))

alpha.diversity.plot

gratia::derivatives(model.richness.gam, term = "s(age_months)") %>% 
  filter(.lower_ci <= 0) %>% 
  summarise(first_month_negative = min(age_months))
# 7.68 months is when lower CI crosses below 0, so roughly 8 months

fig.s1 <- gratia::derivatives(model.richness.gam, term = "s(age_months)") %>%
  ggplot(aes(x = age_months, y = .derivative)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), linewidth = 0.1, alpha = 0.1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment(x = 7.68, xend = 7.68, y = -Inf, yend = 0, linetype = "dashed") +
  ggtitle("First derivative of age_months smooth for richness GAMM") +
  ylab("Derivative") +
  xlab("Age (months)")

fig.s1

