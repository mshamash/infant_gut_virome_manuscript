library(tidyverse)
library(rhmmer)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

# import annotations, descriptions, and AMG list
hmmsearch.kofam <- read_tblout(file.path('annotations', 'kofam-output.txt')) %>%
  select(domain_name, query_name, sequence_score) %>%
  group_by(domain_name) %>%
  arrange(desc(sequence_score)) %>%
  dplyr::slice(1) %>%
  select(-sequence_score)

vibrant.amg.list <- scan("annotations/VIBRANT_AMGs.txt", "")

annotations.kofam <- read_tsv("annotations/kofam.annotations.tsv") %>%
  select(knum, definition) %>%
  `colnames<-`(c("query_name", "description")) %>%
  filter(query_name %in% vibrant.amg.list)

hmmsearch.kofam.annotated <- inner_join(hmmsearch.kofam, annotations.kofam)

# import protein cluster info
protein_clusters <- read_tsv("protein_clusters.tsv", col_names = c("cluster", "protein")) %>%
  mutate(contig = str_remove(protein, "_[^_]+$"))

hmmsearch.kofam.annotated.contigs <- right_join(protein_clusters, hmmsearch.kofam.annotated,
                                               by = c("cluster" = "domain_name"))

# ps.phage.relabund <- readRDS("rds/ps.phage.relabund.rds")

votu_table <- as.data.frame(t(otu_table(ps.phage.relabund)))
class(votu_table) <- "data.frame"
votu_table$sample <- rownames(votu_table)

votu_table.long <- votu_table %>%
  pivot_longer(cols = -sample, names_to = "votu", values_to = "abundance") %>%
  filter(abundance > 0)

# iphop_out.split <- readRDS("rds/iphop_out.split.rds")

amg_table.long <- votu_table.long %>%
  inner_join(hmmsearch.kofam.annotated.contigs, by = c("votu" = "contig")) %>%
  left_join(iphop_out.split %>% select(Virus, family), by = c("votu" = "Virus"))

vibrant.amg.info <- read_tsv("annotations/VIBRANT_KEGG_pathways_summary.tsv") %>%
  separate_rows(KOs, sep = "~") %>%
  rename(KO = KOs) %>%
  filter(KO %in% vibrant.amg.list)

amg_table.long.info <- amg_table.long %>%
  left_join(vibrant.amg.info, by = c("query_name" = "KO"), relationship = "many-to-many")
# saveRDS(amg_table.long.info, "rds/amg_table.long.info.rds")
# amg_table.long.info <- readRDS("rds/amg_table.long.info.rds")

amg_richness <- amg_table.long.info %>%
  group_by(sample) %>%
  summarize(amg_richness = n())

phage_metadata <- read_tsv("metadata.tsv")

amg.matrix <- 
  ### USE PHAGE ABUNDANCE AS AMG ABUNDANCE, NOT CONSIDERING COPY NUMBER ###
  amg_table.long.info %>% 
  select(-Pathway) %>% 
  distinct() %>% 
  group_by(sample, votu, query_name) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(sample, query_name) %>% 
  summarize(amg.rel.abund = sum(abundance)) %>% 
  select(query_name, sample, amg.rel.abund) %>% 
  as.data.frame() %>% 
  left_join(phage_metadata %>% select(run, age_months), by = c("sample" = "run"))


### HEATMAP ###
amg.matrix$age_bin <- ifelse(
  amg.matrix$age_months <= 3,
  0, # months 0,1,2 in bin 1
  ((amg.matrix$age_months-4) %/% 3) + 1
)

heatmap.metabolism <- amg.matrix %>% 
  left_join(vibrant.amg.info %>% 
              select(KO, Metabolism) %>% 
              distinct(), 
            by = c("query_name" = "KO")) %>% 
  filter(!is.na(Metabolism)) %>% 
  group_by(sample, Metabolism) %>% 
  mutate(metabolism.rel.abund = sum(amg.rel.abund)) %>% 
  group_by(age_bin, Metabolism) %>% 
  filter(age_months <= 36) %>% 
  summarise(mean = mean(metabolism.rel.abund)) %>% 
  pivot_wider(values_from = mean, names_from = age_bin, values_fill = 0) %>%
  as.data.frame() %>% 
  column_to_rownames("Metabolism")

heatmap_scaled.amg <- t(apply(heatmap.metabolism, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))

heatmap_colors.amg = viridisLite::viridis(100, option = "viridis")

lgd_abundance.amg = Legend(
  title = "Scaled AMG\nabundance", 
  col_fun = colorRamp2(seq(0, 1, length.out = 100), heatmap_colors.amg), 
  at = c(0, 1), 
  labels = c("min", "max"),
  direction = "vertical"
)

age_labels.amg = c("0-3", "4-6", "7-9", "10-12", "13-15", "16-18", 
               "19-21", "22-24", "25-27", "28-30", "31-33", "34-36")

CombinedPlot <- Heatmap(
  heatmap_scaled.amg,
  col = heatmap_colors.amg,
  show_heatmap_legend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  column_labels = age_labels.amg,
  column_names_side = "bottom",
  column_title_side = "bottom",
  column_names_rot = 45,
  row_names_side = "left",
  column_title = "Age (months)",
  row_names_max_width = unit(10, "cm")
)

amg.heatmap <- ComplexHeatmap::draw(CombinedPlot, 
                     annotation_legend_list = lgd_abundance.amg, 
                     annotation_legend_side = "right")

amg.heatmap


### ENRICHMENT ANALYSIS ###
amg_matrix_votu <- amg_table.long.info %>% 
  select(votu, Pathway) %>% 
  filter(!is.na(Pathway)) %>% 
  distinct() %>% 
  mutate(present = 1) %>% 
  pivot_wider(values_from = present, names_from = votu, values_fill = 0) %>% 
  column_to_rownames("Pathway") %>% 
  t() %>% 
  as.data.frame()

N_total <- 49705 # nrow(amg_matrix_votu)

universe_counts <- data.frame(
  Pathway = colnames(amg_matrix_votu),
  m = colSums(amg_matrix_votu)
)

# iphop_out.split <- readRDS("rds/iphop_out.split.rds")

amg_long <- amg_matrix_votu %>%
  mutate(Virus = rownames(.)) %>%
  left_join(iphop_out.split %>% select(Virus, family)) %>%
  pivot_longer(
    cols = -c(Virus, family),
    names_to = "Pathway",
    values_to = "Present"
  )

phf_sizes <- iphop_out.split %>% 
  ungroup() %>% 
  select(Virus, family) %>% 
  distinct() %>% 
  count(family, name = "k") %>% 
  filter(family != "")
  
phf_hits <- amg_long %>%
  group_by(family, Pathway) %>%
  summarise(q = sum(Present), .groups = "drop") %>%
  filter(q > 0) %>% 
  filter(!is.na(family))

enrichment_results <- phf_hits %>%
  left_join(phf_sizes, by = "family") %>% 
  left_join(universe_counts, by = "Pathway") %>% 
  mutate(n = N_total - m) %>% 
  filter(k >= 5) %>% 
  mutate(
    # q = number of vOTUs in a specific PHF that carry the specific AMG pathway
    # m = total number of vOTUs in entire dataset which DO carry the specific AMG pathway
    # n = total number of vOTUs in entire dataset which DO NOT carry the specific AMG pathway
    # k = number of vOTUs in a specific PHF
    p_val = phyper(q - 1, m, n, k, lower.tail = FALSE),
    # enrichment ratio calculation
    # ratio = (AMG rate within PHF) / (AMG rate in entire dataset)
    enrichment_ratio = (q / k) / (m / N_total)
  ) %>%
  group_by(Pathway) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH")) %>% # FDR correction
  ungroup() %>%
  filter(p_adj < 0.05) # filter for significance

enrichment.plot.data <- enrichment_results %>%
  left_join(vibrant.amg.info %>% select(Metabolism, Pathway) %>% unique(), by = c("Pathway" = "Pathway")) %>% 
  left_join(iphop_out.split %>% ungroup() %>% select(family, phylum) %>% unique(), by = c("family" = "family")) %>% 
  arrange(Metabolism, Pathway) %>%
  mutate(Pathway = factor(Pathway, levels = unique(Pathway))) %>% 
  filter(family %in% (top_phages.phf %>% head(10) %>% pull(family)))

enrichment.plot.data %>% group_by(family) %>% count() %>% arrange(desc(n))

enrichment.plot <- enrichment.plot.data %>% 
  ggplot(aes(x = family, y = Pathway)) +
  geom_point(aes(size = enrichment_ratio, color = -log10(p_adj))) +
  scale_color_viridis_c(option = "plasma", name = "-log10(FDR)", direction = 1) +
  scale_size_continuous(name = "Enrichment\nratio", range = c(2, 6)) +
  facet_grid(Metabolism ~ phylum, scales = "free", space = "free", switch = "both") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        strip.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        panel.spacing = unit(0.1, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
        ) +
  labs(title = "AMG enrichment among prevalent PHFs")

enrichment.plot

