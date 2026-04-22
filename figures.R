library(tidyverse)
library(cowplot)
library(sf)
library(rnaturalearth)
library(magick)
library(pdftools)
library(flextable)

#### FIGURE 1 ####
phage_metadata <- read_tsv("metadata.tsv")

donut.data <- phage_metadata %>% 
  group_by(study) %>% 
  summarise(participants = n_distinct(subjectID), samples = n()) %>% 
  mutate(
    participants.fraction = participants / sum(participants),
    participants.ymax = cumsum(participants.fraction),
    participants.ymin = c(0, head(participants.ymax, n = -1)),
    participants.labelposition.y = (participants.ymax + participants.ymin) / 2,
    participants.labelposition.x = 4.4,
    samples.fraction = samples / sum(samples),
    samples.ymax = cumsum(samples.fraction),
    samples.ymin = c(0, head(samples.ymax, n = -1)),
    samples.labelposition.y = (samples.ymax + samples.ymin) / 2,
    samples.labelposition.x = 4.4,
  ) %>% 
  as.data.frame() %>% 
  `rownames<-`(.$study)

donut.data["Beller2022", "participants.labelposition.x"] <- 2.6
donut.data["Beller2022", "participants.labelposition.y"] <- 0.015
donut.data["Yan2021", "participants.labelposition.x"] <- 2.6
donut.data["Yan2021", "participants.labelposition.y"] <- 0.965
donut.data["Liang2020b", "participants.labelposition.y"] <- 0.17
donut.data["Pannaraj2018", "participants.labelposition.y"] <- 0.285

plot.donut.participants <- 
  donut.data %>% 
  ggplot(aes(ymax = participants.ymax, ymin = participants.ymin, xmax = 4, xmin = 3, fill = study)) +
  geom_rect(color = "white", linewidth = 0.5) +
  geom_text(aes(x = participants.labelposition.x, y = participants.labelposition.y, label = participants, color = study, fontface="bold"), size = 4) +
  annotate("text", x = 0, y = 0, label = "Participants\n(n=999)", size = 5, fontface = "bold") +
  paletteer::scale_colour_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  coord_polar(theta = "y") +
  #xlim(c(1, 4.4)) +
  theme_void() +
  theme(legend.position = "none")
plot.donut.participants


donut.data["Zhao2017", "samples.labelposition.y"] <- 0.98
donut.data["Liang2020b", "samples.labelposition.y"] <- 0.295
donut.data["McCann2018", "samples.labelposition.x"] <- 2.6
donut.data["McCann2018", "samples.labelposition.y"] <- 0.36
donut.data["Pannaraj2018", "samples.labelposition.y"] <- 0.375

plot.donut.samples <- 
  donut.data %>% 
  ggplot(aes(ymax = samples.ymax, ymin = samples.ymin, xmax = 4, xmin = 3, fill = study)) +
  geom_rect(color = "white", linewidth = 0.5) +
  geom_text(aes(x = samples.labelposition.x, y = samples.labelposition.y, label = samples, color = study, fontface = "bold"), size = 4) +
  annotate("text", x = 0, y = 0, label = "Samples\n(n=1,893)", size = 5, fontface = "bold") +
  paletteer::scale_colour_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  coord_polar(theta = "y") +
  #xlim(c(1, 4.4)) +
  theme_void() +
  theme(legend.position = "none")
plot.donut.samples

plot.hist.sampling <- phage_metadata %>% 
  ggplot(aes(x = age_months, fill = study)) +
  geom_histogram() +
  ggbreak::scale_y_break(c(250,600), scale = 0.25, expand = expansion(add = c(0, 0)), space = 0.3) +
  xlab("Age (months)") +
  ylab("Number of samples") +
  scale_x_continuous(breaks = seq(0, 38, by = 1), expand = expansion(add = c(0, 0))) +
  theme_classic() +
  paletteer::scale_colour_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_Green_Orange_12") +
  theme(axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2))

plot.hist.sampling


rm_legend <- function(p){p + theme(legend.position = "none")}
plot.overview <- ggdraw() +
  draw_plot(print(rm_legend(plot.hist.sampling))) +
  draw_plot(plot.donut.samples, x = 0.27, y = 0.35, width = 0.55, height = 0.55) +
  draw_plot(plot.donut.participants, x = 0.55, y = 0.35, width = 0.55, height = 0.55)
plot.overview

world_map <- ne_countries(scale = "small", returnclass = "sf") %>% filter(iso_a3 != "ATA")

map.locations <- data.frame(
  study = c("Beller2022", "Garmaeva2024", "Liang2020a", "Liang2020b", "Lim2015", "Maqsood2019", "McCann2018", 
            "Pannaraj2018", "Shah2023", "Walters2023", "Yan2021", "Zhao2017", "Zhao2017"),
  lat = c(50.841, 52.363, 39.952, 39.989, 38.639, 38.631, 51.882, 
          34.098, 55.676, 37.424, 38.808, 60.170, 58.378),
  lon = c(4.360, 4.904, -75.161, -75.215, -90.279, -90.392, -8.511, 
          -118.290, 12.570, -122.166, 121.315, 24.939, 26.727)
)

set.seed(12)
plot.map <- ggplot() +
  geom_sf(data = world_map, fill = "grey80", color = "grey80") +
  geom_sf(data = st_jitter(st_as_sf(map.locations, 
                          coords = c("lon", "lat"), 
                          crs = 4326), 1.5), 
             aes(fill = study),
             shape = 21,
             color = "black",
             size = 4) +
  coord_sf(crs = "+proj=robin", expand = F) + # Robinson projection
  theme_void() + 
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_Green_Orange_12") + #ggthemes::Classic_Cyclic
  theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(t = 5, r = 0, b = 0, l = 0)
  )
plot.map


legend.plot.overview <- get_legend(plot.hist.sampling)

fig.1ab <- plot_grid(
  plot.map,
  legend.plot.overview,
  plot.overview,
  nrow = 3,
  rel_heights = c(1, 0.2, 1),
  labels = c("A", "", "B"),
  label_size = 22
)

fig.1ab

fig.1c <- ggdraw() + draw_image(image_read_pdf("~/Desktop/Figure1C.pdf", density = 600))


fig.1 <- plot_grid(
  fig.1ab, fig.1c,
  nrow = 1, 
  ncol = 2,
  labels = c("", "C"),
  label_size = 22,
  rel_widths = c(1, 1)
)

fig.1

save_plot("~/Desktop/Figure1.pdf", fig.1, dpi = 600, base_height = 10, base_width = 15)


#### FIGURE 2 ####
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 10, spaceLegend = 0.3) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

fig.2a <- ggdraw() +
  draw_plot(prevalence.rank.plot) +
  draw_plot(addSmallLegend(top.prevalent.plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.7)), pointSize = 2, textSize = 5, spaceLegend = 0.3), 
            x = 0.2, y = 0.375, width = 0.78, height = 0.55)

fig.2b <- ggdraw() +
  draw_plot(prevalence.rank.phf.plot) +
  draw_plot(addSmallLegend(top.prevalent.phf.plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.7)), pointSize = 2, textSize = 10, spaceLegend = 0.3), 
            x = 0.28, y = 0.4, width = 0.7, height = 0.525)

fig.2c <- alpha.diversity.plot

fig.2d <- ordination.plots

fig.2d.table <- permanova.flextable.grob

fig.2 <- plot_grid(
  plot_grid(
    fig.2a, fig.2b, fig.2c,
    nrow = 1, 
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 22
  ),
  plot_grid(
    fig.2d, fig.2d.table, 
    ncol = 3, 
    nrow = 1, 
    rel_widths = c(1, 0.5, 0.5),
    labels = c("D", ""),
    label_size = 22
  ),
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 1)
)

fig.2

save_plot("~/Desktop/Figure2.pdf", fig.2, dpi = 600, base_height = 10, base_width = 20)


#### FIGURE 3 ####
fig.3a <- velocity.plot

fig.3b <- dispersion.plot

fig.3c <- rf.plot

fig.3d <- ggdraw() + draw_image(image_read_pdf("~/Desktop/Figure3D.pdf", density = 600))

fig.3 <- plot_grid(
  plot_grid(
    fig.3a, fig.3b,
    nrow = 1, 
    ncol = 2,
    labels = c("A", "B"),
    label_size = 22,
    rel_widths = c(0.97, 1)
  ),
  plot_grid(
    fig.3c, fig.3d, 
    ncol = 2, 
    nrow = 1, 
    rel_widths = c(0.5, 1),
    labels = c("C", "D"),
    label_size = 22
  ),
  ncol = 1,
  nrow = 2,
  rel_heights = c(0.8, 1)
)

fig.3

save_plot("~/Desktop/Figure3.pdf", fig.3, dpi = 600, base_height = 7.5, base_width = 10)

#### FIGURE 4 ####
fig.4a <- temperate.plot

fig.4b <- temperate.phf.plot

fig.4 <- plot_grid(
  fig.4a, fig.4b,
  nrow = 1, 
  ncol = 2,
  labels = c("A", "B"),
  label_size = 22
)

fig.4

save_plot("~/Desktop/Figure4.pdf", fig.4, dpi = 600, base_height = 4, base_width = 12)
  

#### FIGURE 5 ####
ggsave("~/Desktop/Figure5.pdf", plot = enrichment.plot, width = 15, height = 20, limitsize = FALSE)


#### SUPP FIGURE 1 ####
ggsave("~/Desktop/SuppFigure1.pdf", plot = diet.plot, width = 7, height = 4, limitsize = FALSE)

#### SUPP FIGURE 2 ####
save_plot("~/Desktop/SuppFigure2.pdf", richness.derivative.plot, dpi = 600, base_height = 4, base_width = 7)

#### SUPP FIGURE 4 ####
fig.s4a <- rf.sliding.window.fig

fig.s4b <- losocv.rmse.fig

fig.s4c <- losocv.preds.fig

fig.s4 <- plot_grid(
  plot_grid(
    fig.s4a, fig.s4b,
    nrow = 2, 
    ncol = 1,
    labels = c("A", "B"),
    label_size = 22,
    rel_heights = c(1, 1)
  ),
  fig.s4c,
  ncol = 2,
  nrow = 1,
  labels = c("", "C"),
  label_size = 22,
  rel_widths = c(0.7, 1)
)

fig.s4

save_plot("~/Desktop/SuppFigure4.pdf", fig.s2, dpi = 600, base_height = 6, base_width = 13)

#### SUPP FIGURE 5 ####
pdf("~/Desktop/SuppFigure5.pdf", height = 4, width = 8)
draw(amg.heatmap)
dev.off()

#### SUPP FIGURE 6 ####
fig.s6a <- ggdraw() + draw_image(image_read_pdf("~/Desktop/SuppFigure6A.pdf", density = 400))

fig.s6b <- other.virus.abundance.plot

fig.s6c <- animal_virus.delivery_mode.plot

fig.s6 <- 
  plot_grid(
    fig.s6a,
    plot_grid(
      fig.s6b,
      fig.s6c,
      ncol = 2, 
      nrow = 1, 
      rel_widths = c(1, 1),
      labels = c("B", "C"),
      label_size = 22
    ),
    labels = c("A", ""),
    label_size = 22,
    ncol = 1,
    nrow = 2,
    rel_heights = c(1, 0.6)
  )

fig.s6

ggsave("~/Desktop/SuppFigure6.pdf", plot = fig.s6, width = 10, height = 7, limitsize = FALSE)
