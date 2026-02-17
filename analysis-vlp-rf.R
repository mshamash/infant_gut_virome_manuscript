library(tidyverse)
library(randomForest)
library(reshape2)
library(caret)
library(purrr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

phage_metadata <- phage_metadata <- read_tsv("metadata.tsv")
rownames(phage_metadata) <- phage_metadata$run

rf.otutable.phage <- as.data.frame(t(otu_table(ps.phage.relabund)))
class(rf.otutable.phage) <- "data.frame"

rf.otutable.phage$SampleName <- rownames(rf.otutable.phage)
rf.otutable.phage <- left_join(rf.otutable.phage, 
                               phage_metadata %>% select(run, age_months), 
                               by = c("SampleName" = "run"))

rownames(rf.otutable.phage) <- rf.otutable.phage$SampleName
rf.otutable.phage <- rf.otutable.phage %>% subset(select = -c(SampleName))
names(rf.otutable.phage) <- gsub(x = names(rf.otutable.phage), pattern = "\\|", replacement = ".")

#saveRDS(rf.otutable.phage, "rds/rf.otutable.phage.phf.rds")
# rf.otutable.phage <- readRDS("rds/rf.otutable.phage.phf.rds")

set.seed(123)
trainIndex <- createDataPartition(rf.otutable.phage$age_months, p = 0.7,
                                  list = F, times = 1)

trainData <- rf.otutable.phage[trainIndex, ]
testData <- rf.otutable.phage[-trainIndex, ]

train_x <- as.matrix(trainData %>% select(-age_months))
train_y <- trainData$age_months
test_x <- as.matrix(testData %>% select(-age_months))
test_y <- testData$age_months

set.seed(123)
rf.phage <- randomForest(x = train_x, y = train_y)
# saveRDS(rf.phage, "rds/rf.phage.phf.rds")
# rf.phage <- readRDS("rds/rf.phage.phf.rds")

pred.test <- predict(rf.phage, test_x)
pred.train <- predict(rf.phage, train_x)


### MODEL EVALUATION ###

rf.plot <- 
  rbind(cbind(actual = test_y, pred = pred.test, dataset = "test", run=rownames(test_x)),
      cbind(actual = train_y, pred = pred.train, dataset = "train", run=rownames(train_x))) %>% 
  as.data.frame() %>% 
  mutate(actual = as.numeric(actual), pred = as.numeric(pred)) %>% 
  left_join(phage_metadata %>% select(study, run)) %>% 
  mutate(dataset = factor(dataset, levels = c("train", "test"))) %>% 
  ggplot(aes(x = actual, y = pred, color = dataset)) + 
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", se = F, formula = y ~ splines::ns(x, 3)) +
  geom_abline() +
  scale_x_continuous(breaks = seq(0, 42, by = 6), limits = c(0, 39)) +
  scale_y_continuous(breaks = seq(0, 42, by = 6), limits = c(0, 39)) +
  paletteer::scale_colour_paletteer_d("MetBrewer::Egypt", name = "Dataset", labels = c("Train", "Test")) +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = "bottom") +
  ylab("Predicted age (months)") +
  xlab("Actual age (months)")

rf.plot

# Calculate performance metrics
test_rmse <- sqrt(mean((test_y - pred.test)^2))
test_r2 <- cor(test_y, pred.test)^2
train_rmse <- sqrt(mean((train_y - pred.train)^2))
train_r2 <- cor(train_y, pred.train)^2


### MODEL EVALUATION SLIDING WINDOW ###
rf.results.test <- data.frame(actual = test_y, pred = pred.test)

window_width <- 3
step_size <- 1
min_age <- 1
max_age <- 37

# calculate center points for each window
half_width <- (window_width - 1) / 2
center_ages <- seq(from = (min_age + half_width), to = (max_age - half_width), by = step_size)

calculate_window_error <- function(center) {
  start_age <- center - half_width
  end_age <- center + half_width
  
  window_data <- rf.results.test %>% filter(actual >= start_age, actual <= end_age) # get all samples within this window of actual ages
  
  if(nrow(window_data) == 0) {
    return(data.frame(age_center = center, mae = NA, n_samples = 0)) # return NA if no samples in this window
  }
  
  error <- sqrt(mean((window_data$actual - window_data$pred)^2)) # calculate RMSE
  
  return(data.frame(age_center = center, mae = error, n_samples = nrow(window_data)))
}

# map_dfr() runs the function for each 'center_age' and
# combines the results into a single data frame
sliding_error_df <- map_dfr(center_ages, calculate_window_error)


rf.sliding.window.fig <- ggplot(sliding_error_df, aes(x = age_center, y = mae)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  geom_bar(aes(y = n_samples / 10), stat = "identity", fill = "grey", alpha = 0.3) +
  scale_y_continuous(
    name = "Model error (RMSE, in months)",
    sec.axis = sec_axis(~ . * 10, name = "Number of samples in window")
  ) +
  labs(
    title = "Model performance (3-month sliding window)",
    x = "Center of age window (months)"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 42, by = 6), limits = c(0, 39))


### IMPORTANT TAXA/PHF ###

important.taxa <- randomForest::importance(rf.phage) %>% as.data.frame()

important.taxa$`Host genome` <- rownames(important.taxa)
rownames(important.taxa) <- NULL

# taxtable <- readRDS("rds/taxtable.rds")

important.taxa <- left_join(important.taxa, taxtable %>% select(phylum, class, order, family, `Host genome`))

### calculate R2 for PHFs
important.taxa.r2 <- trainData %>% 
  pivot_longer(cols = -age_months, names_to = "phf", values_to = "relabund") %>% 
  select(phf, relabund, age_months) %>% 
  group_by(phf) %>% 
  summarize(cor = sign(cor(relabund, age_months)),
            r2 = cor(relabund, age_months)^2,
            sign_r2 = r2*cor)

important.taxa <- important.taxa %>% 
  left_join(important.taxa.r2, by = c("Host genome" = "phf")) %>% 
  top_n(20, IncNodePurity)

max_abs_R2 <- max(abs(important.taxa$sign_r2), na.rm = TRUE)
symm_limits <- c(-max_abs_R2, max_abs_R2)

neg_colors <- viridisLite::mako(100,begin = 0.3)
pos_colors <- viridisLite::magma(100,begin = 0.6)
all_colors <- c(neg_colors, rev(pos_colors))

### HEATMAP ###
# ps.phage.relabund.phf <- readRDS("rds/ps.phage.relabund.phf.rds")
ps.phage.relabund.phf.summary <- ps.phage.relabund.phf %>% 
  psmelt() %>%              
  filter(Abundance > 0) %>% 
  filter(family %in% (important.taxa %>%
                        top_n(20, IncNodePurity))$family) %>%
  group_by(Sample, family) %>% 
  summarise(study = study, Abundance = sum(Abundance), age_months = age_months) %>% 
  distinct() %>% 
  ungroup()

ps.phage.relabund.phf.summary$age_bin <- ifelse(
  ps.phage.relabund.phf.summary$age_months <= 3,
  0, # months 0,1,2 in bin 1
  ((ps.phage.relabund.phf.summary$age_months-4) %/% 3) + 1
)

heatmap <- ps.phage.relabund.phf.summary %>% 
  group_by(age_bin, family) %>% 
  filter(age_months <= 36) %>% 
  summarise(mean = mean(Abundance), age_bin = age_bin) %>% 
  distinct() %>% 
  pivot_wider(values_from = mean, names_from = age_bin) %>% 
  as.data.frame()
rownames(heatmap) <- heatmap[, 1]
heatmap <- heatmap[, -1]

heatmap[is.na(heatmap)] <- 0


#### combine heatmap with importance, r2 ####
col_r2 = colorRamp2(
  seq(symm_limits[1], symm_limits[2], length.out = length(all_colors)), 
  all_colors
)

lgd_r2 = Legend(
  title = bquote(bold("sign * R"^2)),
  col_fun = col_r2, 
  at = c(symm_limits[1], 0, symm_limits[2]), 
  labels = c(round(symm_limits[1], 2), "0", round(symm_limits[2], 2))
)

heatmap_colors = viridisLite::viridis(100, option = "viridis")

lgd_abundance = Legend(
  title = "Scaled relative\nabundance", 
  col_fun = colorRamp2(seq(0, 1, length.out = 100), heatmap_colors), 
  at = c(0, 1), 
  labels = c("min", "max"),
  direction = "vertical"
)

heatmap_scaled <- t(apply(heatmap, MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))
row_data <- important.taxa[match(rownames(heatmap_scaled), important.taxa$family), ]

row_anno = rowAnnotation(
  "sign * R2" = anno_simple(row_data$sign_r2, col = col_r2, width = unit(5, "mm"), border = T),
  
  "spacer" = anno_empty(border = FALSE, width = unit(0.5, "mm")),

  "Model importance" = anno_barplot(
    row_data$IncNodePurity, 
    baseline = 0,
    axis_param = list(side = "bottom"),
    gp = gpar(fill = "#D1D1D1", col = "white"),
    width = unit(4, "cm")
  ),
  annotation_label = list(
    "sign * R2" = expression("sign * R"^2)
  ),
  show_annotation_name = c("sign * R2" = TRUE, "spacer" = FALSE, "Model importance" = TRUE)
)

age_labels = c("0-3", "4-6", "7-9", "10-12", "13-15", "16-18", 
               "19-21", "22-24", "25-27", "28-30", "31-33", "34-36")

CombinedPlot <- Heatmap(
  heatmap_scaled,
  col = heatmap_colors,
  show_heatmap_legend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  right_annotation = row_anno,
  column_labels = age_labels,
  column_names_side = "bottom",
  column_title_side = "bottom",
  column_names_rot = 45,
  row_names_side = "left",
  row_names_gp = gpar(
    fontface = "italic"
  ),
  column_title = "Age (months)"
)

rf.heatmap <- ComplexHeatmap::draw(CombinedPlot, 
     annotation_legend_list = packLegend(lgd_abundance, lgd_r2), 
     annotation_legend_side = "right")

rf.heatmap

pdf("~/Desktop/Figure3D.pdf", width = 9, height = 5)
draw(rf.heatmap)
dev.off()
