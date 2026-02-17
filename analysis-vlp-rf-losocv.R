library(tidyverse)
library(caret)
library(ggtext)
library(reshape2)

phage_metadata <- read_tsv("metadata.tsv")

# rf.phage <- readRDS("rds/rf.phage.phf.rds")
# rf.otutable.phage <- readRDS("rds/rf.otutable.phage.phf.rds")
# rf.otutable.phage$run <- rownames(rf.otutable.phage)

rf.otutable.phage <- rf.otutable.phage %>% left_join(phage_metadata %>% select(run, study)) %>% select(-run)

studies <- unique(phage_metadata$study)

# LOSOCV, train a model leaving out a different study eachtime
results <- lapply(studies, function(s) {
  train_df <- rf.otutable.phage %>% filter(study != s)
  test_df  <- rf.otutable.phage %>% filter(study == s)

  set.seed(123)
  rf_mod <- randomForest(x = train_df %>% select(-c(study, age_months)), y = train_df$age_months)
  preds.losocv <- predict(rf_mod, test_df %>% select(-study))
  preds.fullmodel <- predict(rf.phage, test_df %>% select(-study))

  data.frame(
    study = s,
    obs   = test_df$age_months,
    pred.losocv  = preds.losocv,
    pred.full = preds.fullmodel
  )
})

all_preds <- bind_rows(results)
# saveRDS(all_preds, "rds/rf.phage.phf.losocv.results.rds")
# all_preds <- readRDS("rds/rf.phage.phf.losocv.results.rds")


# RMSE
all_preds.stats <- all_preds %>% 
  group_by(study) %>% 
  summarize(
    RMSE.losocv = sqrt(mean((pred.losocv-obs)^2)),
    RMSE.full = sqrt(mean((pred.full-obs)^2))
  )

losocv.rmse.fig <- all_preds.stats %>% 
  melt(id.vars = "study", variable.name = "metric") %>% 
  mutate(metric.class = ifelse(metric %in% c("RMSE.full", "RMSE.losocv"), "RMSE", "")) %>% 
  mutate(model = ifelse(metric %in% c("RMSE.full"), "full", "losocv")) %>% 
  filter(metric %in% c("RMSE.full", "RMSE.losocv")) %>% 
  mutate(model = factor(model, levels = c("full", "losocv"))) %>% 
  ggplot(aes(x = study, y = value, fill = model)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~ factor(metric.class, c("RMSE")), scales = "free", strip.position = "left",
             labeller = as_labeller(c(RMSE = "RMSE (months)") ) ) +
  theme_bw() +
  ylab(NULL) +
  xlab("Study") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(strip.placement = "outside", strip.background = element_blank()) +
  scale_fill_brewer(palette = "Set1", labels = c("Full model", "LOSOCV")) +
  labs(fill = "Model") +
  theme(strip.text = ggtext::element_markdown(),
        strip.text.y.left = ggtext::element_markdown()) +
  theme(aspect.ratio = 0.75)

losocv.rmse.fig

# plot pred vs obs for each model
losocv.preds.fig <- all_preds %>% 
  melt(id.vars = c("study", "obs"), variable.name = "pred") %>% 
  mutate(pred = factor(pred, levels = c("pred.full", "pred.losocv"))) %>% 
  ggplot(aes(x = obs, y = value, color = pred)) +
  geom_jitter(alpha = 0.5) +
  facet_wrap(~ study) +
  scale_x_continuous(breaks = seq(0, 42, by = 6), limits = c(0, 39)) +
  scale_y_continuous(breaks = seq(0, 42, by = 6), limits = c(0, 39)) +
  scale_color_brewer(palette = "Set1", labels = c("Full model", "LOSOCV")) +
  ylab("Predicted age (months)") +
  xlab("Actual age (months)") +
  theme_bw() +
  theme(aspect.ratio = 1) +
  labs(color = "Model")

losocv.preds.fig


