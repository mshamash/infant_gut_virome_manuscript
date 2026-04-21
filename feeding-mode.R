library(tidyverse)


metadata <- read_tsv("metadata.tsv")

metadata.feeding.mode.only <- metadata %>% 
  filter(study != "Garmaeva2024") %>% 
  select(run, study, age_months, diet) %>% 
  group_by(study, age_months, diet) %>% 
  summarize(n_samples = n()) %>% 
  filter(!is.na(diet))

metadata %>% 
  filter(study != "Garmaeva2024") %>% 
  #filter(!is.na(sex)) %>% 
  #filter(!is.na(delivery_mode)) %>% 
  filter(!is.na(diet)) %>% 
  pull(study) %>% 
  table()

diet.plot <- metadata.feeding.mode.only %>% 
  ggplot(aes(x = age_months, y = n_samples, fill = diet)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ study, scales = "free_y") +
  theme_bw() +
  labs(fill = "Diet at time of\nsampling") +
  xlab("Age (months)") +
  ylab("# of samples") +
  theme(aspect.ratio = 0.8)

diet.plot

ggsave("~/Desktop/SuppFigure1.pdf", plot = diet.plot, width = 7, height = 4, limitsize = FALSE)


