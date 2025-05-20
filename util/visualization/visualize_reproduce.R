library(dplyr)
library(readr)
library(ggplot2)

files <- list.files(
  "/nfshomes/lmatalop/DolloLL-CDP-study/simulation/condor_eval", 
  pattern = "\\.csv$", 
  full.names = TRUE
)

df <- bind_rows(lapply(files, read_csv))

df <- df %>%
  mutate(RF_ratio = RF / RF_max)

df$ncells <- as.factor(df$ncells)

custom_colors <- c("condor" = "#2A9D8F")

# Create the boxplot
  ggplot(df, aes(x = factor(ncells), y = relation_accuracy, fill = method, color = method)) +
  geom_boxplot(color='black', width = 0.2, alpha = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(
    color = 'black',
    position = position_dodge(width = 0.6),
    size = 1.5,
    alpha = 0.3
  ) +
  scale_fill_manual(values = custom_colors) +  # Apply custom color palette for box fill
  labs(
    x = "Cells and Mutations",
    y = "Pairwise Ancestral Relation Accuracy",
    fill = "Method",
  ) +
  theme_minimal() + 
  ylim(0, 1)  # Set y-axis limits to 0-1

