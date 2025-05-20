library(dplyr)
library(readr)
library(ggplot2)

files <- list.files(
  path = "/fs/cbcb-lab/ekmolloy/lmatalop/dolloLL-study/simulation", 
  pattern = "_evaluation_contracted\\.csv$", 
  recursive = TRUE,
  full.names = TRUE
)

# Read and combine
df <- bind_rows(lapply(files, read_csv))

# Make ncells a factor so it's treated as categorical on x-axis
df$ncells <- as.factor(df$ncells)

custom_colors <- c("dollocdp" = "#2A9D8F", "condor" = "#E76F51")

#Box Plots

metrics <- c("FP","FN","FPR","FNR","Precision","Recall")

plots <- lapply(metrics, function(metric) {
  ggplot(df, aes(x = ncells, y = .data[[metric]], fill = method)) +
    geom_boxplot(
      aes(group = interaction(ncells,method)),
      position = position_dodge(width = 0.6),
      alpha = 0.6,
      width = 0.4
    ) +
    geom_point(
      aes(group = interaction(ncells,method)),
      position = position_dodge(width = 0.6),
      size = 1.5,
      alpha = 0.3,
      color = "black" # Only points are black outline
    ) +
    scale_fill_manual(values = custom_colors) +
    labs(
      x = "Cells and Mutations",
      y = metric,
      fill = "Method"
    ) +
    theme_minimal()
})

for (p in plots) {
  print(p)
}
