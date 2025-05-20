library(dplyr)
library(readr)
library(ggplot2)

convert_to_seconds <- function(time_str) {
  minutes <- as.numeric(sub("m.*", "", time_str))
  seconds <- as.numeric(sub(".*m|s", "", time_str))
  total_seconds <- (minutes * 60) + seconds
  return(total_seconds)
}

dirs <- c("/nfshomes/lmatalop/DolloLL-CDP-study/simulation
/condor", "/nfshomes/lmatalop/DolloLL-CDP-study/simulation
/dollocdp", "/nfshomes/lmatalop/DolloLL-CDP-study/simulation
/dollocdpLL")

df <- do.call(rbind, lapply(dirs, function(dir) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  do.call(rbind, lapply(files, read.csv))
}))

df$ncells <- as.factor(df$ncells)

custom_colors <- c("dollocdp" = "#2A9D8F", "dollocdpLL" = "#E9C46A", "condor" = "#E76F51")

  ggplot(df, aes(x = factor(ncells), y = seconds, fill = method, color = method)) +
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
    y = "Runtime in seconds",
    fill = "Method",
  ) +
  theme_minimal() 

