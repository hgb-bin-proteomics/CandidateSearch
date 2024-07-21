library(tidyverse)
library(readxl)

b = "1B"
fn <- paste0("benchmark_hpc_", b, ".xlsx")

data <- read_xlsx(fn)

data$Method <- factor(data$Method, levels = c("f32CPU_DV",
                                              "i32CPU_DV",
                                              "f32CPU_SM",
                                              "i32CPU_SM",
                                              "f32CPU_DM",
                                              "i32CPU_DM",
                                              "f32GPU_DV",
                                              "f32GPU_DM"))

grp08_palette = c("#44485D", "#5C6378", "#BDC6D1", "#D5DCE3",
                  "#E1E7EC", "#EDF2F4", "#EF233C", "#E41433")

text_col <- c(rep("white", 2), rep("black", 4), rep("white", 2))

# export as 1200 x 800

ggplot(data, aes(x = Method, y = Mean, fill = Method)) +
  geom_bar(stat="identity", color = "black") +
  geom_errorbar(aes(x = Method, ymin = Mean-SD, ymax = Mean+SD), colour = "black", width = 0.3, linewidth = 1.0) +
  scale_fill_manual(values = grp08_palette) +
  ggtitle(paste0("Execution time of the different methods for benchmark ", b, ".")) +
  xlab("Method") +
  ylab("Time (s)") +
  geom_text(aes(label=round(Mean, 2)), size = 5.0, color = text_col, position = position_stack(vjust = 0.5)) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
