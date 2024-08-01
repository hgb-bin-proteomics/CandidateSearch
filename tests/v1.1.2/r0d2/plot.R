library(tidyverse)
library(readxl)

r = "r0d2"
fn <- paste0(r, ".xlsx")

data <- read_xlsx(fn)

# export to 1000 x 800
ggplot(data, aes(x=as.factor(top_n), y=coverage_perc, fill=as.factor(top_n))) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("#2b2d42", "#8d99ae", "#ef233c", "#d90429")) +
  ggtitle(paste0("Result for run '", r, "'.")) +
  xlab("Top N") +
  ylab("% of psms covered by top n candidates") +
  ylim(c(0, 100)) +
  labs(fill = "Top N") +
  geom_text(aes(label=round(coverage_perc, 2)), size=5.0, color = "white", position = position_stack(vjust = 0.5)) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) 
