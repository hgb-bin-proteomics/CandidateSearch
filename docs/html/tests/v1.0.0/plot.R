library(tidyverse)
library(readxl)

data <- read_xlsx("candidates_psms.xlsx")

facet_labels <- list("peptide"="Identifying peptides", "peptidoform"="Identifying peptidoforms")

facet_labeller <- function(variable, value){
  return(facet_labels[value])
}

ggplot(data, aes(x=as.factor(top_n), y=coverage_perc, fill=as.factor(top_n))) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("#2b2d42", "#8d99ae", "#ef233c", "#d90429")) +
  ggtitle("Dataset: HeLa_1ug_isow2_1hgradient_datasetA_rep2") +
  facet_grid(~ `method`, labeller = facet_labeller) +
  xlab("Top N") +
  ylab("% of psms covered by top n candidates") +
  ylim(c(0, 100)) +
  labs(fill = "Top N") +
  geom_text(aes(label=round(coverage_perc, 2)), size=5.0, color = "white", position = position_stack(vjust = 0.5)) +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) 
