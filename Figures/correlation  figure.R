library(ggplot2)
library(dplyr)
library(reshape2)

setwd("~/Documents/Documents/PhD work/Multispeq project/multispeq paper data")
# Load and clean the dataset
multispeq_paper_coreelation <- read.csv("filtered_multispeq_cutoff.csv")
multispeq_paper_coreelation <- multispeq_paper_coreelation[-1]  # Assuming first column is to be removed

# Remove any non-numeric columns which are not suitable for correlation
df <- multispeq_paper_coreelation[, sapply(multispeq_paper_coreelation, is.numeric)]

# Calculate the correlation matrix
cor_matrix <- cor(df, use = "pairwise.complete.obs")

# Fill the p-value matrix by looping over pairs of variables
variable_names <- names(df)
p_matrix <- matrix(NA, nrow = length(variable_names), ncol = length(variable_names),
                   dimnames = list(variable_names, variable_names))
for (i in 1:(length(variable_names) - 1)) {
  for (j in (i + 1):length(variable_names)) {
    test_result <- cor.test(df[[variable_names[i]]], df[[variable_names[j]]], method = "pearson")
    p_matrix[i, j] <- test_result$p.value
    p_matrix[j, i] <- test_result$p.value  # Symmetric filling
  }
}

# Create a significance matrix
significance_matrix <- p_matrix < 0.05

# Melt the correlation and significance matrices
cor_melted <- melt(cor_matrix)
sig_melted <- melt(significance_matrix)

# Combine the melted data
cor_melted$significance <- sig_melted$value

# Generate a key from the row and column names of cor_matrix to determine positioning
names_list <- rownames(cor_matrix)
name_positions <- setNames(seq_along(names_list), names_list)
cor_melted$triangle <- ifelse(name_positions[cor_melted$Var1] < name_positions[cor_melted$Var2], "upper",
                              ifelse(name_positions[cor_melted$Var1] > name_positions[cor_melted$Var2], "lower", "diag"))

# Define the custom labels with expressions for the traits using uppercase Phi symbol
trait_labels <- c(
  "Relative Chlorophyll", "qL", "PSI-ac",
  bquote(vH^"+"), 
  bquote(NPQ[T]), 
  "Fv'/Fm'",
  bquote(Phi*NO), 
  bquote(gH^"+"), 
  bquote(Phi*PSII), 
  bquote(Phi*NPQ), 
  bquote(ECS[T]), 
  "PSI-oxc", "PSI-orc", "PSI-opc"
)

# Define a common theme with increased tick length and legend text size
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22, color = "black", hjust = 1),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),  # Increase tick length
    axis.ticks = element_line(size = 1.2),  # Increase tick width
    legend.position = "top",
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.background = element_blank(),  # Remove the box around the legend
    legend.key.size = unit(4, "lines"),  # Increase the size of the legend key
    legend.spacing.x = unit(1, "cm"),  # Increase horizontal spacing in the legend
    legend.spacing.y = unit(1, "cm"),  # Increase vertical spacing in the legend
    legend.text = element_text(size = 22),  # Increase legend text size
    plot.title = element_text(hjust = 0.5, color = "black")
  )

# Plotting setup with the common theme applied
p <- ggplot(data = cor_melted, aes(x = Var1, y = Var2)) +
  geom_tile(data = subset(cor_melted, triangle == "diag"), fill = NA, color = NA) +  # Removing display for the diagonal
  # Upper triangle significant correlations with circle size representing the magnitude
  geom_point(data = subset(cor_melted, triangle == "upper" & significance),
             aes(size = abs(value), fill = value), shape = 21, stroke = 0) +
  # Lower triangle significant correlations with text color representing the correlation strength
  geom_text(data = subset(cor_melted, triangle == "lower" & significance),
            aes(label = sprintf("%.2f", value), color = value), size = 8) +
  geom_text(data = subset(cor_melted, triangle == "upper" & !significance), aes(label = "X"), color = "black", size = 8) +
  geom_text(data = subset(cor_melted, triangle == "lower" & !significance), aes(label = "X"), color = "black", size = 8) +
  scale_fill_gradient2(low = "red", mid = "lightgrey", high = "blue", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "") +
  scale_color_gradient2(low = "red", mid = "lightgrey", high = "blue", midpoint = 0, limit = c(-1, 1), space = "Lab", guide = "none") +  # Disable the guide for color gradient
  scale_size_continuous(range = c(3, 12), guide = "none") +  # Adjust the size range and remove legend for size
  common_theme +  # Apply the common theme
  scale_x_discrete(labels = trait_labels) +
  scale_y_discrete(labels = trait_labels) +
  labs(x = NULL, y = NULL)

# Print the plot
print(p)

# Optionally, save the combined plot
ggsave("correlation.svg", plot = p, width = 18, height = 12, dpi = 300)
