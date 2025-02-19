library(dplyr)
library(ggplot2)
library(patchwork)


# Load the datasets
chlorophyll_data <- read.csv('Michael cholorphyll.csv')
genotype_data <- read.csv('transposed_gene_genotype_data.csv')

# Merge datasets on 'Genotype_ID'
merged_datav1 <- merge(chlorophyll_data, genotype_data, by = "GenotypeID") %>%
  filter(!is.na(CHL))  # Remove missing values

# Load the datasets for the second plot
chlorophyll_data2 <- read.csv('filterted_allele_file.csv')
merged_datav2 <- merge(chlorophyll_data2, genotype_data, by = "GenotypeID")

# Load the datasets for the third plot
chlorophyll_data3 <- read.csv('merged_multispeq.csv')
merged_datav3 <- merge(chlorophyll_data3, genotype_data, by = "GenotypeID")

# Custom colors
custom_colors <- c("A|A" = "#FFA500", "G|G" = "steelblue")

# Define a common theme
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.2),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(1, "lines"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.text = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, color = "black")
  )

# Function to format p-values with italic text using bquote
format_p_value <- function(p) {
  if (p < 0.00001) {
    return(bquote(italic(p) < 0.00001))
  } else {
    return(bquote(italic(p) == .(format(p, digits = 4))))
  }
}

# First plot (Chlorophyll)
t_test_result1 <- t.test(CHL ~ Genotype, data = merged_datav1, var.equal = TRUE)
p_value1 <- t_test_result1$p.value
chloro1 <- ggplot(merged_datav1, aes(x = Genotype, y = CHL, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 1.5) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(breaks = c(300, 400, 500, 600)) + 
  labs(title = "", x = "", y = "Chlorophyll") +
  common_theme +
  annotate("text", x = 1.5, y = max(merged_datav1$CHL, na.rm = TRUE) + 0.5, 
           label = format_p_value(p_value1), size = 9, hjust = 0.5)

# Second plot (TPM)
t_test_result2 <- t.test(Zm00001eb377130 ~ Genotype.y, data = merged_datav2, var.equal = TRUE)
p_value2 <- t_test_result2$p.value
chloro2 <- ggplot(merged_datav2, aes(x = Genotype.y, y = Zm00001eb377130, fill = Genotype.y)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 1.5) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(breaks = c(-20, 0, 25, 50)) + 
  labs(title = "", x = "", y = "TPM") +
  common_theme +
  annotate("text", x = 1.5, y = max(merged_datav2$Zm00001eb377130, na.rm = TRUE) + 0.5, 
           label = format_p_value(p_value2), size = 9, hjust = 0.5)

# Third plot (Relative Chlorophyll)
t_test_result3 <- t.test(Relative_Chlorophyll ~ Genotype, data = merged_datav3, var.equal = TRUE)
p_value3 <- t_test_result3$p.value
chloro3 <- ggplot(merged_datav3, aes(x = Genotype, y = Relative_Chlorophyll, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 1.5) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(breaks = c(-20, -10, 0, 10)) + 
  labs(title = "", x = "", y = "Relative Chlorophyll") +
  common_theme +
  annotate("text", x = 1.5, y = max(merged_datav3$Relative_Chlorophyll, na.rm = TRUE) + 0.5, 
           label = format_p_value(p_value3), size = 9, hjust = 0.5)

# Combine plots with custom tags B, C, D
combined_plot <- (chloro3 + ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", vjust = 1.5, hjust = -0.5, size = 12, fontface = "bold") + 
                    chloro1 + ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", vjust = 1.5, hjust = -0.5, size = 12, fontface = "bold") + 
                    chloro2 + ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", vjust = 1.5, hjust = -0.5, size = 12, fontface = "bold")) & 
  theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15))

# Print the combined plot
print(combined_plot)

# Optionally, save the combined plot
ggsave("boxplots.svg", plot = combined_plot, width = 16, height = 8, dpi = 300)
