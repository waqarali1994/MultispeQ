library(dplyr)
library(ggplot2)
library(patchwork)

setwd("~/Documents/Documents/PhD work/Multispeq project/multispeq paper data")
# Load the datasets
chlorophyll_data <- read.csv('Michael cholorphyll.csv')
genotype_data <- read.csv('transposed_gene_genotype_data.csv')

# Merge datasets on 'Genotype_ID'
merged_datav1 <- merge(chlorophyll_data, genotype_data, by = "GenotypeID")
# Filter out rows with missing 'CHL' values
merged_datav1 <- merged_datav1 %>% filter(!is.na(CHL))

# Load the datasets for the second plot
chlorophyll_data2 <- read.csv('filterted_allele_file.csv')

# Merge datasets on 'Genotype_ID'
merged_datav2 <- merge(chlorophyll_data2, genotype_data, by = "GenotypeID")

# Load the datasets for the third plot
chlorophyll_data3 <- read.csv('merged_multispeq.csv')

# Merge datasets on 'Genotype_ID'
merged_datav3 <- merge(chlorophyll_data3, genotype_data, by = "GenotypeID")

# Custom colors
custom_colors <- c("A|A" = "#FFA500", "G|G" = "steelblue")

# Define a common theme with increased tick length and legend text size
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),  # Increase tick length
    axis.ticks = element_line(size = 1.2),  # Increase tick width
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(1, "lines"),
    legend.spacing.x = unit(0.5, "cm"),
    legend.text = element_text(size = 22),  # Set legend text size to 16
    plot.title = element_text(hjust = 0.5, color = "black")
  )


# First plot
t_test_result1 <- t.test(CHL ~ Genotype, data = merged_datav1)
p_value1 <- t_test_result1$p.value
chloro1 <- ggplot(merged_datav1, aes(x = Genotype, y = CHL, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 1.5) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(breaks = c(300,400,500,600)) + 
  labs(title = "", x = "", y = "Chlorophyll") +
  common_theme +
  annotate("text", x = 1.5, y = max(merged_datav1$CHL, na.rm = TRUE) + 0.5, 
           label = ifelse(p_value1 < 0.00001, "p < 0.00001", paste("p =", round(p_value1, 4))), 
           size = 9, hjust = 0.5)

# Second plot
t_test_result2 <- t.test(Zm00001eb377130 ~ Genotype.y, data = merged_datav2)
p_value2 <- t_test_result2$p.value

chloro2 <- ggplot(merged_datav2, aes(x = Genotype.y, y = Zm00001eb377130, fill = Genotype.y)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 1.5) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(breaks = c(-20, 0, 25, 50)) + 
  labs(title = "", x = "", y = "TPM") +
  common_theme +
  annotate("text", x = 1.5, y = max(merged_datav2$Zm00001eb377130, na.rm = TRUE) + 0.5, 
           label = ifelse(p_value2 < 0.00001, "p < 0.00001", paste("p =", round(p_value2, 4))), 
           size = 9, hjust = 0.5)

# Third plot
t_test_result3 <- t.test(Relative_Chlorophyll ~ Genotype, data = merged_datav3)
p_value3 <- t_test_result3$p.value

chloro3 <- ggplot(merged_datav3, aes(x = Genotype, y = Relative_Chlorophyll, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, lwd = 1.5) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.6) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(breaks = c(-20,-10, 0, 10)) + 
  labs(title = "", x = "", y = "Relative Chlorophyll") +
  common_theme +
  annotate("text", x = 1.5, y = max(merged_datav3$Relative_Chlorophyll, na.rm = TRUE) + 0.5, 
           label = ifelse(p_value3 < 0.00001, "p < 0.00001", paste("p =", round(p_value3, 4))), 
           size = 9, hjust = 0.5)

# Combine plots with custom tags B, C, D
combined_plot <- (chloro3 + ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", vjust = 1.5, hjust = -0.5, size = 12, fontface = "bold") + 
                    chloro1 + ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", vjust = 1.5, hjust = -0.5, size = 12, fontface = "bold") + 
                    chloro2 + ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", vjust = 1.5, hjust = -0.5, size = 12, fontface = "bold")) & 
  theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15))


print(combined_plot)

# Optionally, save the combined plot
ggsave("boxplots.svg", plot = combined_plot, width = 16, height = 8, dpi = 300)
