# Load necessary libraries
library(ggplot2)
library(data.table)
library(dplyr)
library(patchwork)

setwd("~/Documents/PhD work/Multispeq project/multispeq/Waqar_MultiSpeQ")
# Read the heritability results data
simple_heritability <- fread("simple_model_heritability.csv")
spats_heritability <- fread("spats_hereitability.csv")

# Extract the order of traits from the simple heritability data
simple_trait_order <- simple_heritability$Trait

# Set the factor levels for the Trait column in both datasets to enforce the desired order
spats_heritability$Trait <- factor(spats_heritability$Trait, levels = simple_trait_order)
simple_heritability$Trait <- factor(simple_heritability$Trait, levels = simple_trait_order)

# Define the raw trait names and ensure they match the data
raw_trait_names <- c('Relative_Chlorophyll', 'qL', 'PS1_Active.Centers', 'vH.', 'NPQt',
                     'FvP_over_FmP', 'PhiNO', 'gH.', 'Phi2', 'PhiNPQ', 'ECS_tau',
                     'PS1_Oxidized.Centers', 'PS1_Over.Reduced.Centers', 'PS1_Open.Centers')

# Define the color mapping based on the Simple model plot
hex_colors <- c("#2CA02C", "#FF7F0E", "brown", "#D62728", "#9467BD", "#8C564B",
                "#E377C2", "#7F7F7F", "#1444BD", "#17BECF", "#321b16", "#d6a851",
                "#BCBD22", "#437f85")
colors_for_traits <- setNames(hex_colors, raw_trait_names)

# Define the labels with expressions for the traits
trait_labels <- c("Relative Chlorophyll", "qL", "PSI-ac", expression(vH^"+"), expression(NPQ[T]), "Fv'/Fm'",
                  expression(Phi*NO), expression(gH^"+"), expression(Phi*PSII),
                  expression(Phi*NPQ), expression(ECS[T]), "PSI-oxc", "PSI-orc",
                  "PSI-opc")

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


# Create the Simple heritability plot
heritability_plot_before <- ggplot(simple_heritability, aes(x = Trait, y = Heritability, fill = Trait)) +
  geom_bar(stat = "identity", color = "NA", width = 0.9) +
  geom_text(aes(label = ifelse(!is.na(Heritability), round(Heritability, 3), ""), y = Heritability + 0.02), 
            size = 8, vjust = 0, color = "black", check_overlap = TRUE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.50),
    expand = expansion(mult = c(0.06, 0.1)),
    name = "Heritability"
  ) +
  scale_fill_manual(values = colors_for_traits) +
  guides(fill = "none") +
  common_theme +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.title.x = element_blank(), # Remove x-axis title
    axis.ticks.x = element_blank()  # Remove x-axis ticks
  )

# Read data from the provided CSV file for variance
variance_multispeq <- fread("variance_components_lme4_rawdata_percentage.csv")

# Ensure the traits are ordered in the same way as the heritability plot
variance_multispeq$Trait <- factor(variance_multispeq$Trait, levels = levels(simple_heritability$Trait))

# Reorder the Factor column to have "Genotype" at the bottom and "Residual" at the top
variance_multispeq$Factor <- factor(variance_multispeq$Factor, levels = c("Residual", "Column", "Row", "Day", "Light Intensity", "Genotype"))

# Define hex color values for variance factors
hex_colors_variance <- c("gray", "orange", "black", "steelblue", "#FFDB58", "darkgreen")

# Create a named vector for colors matching the Factor levels
colors_for_factors <- setNames(hex_colors_variance, c("Residual", "Column", "Row", "Day", "Light Intensity", "Genotype"))

# Variance plot with legend items arranged in one row and x-axis labels removed
variance_plot <- ggplot(variance_multispeq, aes(x = Trait, y = Percentage, fill = Factor)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_fill_manual(
    values = colors_for_factors,
    breaks = c("Genotype", "Light Intensity", "Day", "Row", "Column", "Residual"),
    labels = c("Genotype", "Light Intensity", "Day", "Row", "Column", "Residual"),
    guide = guide_legend(nrow = 1)  # Arrange legend items in one row
  ) +
  labs(x = "", y = "Variance explained (%)", fill = "") +
  common_theme +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = "top",
    legend.spacing.x = unit(0.2, "cm"),  # Adjust spacing between legend items
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.key.size = unit(0.6, "lines"),
    legend.background = element_rect(fill = "white", color = "black", size = 0.5)  # Rectangle border
  )

heritability_plot_after <- ggplot(spats_heritability, aes(x = Trait, y = Heritability, fill = Trait)) +
  geom_bar(stat = "identity", color = "NA", width = 0.9) +
  geom_text(aes(label = ifelse(!is.na(Heritability), round(Heritability, 3), ""), y = Heritability + 0.02), 
            size = 8, vjust = 0, color = "black", check_overlap = TRUE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.50),
    expand = expansion(mult = c(0.06, 0.1)),
    name = "Heritability"
  ) +
  scale_fill_manual(values = colors_for_traits) +
  guides(fill = "none") +
  common_theme +
  scale_x_discrete(labels = trait_labels) +  # Apply custom labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
        axis.title.x = element_blank())  # Remove x-axis title



heritability_plot_before <- heritability_plot_before + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "A", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

variance_plot <- variance_plot + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

heritability_plot_after <- heritability_plot_after + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")


# Combine the plots into a single figure with bold annotations
combined_plot <- heritability_plot_before / variance_plot / heritability_plot_after +
  plot_layout(heights = c(1, 1, 1))  # Bold the plot annotations

# Display the combined plot
print(combined_plot)


# Optionally, save the combined plot
ggsave("Heritability.svg", plot = combined_plot, width = 20, height = 13, dpi = 300)

