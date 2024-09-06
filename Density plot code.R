library(tidyverse)
library(svglite)
library(data.table)
library(patchwork)

setwd("~/Documents/Documents/PhD work/Multispeq project/multispeq paper data")

# Read and prepare the dataset
multispq_densityplot <- fread("filtered_multispeq_cutoff.csv")
multispq_densityplot <- multispq_densityplot[-1]
if ("Genotype" %in% colnames(multispq_densityplot)) {
  multispq_densityplot[, Genotype := NULL]
}

# Define raw trait names and color mapping
raw_trait_names <- c('Relative_Chlorophyll', 'qL', 'PS1_Active.Centers', 'vH.', 'NPQt',
                     'FvP_over_FmP', 'PhiNO', 'gH.', 'Phi2', 'PhiNPQ', 'ECS_tau',
                     'PS1_Oxidized.Centers', 'PS1_Over.Reduced.Centers', 'PS1_Open.Centers')

hex_colors <- c("#2CA02C", "#FF7F0E", "brown", "#D62728", "#9467BD", "#8C564B",
                "#E377C2", "#7F7F7F", "#1444BD", "#17BECF", "#321b16", "#d6a851",
                "#BCBD22", "#437f85")
colors_for_traits <- setNames(hex_colors, raw_trait_names)

# Define the labels with expressions for the traits
trait_labels <- c("Relative Chlorophyll", "qL", "PSI-ac", expression(vH^"+"), expression(NPQ[T]), "Fv'/Fm'",
                  expression(Phi*NO), expression(gH^"+"), expression(Phi*PSII),
                  expression(Phi*NPQ), expression(ECS[T]), "PSI-oxc", "PSI-orc",
                  "PSI-opc")
names(trait_labels) <- raw_trait_names  # Associate labels with their respective trait names

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

# Initialize an empty list to store plots
plot_list <- list()

# Generate plots and store them in the list
for (trait in raw_trait_names) {  # Loop over raw_trait_names
  x_label <- trait_labels[[trait]]  # Use the custom label
  
  # Remove NA values before calculating density
  trait_data <- multispq_densityplot[[trait]]
  trait_data <- trait_data[!is.na(trait_data)]
  
  # Calculate density to determine appropriate y-axis limits and breaks
  density_data <- density(trait_data, na.rm = TRUE)
  y_max <- max(density_data$y)
  y_breaks <- seq(0, y_max, length.out = 3)  # Set three equally spaced breaks on the y-axis
  
  # Set three equally spaced breaks on the x-axis
  x_min <- min(trait_data, na.rm = TRUE)
  x_max <- max(trait_data, na.rm = TRUE)
  x_breaks <- seq(x_min, x_max, length.out = 3)
  
  plot <- ggplot(multispq_densityplot, aes(x = .data[[trait]])) +
    geom_density(fill = colors_for_traits[trait], alpha = 0.7) +
    labs(x = x_label, y = "Density", title = NULL) +
    scale_y_continuous(breaks = y_breaks, labels = scales::number_format(accuracy = 0.01)) +  # Apply custom y-axis breaks with 2 decimal places
    scale_x_continuous(breaks = x_breaks, labels = scales::number_format(accuracy = 0.01)) +  # Apply custom x-axis breaks with 2 decimal places
    common_theme  # Apply the common theme
  
  plot_list[[trait]] <- plot
}

# Combine all plots into a single panel using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 3)

# Print the combined plot
print(combined_plot)

  # Optionally, save the combined plot
  ggsave("combined_density_plots.svg", plot = combined_plot, width = 18, height = 12, dpi = 300)
