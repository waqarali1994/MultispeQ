library(data.table)
library(ggplot2)
library(patchwork)
library(hms)

# Load the dataset
data <- fread("Nebraska2020MultiSpeQPostFlowering.csv")

# Ensure 'Day' is properly formatted
data$Day <- as.Date(data$Day, "%m/%d/%y")  # Convert if necessary

# Define colors for each day
day_colors <- c("July-23" = "darkred", "July-24" = "darkgreen", "July-25" = "navyblue", "July-28" = "darkorange4")

# Define a common theme
common_theme_with_ticks <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.2),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(1, "lines"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.title = element_text(hjust = 0.5, color = "black")
  )

# Function to create a plot for a given day
create_plot <- function(data, day, formatted_day, show_y_label = TRUE, show_x_label = FALSE, custom_y_breaks = c(0, 1000, 2000)) {
  
  # Filter data for the specific day
  data_filtered <- data[data$Day == as.Date(day), ]
  
  # Convert 'Day' to formatted factor for color mapping
  data_filtered$FormattedDay <- factor(formatted_day, levels = names(day_colors))
  
  # Convert Time from decimal hours to hms format
  data_filtered$Time <- hms::as_hms(data_filtered$Time * 3600)
  
  # Define time breaks
  time_breaks <- hms::as_hms(c("09:00:00", "11:00:00", "13:00:00", "14:30:00"))
  
  p <- ggplot(data_filtered, aes(x = Time, y = Light_Intensity, color = FormattedDay)) +
    geom_line(size = 1) +
    scale_color_manual(values = day_colors) +  # Correct color mapping
    scale_y_continuous(
      breaks = custom_y_breaks, 
      limits = c(0, 2000),
      labels = as.character(custom_y_breaks)
    ) +
    scale_x_time(
      limits = c(hms::as_hms("09:00:00"), hms::as_hms("14:30:00")),
      breaks = time_breaks,
      labels = scales::time_format("%I:%M %p")
    ) +
    ylab("Light Intensity") +
    common_theme_with_ticks
  
  if (!show_y_label) {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (!show_x_label) {
    p <- p + theme(axis.title.x = element_blank())
  } else {
    p <- p + xlab("Time (hr)")
  }
  
  return(p)
}

# Create plots for different days
plot_23 <- create_plot(data, "2020-07-23", "July-23", show_y_label = TRUE, show_x_label = FALSE)
plot_24 <- create_plot(data, "2020-07-24", "July-24", show_y_label = TRUE, show_x_label = FALSE)
plot_25 <- create_plot(data, "2020-07-25", "July-25", show_y_label = TRUE, show_x_label = FALSE, custom_y_breaks = c(0, 1000, 2000))
plot_28 <- create_plot(data, "2020-07-28", "July-28", show_y_label = TRUE, show_x_label = TRUE)

# Annotate plots with labels A, B, C, D
plot_23 <- plot_23 + annotate("text", x = -Inf, y = Inf, label = "A", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")
plot_24 <- plot_24 + annotate("text", x = -Inf, y = Inf, label = "B", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")
plot_25 <- plot_25 + annotate("text", x = -Inf, y = Inf, label = "C", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")
plot_28 <- plot_28 + annotate("text", x = -Inf, y = Inf, label = "D", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

# Combine the plots vertically
combined_plot <- plot_23 + plot_24 + plot_25 + plot_28 +
  plot_layout(ncol = 1, guides = 'collect') &
  common_theme_with_ticks &
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 22, color = "black")
  )

# Display the plot
print(combined_plot)

# Save the plot
ggsave("LI.svg", plot = combined_plot, width = 18, height = 12, dpi = 300)

