library(ggplot2)
library(data.table)
library(patchwork)

setwd("~/Documents/PhD work/Multispeq project/multispeq/Waqar_MultiSpeQ")

# Define a common theme with increased tick length
common_theme_with_ticks <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"),# Increase tick length
    axis.ticks = element_line(size = 1.2),  # Increase tick width
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),  # Remove square boxes
    legend.key.size = unit(1, "lines"),
    legend.spacing.x = unit(0.5, "cm"),
    plot.title = element_text(hjust = 0.5, color = "black")
  )

# Function to create a plot for a given date and color
create_plot <- function(data, date, color, formatted_date, show_y_label = TRUE, show_x_label = FALSE, detailed_x_axis = FALSE, custom_y_breaks = c(0, 1000, 2000)) {
  data_filtered <- data[data$Date == as.Date(date), ]
  data_filtered$FormattedDate <- factor(format(data_filtered$Date, "%B-%d"), levels = formatted_date)
  
  p <- ggplot(data_filtered, aes(x = Time, y = Light_Intensity, color = FormattedDate)) +
    geom_line(size = 1) +
    scale_color_manual(values = setNames(color, formatted_date)) +
    scale_y_continuous(breaks = custom_y_breaks, labels = as.character(custom_y_breaks)) +  # Custom Y-axis with specified breaks
    ylab('Light Intensity') +
    common_theme_with_ticks
  
  if (detailed_x_axis) {
    # Detailed x-axis for plot on July 25
    p <- p + scale_x_continuous(
      breaks = c(9.5, 9.75, 10.20),  # 9:30, 9:45, 10:00 in decimal hours
      labels = c("9:30 AM", "9:45 AM", "10:00 AM")
    )
  } else {
    # Normal x-axis with times 10:00 AM, 11:00 AM, 12:00 PM, and 1:30 PM
    p <- p + scale_x_continuous(
      breaks = c(10, 11.5, 13.5),  # 10:00 AM, 11:00 AM, 12:00 PM, 1:30 PM in decimal hours
      labels = c("10:00 AM", "12:00 PM","2:00 PM")
    )
  }
  
  if (!show_y_label) {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (!show_x_label) {
    p <- p + theme(axis.title.x = element_blank())  # Remove the x-axis title
  } else {
    p <- p + xlab("Time (hr)")
  }
  
  return(p)
}

# Read the CSV file
data <- fread("Nebraska2020MultiSpeQPostFlowering.csv")
data$Date <- as.Date(data$Date, "%m/%d/%y")

# Create plots for different dates
plot_23 <- create_plot(data, "2020-07-23", "darkred", "July-23", show_y_label = TRUE, show_x_label = FALSE)
plot_24 <- create_plot(data, "2020-07-24", "darkgreen", "July-24", show_y_label = TRUE, show_x_label = FALSE)  # Normal x-axis
plot_25 <- create_plot(data, "2020-07-25", "navyblue", "July-25", show_y_label = TRUE, show_x_label = FALSE, detailed_x_axis = TRUE, custom_y_breaks = c(0, 750, 1500))  # Detailed x-axis and custom y-axis for July 25
plot_28 <- create_plot(data, "2020-07-28", "darkorange4", "July-28", show_y_label = TRUE, show_x_label = TRUE)  # Normal x-axis, keep the x-axis title

# Manually annotate the second set of plots
plot_23 <- plot_23 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "A", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot_24 <- plot_24 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot_25 <- plot_25 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot_28 <- plot_28 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

# Combine the plots vertically in a single column
combined_plot <- plot_23 + plot_24 + plot_25 + plot_28 +
  plot_layout(ncol = 1, guides = 'collect') &
  common_theme_with_ticks &
  theme(
    legend.position = 'top',
    legend.text = element_text(size = 22, color = "black"),
    plot.tag = element_text(face = "bold"),
  )

# Display the combined plot with a single legend
print(combined_plot)

# Optionally, save the combined plot
ggsave("LI.svg", plot = combined_plot, width = 18, height = 12, dpi = 300)
