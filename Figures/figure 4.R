library(ggplot2)
library(data.table)
library(patchwork)
library(hms)  # For handling times effectively

#setwd("~/Documents/PhD work/Multispeq project/multispeq/Waqar_MultiSpeQ")

# Define a common theme with increased tick length
common_theme_with_ticks <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(size = 1.4, color = "black"),
    axis.ticks.length = unit(0.4, "cm"), # Increase tick length
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
create_plot <- function(data, date, color, formatted_date, show_y_label = TRUE, show_x_label = FALSE, custom_y_breaks = c(0, 1000, 2000)) {
  data_filtered <- data[data$Date == as.Date(date), ]
  data_filtered$FormattedDate <- factor(format(data_filtered$Date, "%B-%d"), levels = formatted_date)
  
  # Convert Time from decimal hours to hms format
  data_filtered$Time <- hms::as_hms(data_filtered$Time * 3600)  # Convert decimal hours to seconds and then to hms
  
  # Set the start and end times in seconds
  start_time <- as.numeric(hms::as_hms("09:00:00"))
  end_time <- as.numeric(hms::as_hms("14:30:00"))  # Set the end time to 2:30 PM (14:30:00)
  
  # Define breaks manually to ensure 2:30 PM is included
  time_breaks <- c(
    hms::as_hms("09:00:00"),
    hms::as_hms("11:00:00"),
    hms::as_hms("13:00:00"),
    hms::as_hms("14:30:00")  # Ensure 2:30 PM is shown
  )
  
  p <- ggplot(data_filtered, aes(x = Time, y = Light_Intensity, color = FormattedDate)) +
    geom_line(size = 1) +
    scale_color_manual(values = setNames(color, formatted_date)) +
    scale_y_continuous(
      breaks = custom_y_breaks, 
      limits = c(0, 2000),  # Set fixed limits for the y-axis
      labels = as.character(custom_y_breaks)  # Ensure consistent labels
    ) +
    ylab('Light Intensity') +
    common_theme_with_ticks +
    scale_x_time(
      limits = c(hms::as_hms(start_time), hms::as_hms(end_time)),  # Set the time range from 9:30 AM to 2:30 PM
      breaks = time_breaks,  # Manually set the breaks to ensure 2:30 PM is included
      labels = scales::time_format("%I:%M %p")  # Format time with AM/PM
    )
  
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
plot_24 <- create_plot(data, "2020-07-24", "darkgreen", "July-24", show_y_label = TRUE, show_x_label = FALSE)
plot_25 <- create_plot(data, "2020-07-25", "navyblue", "July-25", show_y_label = TRUE, show_x_label = FALSE, custom_y_breaks = c(0, 1000, 2000))
plot_28 <- create_plot(data, "2020-07-28", "darkorange4", "July-28", show_y_label = TRUE, show_x_label = TRUE)

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
