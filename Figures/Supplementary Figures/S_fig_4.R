library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)

df=read.csv('extremevalueremoved.csv')
colnames(df)

df=df[,c("Light_Intensity","Day","Phi2","PhiNPQ","Relative_Chlorophyll","qL", "PhiNO")]

light=df[,-2]
day=df[,-1]
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(linewidth = 1.4, color = "black"),
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

# light=light %>%
#   pivot_longer(!Light_Intensity, names_to = 'Trait', values_to = 'value')
# 
# colnames(light)
p1=ggplot(data=light, aes(x=Light_Intensity, y=Phi2)) +
  geom_point(color='#1444BD')+
  theme_classic()+
  ylab(bquote(Phi*PSII))+
  geom_smooth(method = "lm", color='#000000', se=F)+
  xlab('Light Intensity (µmol m-2 s-1)')+
  common_theme
p1
p2=ggplot(data=light, aes(x=Light_Intensity, y=PhiNPQ)) +
  geom_point(color='#17BECF')+
  theme_classic()+
  ylab(bquote(Phi*NPQ))+
  geom_smooth(method = "lm", color='#000000', se=F)+
  xlab('Light Intensity (µmol m-2 s-1)')+
  common_theme

p3=ggplot(data=light, aes(x=Light_Intensity, y=Relative_Chlorophyll)) +
  geom_point(color='#2CA02C')+
  theme_classic()+
  geom_smooth(method = "lm", color='#000000', se=F)+
  ylab('Relative Chlorophyll')+
  xlab('Light Intensity (µmol m-2 s-1)')+
  common_theme

p4=ggplot(data=light, aes(x=Light_Intensity, y=qL)) +
  geom_point(color='#FF7F0E')+
  theme_classic()+
  geom_smooth(method = "lm", color='#000000', se=F)+
  ylab('qL')+
  xlab('Light Intensity (µmol m-2 s-1)')+
  common_theme

p5=ggplot(data=light, aes(x=Light_Intensity, y=PhiNO)) +
  geom_point(color='#E377C2')+
  theme_classic()+
  geom_smooth(method = "lm", color='#000000', se=F)+
  ylab(bquote(Phi*NO))+
  xlab('Light Intensity (µmol m-2 s-1)')+
  common_theme

(p3+p4+p5)/(p1+p2)







# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)

# Read data
df = read.csv('extremevalueremoved.csv')
colnames(df)

# Select specific columns
df = df[, c("Light_Intensity", "Day", "Phi2", "PhiNPQ", "Relative_Chlorophyll", "qL", "PhiNO")]

light = df[, -2]
day = df[, -1]

# Common theme for plots
common_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(size = 22, color = "black"),
    axis.text.y = element_text(size = 22, color = "black"),
    axis.title = element_text(size = 22, color = "black"),
    axis.line = element_line(linewidth = 1.4, color = "black"),
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

# Function to calculate R-squared
calculate_r2 <- function(data, x, y) {
  model <- lm(as.formula(paste(y, "~", x)), data = data)
  summary(model)$r.squared
}

# Calculate R-squared values for each trait
r2_values <- list(
  Phi2 = calculate_r2(light, "Light_Intensity", "Phi2"),
  PhiNPQ = calculate_r2(light, "Light_Intensity", "PhiNPQ"),
  Relative_Chlorophyll = calculate_r2(light, "Light_Intensity", "Relative_Chlorophyll"),
  qL = calculate_r2(light, "Light_Intensity", "qL"),
  PhiNO = calculate_r2(light, "Light_Intensity", "PhiNO")
)

# Plot 1: Phi2
p1 <- ggplot(data = light, aes(x = Light_Intensity, y = Phi2)) +
  geom_point(color = '#1444BD') +
  theme_classic() +
  ylab(bquote(Phi*PSII)) +
  geom_smooth(method = "lm", color = '#000000', se = FALSE) +
  xlab('Light Intensity (µmol m-2 s-1)') +
  annotate("text", x = max(light$Light_Intensity) * 0.8, y = max(light$Phi2) * 0.9, 
           label = paste("R² =", round(r2_values$Phi2, 2)), size = 6, color = "black") +
  common_theme
p1
# Plot 2: PhiNPQ
p2 <- ggplot(data = light, aes(x = Light_Intensity, y = PhiNPQ)) +
  geom_point(color = '#17BECF') +
  theme_classic() +
  ylab(bquote(Phi*NPQ)) +
  geom_smooth(method = "lm", color = '#000000', se = FALSE) +
  xlab('Light Intensity (µmol m-2 s-1)') +
  annotate("text", x = min(light$Light_Intensity+500) * 0.8, y = max(light$PhiNPQ) * 0.9, 
           label = paste("R² =", round(r2_values$PhiNPQ, 2)), size = 6, color = "black") +
  common_theme

p2
# Plot 3: Relative Chlorophyll
p3 <- ggplot(data = light, aes(x = Light_Intensity, y = Relative_Chlorophyll)) +
  geom_point(color = '#2CA02C') +
  theme_classic() +
  geom_smooth(method = "lm", color = '#000000', se = FALSE) +
  ylab('Relative Chlorophyll') +
  xlab('Light Intensity (µmol m-2 s-1)') +
  annotate("text", x = max(light$Light_Intensity) * 0.8, y = min(light$Relative_Chlorophyll, na.rm=T) * 0.9, 
           label = paste("R² =", round(r2_values$Relative_Chlorophyll, 3)), size = 6, color = "black") +
  common_theme
p3
# Plot 4: qL
p4 <- ggplot(data = light, aes(x = Light_Intensity, y = qL)) +
  geom_point(color = '#FF7F0E') +
  theme_classic() +
  geom_smooth(method = "lm", color = '#000000', se = FALSE) +
  ylab('qL') +
  xlab('Light Intensity (µmol m-2 s-1)') +
  annotate("text", x = max(light$Light_Intensity) * 0.8, y = max(light$qL, na.rm = T) * 0.9, 
           label = paste("R² =", round(r2_values$qL, 2)), size = 6, color = "black") +
  common_theme
p4
# Plot 5: PhiNO
p5 <- ggplot(data = light, aes(x = Light_Intensity, y = PhiNO)) +
  geom_point(color = '#E377C2') +
  theme_classic() +
  geom_smooth(method = "lm", color = '#000000', se = FALSE) +
  ylab(bquote(Phi*NO)) +
  xlab('Light Intensity (µmol m-2 s-1)') +
  annotate("text", x = min(light$Light_Intensity+500) * 0.8, y = min(light$PhiNO, na.rm=T) * 0.9, 
           label = paste("R² =", round(r2_values$PhiNO, 2)), size = 6, color = "black") +
  common_theme
p5
# Combine plots
combined_plot <- (p3 + p4 + p5) / (p1 + p2)

print(combined_plot)

# Optionally, save the combined plot
ggsave("linear_unlinear.svg", plot = combined_plot, width = 20, height = 13, dpi = 300)
