library(stringr)
library(dplyr)
library(ggplot2)
library(wesanderson)

setwd("~/Documents/PhD work/Multispeq project/multispeq/Waqar_MultiSpeQ")
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

# Load data
df <- read.csv("multispeqcutoff_final.csv")
df$CHROM <- as.numeric(gsub('chr', '', df$CHROM))
df <- df[sample(1:nrow(df)), ]

# Calculate cumulative base pair positions for each chromosome
data_cum <- df %>%
  group_by(CHROM) %>%
  summarise(max_bp = max(POS)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(CHROM, bp_add)

# Combine cumulative data with GWAS data
gwas_data <- df %>%
  inner_join(data_cum, by = "CHROM") %>%
  mutate(bp_cum = POS + bp_add)

# Check if the 'Trait' column exists in gwas_data
if (!"Trait" %in% colnames(gwas_data)) {
  stop("The 'Trait' column does not exist in gwas_data.")
}

# Print summary of gwas_data to inspect its structure
print(summary(gwas_data))

# Define the raw and formatted trait names in the correct order
raw_trait_names <- c("Relative_Chlorophyll", "qL", "PS1_Active.Centers", "vH.", "NPQt", "FvP_over_FmP",
                     "PhiNO", "gH.", "Phi2", "PhiNPQ", "ECS_tau", "PS1_Oxidized.Centers",
                     "PS1_Over.Reduced.Centers", "PS1_Open.Centers")

formatted_trait_names <- c("Relative Chlorophyll", "qL", "PS1 Active Centers", "vH+", "NPQt", "FvP/FmP",
                           "PhiNO", "gH+", "Phi2", "PhiNPQ", "ECS tAU", "PS1 Oxidized Centers",
                           "PS1 Over Reduced Centers", "PS1 Open Centers")

# Create a named vector for color mapping
hex_colors <- c("#2CA02C", "#FF7F0E", "brown", "#D62728", "#9467BD", "#8C564B",
                "#E377C2", "#7F7F7F", "#1444BD", "#17BECF", "#321b16", "#d6a851",
                "#BCBD22", "#437f85")

colors_for_traits <- setNames(hex_colors, formatted_trait_names)

# Create a mapping from raw to formatted trait names
trait_name_mapping <- setNames(formatted_trait_names, raw_trait_names)

# Check unique values in gwas_data$Trait
unique_phenotypes <- unique(gwas_data$Trait)
print(unique_phenotypes)

# Check which phenotypes are not in the mapping
missing_traits <- setdiff(unique_phenotypes, names(trait_name_mapping))
print(missing_traits)

# Ensure all phenotypes are included in the mapping
if (length(missing_traits) > 0) {
  stop("The following traits are missing in the trait_name_mapping: ", paste(missing_traits, collapse = ", "))
}

# Map the raw trait names in gwas_data to the formatted trait names
gwas_data$Formatted_Trait <- trait_name_mapping[gwas_data$Trait]

# Ensure 'Formatted_Trait' is a factor and set the levels as they should appear in the plot
gwas_data$Formatted_Trait <- factor(gwas_data$Formatted_Trait, levels = formatted_trait_names)

# Check to ensure the number of colors matches the number of unique traits in your data
if (length(hex_colors) != length(unique(gwas_data$Formatted_Trait))) {
  stop("The number of hex colors must match the number of unique traits.")
}

# Define the labels with expressions for the traits using uppercase Phi symbol
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

# Calculate the cumulative positions for chromosome rectangles
chr_limits <- gwas_data %>%
  group_by(CHROM) %>%
  summarise(start = min(bp_cum), end = max(bp_cum))

# Calculate axis settings for the plot
axis_set <- gwas_data %>%
  group_by(CHROM) %>%
  summarize(center = mean(bp_cum))

# Initialize the ggplot object
p <- ggplot() +
  scale_x_continuous(breaks = axis_set$center, labels = paste0("Chr", sprintf("%02d", axis_set$CHROM))) +
  scale_y_continuous(expand = c(0, 0.02), limits = c(0, 0.4)) +
  common_theme  # Apply the common theme with ticks here

# Loop over the traits to add points to the plot
for (i in seq_along(formatted_trait_names)) {
  trait_data <- gwas_data %>%
    filter(Formatted_Trait == formatted_trait_names[i])
  if (nrow(trait_data) > 0) {
    p <- p + geom_point(data = trait_data,
                        aes(x = bp_cum, y = support, color = Formatted_Trait),
                        alpha = 0.9, size = 4)
  }
}

# Add horizontal lines for threshold levels
threshold_value <- 0.2
p <- p + geom_hline(yintercept = threshold_value, color = "black", linetype = "dashed", lwd = 0.5)

# Add the color scale for 'Trait' outside the loop
p <- p + scale_color_manual(values = colors_for_traits, name = "Trait",
                            breaks = formatted_trait_names, labels = trait_labels)

# Label specific hits based on RMIP significance
label_traits <- c("qL", "Relative Chlorophyll", "Phi2", "PhiNPQ")
trait_numbers <- c("Phi2" = "1", "PhiNPQ" = "2", "qL" = "3", "Relative Chlorophyll" = "4")
rmip_threshold <- 0.2

# Filter the data for the specific traits and RMIP above the threshold, prioritize by RMIP
label_data <- gwas_data %>%
  filter(Formatted_Trait %in% names(trait_numbers) & support > rmip_threshold) %>%
  arrange(desc(support)) %>%
  group_by(bp_cum) %>%
  mutate(Label_Number = case_when(
    Formatted_Trait == "Phi2" & CHROM == 8 ~ "E",
    Formatted_Trait == "Phi2" & CHROM == 9 ~ "B",
    Formatted_Trait == "PhiNPQ" & CHROM == 3 ~ "F",
    Formatted_Trait == "PhiNPQ" & CHROM == 9 ~ "C",
    Formatted_Trait == "qL" ~ "D",
    Formatted_Trait == "Relative Chlorophyll" ~ "G",
    TRUE ~ NA_character_
  ))

# Add labels to the plot using the numbering, with increased font size and bold text
p <- p +
  geom_text(data = label_data, aes(x = bp_cum, y = support, label = Label_Number),
            vjust = -1, hjust = 0.5, size = 6, fontface = "bold", color = "black")

# Final touches to the plot with increased legend text size
p <- p + 
  labs(x = "", y = "RMIP") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  common_theme  # Use the common theme here, ensuring consistency

# Print the plot with numbered labels
print(p)




setwd("~/Documents/PhD work/Multispeq project/multispeq/Waqar_MultiSpeQ")

# Load necessary libraries
library(data.table)
library(ggplot2)
library(patchwork)


# Define a common theme with increased tick length
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
    plot.title = element_text(hjust = 0.5, color = "black")
  )


chr9_94775951 <- fread("LD_values/LD_chr9_94775951.ld", data.table = F)
#chr9.1 <- chr9.1[-which(chr9.1$SNP_B == "chr9_21755970"),]

#head(chr9_94775951)

#summary(chr9_94775951)

hitPOS.94775951 <- 94775951


minW.1 <- min(chr9_94775951$BP_B[which(chr9_94775951$R2 >= 0.6)])
maxw.1 <- max(chr9_94775951$BP_B[which(chr9_94775951$R2 >= 0.6)])


### SNPs track
SNPS.1 <- chr9_94775951 %>%
  filter(BP_B > minW.1-50000 & BP_B < maxw.1+ 50000)
head(SNPS.1)

### Gene track
genes.desc <- fread("B73_geneModels_v5.csv", data.table = F)
# colnames(genes.desc)

genes.desc.1 <- genes.desc %>%
  filter(chr==9, start > minW.1 &  end < maxw.1)

genes.desc.1.inside <- genes.desc %>%
  filter(chr==9, start > minW.1-width &  end < maxw.1)


genes.desc.1.inside <- genes.desc.1.inside[,c(1:4,6,8,9,14:23)]

#third.hitPOS.window2$nothing <- NA
#third.hitPOS.window2$nothing[grepl("Zm00001eb300700", third.hitPOS.window2$B73_V5)] <- "blue4"
#third.hitPOS.window2$nothing[!grepl("Zm00001eb300700", third.hitPOS.window2$B73_V5)] <- "grey"

#summary(SNPS.1)

SNPS.1$R2comp <- NA
SNPS.1$R2comp[which(SNPS.1$R2 >= 0.6)] <- "High LD"
SNPS.1$R2comp[which(SNPS.1$R2 < 0.6)] <- "low LD"
SNPS.1$R2comp[which(SNPS.1$R2 == 1)] <- "associated SNP"

plot1 <- ggplot() + 
  geom_point(data=SNPS.1[SNPS.1$R2comp == "associated SNP", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, stroke = 0, size = 5, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.1[SNPS.1$R2comp == "High LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.1[SNPS.1$R2comp == "low LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = c("blue", "red", "grey"),
                    breaks = c("associated SNP", "High LD", "low LD")) +
  geom_point(data=genes.desc.1.inside[genes.desc.1.inside$strand=="+",],
             aes(x=((start+end)/2)/1000000, y=1.1), colour="black" , shape="\u25BA", size=6) +
  geom_point(data=genes.desc.1.inside[genes.desc.1.inside$strand=="-",], 
             aes(x=((start+end)/2/1000000), y=1.05), colour="black", shape="\u25C4", size=6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #scale_color_manual(values = c("blue4", "grey")) +
  #geom_hline(yintercept = -log10(0.05/nrow(gwas.dat)), linetype=2) + 
  #geom_vline(xintercept = minW.1/1000000, linetype="dotted", color = "black", linewidth=1) +
  #geom_vline(xintercept = maxw.1/1000000, linetype="dotted", color = "black", linewidth=1) +
  common_theme +
  ylab(expression(R^2)) + 
  xlab(paste0("Chromosome 9 (Mb)")) +
  theme(legend.position = "none")

plot1

#ggplot_build(plot1)$layout$panel_params[[1]]$x.range
#summary(SNPS.1)

# Gene top 2
chr9_25683346 <- fread("LD_values/LD_chr9_25683346.ld", data.table = F)
hitPOS.2 <- 25683346

minW.2 <- min(chr9_25683346$BP_B[which(chr9_25683346$R2 >= 0.6)])
maxw.2 <- max(chr9_25683346$BP_B[which(chr9_25683346$R2 >= 0.6)])

SNPS.2 <- chr9_25683346 %>%
  filter(BP_B > minW.2 - 50000 & BP_B < maxw.2 + 50000)

genes.desc.2  <- genes.desc %>%
  filter(chr == 9, start > minW.2 - 50000 & end < maxw.2 + 50000)
genes.desc.2.inside <- genes.desc %>%
  filter(chr == 9, start > minW.2 - 1000 & end < maxw.2)

SNPS.2$R2comp <- NA
SNPS.2$R2comp[which(SNPS.2$R2 >= 0.6)] <- "High LD"
SNPS.2$R2comp[which(SNPS.2$R2 < 0.6)] <- "low LD"
SNPS.2$R2comp[which(SNPS.2$R2 == 1)] <- "associated SNP"

plot2 <- ggplot() + 
geom_point(data=SNPS.2[SNPS.2$R2comp == "associated SNP", ],
           aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
           shape = 21, stroke = 0, size = 5, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.2[SNPS.2$R2comp == "High LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.2[SNPS.2$R2comp == "low LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = c("#17BECF", "red", "grey"),
                    breaks = c("associated SNP", "High LD", "low LD")) +
  geom_point(data=genes.desc.2[genes.desc.2$strand=="+",],
             aes(x = ((start + end) / 2) / 1000000, y = 1.1), 
             colour = "black", shape = "\u25BA", size = 6) +
  geom_point(data=genes.desc.2[genes.desc.2$strand=="-",], 
             aes(x = ((start + end) / 2) / 1000000, y = 1.05), 
             colour = "black", shape = "\u25C4", size = 6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #geom_vline(xintercept = minW.2 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  #geom_vline(xintercept = maxw.2 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  common_theme +
  ylab(expression(R^2)) + 
  xlab(paste0("Chromosome 9 (Mb)")) +
  theme(legend.position = "none")

plot2

ggplot_build(plot2)$layout$panel_params[[1]]$x.range
summary(SNPS.2)


# Gene top 3
chr10_83741616 <- fread("LD_values/LD_chr10_83741616.ld", data.table = F)
hitPOS.3 <- 83741616

minW.3 <- min(chr10_83741616$BP_B[which(chr10_83741616$R2 >= 0.6)])
maxw.3 <- max(chr10_83741616$BP_B[which(chr10_83741616$R2 >= 0.6)])

SNPS.3 <- chr10_83741616 %>%
  filter(BP_B > minW.3 - 50000 & BP_B < maxw.3 + 50000)

genes.desc.3 <- genes.desc %>%
  filter(chr == 10, start > minW.3 - 50000 & end < maxw.3 + 50000)
genes.desc.3.inside <- genes.desc %>%
  filter(chr == 10, start > minW.3 - 1000 & end < maxw.3)

SNPS.3$R2comp <- NA
SNPS.3$R2comp[which(SNPS.3$R2 >= 0.6)] <- "High LD"
SNPS.3$R2comp[which(SNPS.3$R2 < 0.6)] <- "low LD"
SNPS.3$R2comp[which(SNPS.3$R2 == 1)] <- "associated SNP"

plot3 <- ggplot() + 
  geom_point(data=SNPS.3[SNPS.3$R2comp == "associated SNP", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, stroke = 0, size = 5, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.3[SNPS.3$R2comp == "High LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.3[SNPS.3$R2comp == "low LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = c("#FF7F0E", "red", "grey"),
                    breaks = c("associated SNP", "High LD", "low LD")) +
  geom_point(data=genes.desc.3[genes.desc.3$strand=="+",],
             aes(x = ((start + end) / 2) / 1000000, y = 1.1), 
             colour = "black", shape = "\u25BA", size = 6) +
  geom_point(data=genes.desc.3[genes.desc.3$strand=="-",], 
             aes(x = ((start + end) / 2) / 1000000, y = 1.05), 
             colour = "black", shape = "\u25C4", size = 6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #geom_vline(xintercept = minW.3 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  #geom_vline(xintercept = maxw.3 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  common_theme +
  ylab(expression(R^2)) + 
  xlab(paste0("Chromosome 10 (Mb)")) +
  theme(legend.position = "none")

plot3

# Gene top 4
chr8_102200025 <- fread("LD_values/LD_chr8_102200025.ld", data.table = F)
hitPOS.4 <- 102200025

minW.4 <- min(chr8_102200025$BP_B[which(chr8_102200025$R2 >= 0.6)])
maxw.4 <- max(chr8_102200025$BP_B[which(chr8_102200025$R2 >= 0.6)])

SNPS.4 <- chr8_102200025 %>%
  filter(BP_B > minW.4 - 100000 & BP_B < maxw.4 + 100000)

genes.desc.4 <- genes.desc %>%
  filter(chr == 8, start > minW.4 - 1000 & end < maxw.4 + 1000)
genes.desc.4.inside <- genes.desc %>%
  filter(chr == 8, start > minW.4 - 1000 & end < maxw.4)

SNPS.4$R2comp <- NA
SNPS.4$R2comp[which(SNPS.4$R2 >= 0.6)] <- "High LD"
SNPS.4$R2comp[which(SNPS.4$R2 < 0.6)] <- "low LD"
SNPS.4$R2comp[which(SNPS.4$R2 == 1)] <- "associated SNP"

plot4 <- ggplot() + 
  geom_point(data=SNPS.4[SNPS.4$R2comp == "associated SNP", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, stroke = 0.5, size = 5, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.4[SNPS.4$R2comp == "High LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.4[SNPS.4$R2comp == "low LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = c("blue", "red", "grey"),
                    breaks = c("associated SNP", "High LD", "low LD")) +
  geom_point(data=genes.desc.4[genes.desc.4$strand=="+",],
             aes(x = ((start + end) / 2) / 1000000, y = 1.1), 
             colour = "black", shape = "\u25BA", size = 6) +
  geom_point(data=genes.desc.4.inside[genes.desc.4$strand=="-",], 
             aes(x = ((start + end) / 2) / 1000000, y = 1.05), 
             colour = "black", shape = "\u25C4", size = 6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  #geom_vline(xintercept = minW.4 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
 # geom_vline(xintercept = maxw.4 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  common_theme +
  ylab(expression(R^2)) + 
  xlab(paste0("Chromosome 8 (Mb)")) +
  theme(legend.position = "none")
plot4

# Gene top 5
chr3_149998190 <- fread("LD_values/LD_chr3_149998190.ld", data.table = F)
hitPOS.5 <- 149998190

minW.5 <- min(chr3_149998190$BP_B[which(chr3_149998190$R2 >= 0.6)])
maxw.5 <- max(chr3_149998190$BP_B[which(chr3_149998190$R2 >= 0.6)])

SNPS.5 <- chr3_149998190 %>%
  filter(BP_B > minW.5 - 50000 & BP_B < maxw.5 + 50000)

genes.desc.5 <- genes.desc %>%
  filter(chr == 3, start > minW.5 - 1000 & end < maxw.5 + 1000)
genes.desc.5.inside <- genes.desc %>%
  filter(chr == 3, start > minW.5 - 1000 & end < maxw.5)

SNPS.5$R2comp <- NA
SNPS.5$R2comp[which(SNPS.5$R2 >= 0.6)] <- "High LD"
SNPS.5$R2comp[which(SNPS.5$R2 < 0.6)] <- "low LD"
SNPS.5$R2comp[which(SNPS.5$R2 == 1)] <- "associated SNP"

plot5 <- ggplot() + 
  geom_point(data=SNPS.5[SNPS.5$R2comp == "associated SNP", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, stroke = 0, size = 5, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.5[SNPS.5$R2comp == "High LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.5[SNPS.5$R2comp == "low LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = c("#17BECF", "red", "grey"),
                    breaks = c("associated SNP", "High LD", "low LD")) +
  geom_point(data=genes.desc.5.inside[genes.desc.5.inside$strand=="+",],
             aes(x = ((start + end) / 2) / 1000000, y = 1.1), 
             colour = "black", shape = "\u25BA", size = 6) +
  geom_point(data=genes.desc.5.inside[genes.desc.5.inside$strand=="-",], 
             aes(x = ((start + end) / 2) / 1000000, y = 1.05), 
             colour = "black", shape = "\u25C4", size = 6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
 # geom_vline(xintercept = minW.5 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  #geom_vline(xintercept = maxw.5 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  common_theme +
  ylab(expression(R^2)) + 
  xlab(paste0("Chromosome 3 (Mb)")) +
  theme(legend.position = "none")
plot5

# Gene top 6
chr9_21755970 <- fread("LD_values/LD_chr9_21755970.ld", data.table = F)
hitPOS.6 <- 21755970

minW.6 <- min(chr9_21755970$BP_B[which(chr9_21755970$R2 >= 0.6)])
maxw.6 <- max(chr9_21755970$BP_B[which(chr9_21755970$R2 >= 0.6)])

SNPS.6 <- chr9_21755970 %>%
  filter(BP_B > minW.6 - 70000 & BP_B < maxw.6 + 70000)

genes.desc.6 <- genes.desc %>%
  filter(chr == 9, start > minW.6 - 1000 & end < maxw.6)
genes.desc.6.inside <- genes.desc %>%
  filter(chr == 9, start > minW.6 - 1000 & end < maxw.6)

SNPS.6$R2comp <- NA
SNPS.6$R2comp[which(SNPS.6$R2 >= 0.6)] <- "High LD"
SNPS.6$R2comp[which(SNPS.6$R2 < 0.6)] <- "low LD"
SNPS.6$R2comp[which(SNPS.6$R2 == 1)] <- "associated SNP"

plot6 <- ggplot() + 
  geom_point(data=SNPS.6[SNPS.6$R2comp == "associated SNP", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21,stroke = 0, size = 5, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.6[SNPS.6$R2comp == "High LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  geom_point(data=SNPS.6[SNPS.6$R2comp == "low LD", ],
             aes(x = BP_B / 1000000, y = R2, fill = R2comp), 
             shape = 21, color = "black", stroke = 0.5, size = 3, alpha = 0.9) +  # Adjust size here
  scale_y_continuous(breaks = c(0, 0.5, 1)) + 
  scale_fill_manual(values = c("#2CA02C", "red", "grey"),
                    breaks = c("associated SNP", "High LD", "low LD")) +
  geom_point(data=genes.desc.6.inside[genes.desc.6.inside$strand=="+",],
             aes(x = ((start + end) / 2) / 1000000, y = 1.1), 
             colour = "black", shape = "\u25BA", size = 6) +
  geom_point(data=genes.desc.6.inside[genes.desc.6.inside$strand=="-",], 
             aes(x = ((start + end) / 2) / 1000000, y = 1.05), 
             colour = "black", shape = "\u25C4", size = 6) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
 # geom_vline(xintercept = minW.6 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
 # geom_vline(xintercept = maxw.6 / 1000000, linetype = "dotted", color = "black", linewidth = 1) +
  common_theme +
  ylab(expression(R^2)) + 
  xlab(paste0("Chromosome 9 (Mb)")) +
  theme(legend.position = "none")

plot6


# Assuming p, plot1, plot2, plot3, plot4, plot5, plot6 are already defined
# Manually annotate the first set of plots

p <- p + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "A", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot1 <- plot1 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot2 <- plot2 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot3 <- plot3 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

# Manually annotate the second set of plots
plot4 <- plot4 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "E", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot5 <- plot5 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "F", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

plot6 <- plot6 + 
  ggplot2::annotate("text", x = -Inf, y = Inf, label = "G", vjust = 2, hjust = -0.1, size = 12, fontface = "bold")

# Combine the plots into a grid
combined_plot <- (p) / 
  (plot1 + plot2 + plot3) / 
  (plot4 + plot5 + plot6) +
  plot_layout(ncol = 1)  # Ensures the tags are bold

# Display the combined plot
print(combined_plot)


# Optionally, save the combined plot
ggsave("GWAS.svg", plot = combined_plot, width = 18, height = 12, dpi = 300)
