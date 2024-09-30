library(SpATS)
library(tidyverse)
library(data.table) # For fread function

# Load the data
fiberblues <- fread('Bugeater2021_merged_v6.2_Raw.csv')

# Define the list of phenotypes to analyze
phenotypes <- c("ProteinPCT", "OilPCT", "FiberPCT", "AshPCT", "StarchPCT", "MoisturePCT", "KernelColor") # Ensure these column names match your data

# Define column and row knots
columnKnots <- 22 # Adjust as necessary
rowKnots <- 21
numGenotypes <- length(unique(fiberblues$Genotype)) # Ensure this column matches your data


# Perform analysis for the first phenotype
model <- SpATS(response = phenotypes[1], genotype = 'Genotype', genotype.as.random = FALSE, 
               spatial = ~ SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
               data = fiberblues)

# Extract BLUEs and include Genotype IDs
spatiallyCorrectedBLUEs <- as.data.frame(model$coeff[1:numGenotypes])
colnames(spatiallyCorrectedBLUEs) <- phenotypes[1]
spatiallyCorrectedBLUEs <- as_tibble(spatiallyCorrectedBLUEs, rownames = 'Genotype')

# Plot the fitted spatial trend and title the plot with the phenotype
plot.SpATS(model, main = phenotypes[1])
boxplot(fiberblues$ProteinPCT)

# Loop through the remaining phenotypes
for(i in 2:length(phenotypes)) {
  model <- SpATS(response = phenotypes[i], genotype = 'Genotype', genotype.as.random = FALSE, 
                 spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
                 data = fiberblues)
  
  # Extract BLUEs for the current phenotype
  df <- as.data.frame(model$coeff[1:numGenotypes])
  colnames(df) <- phenotypes[i]
  df <- as_tibble(df, rownames = 'Genotype')
  
  # Add it to the main dataframe
  spatiallyCorrectedBLUEs <- full_join(spatiallyCorrectedBLUEs, df, by = 'Genotype')
  
  # Plot the fitted spatial trend and title the plot with the phenotype
  plot.SpATS(model, main = phenotypes[i])
}

# Save the BLUEs to a CSV file
write.csv(spatiallyCorrectedBLUEs, "spatially_corrected_fiberBLUEs.csv", row.names = FALSE)
