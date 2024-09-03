library(SpATS)
library(tidyverse)
library(data.table) # For fread function

# Load the data
multispeq <- fread('extremevalueremoved.csv')

# Define the list of phenotypes to analyze
phenotypes <- c('ECS_.mAU', 'ECS_tau', 'gH.', 'vH.', 'Phi2', 'PhiNPQ', 'PhiNO', 
                'qL', 'NPQt', 'FvP_over_FmP', 'PS1_Active.Centers', 'PS1_Open.Centers', 
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers', 'Light_Intensity', 
                'Relative_Chlorophyll') # Ensure these column names match your data

# Define column and row knots
columnKnots <- 31 # Adjust as necessary
rowKnots <- 15
numGenotypes <- length(unique(multispeq$Genotypes)) # Ensure this column matches your data

# Convert Day to a factor
multispeq$Day <- as.factor(multispeq$Day)

# Perform analysis for the first phenotype
model <- SpATS(response = phenotypes[1], genotype = 'Genotypes', genotype.as.random = FALSE, 
               spatial = ~ SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
               fixed = ~ Light_Intensity,
               random = ~ Day, 
               data = multispeq)

# Extract BLUEs and include Genotype IDs
spatiallyCorrectedBLUEs <- as.data.frame(model$coeff[1:numGenotypes])
colnames(spatiallyCorrectedBLUEs) <- phenotypes[1]
spatiallyCorrectedBLUEs <- as_tibble(spatiallyCorrectedBLUEs, rownames = 'Genotypes')

# Plot the fitted spatial trend and title the plot with the phenotype
plot.SpATS(model, main = phenotypes[1])
boxplot(multispeq$ECS_.mAU)

# Loop through the remaining phenotypes
for(i in 2:length(phenotypes)) {
  model <- SpATS(response = phenotypes[i], genotype = 'Genotypes', genotype.as.random = FALSE, 
                 spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), 
                 fixed = Light_Intensity,
                 random = ~ Day, 
                 data = multispeq)
  
  # Extract BLUEs for the current phenotype
  df <- as.data.frame(model$coeff[1:numGenotypes])
  colnames(df) <- phenotypes[i]
  df <- as_tibble(df, rownames = 'Genotypes')
  
  # Add it to the main dataframe
  spatiallyCorrectedBLUEs <- full_join(spatiallyCorrectedBLUEs, df, by = 'Genotypes')
  
  # Plot the fitted spatial trend and title the plot with the phenotype
  plot.SpATS(model, main = phenotypes[i])
}

# Save the BLUEs to a CSV file
write.csv(spatiallyCorrectedBLUEs, "spatially_corrected_BLUEs.csv", row.names = FALSE)
