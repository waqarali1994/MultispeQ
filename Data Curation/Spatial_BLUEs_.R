# Load necessary libraries
library(SpATS)
library(tidyverse)
library(data.table)

setwd("~/Documents")
# Load the data
fiberblues <- fread('extremevalueremoved.csv')

# Define the list of phenotypes to analyze
phenotypes <- c("Relative_Chlorophyll", "qL", "PS1_Active.Centers", "vH.", "NPQt", "FvP_over_FmP",
                "PhiNO", "gH.", "Phi2", "PhiNPQ", "ECS_tau", "PS1_Oxidized.Centers",
                "PS1_Over.Reduced.Centers", "PS1_Open.Centers") # Ensure these column names match your data

# Define column and row knots for spatial modeling
columnKnots <- 15 # Adjust as necessary
rowKnots <- 31
numGenotypes <- length(unique(fiberblues$Genotypes)) # Ensure 'Genotypes' column matches your data

fiberblues$Day <- as.factor(fiberblues$Day)
#fiberblues$Device_ID <- as.factor(fiberblues$Device_ID)


# Initialize an empty list to store BLUEs
spatiallyCorrectedBLUEs <- data.frame(Genotypes = unique(fiberblues$Genotypes))

# Perform analysis for each phenotype
for (phenotype in phenotypes) {
  cat("Processing phenotype:", phenotype, "\n")

  # Fit the SpATS model
  model <- SpATS(
    response = phenotype,
    genotype = 'Genotypes',
    genotype.as.random = FALSE,
    spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
    fixed = ~ Light_Intensity + Ambient_Temperature,
    random = ~ Day,
    data = fiberblues)
  pred2.model <- predict(model, which = "Genotypes")
  preds = pred2.model[,c("Genotypes","predicted.values")]
  colnames(preds) <- c("Genotypes",phenotype)
  spatiallyCorrectedBLUEs <- merge(spatiallyCorrectedBLUEs,preds,by="Genotypes")
  plot.SpATS(model, main = paste("Fitted Spatial Trend for", phenotype))
}


# Save the combined BLUEs to a CSV file
write.csv(spatiallyCorrectedBLUEs, "spatially_corrected_multispeq_blues_with_deviceid.csv", row.names = FALSE)

cat("Analysis complete. Results saved to 'spatially_corrected_multispeq_BLUEs_v1.csv'.\n")
