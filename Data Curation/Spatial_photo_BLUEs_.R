# Load necessary libraries
library(SpATS)
library(tidyver)
library(data.table)

setwd("~/Documents/PhD work/Photosyntheiss_project/Photosynthsis_paper_data/Anlysis_data/")
# Load the data
fiberblues <- fread('extremevalueremoved.csv')

# Define the list of phenotypes to analyze
phenotypes <- c("Relative_Chlorophyll", "qL", "PS1_Active.Centers", "vH.", "NPQt", "FvP_over_FmP",
                "PhiNO", "gH.", "Phi2", "PhiNPQ", "ECS_tau", "PS1_Oxidized.Centers",
                "PS1_Over.Reduced.Centers", "PS1_Open.Centers") # Ensure these column names match your data

# Define column and row knots for spatial modeling
columnKnots <- 15 # Adjust as necessary
rowKnots <- 31 
numGenotypes <- length(unique(fiberblues$Genotype)) # Ensure 'Genotypes' column matches your data

fiberblues$Day <- as.factor(fiberblues$Day)
#fiberblues$Device_ID <- as.factor(fiberblues$Device_ID)


# Initialize an empty list to store BLUEs
spatiallyCorrectedBLUEs <- data.frame(Genotype = unique(fiberblues$Genotype))

# Perform analysis for each phenotype
for (phenotype in phenotypes) {
  cat("Processing phenotype:", phenotype, "\n")

  # Fit the SpATS model
  model <- SpATS(
    response = phenotype,
    genotype = 'Genotype',
    genotype.as.random = FALSE,
    spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
    fixed = ~ Light_Intensity + Ambient_Temperature,
    random = ~ Day,
    data = fiberblues)
  pred2.model <- predict(model, which = "Genotype")
  preds = pred2.model[,c("Genotype","predicted.values")]
  colnames(preds) <- c("Genotype",phenotype)
  spatiallyCorrectedBLUEs <- merge(spatiallyCorrectedBLUEs,preds,by="Genotype")
  plot.SpATS(model, main = paste("Fitted Spatial Trend for", phenotype))
}


# Save the combined BLUEs to a CSV file
write.csv(spatiallyCorrectedBLUEs, "jan29_photo_blues.csv", row.names = FALSE)

cat("Analysis complete. Results saved to 'spatially_corrected_multispeq_BLUEs_v1.csv'.\n")








setwd('~/Documents/PhD work/Photosyntheiss_project/Photosynthsis_paper_data/Anlysis_data/')
df1 <- fread('filtered_multispeq_cutoff.csv')








df <- fread('Feb1_excluding_ambient_extreme_removed.csv')

hist(df$Relative_Chlorophyll, breaks = 40)
hist(df$ECS_tau, breaks = 40)
hist(df$gH., breaks = 40)
hist(df$vH., breaks = 40)
hist(df$Phi2, breaks = 40)
hist(df$PhiNPQ, breaks = 40)
hist(df$PhiNO, breaks = 40)
hist(df$qL, breaks = 40)
hist(df$NPQt, breaks = 40)
hist(df$FvP_over_FmP, breaks = 40)
hist(df$PS1_Active.Centers, breaks = 40)
hist(df$PS1_Open.Centers, breaks = 40)
hist(df$PS1_Over.Reduced.Centers, breaks = 40)
hist(df$PS1_Oxidized.Centers, breaks = 40)


data <- fread('Feb1_blues_.csv')

# Define the cutoff values for each trait
cutoffs <- list(
  `Relative_Chlorophyll` = c(-6, 15),
  `ECS_tau` = c(-0.002, 0.0050),
  `gH.` = c(-80, 80),
  `vH.` = c(-0.07, 0.01),
  `Phi2` = c(-0.1, 0.06),
  `PhiNPQ` = c(NA, 0.10),
  `PhiNO` = c(-0.05, 0.02),
  `qL` = c(NA, 0.16),
  `NPQt` = c(-1.0, 1.4),
  `FvP_over_FmP` = c(-0.1, 0.1), # Corrected from Fvp_over_FmP
  `PS1_Active.Centers` = c(-2.5, 1.5),
  `PS1_Open.Centers` = c(-0.25, 0.55),
   `PS1_Over.Reduced.Centers` = c(-0.55, 0.4),
  `PS1_Oxidized.Centers` = c(-0.4, 0.6)
)

# Apply the cutoff filters
filtered_data <- data

for (trait in names(cutoffs)) {
  if (trait %in% colnames(data)) { # Ensure trait exists in the dataset
    lower <- cutoffs[[trait]][1]
    upper <- cutoffs[[trait]][2]
    
    if (!is.na(lower) & !is.na(upper)) {
      filtered_data[[trait]] <- ifelse(data[[trait]] >= lower & data[[trait]] <= upper, data[[trait]], NA)
    } else if (!is.na(lower)) {
      filtered_data[[trait]] <- ifelse(data[[trait]] >= lower, data[[trait]], NA)
    } else if (!is.na(upper)) {
      filtered_data[[trait]] <- ifelse(data[[trait]] <= upper, data[[trait]], NA)
    }
  } else {
    warning(paste("Trait", trait, "not found in data. Skipping."))
  }
}



# Save the filtered data to a CSV file
write.csv(filtered_data, "Feb1_excluding_ambient_extreme_removed.csv", row.names = FALSE)

cat("Filtered data saved to 'filtered_cutoff_data.csv'.\n")









# Load necessary libraries
library(SpATS)
library(tidyverse)
library(data.table)

# Load data
multispeq <- fread('extremevalueremoved.csv')

# Define phenotype column names
phenotypes <- c("Relative_Chlorophyll", "qL", "PS1_Active.Centers", "vH.", "NPQt", "FvP_over_FmP",
                "PhiNO", "gH.", "Phi2", "PhiNPQ", "ECS_tau", "PS1_Oxidized.Centers",
                "PS1_Over.Reduced.Centers", "PS1_Open.Centers")

# Define number of knots
columnKnots <- 15
rowKnots <- 31

# Ensure correct data types
multispeq$Day <- as.factor(multispeq$Day)
#multispeq$Ambient_Temperature <- as.factor(multispeq$Ambient_Temperature)
#multispeq$Light_Intensity <- as.factor(multispeq$Light_Intensity)
# Get the number of genotypes
numGenotypes <- length(unique(multispeq$Genotype))

# Initialize a dataframe for spatially corrected BLUEs
spatiallyCorrectedBLUEs <- tibble(Genotype = unique(multispeq$Genotype))

# Function to process phenotypes and update BLUEs dataframe
process_phenotypes <- function(phenotype, data, base_df) {
  # Fit the SpATS model
  model <- SpATS(response = phenotype, 
                 genotype = 'Genotype', 
                 genotype.as.random = FALSE, 
                 spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), 
                 fixed = ~ Light_Intensity,
                 random = ~ Day, 
                 data = data)
  
  # Extract coefficients and add to the base dataframe
  coeff_df <- as.data.frame(model$coeff[1:numGenotypes])
  colnames(coeff_df) <- phenotype
  coeff_df <- as_tibble(coeff_df, rownames = 'Genotype')
  
  # Merge with the main dataframe
  updated_df <- full_join(base_df, coeff_df, by = 'Genotype')
  
  # Plot the fitted spatial trend
  plot.SpATS(model, main = phenotype)
  
  return(updated_df)
}

# Loop through all phenotypes
for (phenotype in phenotypes) {
  spatiallyCorrectedBLUEs <- process_phenotypes(phenotype, multispeq, spatiallyCorrectedBLUEs)
}


write.csv(spatiallyCorrectedBLUEs, "Feb1_blues_.csv")