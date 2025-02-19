library(SpATS)
library(tidyverse)

setwd("~/Documents")
# Load the data
multispeq <- read.csv("averaged_dataset_by_PlotID.csv")

# Convert Day to a factor
multispeq$Day <- as.factor(multispeq$Day)

phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers', 
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')

columnKnots <- 15
rowKnots <- 31
numPlots <- length(unique(multispeq$PlotID))

# Get unique Plot IDs and their corresponding Genotype IDs, converting PlotID to character
plot_genotype <- multispeq %>% 
  select(PlotID, Genotypes) %>% 
  distinct() %>%
  mutate(PlotID = as.character(PlotID))

# Perform analysis for the first phenotype
model <- SpATS(response = phenotypes[1], genotype = 'PLotID',
               genotype.as.random = FALSE, 
               spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
               fixed = ~ Light_Intensity + Ambient_Temperature,
               random = ~ Day,
               data = multispeq)

# Extract BLUPs and include Genotype IDs, ensuring PlotID is character
spatiallyCorrectedBLUPs <- as.data.frame(model$coeff[1:numPlots])
colnames(spatiallyCorrectedBLUPs) <- phenotypes[1]
spatiallyCorrectedBLUPs <- as_tibble(spatiallyCorrectedBLUPs, rownames = 'PlotID') %>%
  mutate(PlotID = as.character(PlotID)) %>%  # Ensure PlotID is character
  left_join(plot_genotype, by = 'PlotID')  # Join with genotype IDs

plot.SpATS(model, main = phenotypes[1])

# Loop through the remaining phenotypes
for(i in 2:length(phenotypes)) {
  model <- SpATS(response = phenotypes[i], genotype = 'PlotID', 
                 genotype.as.random = FALSE,
                 spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), 
                 fixed = ~ Light_Intensity + Ambient_Temperature,
                 random = ~ Day,  # Ensure Day is included as a random factor
                 data = multispeq)
  
  # Extract BLUPs and include Genotype IDs
  df <- as.data.frame(model$coeff[1:numPlots])
  colnames(df) <- phenotypes[i]
  df <- as_tibble(df, rownames = 'PlotID') %>%
    mutate(PlotID = as.character(PlotID)) %>%  # Ensure PlotID is character
    left_join(plot_genotype, by = 'PlotID')  # Join with genotype IDs
  
  # Merge with the main BLUPs data frame
  spatiallyCorrectedBLUPs <- full_join(spatiallyCorrectedBLUPs, df, by = c('PlotID', 'Genotypes'))
  
  plot.SpATS(model, main = phenotypes[i])
}

# Save the BLUPs as a CSV
write.csv(as.data.frame(spatiallyCorrectedBLUPs), "cross_validating_BLUPs.csv", row.names = FALSE)












library(SpATS)
library(tidyverse)

setwd("~/Documents")

# Load the data
multispeq <- read.csv("extremevalueremoved.csv")

# Convert Day to a factor
multispeq$Day <- as.factor(multispeq$Day)

phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers', 
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')

columnKnots <- 15
rowKnots <- 31

# Get unique Plot IDs and their corresponding Genotype IDs
plot_genotype <- multispeq %>% 
  select(PlotID, Genotypes) %>% 
  distinct() %>%
  mutate(PlotID = as.character(PlotID))

# Initialize an empty data frame for storing results
spatiallyCorrectedBLUEs <- NULL

# Loop through all phenotypes
for (phenotype in phenotypes) {
  model <- SpATS(
    response = phenotype, 
    genotype = 'PlotID', 
    genotype.as.random = FALSE, # Treat genotype as fixed for BLUEs
    spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), 
    fixed = ~ Light_Intensity + Ambient_Temperature, 
    random = ~ Day,  # Include Day as a random factor
    data = multispeq
  )
  
  # Extract BLUEs for the current phenotype
  df <- as.data.frame(model$coeff[grep("PlotID", rownames(model$coeff))])
  colnames(df) <- phenotype
  df <- as_tibble(df, rownames = 'PlotID') %>%
    mutate(PlotID = sub("PlotID_", "", PlotID)) %>% # Remove prefix if present in PlotID
    left_join(plot_genotype, by = 'PlotID')
  
  # Merge with the main BLUEs data frame
  if (is.null(spatiallyCorrectedBLUEs)) {
    spatiallyCorrectedBLUEs <- df
  } else {
    spatiallyCorrectedBLUEs <- full_join(spatiallyCorrectedBLUEs, df, by = c('PlotID', 'Genotypes'))
  }
  
  # Plot the model results
  plot.SpATS(model, main = phenotype)
}

# Save the BLUEs as a CSV
write.csv(as.data.frame(spatiallyCorrectedBLUEs), "cross_validating_BLUEs.csv", row.names = FALSE)



























#only plot level values through Coef


library(SpATS)
library(tidyverse)

setwd("~/Documents")
# Load the data
multispeq <- read.csv("extremevalueremoved.csv")

# Convert Day to a factor
multispeq$Day <- as.factor(multispeq$Day)

phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers', 
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')

columnKnots <- 15
rowKnots <- 31

# Initialize an empty data frame to store fitted values
spatiallyCorrectedValues <- multispeq

# Perform analysis for each phenotype
for (i in seq_along(phenotypes)) {
  # Fit the SpATS model
  model <- SpATS(response = phenotypes[i], genotype = 'PlotID', 
                 genotype.as.random = TRUE,
                 spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), 
                 fixed = ~ Light_Intensity + Ambient_Temperature,
                 random = ~ Day,  # Ensure Day is included as a random factor
                 data = multispeq)
  
  # Extract spatially corrected fitted values for the phenotype
  spatiallyCorrectedValues[[paste0("Fitted_", phenotypes[i])]] <- coef(model)
}

# Save the updated dataset with fitted values as new columns
write.csv(spatiallyCorrectedValues, "spatially_corrected_with_traits_coef.csv", row.names = FALSE)


setwd('~/Documents/PhD work/Photosyntheiss_project/Photosynthsis_paper_data/Anlysis_data/')

data <- fread('spatially_corrected_with_traits_coef.csv')



library(dplyr)

# Assuming your data frame is named 'data'
cleaned_averaged_data <- data %>%
  # Group by PlotID and Genotypes
  group_by(PlotID, Genotypes) %>%
  # Calculate the mean for all numeric columns
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Save the cleaned and averaged data to a CSV file
write.csv(cleaned_averaged_data, "cleaned_averaged_data.csv", row.names = FALSE)












#predict only plot level


# Load necessary libraries
library(SpATS)

# Read the input data
multispeq <- read.csv("extremevalueremoved.csv")

# Convert Day to a factor
multispeq$Day <- as.factor(multispeq$Day)

# Define phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers', 
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')

# Define spatial knot settings
columnKnots <- 15
rowKnots <- 31

# Initialize an empty data frame to store results
spatiallyCorrectedBLUEs <- unique(multispeq[, c("PlotID", "Genotypes")])

# Perform analysis for each phenotype
for (phenotype in phenotypes) {
  # Fit the SpATS model
  model <- SpATS(
    response = phenotype,
    genotype = 'PlotID',  # Predicting for PlotID
    genotype.as.random = FALSE,
    spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)),
    fixed = ~ Light_Intensity + Ambient_Temperature,
    random = ~ Day,
    data = multispeq
  )
  
  # Predict values for PlotID
  pred2.model <- predict(model, which = "PlotID")  # Predict for PlotID directly
  preds <- pred2.model[, c("PlotID", "predicted.values")]  # Predictions for PlotID
  
  # Rename columns
  colnames(preds) <- c("PlotID", phenotype)
  
  # Merge predictions with the main data frame
  spatiallyCorrectedBLUEs <- merge(spatiallyCorrectedBLUEs, preds, by = "PlotID", all.x = TRUE)
  
  # Optional: Plot the spatial trend for each phenotype
  plot.SpATS(model, main = paste("Fitted Spatial Trend for", phenotype))
}

# Save the updated dataset with predicted values and Genotypes as new columns
write.csv(spatiallyCorrectedBLUEs, "spatially_corrected_with_traits_predict_v1.csv", row.names = FALSE)