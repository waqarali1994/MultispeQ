

#1 simple model heritability

library(lme4)
library(data.table)

setwd("~/Documents/PhD work/Multispeq project/multispeq")

data <- fread('averaged_phenotypic_data.csv')

# Define the list of phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers',
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')
# Define the function to calculate heritability
calculate_heritability <- function(data, trait, genotype_col = "Genotypes") {
  # Fit the mixed model for the given trait
  formula <- as.formula(paste(trait, "~ (1|", genotype_col, ")"))
  model <- lmer(formula, data = data)
  # Extract variance components
  var_comp <- as.data.frame(VarCorr(model))
  # Extract genotype variance
  var_genotype <- var_comp[var_comp$grp == genotype_col, "vcov"]
  # Extract residual variance
  var_residual <- attr(VarCorr(model), "sc")^2
  # Calculate heritability
  heritability <- var_genotype / (var_genotype + var_residual / 2)
  return(heritability)
}
# Assuming your dataset is named 'data' and has columns 'Genotypes' and traits
heritability_results <- sapply(phenotypes, function(trait) {
  calculate_heritability(data, trait)
})
# Convert the results to a data frame for better readability
heritability_df <- data.frame(Trait = phenotypes, Heritability = heritability_results)
# Print the heritability results
print(heritability_df)
write.csv(heritability_df,"simple_model_heritability.csv")





#2 SpATS BLUPs based heritability in lme4

library(lme4)
library(data.table)

setwd("~/Documents/PhD work/Multispeq project/multispeq")

data <- fread('cross_validating_BLUPs.csv')

# Define the list of phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers',
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')
# Define the function to calculate heritability
calculate_heritability <- function(data, trait, genotype_col = "Genotypes") {
  # Fit the mixed model for the given trait
  formula <- as.formula(paste(trait, "~ (1|", genotype_col, ")"))
  model <- lmer(formula, data = data)
  # Extract variance components
  var_comp <- as.data.frame(VarCorr(model))
  # Extract genotype variance
  var_genotype <- var_comp[var_comp$grp == genotype_col, "vcov"]
  # Extract residual variance
  var_residual <- attr(VarCorr(model), "sc")^2
  # Calculate heritability
  heritability <- var_genotype / (var_genotype + var_residual / 2)
  return(heritability)
}
# Assuming your dataset is named 'data' and has columns 'Genotypes' and traits
heritability_results <- sapply(phenotypes, function(trait) {
  calculate_heritability(data, trait)
})


# Convert the results to a data frame for better readability
heritability_df <- data.frame(Trait = phenotypes, Heritability = heritability_results)
# Print the heritability results
print(heritability_df)
write.csv(heritability_df,"spats_blups_based_heritability.csv")



#3 lme4 heriatbility including additional factors



setwd("~/Documents/PhD work/Photosyntheiss_project/Photosynthsis_paper_data/Anlysis_data/")
# Load required libraries
library(lme4)
library(data.table)



library(data.table)

# Load the dataset
data <- fread('extremevalueremoved.csv')

# Average the dataset by PlotID
averaged_data <- data[, lapply(.SD, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1]),
                      by = PlotID]



# Save the averaged dataset
fwrite(averaged_data, "averaged_dataset_by_PlotID.csv")



library(lme4)
library(data.table)

# Load the dataset
data <- fread('extremevalueremoved.csv')

# Ensure relevant columns are factors
data$Day <- as.factor(data$Day)
data$Row <- as.factor(data$Row)
data$Column <- as.factor(data$Column)

# Define the list of phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers',
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')

# Define the function to calculate heritability
calculate_heritability <- function(data, trait) {
  tryCatch({
    # Construct the formula with Light_Intensity as a fixed effect
    model <- lmer(
      as.formula(paste0(trait, " ~ Light_Intensity + Ambient_Temperature + (1|Genotypes) + 
                        (1|Row) + (1|Day)")),
      data = data
    )
    
    class(model)
    
    # Extract variance components
    var_comp <- as.data.frame(VarCorr(model))
    
    # Extract genotype variance
    var_genotype <- var_comp[var_comp$grp == "Genotypes", "vcov"]
    
    # Extract residual variance
    var_residual <- attr(VarCorr(model), "sc")^2
    
    # Calculate heritability
    heritability <- var_genotype / (var_genotype + var_residual / 2)
    
    return(heritability)
  }, error = function(e) {
    # Return NA if an error occurs
    message(paste("Error with trait:", trait, "-", e$message))
    return(NA)
  })
}

# Calculate heritability for each phenotype
heritability_results <- sapply(phenotypes, function(trait) {
  calculate_heritability(data, trait)
})

# Convert the results to a data frame for better readability
heritability_df <- data.frame(Trait = phenotypes, Heritability = heritability_results)

# Print the heritability results
print(heritability_df)

# Save the results to a CSV file
write.csv(heritability_df, "lme4_heritability_additional_.csv", row.names = FALSE)






#4 SpATS only spatail effects heriatbility

library(SpATS)
library(tidyverse)
library(data.table)


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


# Initialize an empty data frame to store BLUEs and heritability
spatiallyCorrectedBLUEs <- data.frame(Genotypes = unique(fiberblues$Genotypes))
heritabilityResults <- data.frame(Phenotype = character(), Heritability = numeric(), stringsAsFactors = FALSE)

# Perform analysis for each phenotype
for (phenotype in phenotypes) {
  cat("Processing phenotype:", phenotype, "\n")
  
  # Fit the SpATS model
  model <- SpATS(
    response = phenotype,
    genotype = 'Genotypes',
    genotype.as.random = TRUE,
    spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), data = fiberblues)
  
  # Extract predicted values for genotypes
  pred2.model <- predict(model, which = "Genotypes")
  preds = pred2.model[,c("Genotypes","predicted.values")]
  colnames(preds) <- c("Genotypes",phenotype)
  spatiallyCorrectedBLUEs <- merge(spatiallyCorrectedBLUEs, preds, by="Genotypes")
  
  # Calculate heritability using getHeritability
  h2 <- getHeritability(model)
  heritabilityResults <- rbind(heritabilityResults, data.frame(Phenotype = phenotype, Heritability = h2))
  
  # Plot the spatial trend
  plot.SpATS(model, main = paste("Fitted Spatial Trend for", phenotype))
}

# Save the combined BLUEs and heritability results to CSV files
write.csv(spatiallyCorrectedBLUEs, "SpATS_spatially_corrected_multispeq_blups.csv", row.names = FALSE)
write.csv(heritabilityResults, "SpATS_heritability_results_only_sptial.csv", row.names = FALSE)

cat("Analysis complete. Results saved to 'spatially_corrected_multispeq_blues.csv' and 'heritability_results.csv'.\n")












#5 SpATS including additional factors heriatbility

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

# Initialize an empty data frame to store BLUEs and heritability
spatiallyCorrectedBLUEs <- data.frame(Genotypes = unique(fiberblues$Genotypes))
heritabilityResults <- data.frame(Phenotype = character(), Heritability = numeric(), stringsAsFactors = FALSE)

# Perform analysis for each phenotype
for (phenotype in phenotypes) {
  cat("Processing phenotype:", phenotype, "\n")
  
  # Fit the SpATS model
  model <- SpATS(
    response = phenotype,
    genotype = 'Genotypes',
    genotype.as.random = TRUE,
    spatial = ~SAP(Column, Row, nseg = c(columnKnots, rowKnots)), 
    fixed = ~ Light_Intensity + Ambient_Temperature,
    random = ~ Day,
    data = fiberblues)
  
  # Extract predicted values for genotypes
  pred2.model <- predict(model, which = "Genotypes")
  preds = pred2.model[,c("Genotypes","predicted.values")]
  colnames(preds) <- c("Genotypes",phenotype)
  spatiallyCorrectedBLUEs <- merge(spatiallyCorrectedBLUEs, preds, by="Genotypes")
  
  # Calculate heritability using getHeritability
  h2 <- getHeritability(model)
  heritabilityResults <- rbind(heritabilityResults, data.frame(Phenotype = phenotype, Heritability = h2))
  
  # Plot the spatial trend
  plot.SpATS(model, main = paste("Fitted Spatial Trend for", phenotype))
}

# Save the combined BLUEs and heritability results to CSV files
write.csv(spatiallyCorrectedBLUEs, "spatially_corrected_multispeq_blups.csv", row.names = FALSE)
write.csv(heritabilityResults, "SpATS_heritability_results_additional_factors.csv", row.names = FALSE)

cat("Analysis complete. Results saved to 'spatially_corrected_multispeq_blues.csv' and 'heritability_results.csv'.\n")
