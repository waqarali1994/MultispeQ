# Load necessary libraries
library(lme4)
library(dplyr)

# Load the data
mergedblups <- read.csv("extremevalueremoved.csv")

# Modify Row and Column values by appending A and B
mergedblups$Row <- paste0(mergedblups$Row, "A")
mergedblups$Column <- paste0(mergedblups$Column, "B")

# Convert Row and Column to factors if they are not already
mergedblups$Row <- as.factor(mergedblups$Row)
mergedblups$Column <- as.factor(mergedblups$Column)

# Define the list of phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_.mAU', 'ECS_tau', 'gH.', 'vH.',
                'Phi2', 'PhiNPQ', 'PhiNO', 'qL', 'NPQt', 'FvP_over_FmP',
                'PS1_Active.Centers', 'PS1_Open.Centers', 
                'PS1_Over.Reduced.Centers', 'PS1_Oxidized.Centers')

# Initialize an empty list to store variance components for all traits
varcomp_list <- list()

# Loop through each phenotype to fit a model and extract variance components
for (trait in phenotypes) {
  # Fit the mixed-effects model
  formula <- as.formula(paste(trait, "~ (1|Genotypes) + (1|Light_Intensity) + (1|Day) + (1|Row) + (1|Column)"))
  model_random <- lmer(formula, data = mergedblups)
  
  # Extract variance components
  var_comp <- as.data.frame(VarCorr(model_random))
  
  # Add trait name to the results for easier identification
  var_comp$Trait <- trait
  
  # Store variance components in the list
  varcomp_list[[trait]] <- var_comp
}

# Combine all variance components into a single data frame
combined_varcomp <- do.call(rbind, varcomp_list)

# Convert variance components to percentage of total variance
combined_varcomp <- combined_varcomp %>%
  group_by(Trait) %>%
  mutate(Variance_Percentage = (vcov / sum(vcov)) * 100)

# Print the variance components
print(combined_varcomp)

# Write variance components to a CSV file
write.csv(combined_varcomp, "variance_components.csv", row.names = FALSE)
