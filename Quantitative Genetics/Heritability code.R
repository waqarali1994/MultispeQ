library(lme4)
library(data.table)

setwd("~/Documents/PhD work/Multispeq project/multispeq")

data <- fread('averaged_phenotypic_data.csv')

# Define the list of phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_.mAU', 'ECS_tau', 'gH.', 'vH.',
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
write.csv(heritability_df,"simple model heritabilioty.csv")







library(lme4)
library(data.table)

setwd("~/Documents/PhD work/Multispeq project/multispeq")

data <- fread('cross_validating_BLUPs.csv')

# Define the list of phenotypes to analyze
phenotypes <- c('Relative_Chlorophyll', 'ECS_.mAU', 'ECS_tau', 'gH.', 'vH.',
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
write.csv(heritability_df,"spats_heritabilioty.csv")