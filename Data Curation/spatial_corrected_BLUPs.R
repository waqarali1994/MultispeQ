library(SpATS)
library(tidyverse)

setwd("~/Documents")
# Load the data
multispeq <- read.csv("extremevalueremoved.csv")

# Convert Day to a factor
multispeq$Day <- as.factor(multispeq$Day)

phenotypes <- c('Relative_Chlorophyll', 'ECS_.mAU', 'ECS_tau', 'gH.', 'vH.',
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
model <- SpATS(response = phenotypes[1], genotype = 'PlotID',
               genotype.as.random = TRUE, 
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
                 genotype.as.random = TRUE,
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

