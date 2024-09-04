# Install and load necessary packages
# install.packages("rMVP")
# install.packages("data.table")
library(rMVP)
library(data.table)

# Read in the command line arguments
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# Ensure command line arguments are correctly passed
if (length(args) < 1) {
  stop("Insufficient arguments provided. Please provide the trait name.")
}

# Extract arguments
trait_arg <- args[1]

# Set the working directory
setwd("/work/schnablelab/waqarali/Maizefiber/imputedvcf/eGWAS2020multispeq")

# Load phenotype data
phenotype <- read.table("mvp.plink.phe", head = TRUE)

# Load raw phenotype data
phenotypeRAW <- fread("finalBLUP-adjusted_eqtlphenotypefile.csv", data.table = F)

# Verify the structure of phenotypeRAW
str(phenotypeRAW)

# Check the column names
colnames(phenotypeRAW)

# Rename columns for clarity if needed
colnames(phenotypeRAW) <- c("Genotype", trait_arg)

# Merge phenotype data with existing phenotype information
phenotype2 <- merge(phenotypeRAW, phenotype, by.x = "Genotype", by.y = 1)

# Extract the trait of interest
trait <- phenotype2[, c("Genotype", trait_arg)]

# Verify the structure of the trait data
str(trait)
print(head(trait))

# Load genotype and kinship data
pc <- attach.big.matrix("mvp.plink.pc.desc")[]
genotype <- attach.big.matrix("mvp.plink.geno.desc")
map <- data.table::fread("mvp.plink.geno.map", data.table = F)
Kin <- attach.big.matrix("mvp.plink.kin.desc")

# Perform GWAS using MLM model
imMVP <- MVP(
  phe = trait,
  geno = genotype,
  map = map,
  # CV.MLM = pc,  # Uncomment if you want to add PCs as covariates
  nPC.MLM = 3,  # Set to 0 if not using PCs
  K = Kin,
  vc.method = "EMMA",  # Only works for MLM
  method = c("MLM"),  # Specify the MLM method
  ncpus = 16,
  maxLoop = 10,
  file.output = F,
  p.threshold = 0.05 / nrow(map)
)

# Check the structure of the imMVP object
str(imMVP)

# Extract and combine results
if (is.null(imMVP$map) || is.null(imMVP$mlm.results)) {
  stop("GWAS analysis did not produce expected results. Check the inputs and parameters.")
}

# Combine results
imMVP_results <- cbind(imMVP$map, imMVP$mlm.results)

# Verify the structure of the combined results
str(imMVP_results)
print(head(imMVP_results))

# Define output file name
output_file <- "results/BLUPs_mlm_adjusted_analysis_resultschlorophyll.csv.gz"

# Create the results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Write results to a gzipped CSV file
data.table::fwrite(imMVP_results, output_file)
