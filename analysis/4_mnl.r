# ===========================================================================
# Title: 4_mnl.r
# Date: 2026-07-20
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library & define core setting
# ------------------------------------------------------------------------

# Clear memory
rm(list = ls())

# Load packages
library(here)
library(apollo)
library(ggplot2)

# Source universal functions
source(here("utils.r"))

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = "MNL",
    modelName = "4_mnl",
    indivID = "id",
    outputDirectory = here("output")
)

# ------------------------------------------------------------------------
# SECTION 2: Load data & apply transformations
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset.csv"), header = TRUE)

# Exclude repeated task
database <- data[data$rptTask != 1, ]

# # Exclude respondents with short completion times
# database <- database[database$time >= 4, ]

# ------------------------------------------------------------------------
# SECTION 3: Define model parameters
# ------------------------------------------------------------------------

# List all parameters & starting values
apollo_beta <- c(
    ascPdt = 0,
    ascOpt = 0, # Alternative-specific constant for opt-out
    bInf_mab = 0, # Alternative-specific coefficient
    bPreg_mv = 0, # Alternative-specific coefficient
    bNoAe = 0,
    bRedness = 0,
    bPreterm_mv = 0, # Alternative-specific coefficient
    bHosp1.5 = 0,
    bHosp1.2 = 0,
    bHosp0.9 = 0,
    bHosp0.3 = 0,
    bCost0 = 0,
    bCost0.05 = 0,
    bCost0.1 = 0,
    bCost0.3 = 0,
    bCost1 = 0,
    bOptIntRisk = 0,
    bHospIntRisk = 0
)
# List parameters to be fixed at their starting values
apollo_fixed <- c("ascPdt", "bInf_mab", "bNoAe", "bHosp1.5", "bCost0")

# ------------------------------------------------------------------------
# SECTION 4: Group & validate inputs
# ------------------------------------------------------------------------

apollo_inputs <- apollo_validateInputs()

# ------------------------------------------------------------------------
# SECTION 5: Define model & likelihood function
# ------------------------------------------------------------------------

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    # Attach inputs & detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    # Create list of probabilities
    P <- list()
    
    # List utilities
    V <- list()
    V[["mv"]] <- (
        bPreg_mv * mv_preg +
        bRedness * mv_redness +
        bPreterm_mv * mv_preterm +
        (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
        bCost0.05 * mv_cost0.05 +
        bCost0.1 * mv_cost0.1 +
        bCost0.3 * mv_cost0.3 +
        bCost1 * mv_cost1
    )
    V[["mab"]] <- (
        bRedness * mab_redness +
        (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
        bCost0.05 * mab_cost0.05 +
        bCost0.1 * mab_cost0.1 +
        bCost0.3 * mab_cost0.3 +
        bCost1 * mab_cost1
    )
    V[["opt"]] <- ascOpt + bOptIntRisk * risk

    # Define MNL model settings
    mnl_settings <- list(
        alternatives = c(mv = 1, mab = 2, opt = 3),
        avail = 1, # All alternatives available
        choiceVar = choice_stage2,
        utilities = V
    )

    # Compute probabilities using MNL model
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)

    # Take product across observation for same individual
    P <- apollo_panelProd(P, apollo_inputs, functionality)

    # Prepare & return outputs of function
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
}

# ------------------------------------------------------------------------
# SECTION 6: Estimate model
# ------------------------------------------------------------------------

fit <- apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
)

# ------------------------------------------------------------------------
# SECTION 7: Output model results
# ------------------------------------------------------------------------

# Show 2-sided p values for coefficients
modelOutput_settings <- list(printPVal = 2)
saveOutput_settings <- list(printPVal = 2)

# Output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# Save model output to WD
apollo_saveOutput(fit, saveOutput_settings)

# ------------------------------------------------------------------------
# SECTION 8: Plot model coefficients
# ------------------------------------------------------------------------

# Read CSV
result <- read.csv(here("output", "4_mnl_estimates.csv"))

# Prepare & plot line plot
mnl <- linePrep(result)
linePlot(mnl)
