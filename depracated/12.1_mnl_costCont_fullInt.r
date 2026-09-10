# ===========================================================================
# Title: 12.1_mnl_costCont_fullInt.r
# Date: 2026-07-20
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library & define core setting
# ------------------------------------------------------------------------

# Clear memory
rm(list = ls())

# Load packages
library(here)
library(dplyr)
library(apollo)

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = ("MNL; continuous cost & full interactions (complete)"),
    modelName = "12.1_mnl_costCont_fullInt",
    indivID = "id",
    outputDirectory = here("output")
)

# ------------------------------------------------------------------------
# SECTION 2: Load data & apply transformations
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset_complete.csv"), header = TRUE)

# Exclude repeated task
database <- data[data$rptTask != 1, ]

# # Exclude respondents with short completion times
# database <- database[database$time >= 4, ]

# ------------------------------------------------------------------------
# SECTION 3: Define model parameters
# ------------------------------------------------------------------------

# List all parameters & starting values
apollo_beta <- c(
    ascOpt = 0,
    ascOptIntAge = 0, # Interaction between opt & ageStd
    ascOptIntEdu = 0, # Interaction between opt & eduGrp
    ascOptIntGov = 0, # Interaction between opt & govTrust
    ascOptIntVacConf = 0, # Interaction between opt & vacConf
    bPreg_mv = 0,
    bPregIntEdu_mv = 0, # Interaction between preg & eduGrp
    bPregIntGov_mv = 0, # Interaction between preg & govTrust
    bPregIntVacConf_mv = 0, # Interaction between preg & vacConf
    bRedness = 0 ,
    bPreterm_mv = 0,
    bPretermIntGov_mv = 0, # Interaction between preterm & govTrust
    bPretermIntVacConf_mv = 0, # Interaction between preterm & vacConf
    bHosp1.2 = 0,
    bHosp0.9 = 0,
    bHosp0.3 = 0,
    bCost = 0,
    bCostIntAge = 0, # Interaction between cost & ageStd
    bCostIntSib = 0, # Interaction between cost & sibGrp
    bCostIntEdu = 0, # Interaction between cost & eduGrp
    bCostIntGov = 0, # Interaction between cost & govTrust
    bOptIntRisk = 0,
    bHospIntRisk = 0
)
# List parameters to be fixed at their starting values
apollo_fixed <- c()

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
        (bPreg_mv +
        (bPregIntEdu_mv * eduGrp) +
        (bPregIntGov_mv * govTrust) +
        (bPregIntVacConf_mv * vacConf)) * mv_preg +
        bRedness * mv_redness +
        (bPreterm_mv +
        (bPretermIntGov_mv * govTrust) +
        (bPretermIntVacConf_mv * vacConf)) * mv_preterm +
        (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
        (bCost +
        (bCostIntAge * ageStd) + 
        (bCostIntSib * sibGrp) + 
        (bCostIntEdu * eduGrp) + 
        (bCostIntGov * govTrust)) * mv_cost
    )
    V[["mab"]] <- (
        bRedness * mab_redness +
        (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
        (bCost +
        (bCostIntAge * ageStd) +
        (bCostIntSib * sibGrp) +
        (bCostIntEdu * eduGrp) +
        (bCostIntGov * govTrust)) * mab_cost
    )
    V[["opt"]] <- (
        ascOpt +
        bOptIntRisk * risk +
        ascOptIntAge * ageStd +
        ascOptIntEdu * eduGrp +
        ascOptIntGov * govTrust +
        ascOptIntVacConf * vacConf
    )

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

# Estimate model
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
