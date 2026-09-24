# ===========================================================================
# Title: 9.1_mnl_costCont_pregInt.r
# Date: 2026-09-23
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

# Source universal functions
source(here("utils.r"))

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = "MNL; continuous cost & recipient interaction (complete)",
    modelName = "9.1_mnl_costCont_pregInt",
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
    ascOpt = 0, # Alternative-specific constant for opt-out
    bPreg_mv = 0, # Alternative-specific coefficient
    bRedness = 0,
    bPreterm_mv = 0, # Alternative-specific coefficient
    bHosp1.2 = 0,
    bHosp0.9 = 0,
    bHosp0.3 = 0,
    bCost = 0, # Linear cost coefficient (per 100,000 JPY)
    bOptIntRisk = 0,
    bHospIntRisk = 0,
    bPregIntAge = 0, # Interaction between preg & ageStd
    bPregIntSib = 0, # Interaction between preg & sibGrp
    bPregIntEdu = 0, # Interaction between preg & eduGrp
    bPregIntInc = 0, # Interaction between preg & incGrp
    bPregIntGov = 0, # Interaction between preg & govTrust
    bPregIntVacImpt = 0, # Interaction between preg & vacImpt
    bPregIntVacSafe = 0, # Interaction between preg & vacSafe
    bPregIntVacEff = 0, # Interaction between preg & vacEff
    bPregIntVacBlf = 0 # Interaction between preg & vacBlf
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
        (bPregIntAge * ageStd) +
        (bPregIntSib * sibGrp) +
        (bPregIntEdu * eduGrp) +
        (bPregIntInc * incGrp) +
        (bPregIntGov * govTrust) +
        (bPregIntVacImpt * vacImpt) +
        (bPregIntVacSafe * vacSafe) +
        (bPregIntVacEff * vacEff) +
        (bPregIntVacBlf * vacBlf)) * mv_preg +
        bRedness * mv_redness +
        bPreterm_mv * mv_preterm +
        (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
        bCost * mv_cost
    )
    V[["mab"]] <- (
        bRedness * mab_redness +
        (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
        bCost * mab_cost
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

# Output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# Write custom output to text file
sink(here("output", "9.1_mnl_costCont_pregInt_output.txt"))
cat("===========================================================================\n")
cat("Title: 9.1_mnl_costCont_pregInt_output.txt\n")
cat("Date: 2026-09-23\n")
cat("===========================================================================\n")

format_cca_results(fit)

# Switch off writing
sink()