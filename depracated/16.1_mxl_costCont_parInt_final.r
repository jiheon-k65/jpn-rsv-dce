# ===========================================================================
# Title: 16.1_mxl_costCont_parInt_final.r
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
    modelDescr = ("MXL; continuous cost & parsimonious interactions (complete)"),
    modelName = "16.1_mxl_costCont_parInt_final",
    indivID = "id",
    outputDirectory = here("output"),
    nCores = 12
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
    delta_ascOpt = 0, # Constant
    gamma_ascOptIntGov = 0, # Interaction between opt & govTrust
    gamma_ascOptIntVacConf = 0, # Interaction between opt & vacConf
    sigma_ascOpt = 0,
    mu_bPreg_mv = 0,
    sigma_bPreg_mv = 0,
    bRedness = 0 ,
    mu_bPreterm_mv = 0,
    sigma_bPreterm_mv = 0,
    bHosp1.2 = 0,
    bHosp0.9 = 0,
    bHosp0.3 = 0,
    delta_bCost = 0,
    gamma_bCostIntEdu = 0, # Interaction between cost & eduGrp
    sigma_bCost = 0,
    bOptIntRisk = 0,
    bHospIntRisk = 0
)
# List parameters to be fixed at their starting values
apollo_fixed <- c()

# ------------------------------------------------------------------------
# SECTION 4: Define draws & random parameters
# ------------------------------------------------------------------------

# Set parameters for generating draws
apollo_draws <- list(
    interDrawsType = "halton",
    interNDraws = 100,
    interUnifDraws = c(),
    interNormDraws = c(
        "draws_ascOpt",
        "draws_bPreg",
        "draws_bPreterm",
        "draws_bCost"
    ),
    intraDrawsType = "",
    intraNDraws = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
)

# Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
    randcoeff <- list()
    randcoeff[["ascOpt"]] <- (
        delta_ascOpt +
        gamma_ascOptIntGov * govTrust +
        gamma_ascOptIntVacConf * vacConf +
        sigma_ascOpt * draws_ascOpt
    )
    randcoeff[["bPreg_mv"]] <- (
        mu_bPreg_mv +
        sigma_bPreg_mv * draws_bPreg
    )
    randcoeff[["bPreterm_mv"]] <- (
        mu_bPreterm_mv +
        sigma_bPreterm_mv * draws_bPreterm
    )
    randcoeff[["bCost"]] <- (
        delta_bCost +
        gamma_bCostIntEdu * eduGrp +
        sigma_bCost * draws_bCost
    )
    return(randcoeff)
}

# ------------------------------------------------------------------------
# SECTION 5: Group & validate inputs
# ------------------------------------------------------------------------

apollo_inputs <- apollo_validateInputs()

# ------------------------------------------------------------------------
# SECTION 6: Define model & likelihood function
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

    # Average across inter-individual draws
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)

    # Prepare & return outputs of function
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
}

# ------------------------------------------------------------------------
# SECTION 7: Estimate model
# ------------------------------------------------------------------------

# Search for starting values
apollo_beta <- apollo_searchStart(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities,
    apollo_inputs
)
# Estimate model
fit <- apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities, 
    apollo_inputs
)

# ------------------------------------------------------------------------
# SECTION 8: Output model results
# ------------------------------------------------------------------------

# Show 2-sided p values for coefficients
modelOutput_settings <- list(printPVal = 2)
saveOutput_settings <- list(printPVal = 2)

# Output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# Save model output to WD
apollo_saveOutput(fit, saveOutput_settings)

# Save estimates & variance-covariance matrix to WD
saveRDS(list(
    est = fit$estimate,
    varcov = fit$varcov
), here("output", "16.1_mxl_costCont_parInt_final_results.rds"))