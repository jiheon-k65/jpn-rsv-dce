# ===========================================================================
# Title: 7_mxl_costCont
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

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = "MXL; continuous cost",
    modelName = "7_mxl_costCont",
    indivID = "id",
    outputDirectory = here("output"),
    nCores = 12
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
    mu_ascOpt = 0,
    sigma_ascOpt = 0,
    bInf_mab = 0,
    mu_bPreg_mv = 0,
    sigma_bPreg_mv = 0,
    bNoAe = 0,
    mu_bRedness = 0,
    sigma_bRedness = 0,
    mu_bPreterm_mv = 0,
    sigma_bPreterm_mv = 0,
    bHosp15 = 0,
    mu_bHosp1.2 = 0,
    sigma_bHosp1.2 = 0,
    mu_bHosp0.9 = 0,
    sigma_bHosp0.9 = 0,
    mu_bHosp0.3 = 0,
    sigma_bHosp0.3 = 0,
    mu_bCost = 0,
    sigma_bCost = 0,
    bOptIntRisk = 0, # Assumed fixed effect
    bHospIntRisk = 0 # Assumed fixed effect
)
# List parameters to be fixed at their starting values
apollo_fixed <- c("ascPdt", "bInf_mab", "bNoAe", "bHosp15")

# ------------------------------------------------------------------------
# SECTION 4: Define draws & random parameters
# ------------------------------------------------------------------------

# To set parameters for generating draws
apollo_draws <- list(
    interDrawsType = "sobol", # Halton not recommended for >5
    interNDraws = 100,
    interUnifDraws = c(),
    interNormDraws = c(
        "draws_ascOpt",
        "draws_bPreg",
        "draws_bRedness",
        "draws_bPreterm",
        "draws_bHosp1.2",
        "draws_bHosp0.9",
        "draws_bHosp0.3",
        "draws_bCost"
    ),
    intraDrawsType = "",
    intraNDraws = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
)

# To create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
    randcoeff <- list()
    randcoeff[["ascOpt"]] <- mu_ascOpt + sigma_ascOpt * draws_ascOpt
    randcoeff[["bPreg_mv"]] <- mu_bPreg_mv + sigma_bPreg_mv * draws_bPreg
    randcoeff[["bRedness"]] <- mu_bRedness + sigma_bRedness * draws_bRedness
    randcoeff[["bPreterm_mv"]] <- mu_bPreterm_mv + sigma_bPreterm_mv * draws_bPreterm
    randcoeff[["bHosp1.2"]] <- mu_bHosp1.2 + sigma_bHosp1.2 * draws_bHosp1.2
    randcoeff[["bHosp0.9"]] <- mu_bHosp0.9 + sigma_bHosp0.9 * draws_bHosp0.9
    randcoeff[["bHosp0.3"]] <- mu_bHosp0.3 + sigma_bHosp0.3 * draws_bHosp0.3
    randcoeff[["bCost"]] <- mu_bCost + sigma_bCost * draws_bCost
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
    # To attach inputs & detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    # To create list of probabilities
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

    # To compute probabilities using MNL model
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)

    # To take product across observation for same individual
    P <- apollo_panelProd(P, apollo_inputs, functionality)

    # To average across inter-individual draws
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)

    # To prepare & return outputs of function
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

# To estimate model
fit <- apollo_estimate(
    apollo_beta,
    apollo_fixed,
    apollo_probabilities, 
    apollo_inputs
)

# ------------------------------------------------------------------------
# SECTION 8: Output model results
# ------------------------------------------------------------------------

# To show 2-sided p values for coefficients
modelOutput_settings <- list(printPVal = 2)
saveOutput_settings <- list(printPVal = 2)

# To output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# To save model output to WD
apollo_saveOutput(fit, saveOutput_settings)
