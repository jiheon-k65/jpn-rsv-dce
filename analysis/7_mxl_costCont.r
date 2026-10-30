# ===========================================================================
# Title: 7_mxl_costCont
# Date: 2026-10-27
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library & define core setting
# ------------------------------------------------------------------------

# Clear memory
rm(list = ls())

# Load packages
library(here)
library(apollo)

# Source universal functions
source(here("analysis", "utils.r"))

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = "MXL; linear cost",
    modelName = "7_mxl_costCont",
    indivID = "id",
    outputDirectory = here("output"),
    nCores = 8
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
    mu_opt = 0,
    sigma_opt = 0,
    mu_preg_mv = 0,
    sigma_preg_mv = 0,
    mu_redness = 0,
    sigma_redness = 0,
    mu_preterm_mv = 0,
    sigma_preterm_mv = 0,
    mu_hosp1.2 = 0,
    sigma_hosp1.2 = 0,
    mu_hosp0.9 = 0,
    sigma_hosp0.9 = 0,
    mu_hosp0.3 = 0,
    sigma_hosp0.3 = 0,
    ln_mu_cost = 0,
    ln_sigma_cost = 0,
    mu_optIntRisk = 0,
    sigma_optIntRisk = 0,
    mu_hospIntRisk = 0,
    sigma_hospIntRisk = 0
)
# List parameters to be fixed at their starting values
apollo_fixed <- c()

# ------------------------------------------------------------------------
# SECTION 4: Define draws & random parameters
# ------------------------------------------------------------------------

# To set parameters for generating draws
apollo_draws <- list(
    interDrawsType = "sobol", # Halton not recommended for >5
    interNDraws = 1000,
    interUnifDraws = c(),
    interNormDraws = c(
        "draws_opt",
        "draws_preg",
        "draws_redness",
        "draws_preterm",
        "draws_hosp1.2",
        "draws_hosp0.9",
        "draws_hosp0.3",
        "draws_cost",
        "draws_optIntRisk",
        "draws_hospIntRisk"
    ),
    intraDrawsType = "",
    intraNDraws = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
)

# To create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
    randcoeff <- list()
    randcoeff[["ascOpt"]] <- mu_opt + sigma_opt * draws_opt
    randcoeff[["bPreg_mv"]] <- mu_preg_mv + sigma_preg_mv * draws_preg
    randcoeff[["bRedness"]] <- mu_redness + sigma_redness * draws_redness
    randcoeff[["bPreterm_mv"]] <- mu_preterm_mv + sigma_preterm_mv * draws_preterm
    randcoeff[["bHosp1.2"]] <- mu_hosp1.2 + sigma_hosp1.2 * draws_hosp1.2
    randcoeff[["bHosp0.9"]] <- mu_hosp0.9 + sigma_hosp0.9 * draws_hosp0.9
    randcoeff[["bHosp0.3"]] <- mu_hosp0.3 + sigma_hosp0.3 * draws_hosp0.3
    randcoeff[["bCost"]] <- - exp(ln_mu_cost + ln_sigma_cost * draws_cost) # Lognormal distribution (ensures -ve coefficient)
    randcoeff[["bOptIntRisk"]] <- mu_optIntRisk + sigma_optIntRisk * draws_optIntRisk
    randcoeff[["bHospIntRisk"]] <- mu_hospIntRisk + sigma_hospIntRisk * draws_hospIntRisk
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

# To output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# Write custom output to text file
sink(here("output", "7_mxl_costCont_output.txt"))
cat("===========================================================================\n")
cat("Title: 7_mxl_costCont_output.txt\n")
cat("Date: 2026-10-27\n")
cat("===========================================================================\n")

format_cca_results(fit)

cat("\nNote:\n")
cat("Negative scale estimates may arise due to non-centered parameterization")

# Switch off writing
sink()