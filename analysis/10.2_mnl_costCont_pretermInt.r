# ===========================================================================
# Title: 10.2_mnl_costCont_pretermInt.r
# Date: 2026-10-18
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
source(here("analysis", "utils.r"))

# Initialize code
apollo_initialise()

# Set core controls
apollo_control <- list(
    modelDescr = "MNL; linear cost & preterm interactions (imputed)",
    modelName = "10.2_mnl_costCont_pretermInt",
    indivID = "id",
    outputDirectory = here("output")
)

# ------------------------------------------------------------------------
# SECTION 2: Load data & apply transformations
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset_imputed.csv"), header = TRUE)

# Exclude repeated task
data <- data[data$rptTask != 1, ]

# # Exclude respondents with short completion times
# data <- data[data$time >= 4, ]

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
    bPretermIntAge = 0, # Interaction between preterm & ageStd
    bPretermIntSib = 0, # Interaction between preterm & sibGrp
    bPretermIntEdu = 0, # Interaction between preterm & eduGrp
    bPretermIntInc = 0, # Interaction between preterm & incGrp
    bPretermIntGov = 0, # Interaction between preterm & govTrust
    bPretermIntVacImpt = 0, # Interaction between preterm & vacImpt
    bPretermIntVacSafe = 0, # Interaction between preterm & vacSafe
    bPretermIntVacEff = 0, # Interaction between preterm & vacEff
    bPretermIntVacBlf = 0 # Interaction between preterm & vacBlf
)
# List parameters to be fixed at their starting values
apollo_fixed <- c()

# ------------------------------------------------------------------------
# SECTION 4: Define model & likelihood function
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
        (bPreterm_mv +
        (bPretermIntAge * ageStd) +
        (bPretermIntSib * sibGrp) +
        (bPretermIntEdu * eduGrp) +
        (bPretermIntInc * incGrp) +
        (bPretermIntGov * govTrust) +
        (bPretermIntVacImpt * vacImpt) +
        (bPretermIntVacSafe * vacSafe) +
        (bPretermIntVacEff * vacEff) +
        (bPretermIntVacBlf * vacBlf)) * mv_preterm +
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
# SECTION 5: Estimate model
# ------------------------------------------------------------------------

# Number of imputations
M <- 25

# Initialize list to store results
est_list <- vector("list", M) # Maximum likelihood estimates
varcov_list <- vector("list", M) # Variance-covariance matrices
robvarcov_list <- vector("list", M) # Robust variance-covariance matrices
aic_list <- vector("list", M)
bic_list <- vector("list", M)

for (i in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nRunning model for imputation", i, "\n\n")
    # Load imputed dataset for current imputation
    database <- data %>% filter(imp == i)
    # Validate Apollo inputs
    apollo_inputs <- apollo_validateInputs()
    # Estimate model
    fit <- apollo_estimate(
        apollo_beta,
        apollo_fixed,
        apollo_probabilities,
        apollo_inputs
    )
    # Store estimates and variance-covariance matrix
    est_list[[i]] <- fit$estimate
    varcov_list[[i]] <- fit$varcov
    robvarcov_list[[i]] <- fit$robvarcov
    aic_list[[i]] <- fit$AIC
    bic_list[[i]] <- fit$BIC
}

# ------------------------------------------------------------------------
# SECTION 6: Pool results using Rubin's rules
# ------------------------------------------------------------------------

# Write output to text file
sink(here("output", "10.2_mnl_costCont_pretermInt_output.txt"))
cat("===========================================================================\n")
cat("Title: 10.2_mnl_costCont_pretermInt_output.txt\n")
cat("Date: 2026-10-18\n")
cat("===========================================================================\n")

pool_coeff(est_list, robvarcov_list, M)
summ_aic_bic(aic_list, bic_list)

# Switch off writing
sink()