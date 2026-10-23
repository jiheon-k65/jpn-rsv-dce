# ===========================================================================
# Title: 17.2_hcm_costCont_ranInt_ranCoeff.r
# Date: 2026-10-23
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library & define core setting
# ------------------------------------------------------------------------

# Activate renv (for HPC)
renv::activate()

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
    modelDescr = ("HCM; linear cost, random intercept & random coefficients (imputed)"),
    modelName = "17.2_hcm_costCont_ranInt_ranCoeff",
    indivID = "id",
    outputDirectory = here("output"),
    nCores = 8
)

# ------------------------------------------------------------------------
# SECTION 2: Load data & apply transformations
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset_imputed.csv"), header = TRUE)

# Recode indicators to start from 1 (instead of 0)
indicator_cols <- c("govTrust", "vacImpt", "vacSafe", "vacEff", "vacBlf")
data[indicator_cols] <- data[indicator_cols] + 1

# Exclude repeated task
data <- data[data$rptTask != 1, ]

# # Exclude respondents with short completion times
# data <- data[data$time >= 4, ]

# ------------------------------------------------------------------------
# SECTION 3: Define model parameters
# ------------------------------------------------------------------------

# List all parameters & starting values
apollo_beta <- c(
    sigma_pdt = 0,
    mu_preg_mv = 0,
    sigma_preg_mv = 0,
    bRedness = 0 ,
    mu_preterm_mv = 0,
    sigma_preterm_mv = 0,
    mu_hosp1.2 = 0,
    sigma_hosp1.2 = 0,
    bHosp0.9 = 0,  # Fixed as SD was not different from zero
    mu_hosp0.3 = 0,
    sigma_hosp0.3 = 0,
    ln_mu_cost = 0,
    ln_sigma_cost = 0,
    bHospIntRisk = 0,
    bOptIntLV = 1,
    bOptIntRisk = 0,
    gamma_age = 0,
    gamma_sib = 0,
    gamma_edu = 0,
    gamma_inc = 0,
    zeta_govTrust = 0,
    zeta_vacImpt = 0,
    zeta_vacSafe = 0,
    zeta_vacEff = 0,
    zeta_vacBlf = 0,
    tau_govTrust1 = -1,
    tau_govTrust2 = 0,
    tau_govTrust3 = 1,
    tau_vacImpt1 = -1,
    tau_vacImpt2 = 0,
    tau_vacImpt3 = 1,
    tau_vacSafe1 = -1,
    tau_vacSafe2 = 0,
    tau_vacSafe3 = 1,
    tau_vacEff1 = -1,
    tau_vacEff2 = 0,
    tau_vacEff3 = 1,
    tau_vacBlf1 = -1,
    tau_vacBlf2 = 0,
    tau_vacBlf3 = 1
)
# List parameters to be fixed at their starting values
apollo_fixed <- c()

# ------------------------------------------------------------------------
# SECTION 4: Define draws & random parameters
# ------------------------------------------------------------------------

# Set parameters for generating draws
apollo_draws <- list(
    interDrawsType = "sobol",
    interNDraws = 1000,
    interUnifDraws = c(),
    interNormDraws = c(
        "eta",
        "draws_pdt",
        "draws_preg",
        "draws_preterm",
        "draws_hosp1.2",
        "draws_hosp0.3",
        "draws_cost"
    ),
    intraDrawsType = "",
    intraNDraws = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
)

# Create random parameters
apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
    randcoeff <- list()
    randcoeff[["LV"]] <- (
        gamma_age * ageStd +
        gamma_sib * sibGrp +
        gamma_edu * eduGrp +
        gamma_inc * incGrp +
        eta
    )
    randcoeff[["alpha"]] <- sigma_pdt * draws_pdt # Random intercept for product (mimics nested logit)
    randcoeff[["bPreg_mv"]] <- mu_preg_mv + sigma_preg_mv * draws_preg
    randcoeff[["bPreterm_mv"]] <- mu_preterm_mv + sigma_preterm_mv * draws_preterm
    randcoeff[["bHosp1.2"]] <- mu_hosp1.2 + sigma_hosp1.2 * draws_hosp1.2
    randcoeff[["bHosp0.3"]] <- mu_hosp0.3 + sigma_hosp0.3 * draws_hosp0.3
    randcoeff[["bCost"]] <- - exp(ln_mu_cost + ln_sigma_cost * draws_cost) # Lognormal (ensures -ve coefficient)
    return(randcoeff)
}

# ------------------------------------------------------------------------
# SECTION 5: Define model & likelihood function
# ------------------------------------------------------------------------

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    # Attach inputs & detach after function exit
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    # Create list of probabilities
    P <- list()

    # Likleihood of indicators
    ol_settings1 = list(
        outcomeOrdered = govTrust,
        utility = zeta_govTrust * LV,
        tau = list(tau_govTrust1, tau_govTrust2, tau_govTrust3),
        rows = (task==1)
    )
    ol_settings2 = list(
        outcomeOrdered = vacImpt,
        utility = zeta_vacImpt * LV,
        tau = list(tau_vacImpt1, tau_vacImpt2, tau_vacImpt3),
        rows = (task==1)
    )
    ol_settings3 = list(
        outcomeOrdered = vacSafe,
        utility = zeta_vacSafe * LV,
        tau = list(tau_vacSafe1, tau_vacSafe2, tau_vacSafe3),
        rows = (task==1)
    )
    ol_settings4 = list(
        outcomeOrdered = vacEff,
        utility = zeta_vacEff * LV,
        tau = list(tau_vacEff1, tau_vacEff2, tau_vacEff3),
        rows = (task==1)
    )
    ol_settings5 = list(
        outcomeOrdered = vacBlf,
        utility = zeta_vacBlf * LV,
        tau = list(tau_vacBlf1, tau_vacBlf2, tau_vacBlf3),
        rows = (task==1)
    )
    P[["indic_govTrust"]] <- apollo_ol(ol_settings1, functionality)
    P[["indic_govTrust"]] <- apollo_panelProd(P[["indic_govTrust"]], apollo_inputs, functionality)
    P[["indic_vacImpt"]] <- apollo_ol(ol_settings2, functionality)
    P[["indic_vacImpt"]] <- apollo_panelProd(P[["indic_vacImpt"]], apollo_inputs, functionality)
    P[["indic_vacSafe"]] <- apollo_ol(ol_settings3, functionality)
    P[["indic_vacSafe"]] <- apollo_panelProd(P[["indic_vacSafe"]], apollo_inputs, functionality)
    P[["indic_vacEff"]] <- apollo_ol(ol_settings4, functionality)
    P[["indic_vacEff"]] <- apollo_panelProd(P[["indic_vacEff"]], apollo_inputs, functionality)
    P[["indic_vacBlf"]] <- apollo_ol(ol_settings5, functionality)
    P[["indic_vacBlf"]] <- apollo_panelProd(P[["indic_vacBlf"]], apollo_inputs, functionality)

    # Likelihood of choices
    V <- list()
    V[["mv"]] <- (
        bPreg_mv * mv_preg +
        bRedness * mv_redness +
        bPreterm_mv * mv_preterm +
        (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
        bCost * mv_cost +
        alpha
    )
    V[["mab"]] <- (
        bRedness * mab_redness +
        (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
        (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
        (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
        bCost * mab_cost +
        alpha
    )
    V[["opt"]] <- (
        bOptIntLV * LV +
        bOptIntRisk * risk
    )

    # Define MNL model settings
    mnl_settings <- list(
        alternatives = c(mv = 1, mab = 2, opt = 3),
        avail = 1, # All alternatives available
        choiceVar = choice_stage2,
        utilities = V
    )

    # Compute probabilities for MNL model component
    P[["choice"]] <- apollo_mnl(mnl_settings, functionality)
    P[["choice"]] <- apollo_panelProd(P[["choice"]], apollo_inputs, functionality)

    # Likelihood of whole model
    P <- apollo_combineModels(P, apollo_inputs, functionality)

    # Average across inter-individual draws
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)

   # Prepare & return outputs of function
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
}

# ------------------------------------------------------------------------
# SECTION 6: Estimate model
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
    if (i == 1) {
        fit <- apollo_estimate(
            apollo_beta,
            apollo_fixed,
            apollo_probabilities,
            apollo_inputs
        )
    } else {
        # Use previous estimates as starting values
        fit <- apollo_estimate(
            fit$estimate,
            apollo_fixed,
            apollo_probabilities,
            apollo_inputs
        )
    }
    # Store estimates and variance-covariance matrix
    est_list[[i]] <- fit$estimate
    varcov_list[[i]] <- fit$varcov
    robvarcov_list[[i]] <- fit$robvarcov
    aic_list[[i]] <- fit$AIC
    bic_list[[i]] <- fit$BIC
}

# Save estimates and variance-covariance matrices to RDS files
saveRDS(
    list(
        est = est_list,
        varcov = varcov_list,
        robvarcov = robvarcov_list
    ),
    here("output", "17.2_hcm_costCont_ranInt_ranCoeff_results.rds")
)

# ------------------------------------------------------------------------
# SECTION 7: Pool results using Rubin's rules
# ------------------------------------------------------------------------

# Write output to text file
sink(here("output", "17.2_hcm_costCont_ranInt_ranCoeff_output.txt"))
cat("===========================================================================\n")
cat("Title: 17.2_hcm_costCont_ranInt_ranCoeff_output.txt\n")
cat("Date: 2026-10-23\n")
cat("===========================================================================\n")

pool_coeff(est_list, robvarcov_list, M)
summ_aic_bic(aic_list, bic_list)

# Switch off writing
sink()