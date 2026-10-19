# ===========================================================================
# Title: 13.1_hcm_costCont_ranInt.r
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
    modelDescr = ("HCM; continuous cost (complete) & random intercept"),
    modelName = "13.1_hcm_costCont_ranInt",
    indivID = "id",
    outputDirectory = here("output"),
    nCores = 8
)

# ------------------------------------------------------------------------
# SECTION 2: Load data & apply transformations
# ------------------------------------------------------------------------

# Load data from WD
data <- read.csv(here("data", "jpn_dataset_subset_complete.csv"), header = TRUE)

# Recode indicators to start from 1 (instead of 0)
indicator_cols <- c("govTrust", "vacImpt", "vacSafe", "vacEff", "vacBlf")
data[indicator_cols] <- data[indicator_cols] + 1

# Exclude repeated task
database <- data[data$rptTask != 1, ]

# # Exclude respondents with short completion times
# database <- database[database$time >= 4, ]

# ------------------------------------------------------------------------
# SECTION 3: Define model parameters
# ------------------------------------------------------------------------

# List all parameters & starting values
apollo_beta <- c(
    sigma_pdt = 0,
    bPreg_mv = 0,
    bRedness = 0 ,
    bPreterm_mv = 0,
    bHosp1.2 = 0,
    bHosp0.9 = 0,
    bHosp0.3 = 0,
    bCost = 0,
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
    interDrawsType = "halton",
    interNDraws = 1000,
    interUnifDraws = c(),
    interNormDraws = c("eta", "draws_pdt"),
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
    randcoeff[["alpha"]] <- sigma_pdt * draws_pdt # Random intercept for product
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

    # Likleihood of choices
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
# SECTION 7: Estimate model
# ------------------------------------------------------------------------

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

# Output model to terminal
apollo_modelOutput(fit, modelOutput_settings)

# Write custom output to text file
sink(here("output", "13.1_hcm_costCont_ranInt_output.txt"))
cat("===========================================================================\n")
cat("Title: 13.1_hcm_costCont_ranInt_output.txt\n")
cat("Date: 2026-10-18\n")
cat("===========================================================================\n")

format_cca_results(fit)

# Switch off writing
sink()