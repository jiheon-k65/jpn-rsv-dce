# ===========================================================================
# Title: 18_hcm_costCont_ranInt_ranCoeff_choiceProb_mvCost.r
# Date: 2026-09-10
# ===========================================================================

# ------------------------------------------------------------------------
# SECTION 1: Load library
# ------------------------------------------------------------------------

# Load packages
library(here)
library(dplyr)
library(MASS)

# ------------------------------------------------------------------------
# SECTION 2: Load data & results
# ------------------------------------------------------------------------

data <- read.csv(here("data", "jpn_dataset_subset_imputed.csv"), header = TRUE)

results <- readRDS(here("output", "14.2_hcm_costCont_ranInt_ranCoeff_results.rds"))
est_list <- results$est
varcov_list <- results$varcov

# ------------------------------------------------------------------------
# SECTION 3: Define inputs for calculations
# ------------------------------------------------------------------------

# Attribute levels for each alternative
risk <- 0
mv_hosp1.2 <- 0
mv_hosp0.9 <- 1 # 70% efficacy for MV
mv_hosp0.3 <- 0
mv_cost_base <- 0.24 # 24,000 JPY for MV
mv_redness <- 0
mv_preterm <- 0
mv_preg <- 1 # Pregnant women as recipient for MV
mab_hosp1.2 <- 0
mab_hosp0.9 <- 0
mab_hosp0.3 <- 1 # 90% efficacy for mAb
mab_cost <- 0.45 # 45,000 JPY for mAb
mab_redness <- 0

# Cost grid for MV
cost_grid <- seq(0, 0.4, by = 0.02)

# Keep only first imputation & unique respondents
df <- data %>% filter(imp == 1) %>% distinct(id, .keep_all = TRUE)
stopifnot(nrow(df) == length(unique(df$id)))

N <- nrow(df) # No. of respondents
M <- max(data$imp) # No. of imputations
S <- 200 # No. of simulations
K <- length(cost_grid) # No. of cost values for MV

# ------------------------------------------------------------------------
# SECTION 4: Simulate parameter values
# ------------------------------------------------------------------------

# Obtain respondent characteristics
ageStd <- df$ageStd
sibGrp <- df$sibGrp
eduGrp <- df$eduGrp
incGrp <- df$incGrp

# Initialize arrays to hold probabilities
prob_mv_array <- array(NA, dim = c(S, M, K))
prob_mab_array <- array(NA, dim = c(S, M, K))
prob_opt_array <- array(NA, dim = c(S, M, K))

for (m in 1:M) {
    cat("\n===========================================================================\n")
    cat("\nCalculating probabilities for imputation", m, "\n")

    # Extract estimates & variance-covariance matrix
    theta_hat <- est_list[[m]]
    sigma <- varcov_list[[m]]

    # Draw parameter values from multivariate normal distribution
    set.seed(123 + m)
    theta <- mvrnorm(n = S, mu = theta_hat, Sigma = sigma)
    colnames(theta) <- names(theta_hat)

    # Simulate random coefficients
    for (s in 1:S) {
        # Isolate fixed coefficients
        bHosp1.2 <- theta[s, "bHosp1.2"]
        bHosp0.9 <- theta[s, "bHosp0.9"]
        bHosp0.3 <- theta[s, "bHosp0.3"]
        bRedness <- theta[s, "bRedness"]
        bHospIntRisk <- theta[s, "bHospIntRisk"]
        bOptIntLV <- theta[s, "bOptIntLV"]
        bOptIntRisk <- theta[s, "bOptIntRisk"]

        # Simulate bCost for each respondent
        ln_mu_bCost <- theta[s, "ln_mu_bCost"]
        ln_sigma_bCost <- theta[s, "ln_sigma_bCost"]
        bCost <- -rlnorm(N, mean = ln_mu_bCost, sd = abs(ln_sigma_bCost)) # Lognormal

        # Simulate bPreg_mv for each respondent
        mu_bPreg_mv <- theta[s, "mu_bPreg_mv"]
        sigma_bPreg_mv <- theta[s, "sigma_bPreg_mv"]
        bPreg_mv <- rnorm(N, mean = mu_bPreg_mv, sd = abs(sigma_bPreg_mv))

        # Simulate bPreterm_mv for each respondent
        mu_bPreterm_mv <- theta[s, "mu_bPreterm_mv"]
        sigma_bPreterm_mv <- theta[s, "sigma_bPreterm_mv"]
        bPreterm_mv <- rnorm(N, mean = mu_bPreterm_mv, sd = abs(sigma_bPreterm_mv))

        # Simulate alpha for each respondent
        sigma_pdt <- theta[s, "sigma_pdt"]
        alpha <- rnorm(N, mean = 0, sd = abs(sigma_pdt))

        # Simulate LV for each respondent
        gamma_age <- theta[s, "gamma_age"]
        gamma_sib <- theta[s, "gamma_sib"]
        gamma_edu <- theta[s, "gamma_edu"]
        gamma_inc <- theta[s, "gamma_inc"]
        eta <- rnorm(N, mean = 0, sd = 1) # Random disturbance
        LV <- (
            gamma_age * ageStd +
            gamma_sib * sibGrp +
            gamma_edu * eduGrp +
            gamma_inc * incGrp +
            eta
        )
        stopifnot(length(LV) == N)

        for (k in seq_along(cost_grid)) {
            # Update mv_cost for each iteration
            mv_cost <- cost_grid[k]

            # Calculate utilities for each alternative
            V_mv <- (
                bPreg_mv * mv_preg +
                bRedness * mv_redness +
                bPreterm_mv * mv_preterm +
                (bHosp1.2 + bHospIntRisk * risk) * mv_hosp1.2 +
                (bHosp0.9 + bHospIntRisk * risk) * mv_hosp0.9 +
                (bHosp0.3 + bHospIntRisk * risk) * mv_hosp0.3 +
                bCost * mv_cost +
                alpha
            )
            V_mab <- (
                bRedness * mab_redness +
                (bHosp1.2 + bHospIntRisk * risk) * mab_hosp1.2 +
                (bHosp0.9 + bHospIntRisk * risk) * mab_hosp0.9 +
                (bHosp0.3 + bHospIntRisk * risk) * mab_hosp0.3 +
                bCost * mab_cost +
                alpha
            )
            V_opt <- (
                bOptIntRisk * risk +
                bOptIntLV* LV
            )
            stopifnot(
                length(V_mv) == N,
                length(V_mab) == N,
                length(V_opt) == N
            )

            # Calculate choice probabilities (Dim: N x S x M x K)
            prob_total <- exp(V_mv) + exp(V_mab) + exp(V_opt)
            prob_mv <- exp(V_mv) / prob_total
            prob_mab <- exp(V_mab) / prob_total
            prob_opt <- exp(V_opt) / prob_total

            # Calculate & store means across respondents (Dim: S x M x K)
            prob_mv_array[s, m, k] <- mean(prob_mv)
            prob_mab_array[s, m, k] <- mean(prob_mab)
            prob_opt_array[s, m, k] <- mean(prob_opt)
        }
    }
}

# Initialize list to hold results for each cost value
summary_list <- vector("list", length = length(cost_grid))

for (k in seq_along(cost_grid)) {
    cat("\n===========================================================================\n")
    cat("\nCalculating mean & 95% UI for cost value", k, "\n")
    
    summary_list[[k]] <- data.frame(
        Cost = cost_grid[k],
        Alt = c("MV", "mAb", "Opt-out"),
        Mean = c(
            mean(prob_mv_array[, , k]),
            mean(prob_mab_array[, , k]),
            mean(prob_opt_array[, , k])
        ),
        Lower_UI = c(
            quantile(as.vector(prob_mv_array[, , k]), 0.025),
            quantile(as.vector(prob_mab_array[, , k]), 0.025),
            quantile(as.vector(prob_opt_array[, , k]), 0.025)
        ),
        Upper_UI = c(
            quantile(as.vector(prob_mv_array[, , k]), 0.975),
            quantile(as.vector(prob_mab_array[, , k]), 0.975),
            quantile(as.vector(prob_opt_array[, , k]), 0.975)
        )
    )
}

# Bind results into a single data frame
summary_combined <- bind_rows(summary_list) %>%
    mutate(across(Mean:Upper_UI, ~ round(.x, 5)))

# Save results as CSV
write.csv(
    summary_combined,
    file = here("data", "choiceProb_mvCost.csv"),
    row.names = FALSE
)
cat("\nFile successfully saved as choiceProb_mvCost.csv\n")